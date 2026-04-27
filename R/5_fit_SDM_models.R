#### SCRIPT INFO ####
# This script fits final MaxEnt models per species using a fixed set of predictors.
#
# Main steps:
# 1. Load cleaned and filtered occurrences.
# 2. Load environmental predictors (WorldClim + ENVIREM), crop and mask to area M.
# 3. For each species:
#    a. Apply spatial and/or environmental filter to occurrences.
#    b. Fit block-partitioned MaxEnt models with multiple feature classes and RM values.
#    c. Save model results, variable importance, and response curves.
#    d. Select best model(s) by multiple criteria (AICc, AUC, omission rate).
#    e. Generate and save null models for statistical significance testing.
#    f. Plot variable marginal responses.
#
# Notes:
# - This script supports both MaxEnt.jar and maxnet as backends.
# - Parallelization is used for model fitting.
# - Null models are used to assess statistical support vs. random expectations.

#### SETTINGS ####
options(java.parameters = "-Xmx16g")  # Set Java memory limit for MaxEnt.jar

#### LIBRARIES ####
library(terra)         # Raster data handling
library(ENMeval)       # MaxEnt tuning and evaluation
library(tidyverse)     # Data wrangling
library(sf)            # Vector data handling
library(doParallel)    # Parallel backend
library(ENMTools)      # Null model generation, marginal plots

#### USER DEFINED FUNCTIONS ####
source('R/udf/funcs.R')
source('R/udf/margin_plots.R')

##### SELECT VARIABLES TO INCLUDE IN FINAL MODEL #####
varsIn <- c('climaticMoistureIndex',
            'PETColdestQuarter',
            'PETWettestQuarter',
            'bio_17',
            'bio_2',
            'bio_6',
            'PETseasonality')

##### LOAD OCCURRENCES #####
# Load all cleaned and flagged occurrences
allOccs  <- read_csv('outputs/occs/all_occs_flagged.csv')

# Example: Split Lasiurus varius occurrences by season
lasiurus <- allOccs %>%
  filter(species == 'Lasiurus varius') %>%
  mutate(month = dmy(date) %>% as_date() %>% month()) %>%
  filter(!is.na(month)) %>%
  mutate(month = ifelse(month %in% c(10,11,12,1,2,3), 'summer', 'winter')) %>%
  mutate(species = str_c(species, '_', month)) %>%
  split(.$species)

# Split all occurrences by species
occsList <- allOccs %>% split(.$species)
rm(allOccs)

# Add Lasiurus seasonal splits back to list
occsList <- c(occsList, lasiurus)

##### LOAD MODELING AREA #####
aream <- read_sf('outputs/aream/new_area_m.geojson')

##### PARALLEL SETTINGS #####
nCores <- 8
setParallel <- ifelse(nCores > 1, TRUE, FALSE)

##### LOAD PREDICTORS #####
# Load and name WorldClim predictors
biopreds <- rast(list.files('inputs/worldclim/2.5m/', '.tif', full.names = TRUE))
names(biopreds) <- gsub('wc2.1_2.5m_', '', names(biopreds))

# Load and name ENVIREM predictors
envirem <- rast(list.files('/mnt/2TB/GIS/Rasters/Clima/ENVIREM/2.5 m/current/', '.bil$', full.names = TRUE))
names(envirem) <- gsub('current_2-5arcmin_', '', names(envirem))

# Combine predictors and crop/mask to Area M
preds <- c(biopreds %>% crop(aream), envirem %>% crop(aream))
preds <- terra::mask(terra::crop(preds[[varsIn]], aream), terra::vect(aream))

# Clean up
rm(biopreds, envirem)
gc()

##### MODELING SETTINGS #####
algorithm   <- 'maxent.jar'    # Options: 'maxnet', 'maxent.jar'
rasterPreds <- TRUE            # Whether to generate raster outputs
RMlists     <- c(seq(1, 4, 0.5), 5, 7, 8)  # Regularization multipliers

##### SPECIES LOOP #####
for(i in seq_along(occsList)){
  
  ##### SET FEATURE CLASSES BASED ON N RECORDS #####
  if(NROW(occsList[[i]]) > 80){
    FClist <- c('L', 'Q', 'LQ', 'H', 'LQH', 'LQHP')
  } else {
    FClist <- c('L', 'LQ', 'Q', 'LQP')
  }
  
  ##### PREPARE SPECIES-SPECIFIC DATA #####
  cat(paste0("Running model for ", names(occsList[i]), " at: ", Sys.time(), '\n'))
  
  vars <- varsIn
  m.mask <- aream
  spM1 <- preds[[vars]] %>% mask(m.mask)
  
  ##### SAMPLE BACKGROUND POINTS #####
  bg <- terra::spatSample(spM1[[1]], 10000, xy = TRUE, na.rm = TRUE)[, 1:2]
  names(bg) <- c('lon', 'lat')
  
  # Force NA in spM1 if any NA in any layer
  spM1 <- terra::app(spM1, fun = \(x) {if(sum(is.na(x)) > 0) x * NA else x})
  
  ##### SET OUTPUT DIRECTORY #####
  saveDir <- paste0('outputs/models/fitting/', names(occsList[i]))
  if (!dir.exists(saveDir)) dir.create(saveDir, recursive = TRUE)
  
  ##### APPLY OCCURRENCE FILTER #####
  occs <- occsList[[i]]
  if(all(occs$use.env.filter)){
    occs <- occs %>% filter(env.unique, spat.unique)
  } else {
    occs <- occs %>% filter(spat.unique)
  }
  
  ##### CHECK IF MODEL ALREADY EXISTS #####
  if(!file.exists(paste0(saveDir, "/block_results.csv"))){
    
    ##### FIT ENMEVAL BLOCK MODELS #####
    enm <- ENMevaluate(occs = occs[, c('lon', 'lat')],
                       envs = raster::stack(spM1),
                       bg = bg,
                       tune.args = list(fc = FClist, rm = RMlists),
                       partitions = "block",
                       algorithm = algorithm,
                       numCores = nCores,
                       doClamp = FALSE,
                       parallel = TRUE)
    
    ##### SAVE MODEL RESULTS #####
    write_csv(enm@results, paste0(saveDir, '/block_results.csv'))
    
    ##### SAVE VARIABLE IMPORTANCE #####
    imp <- eval.variable.importance(enm) %>%
      map(~.x %>% pivot_longer(cols = c('percent.contribution', 'permutation.importance'))) %>%
      enframe() %>%
      unnest(cols = 2, names_repair = 'universal')
    colnames(imp) <- c('model', 'predictor', 'variable', 'value')
    write_csv(imp, paste0(saveDir, '/block_varimportance.csv'))
    
    ##### SELECT BEST MODELS BY MULTIPLE CRITERIA #####
    cases <- c('AICc', 'auc.train', 'or.mtp.avg')
    
    for(c in cases){
      
      ##### SELECT BEST MODEL #####
      best <- enm@results %>% filter(auc.train >= quantile(auc.train)[4])
      if(c == 'auc.train'){
        best <- best %>% filter(!!ensym(c) == max(!!ensym(c)))
      } else {
        best <- best %>% filter(!!ensym(c) == min(!!ensym(c)))
      }
      if(nrow(best) > 1) best <- best %>% sample_n(1)
      
      ##### SAVE BEST MODEL SETTINGS #####
      write_csv(best, paste0(saveDir, '/block_best_', c, 'tune_args.csv'))
      
      ##### PLOT BEST MODEL PREDICTION #####
      r <- enm@predictions[[best$tune.args]]
      png(paste0(saveDir, '/block_best_', c, '.png'), width = 12, height = 15, units = 'cm', res = 300)
      plot(r, col = viridis::rocket(15))
      plot(aream, border = 'gray90', col = alpha('white', 0.1), add = TRUE, reset = FALSE, lwd = 0.5)
      points(occs[, c('lon', 'lat')], pch = '+', col = 'white', cex = 1.5)
      points(occs[, c('lon', 'lat')], pch = '+', col = 'black', cex = 1.2)
      title(main = paste0('Block best ', c, ': ', best$tune.args))
      dev.off()
      
      ##### PLOT RESPONSE CURVES #####
      png(paste0(saveDir, '/block_best_', c, '_var_response.png'), width = 17, height = 15, units = 'cm', res = 300)
      dismo::response(eval.models(enm)[[best$tune.args]])
      dev.off()
      
      ##### GENERATE NULL MODELS #####
      e <- ENMevaluate(occs = occs[, c('lon', 'lat')],
                       envs = raster::stack(spM1),
                       bg = bg,
                       tune.args = list(fc = best$fc, rm = as.double(as.character(best$rm))),
                       partitions = "block",
                       algorithm = algorithm,
                       numCores = nCores,
                       doClamp = FALSE,
                       parallel = FALSE)
      
      mod.null <- ENMnulls(e,
                           mod.settings = list(fc = best$fc, rm = as.double(as.character(best$rm))),
                           numCores = nCores,
                           removeMxTemp = TRUE,
                           no.iter = 100)
      
      ##### PLOT NULL MODEL RESULTS #####
      png(paste0(saveDir, '/block_best_', c, '_null_models.png'), width = 17, height = 15, units = 'cm', res = 300)
      plot(evalplot.nulls(mod.null, stats = c('or.10p', 'auc.val', 'auc.train'), plot.type = 'histogram'))
      dev.off()
      
      ##### SAVE NULL MODEL RESULTS #####
      write_csv(null.emp.results(mod.null), paste0(saveDir, '/block_best_', c, '_null_models.csv'))
      
      ##### PLOT MARGINAL RESPONSE CURVES #####
      dir.create(paste0(saveDir, '/response_plots'), showWarnings = FALSE)
      for(n in names(spM1)){
        png(paste0(saveDir, '/response_plots/', n, '_', c, '.png'), width = 17, height = 15, units = 'cm', res = 300)
        plot(my_margin_plot(enm@models[[best$tune.args[1]]], env = spM1, n, bg, occs[, c('lon', 'lat')]))
        dev.off()
      }
    }
    
  } else {
    cat(paste0(saveDir, "/block_results.csv already exists. Skipping to next model\n"))
  }
  
} # end species loop
