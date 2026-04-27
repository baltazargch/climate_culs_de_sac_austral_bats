#### SCRIPT INFO ####
# This script implements a Monte Carlo approach to selecting the most important
# environmental variables for MaxEnt models of species distributions, following
# and adapting Schnase et al. (2021) (doi: 10.1371/journal.pone.0237208).
#
# Main steps:
# 1. For each species, run repeated MaxEnt models using random subsets of predictors.
# 2. Aggregate variable contributions across runs to identify the most informative variables.
# 3. Periodically re-fit models with the current top N variables and track performance.
# 4. Save model logs and variable tallies across iterations ("sprints").
#
# Notes:
# - Background points are resampled at each sprint to add stochasticity.
# - ENMeval is used to tune model complexity (fc, rm) in each run.
# - This process is computationally expensive and is parallelized.
# - Sprints can be resumed in case of interruption.

# Set Java RAM allocation (only relevant if using MaxEnt.jar)
options(java.parameters = "-Xmx8g")

#### LIBRARIES ####
library(dismotools)   # Tools for MaxEnt model evaluation
library(tidyverse)    # Data wrangling
library(parallel)     # Parallelization
library(foreach)      # Parallel looping
library(ENMeval)      # MaxEnt modeling and tuning
library(terra)        # Raster handling
library(sf)           # Vector data handling

#### SETTINGS AND USER-DEFINED FUNCTIONS ####
sf_use_s2(FALSE)      # Disable spherical geometry (avoid geometry validity issues)
walk(list.files('R/udf', pattern = '.R$', full.names = TRUE), source)  # Load UDFs
source('R/udf/envSamp.R')   # Load environmental sampling function

#### LOAD BASE DATA ####

# Load cleaned and flagged occurrences
allOccs  <- read_csv('outputs/occs/all_occs_flagged.csv')

# Example: Seasonal split of Lasiurus varius occurrences
# --> Example of subsetting occurrences by season for exploratory models
lasiurus <- allOccs %>%
  filter(species == 'Lasiurus varius') %>%
  mutate(month = dmy(date) %>% as_date() %>% month()) %>%
  filter(!is.na(month)) %>%
  mutate(month = ifelse(month %in% c(10,11,12,1,2,3), 'summer', 'winter')) %>%
  mutate(species = str_c(species, '_', month)) %>%
  split(.$species)

# Split occurrences by species
occsList <- allOccs %>% split(.$species)
rm(allOccs)

# Load species log to filter species eligible for modeling
sppLog   <- read.csv('outputs/species_log.csv')

# Load merged modeling areas (Area M)
modArea <- read_sf('outputs/aream/new_area_m.geojson')

# Load predictor layers (WorldClim + ENVIREM)
biopreds <- rast(list.files('inputs/worldclim/2.5m/', '.tif', full.names = TRUE))
names(biopreds) <- gsub('wc2.1_2.5m_', '', names(biopreds))

envirem <- rast(list.files('/mnt/2TB/GIS/Rasters/Clima/ENVIREM/2.5 m/current/', '.bil$', 
                           full.names = TRUE))
names(envirem) <- gsub('current_2-5arcmin_', '', names(envirem))

# Crop predictors to modeling extent (Area M)
preds <- list(
  bio  = biopreds %>% crop(modArea),
  env  = envirem %>% crop(modArea),
  all  = c(biopreds %>% crop(modArea), envirem %>% crop(modArea))
)
rm(biopreds, envirem)
gc()

# Apply environmental filter to Lasiurus varius
# --> Example of pre-filtering occurrences using PCA of climate space
lasiurus <- lapply(lasiurus, \(x) {
  climOccs <- extract(preds$bio, x[, c('lon', 'lat')], cell = TRUE, ID = FALSE)
  climPCA <- prcomp(climOccs)$x
  flg <- envSample(
    x[, c('lon', 'lat')], list(climPCA[,1], climPCA[,2]),
    res = list(diff(range(climPCA[,1]))/40, diff(range(climPCA[,2]))/40),
    do.plot = FALSE
  )
  climOccs <- extract(preds$bio, x[, c('lon', 'lat')], cell = TRUE, ID = FALSE)
  x$spat.unique <- !duplicated(climOccs$cell)
  x$env.unique <- flg
  x$use.env.filter <- sum(flg) > 5
  return(x)
})

# Merge seasonal subset back into main occurrences list
occsList <- c(occsList, lasiurus)
rm(lasiurus)

################################################.
########### SET MONTECARLO VARIABLES ###########
################################################.

# Monte Carlo simulation parameters:
# - totalVariables: number of available predictors per dataset
# - varsRun: number of predictors used per run
# - sprintRuns: number of runs per sprint
# - topVars: number of top predictors selected after each sprint
# - totalSprints: total number of sprints

totalVariables <- map_dbl(preds, nlyr)
varsRun        <- 6
sprintRuns     <- 10
topVars        <- 6
totalSprints   <- 50 

#### MONTE CARLO SIMULATION ####

dir.create(dirMC <- "outputs/montecarlo/", recursive = TRUE)

##### Species loop #####
for(sp in seq_along(occsList)){
  x <- occsList[[sp]]
  
  # Skip species that cannot be modeled (from species log)
  if(x$species[1] %in% sppLog$species[!sppLog$canbemodeled.nobias]) next
  
  # Apply environmental or spatial filter
  use.env <- ifelse(all(x$use.env.filter), TRUE, FALSE)
  x <- if(use.env) x %>% filter(env.unique) else x %>% filter(spat.unique)
  
  ##### Predictor sets loop #####
  for(p in seq_along(preds)){
    ourDir <- paste0(dirMC, x$species[1], '/', names(preds)[p])
    dir.create(paste0(ourDir, '/runs'), recursive = TRUE, showWarnings = FALSE)
    
    # Attempt to resume from previous progress
    advance <- try({read.csv(paste0(ourDir, '/Summary.csv'))})
    sprint_count <- if(inherits(advance, 'try-error')) 1 else max(advance$Sprint) + 1
    
    if(sprint_count > totalSprints) next
    
    # Prepare occurrence data
    occ <- x %>% select(lon, lat) %>% rename(longitude = lon, latitude = lat)
    
    # Tuning parameters for ENMeval
    nrm <- c(0.5, 1, 3, 5)
    fc <- c('L', 'Q', 'LQ', 'LQP')
    
    # Sort predictor names numerically
    lnames <- names(preds[[p]])[str_order(names(preds[[p]]), numeric = TRUE)]
    
    # Extract environmental values for occurrences
    occ <- cbind(occ, terra::extract(preds[[p]][[lnames]], occ)[,-1])
    
    # Initialize variable tally table
    tallyTable <- array(data = 0, dim = c(totalVariables[p], 4))
    dimnames(tallyTable)[[2]] <- c("Layer", "Count", "PI Sum", "PI Avg")
    tallyTable[,1] <- 1:totalVariables[p]
    
    # Mask and crop base raster for background sampling
    basePred <- mask(crop(preds[[p]][[1]], modArea), as(modArea, 'SpatVector'))
    
    ##### Sprint loop #####
    while (sprint_count <= totalSprints) {
      message('Running for species ', x$species[1], ' (', sp, ') at sprint ', sprint_count, 
              '.\n Modeling for predictor set: ', names(totalVariables)[p])
      
      # Resample background points (adds stochasticity across sprints)
      bg <- terra::spatSample(basePred, 10000, xy = TRUE, values = FALSE, na.rm = TRUE) %>%
        as.data.frame() %>%
        rename(longitude = x, latitude = y)
      
      bg <- cbind(bg, terra::extract(preds[[p]][[lnames]], bg)[,-1])
      
      # Setup parallel backend
      cores <- ifelse(10 >= detectCores(), round(detectCores()/2), 10)
      cl <- makeCluster(cores)
      doParallel::registerDoParallel(cl)
      
      ##### Foreach runs loop #####
      run_table <- foreach(run_count = 1:sprintRuns,
                           .packages = c("raster", "terra", "ENMeval", "dismotools", "dplyr")) %dopar% {
                             
                             # Run until valid model found (AICc computed)
                             res_run <- data.frame()
                             while(nrow(res_run) == 0){
                               random_variables <- sample(1:totalVariables[p], varsRun)
                               vars <- c('longitude', 'latitude', lnames[random_variables])
                               occs <- occ[ , vars ]
                               bgs  <- bg [ , vars ]
                               
                               modelRun <- ENMevaluate(occs = occs, bg = bgs,
                                                       tune.args = list(fc = fc, rm = nrm),
                                                       partitions = 'block',
                                                       doClamp = FALSE,
                                                       algorithm = "maxent.jar",
                                                       taxon.name = x$species[1],
                                                       quiet = TRUE)
                               
                               res_run <- modelRun@results %>%
                                 filter(auc.train > 0.5, !is.na(delta.AICc))
                             }
                             
                             # Select best model by AICc
                             xb <- which.min(res_run$delta.AICc)
                             if(length(xb) > 1){
                               xb <- xb[which.min(res_run$auc.diff.avg[xb])]
                             }
                             
                             # Extract permutation importance
                             permImport <- maxent_get_results(modelRun@models[[xb]], "importance")
                             
                             # Update tally table
                             tallyRun <- tallyTable
                             tallyRun[random_variables,2] <- tallyRun[random_variables,2] + 1
                             tallyRun[random_variables,3] <- tallyRun[random_variables,3] + unname(permImport)
                             
                             # Extract run-level metrics
                             aicMod <- eval.results(modelRun)[xb,]
                             outLine <- data.frame(Run = run_count,
                                                   AICc = aicMod$AICc,
                                                   avgAUC = aicMod$auc.val.avg,
                                                   orMTP = aicMod$or.mtp.avg,
                                                   RM = aicMod$rm,
                                                   FC = aicMod$fc,
                                                   Parameters = aicMod$ncoef)
                             
                             list(run_out = outLine, tally_out = tallyRun)
                           } # End foreach
      stopCluster(cl)
      
      ##### Sprint tally aggregation #####
      # (you already have this part perfect — no need to rewrite)
      ...
      
      ##### Sprint summary and saving #####
      # (you already have this part perfect — no need to rewrite)
      ...
      
      # Increment sprint counter
      sprint_count <- sprint_count + 1
    } # end while (sprint loop)
  } # end preds loop
} # end species loop
