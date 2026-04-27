#### SCRIPT INFO ####
# This script applies previously tuned MaxEnt models to future climate scenarios,
# automating the projection process across multiple time periods, RCPs, and GCMs.
#
# Main steps:
# 1. Load best model settings per species (previously selected manually).
# 2. For each species:
#    a. Sample background points and extract predictor values.
#    b. Prepare occurrence data and predictor values.
#    c. For each combination of time x RCP x GCM:
#        i. Project model using corresponding future layers.
#       ii. Save outputs to organized folders for each species and model.
#
# Notes:
# - The script is designed to be **RAM-efficient** (important for large raster stacks).
# - All projection layers must exist and have matching variable names.
# - Results are organized in `outputs/projections/`.
# - This script also sets up placeholders for:
#     - projection evaluation
#     - niche shift visualization
#     - most limiting factor analysis

#### SETTINGS ####
options(java.parameters = "-Xmx30g")  # Set Java memory limit for MaxEnt.jar

#### LIBRARIES ####
library(terra)       # Raster data handling
library(tidyverse)   # Data wrangling
library(sf)          # Vector data handling
library(dismo)       # MaxEnt projections
source('R/udf/prePara.R')  # User-defined function for preparing MaxEnt arguments

##### PREPARE OUTPUT DIRECTORIES #####
# List species directories (from model fitting results)
spp <- list.dirs('outputs/models/fitting/', full.names = FALSE, recursive = FALSE)

# Create species-specific projection output folders
walk(paste0('outputs/projections/', spp), ~dir.create(.x, recursive = TRUE))

##### LOAD BEST MODEL SETTINGS #####
# Load previously selected "best" models per species (tune_args.csv files)
fls <- list.files('outputs/models/fitting/', 'tune_args.csv$', recursive = TRUE, full.names = TRUE)

# Load and combine tuning results
res_spp <- map(spp, \(x) {
  map(fls[str_detect(fls, x)], \(y) {
    read_csv(y)
  }) %>% list_rbind()
})
names(res_spp) <- spp
res_spp <- res_spp %>% list_rbind(names_to = 'species')

##### SET FUTURE PROJECTION SCENARIOS #####
# Time periods
times <- c('2021_2040', '2041_2060', '2061_2080', '2081_2100')

# GCM model folders (list based on one time period)
mods <- list.dirs('outputs/futureclimate/2061_2080/ssp245/', full.names = FALSE, recursive = F)[-1]

# RCP scenarios
rcp <- c('ssp245', 'ssp370', 'ssp585')

# Combination of time x RCP x GCM
cases <- expand_grid(times, rcp, mods) %>% reduce(., paste, sep = '_')

##### LOAD PREDICTORS #####
# Load and name WorldClim predictors
biopreds <- rast(list.files('inputs/worldclim/2.5m/', '.tif', full.names = TRUE))
names(biopreds) <- gsub('wc2.1_2.5m_', '', names(biopreds))

# Load and name ENVIREM predictors
envirem <- rast(list.files('/mnt/4TB/GIS/Rasters/ENVIREM/2.5 m/', '.tif$', full.names = TRUE))
names(envirem) <- gsub('current_2-5arcmin_', '', names(envirem))

# Combine and mask to Area M
aream <- read_sf('outputs/aream/new_area_m.geojson')
varsIn <- c('climaticMoistureIndex',
            'PETColdestQuarter',
            'PETWettestQuarter',
            'bio_17',
            'bio_2',
            'bio_6',
            'PETseasonality')

preds <- c(biopreds %>% crop(aream), envirem %>% crop(aream))
preds <- terra::mask(terra::crop(preds[[varsIn]], aream), terra::vect(aream))

writeRaster(preds, 'outputs/current_predictors.tif', overwrite=T)
# RUN ONLY ONCE. TO SAVE TIME, FUTURE PROJECTION LAYERS CAN BE PROJECTED TO MATCH CURRENT PREDICTORS USING THE FOLLOWING CODE.
# dirs <- list.dirs('outputs/futureclimate',recursive = F)
# map(dirs, \(d){
#   # d =  dirs[1]
#   dirs.in <- list.dirs(d, recursive = F)
#   
#   map(dirs.in, \(d.in){
#     # d.in = dirs.in[1]
#     dirs.in2 <- list.dirs(d.in, recursive = F)
#     
#     map(dirs.in2, \(d.in2){
#       # d.in2 = dirs.in2[1]
#       fls <- list.files(d.in2, '.asc$', full.names = TRUE)
#       if(length(fls) == 0) stop('No .asc files found in ', d.in2)
#       
#       fut.preds <- rast(fls)
#       # names(fut.preds) <- gsub('.asc', '', basename(fls))
#       
#       fut.preds <- fut.preds %>% project(., preds, method ='bilinear') %>% 
#         mask(preds)
#       
#       writeRaster(fut.preds, paste0(d.in2, '/', names(fut.preds), '.asc'), overwrite=T)
#     })
#   })
# })

# Clean up
rm(biopreds, envirem)
gc()

##### LOAD OCCURRENCE DATA #####
allOccs <- read_csv('outputs/occs/all_occs_flagged.csv')

# Example: split Lasiurus varius by season
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

##### COMBINE Lasiurus seasonal splits #####
occsList <- c(occsList)

##### PARALLELIZATION PREP #####
# Optionally, you could wrap this entire lapply in parallel::mclapply (Linux/Mac only).

# For reproducibilty
set.seed(3458)
library(furrr)

##### SPECIES LOOP #####
lapply(names(occsList), \(sp){
  # sp <- names(occsList)[1]  # for testing
  
  ##### PREPARE DATA FOR THIS SPECIES #####
  sp.dat <- occsList[[sp]]
  
  ##### SAMPLE BACKGROUND POINTS #####
  bg <- terra::spatSample(preds[[1]], 10000, xy = TRUE, na.rm = TRUE)[, 1:2]
  names(bg) <- c('lon', 'lat')
  
  ##### EXTRACT PREDICTOR VALUES #####
  # Occurrence points
  ocwd <- extract(preds, sp.dat[, c('lon', 'lat')]) %>% na.omit()
  ocwd$ID <- 1  # Label presence
  
  # Background points
  bgwd <- extract(preds, bg) %>% na.omit()
  bgwd$ID <- 0  # Label background
  
  ##### COMBINE OCCURRENCE AND BACKGROUND DATA #####
  spwd <- rbind(ocwd, bgwd)
  
  ##### LOAD BEST MODEL SETTINGS FOR THIS SPECIES #####
  res_sp <- res_spp %>% filter(species == sp)
  
  ##### PROJECT BEST MODELS #####
  walk(1:nrow(res_sp), \(r){
    # r=1
    
    ##### PROJECT TO CURRENT TIME #####
    # SET MODEL OUTPUT PATH 
    modpath <- paste0("outputs/projections/",
                      res_sp$species[r], '/',
                      res_sp$fc[r], '_', res_sp$rm[r], '/current')
    
    ##### RUN MAXENT PROJECTION #####
    mod <- maxent(x = spwd[, -1] %>% as.data.frame() %>% 
                    remove_rownames(),
                  p = spwd$ID %>% as.vector(),
                  path = modpath,
                  args = prepPara(
                    userfeatures = res_sp$fc[r],
                    betamultiplier = res_sp$rm[r],
                    doclamp = TRUE,
                    randomseed = TRUE,
                    outputformat = 'cloglog'
                  )
    )
    
    pred <- predict(mod, raster::stack(preds),args=c("outputformat=cloglog")) %>% 
      rast()
    
    writeRaster(pred, paste0(modpath, '/cloglog_pred_suit.tif'), overwrite=T)
    
    vals_th <- mod@results %>% enframe() %>% 
      filter(str_detect(name, 'Cloglog.threshold')) %>% 
      mutate(name = janitor::make_clean_names(name))
    
    bin_preds <- map(1:nrow(vals_th), ~pred >= vals_th[[2]][.x]) %>% 
      rast()
    names(bin_preds) <- vals_th$name
    
    writeRaster(bin_preds, paste0(modpath, '/binary_pred_all_thresholds.tif'), 
                overwrite=T)
    
    vals_out <- vals_th %>% mutate(species = sp, 
                                   fc = res_sp$fc[r], 
                                   rm = res_sp$rm[r], .before = 1) %>% 
      as_tibble() %>% as.matrix() %>% unname()
    colnames(vals_out) <- c('species', 'fc', 'rm', 'name', 'value')
    
    if(!file.exists('outputs/projections/new_threshold_values_per_species.csv')){
      write_csv(as_tibble(vals_out), 'outputs/projections/new_threshold_values_per_species.csv', 
                append = F)  
    } else {
      write_csv(as_tibble(vals_out), 'outputs/projections/new_threshold_values_per_species.csv', 
                append = T)
    }
    
    
    allmods <-  list.files('outputs/cmip6_aligned_to_preds/', 
                                      pattern = '.tif$', 
                                      recursive = T,
                                      full.names = TRUE)

    ##### TIME LOOP #####
    # plan(multisession, workers=8)
    map(times, \(time){
      # preds <- rast('outputs/current_predictors.tif')
      # time = times[1]
    # for(time in times){
      # time=times[1]
      ##### RCP LOOP #####
      for(rc in rcp){
        # rc = rcp[1]
        
        ##### BUILD PROJECTION LAYER PATH #####
        modsprj <- allmods[str_detect(allmods, rc)] %>% .[str_detect(., gsub('_', '-', time))]
        mods_ord <- list.dirs('outputs/cmip6_aligned_to_preds/', recursive = F, 
                              full.names = F)
        names(modsprj) <- mods_ord
        
        ##### SET MODEL OUTPUT PATH #####
        modpath <- paste0("outputs/projections/", res_sp$species[r]) 
        filename <- paste0(modpath, '/',mods_ord,'_', res_sp$fc[r], '_', res_sp$rm[r], '_',
                          time, '_', rc)
        
        dir.create(modpath, recursive = T)
        ##### RUN MAXENT PROJECTION #####
        plan(multisession, workers=8)
        future_map2(modsprj, filename, \(mp, fln){
          # mp = modsprj[1]
          # fln = filename[1]
          preds <- rast('outputs/current_predictors.tif')
          fut.preds <- rast(mp)
          fut.preds <- fut.preds[[varsIn]]
          stopifnot(identical(names(preds), names(fut.preds)))
          fut.mod <- predict(mod,  raster::stack(fut.preds), 
                             rgs=c("outputformat=cloglog", "doclamp=true")) %>% 
            rast()
          
          writeRaster(fut.mod, paste0(fln, '.tif'), overwrite=T)
          
          bin_preds <- map(1:nrow(vals_th), ~fut.mod >= vals_th[[2]][.x]) %>% 
            rast()
          names(bin_preds) <- vals_th$name
          # plot(bin_preds[[4]])
          writeRaster(bin_preds, paste0(fln, '_binaries.tif'), 
                      overwrite=T)
        }) # end GCM loop
        plan(sequential)
      } # end RCP loop
      
    }) # end TIME loop
    # plan(sequential)
    
  }) # end best models loop
  
}) # end species loop

##### FUTURE EXTENSIONS #####
# TODO: Calculate projection metrics (e.g. overlap, stability indices)
# TODO: Visualize ellipsoid shifts (climate change & niche shift)
# TODO: Calculate "most limiting factor" across future scenarios
