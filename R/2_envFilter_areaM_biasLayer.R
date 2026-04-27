#### SCRIPT INFO ####
# This script prepares cleaned and spatially/environmentally-filtered occurrence data,
# derives species-specific modeling areas ("area M"), and creates spatial bias layers
# to be used in species distribution modeling (SDM).
#
# The script follows three main steps:
# 1. Environmental filtering of occurrences based on WorldClim layers
#    (Varela et al. 2014: doi: 10.1111/j.1600-0587.2013.00441.x).
# 2. Generation of species-specific modeling areas using ecoregions (Barve et al. 2011).
# 3. Creation of spatial bias layers for each species.
#
# Notes:
# - WorldClim bio8-9, bio18-19 are excluded due to known artifacts.
# - Modeling areas are based on ecoregions intersected with occurrence convex hulls
#   or alpha-hulls if possible.
# - Bias layers are density surfaces of occurrences clipped to area M.

#### LIBRARIES ####
library(tidyverse)    # For data wrangling
library(parallel)     # For parallelization (not used explicitly here)
library(foreach)      # For parallel loops (not used explicitly here)
library(terra)        # For raster data manipulation
library(sf)           # For vector data manipulation

#### SETTINGS AND USER-DEFINED FUNCTIONS ####
sf_use_s2(FALSE)      # Disable spherical geometry to avoid some issues with polygon validity
source('R/udf/funcs.R')
source('R/udf/envSamp.R')
source('R/udf/alphahull_modified.R')
source('R/udf/alphahulls2sf.R')

#### BASE DATA ####
# Load bioclimatic predictors (WorldClim 2.5m resolution)
# Note: Exclude known problematic variables (bio8-9, bio18-19)
bioclim <- rast(list.files('inputs/worldclim/2.5m/', pattern = '.tif', full.names = TRUE))

#### STEP 1: ENVIRONMENTAL FILTERING ####
# Goal: Flag occurrences that are environmentally unique in climatic space,
# based on a PCA of climate variables. 
# If too few unique records remain, fallback to a spatial filter based on pixel uniqueness.

# Load cleaned occurrence records
allOccs <- lapply(
  list.files('inputs/records/', '.csv', full.names = TRUE), 
  \(x) read_delim(file = x, delim = ';')
) %>% do.call(rbind, .) %>% filter(!is.na(lon))

# Split occurrences by species
occsList <- allOccs %>% split(allOccs$species)
rm(allOccs)
gc()

# Apply filtering per species
occsFlag <- lapply(
  occsList, 
  \(x) {
    # Extract climate values at occurrence points
    climOccs <- extract(bioclim, x[, c('lon', 'lat')], cell = TRUE, ID = FALSE)
    
    # Initial spatial filter: unique climate grid cells
    x$spat.unique <- !duplicated(climOccs$cell)
    
    # If too few unique records, skip environmental filter
    if(sum(x$spat.unique) < 5){
      x$env.unique <- FALSE
      x$use.env.filter <- FALSE
      return(x)
    } else {
      # Perform PCA on climate variables (excluding NAs if needed)
      climOccs <- climOccs[, -ncol(climOccs)]  # Remove 'cell' column
      occsNA <- which(is.na(climOccs), arr.ind = TRUE)[, 1] %>% unique()
      occsIn <- which(!is.na(climOccs), arr.ind = TRUE)[, 1] %>% unique()
      
      if(length(occsNA) > 0){
        climPCA <- prcomp(climOccs[-occsNA, ])$x
        flgP <- envSample(
          x[-occsNA, c('lon', 'lat')], 
          list(climPCA[, 1], climPCA[, 2]), 
          res = list(diff(range(climPCA[, 1])) / 40, diff(range(climPCA[, 2])) / 40),
          do.plot = FALSE
        )
        flg <- vector(length = nrow(x))
        flg[occsIn] <- flgP
      } else {
        climPCA <- prcomp(climOccs)$x
        flg <- envSample(
          x[, c('lon', 'lat')], 
          list(climPCA[, 1], climPCA[, 2]), 
          res = list(diff(range(climPCA[, 1])) / 40, diff(range(climPCA[, 2])) / 40),
          do.plot = FALSE
        )
      }
      
      x$env.unique <- flg
      x$use.env.filter <- sum(flg) > 5
      return(x)
    }
  }
)

# Combine filtered occurrences and save
occsFlag <- do.call(rbind, occsFlag)
saveRDS(occsFlag, 'outputs/occs/occs_Cleaned_australBats.rds')
gc()

#### STEP 2: MAKE MODELING AREA (Area M) ####
# Goal: Derive species-specific modeling area using ecoregions + occurrence-based hulls.

# Load ecoregions
ecoregions <- st_read('inputs/bioregions/southerconeecoregions.gpkg') %>% st_make_valid()
stopifnot(st_is_valid(ecoregions))  # Validate geometry

# Create output directory for area M
dir.create('outputs/aream', recursive = TRUE)

# Prepare data per species
occsDt <- split(occsFlag, occsFlag$species)

# Loop over species and derive area M
spMaps <- lapply(occsDt, \(x) {
  x <- x %>% filter(spat.unique)
  
  filSave <- paste0('outputs/aream/', x$species[1], '.geojson')
  
  if(nrow(x) < 3){
    return(data.frame(species = x$species[1], hasmap = FALSE, rows = nrow(x)))
  } else if(file.exists(filSave)) {
    return(data.frame(species = x$species[1], hasmap = TRUE, rows = nrow(x)))
  } else {
    advance_message('Austral bats:', x$species[1], 'estimating area M',
                    appendLF = x$species[1] == names(occsDt)[length(occsDt)])
    
    m <- try({
      make_area_m(x, 'ecoregion', ecoregion = ecoregions,
                  save.area.m.as.shapes = TRUE,
                  loncol = 'lon', latcol = 'lat',
                  file = filSave)
    })
    gc()
    
    # Fallback if alpha-hull fails: use buffered convex hull + ecoregion intersection
    if(class(m)[1] == 'try-error'){
      coordinates(x) <- cbind(x$lon, x$lat)
      proj4string(x) <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
      
      bg.mcp <- as(mcp(as(x, 'SpatialPoints'), percent = 100), 'sf')
      area.m <- bg.mcp %>%
        st_buffer(3) %>%
        st_transform(st_crs(x)) %>%
        st_crop(ecoregions) %>%
        st_intersection(ecoregions) %>%
        st_union()
      
      sf::write_sf(area.m, filSave)
    }
    
    return(data.frame(species = x$species[1], hasmap = TRUE, rows = nrow(x)))
  }
})

# Combine area M results
spMaps %>% do.call(rbind, .)

#### STEP 3: MAKE BIAS LAYER ####
# Goal: Create species-specific bias layers (smoothed kernel density raster),
# constrained to area M, to correct for spatial sampling bias.

# Use first bioclim layer as template raster
baseRaster <- bioclim$wc2.1_2.5m_bio_1

# Load updated area M (single file, e.g. after manual adjustments)
modArea <- st_read('outputs/aream/new_area_m.geojson') %>% st_make_valid()

# Create output directory for bias layers
dir.create('outputs/biaslayers', recursive = TRUE)

# Prepare list
occsDt <- split(occsFlag, occsFlag$species)

# Loop over species to create bias layers
sppBias <- lapply(names(occsDt), \(sp){
  f <- modArea
  
  # Bias layer only created if area M exists
  if(is_empty(f)){
    return(data.frame(species = sp, hasbias = FALSE))
  } else {
    if(file.exists(paste0('outputs/biaslayers/', sp, '.tif'))) {
      return(data.frame(species = sp, hasbias = TRUE))
    } else {
      x <- occsDt[[sp]]
      
      # Determine whether to apply environmental or spatial filter
      use.env <- ifelse(all(occsDt[[sp]]$use.env.filter), TRUE, FALSE)
      
      if(use.env) {
        x <- x %>% filter(env.unique)
      } else {
        x <- x %>% filter(spat.unique)
      }
      
      spArea <- f
      
      # Generate bias layer (see biaslayer() in funcs.R)
      spBias <- try({
        biaslayer(x, 'lon', 'lat',
                  mask(crop(baseRaster, spArea), spArea))
      })
      
      if(class(spBias) == 'try-error') {
        return(data.frame(species = sp, hasbias = FALSE))
      } else {
        writeRaster(spBias, paste0('outputs/biaslayers/', sp, '.tif'))
        return(data.frame(species = sp, hasbias = TRUE))
      }
    }
  }
})

#### FINAL STEP: SUMMARY TABLE ####
# Combine info from area M, record counts, and bias layer into one summary table
finalSpp <- left_join(spMaps %>% do.call(rbind, .), 
                      sppBias %>% do.call(rbind, .), 
                      by = 'species') %>% as_tibble()

# Flag species eligible for modeling
finalSpp <- finalSpp %>%
  mutate(enoughrecords = rows >= 5) %>%
  dplyr::select(-rows)

# Species that can be modeled with and without bias layer
finalSpp$canbemodeled <- rowSums(finalSpp[, -1]) == 3 
finalSpp$canbemodeled.nobias <- rowSums(finalSpp[, c(2, 4)]) == 2
finalSpp$nrecords <- spMaps %>% do.call(rbind, .) %>% .$rows

# Save summary
write_csv(finalSpp, 'outputs/species_log.csv')
