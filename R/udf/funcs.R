# Read all csv records in folder dir
read.occs <- function(dir='records/', ext = '.csv$', rm.dups = TRUE){
  
  OCCS    <- lapply(list.files(dir, pattern = ext, full.names = T), read.csv)
  OCCS <- lapply(OCCS, \(x) {colnames(x)[1:3] <- c('species', 'lon', 'lat'); x})
  
  names(OCCS) <- gsub("\\.csv$", "", list.files(dir, pattern = ext))
  
  if(rm.dups) OCCS <- lapply(OCCS, function(x) { subset(x, rep != TRUE) })
  
  return(OCCS)
}

# Standardized file names and keep only
standardize.filenames <- function(recordDir = NULL, cols.keep = 1:3){
  for (i in list.files(recordDir, full.names = T)) {
    newName <- gsub('_', ' ', i)
    c <- read.csv(i)
    colnames(c)[cols.keep] <- c('species', 'lon', 'lat')
    write_csv(c[,cols.keep], newName)
  }
  
  for (i in list.files(recordDir, pattern = '_', full.names = T)) {
    file.remove(i)
  }
}

# Assess montecarlo convergence
asess_convergence <- function(outdir, advance,
                              total_number_of_sprints, OCCS, sp,
                              n_top_vars){
  if (!rlang::is_empty(advance)){
    sprint_advance <- advance[grep('Sprint', list.files(paste0(out_dir, '/runs')))] %>%
      str_split(., '-') %>%
      sapply(., \(x){x[2]}) %>%
      str_remove(., 'Sprint') %>%
      as.numeric() %>%  max()
    
    sp.dir <- paste0(out_dir, '/MC_01-Summary.csv')
    mc <- read.csv(sp.dir)
    minus <- ifelse(nrow(mc) > 10, 9, nrow(mc))
    
    last <- mc[as.numeric(NROW(mc)-minus):NROW(mc), -c(1:10)] %>%
      apply(., 1, sort, simplify = F) %>% do.call(rbind, .)
    
    if(nrow(mc) < 49){
      converged <- FALSE
    } else {
      converged <- all(sapply(1:n_top_vars, \(x) {all(last[1:9, x] == last[10, x])}))
    }
    
    if(sprint_advance >= total_number_of_sprints) {
      message(paste0(names(OCCS)[sp], ' already complete.\n'))
      # next
    } else if(isTRUE(converged)){
      message(paste0(names(OCCS)[sp], ' already converged.\n'))
      # next
    } else {
      sprint_count <<- sprint_advance + 1
      message(paste0(names(OCCS)[sp], ' already has advance. Starting at #',
                     sprint_count, ' Sprint \n'))
    }
    return(converged)
  } else {
    message(paste0(names(OCCS)[sp], ' has no advance. Starting at first Sprint \n'))
  }
}

# Make modeling areas based on different criteria
make_area_m <- function(sp.occ=NULL, method=c('simple', 'ecoregion'), ecoregion = NULL,
                        save.area.m.as.shapes=T, file=NULL,
                        buffer.eco=0.5, buffer.simple=3, cut.terrestrial=T,
                        cutline=NULL, loncol = "decimallongitude", 
                        latcol="decimallatitude",
                        crs.proj = "+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext  +no_defs",
                        crs.out ="+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0") {
  require(sp)
  require(sf)
  require(adehabitatHR)
  
  # sp.occ <- x[,c(1,4,5)] %>% na.omit()
  # ecoregion <- ecoregions
  # crs.out ="+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
  # crs.proj = "+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext  +no_defs"
  # loncol = "lon"
  # latcol="lat"
  
  coordinates(sp.occ) <-  cbind(sp.occ[[loncol]], sp.occ[[latcol]])
  proj4string(sp.occ) <-  CRS(crs.out)

  if(method == 'simple'){
    
    bg.mcp <- as(mcp(as(sp.occ, 'SpatialPoints'), percent = 100), 'sf')
    area.m <- st_buffer(bg.mcp, buffer.simple * 1e05) %>% st_transform(crs.out)
    
    if(cut.terrestrial & !is.null(cutline)) {
      bbox <- st_bbox(area.m)
      
      cutted <- st_intersection(cutline, st_as_sfc(bbox))
      
      area.m <- area.m %>% st_intersection(cutted) %>% st_union()
    }
    
    if(save.area.m.as.shapes) sf::write_sf(area.m, file)
    return(area.m)
    
  } else if(method == 'ecoregion' & !is.null(ecoregion)) {
    
    # require(rangeBuilder)
    mcp <- try({ahull_modified(sp.occ@coords,
                          fraction = 1,
                          partCount = 1, 
                          buff = 100000,
                          verbose=F)})
    
    if(class(mcp) == 'try-error'){
      mcp <- list()
      bg.mcp <- as(mcp(as(sp.occ, 'SpatialPoints'), percent = 100), 'sf')
      area.m <- st_buffer(bg.mcp, 10 * 1e05) %>% st_transform(crs.out)
      mcp[[1]] <- area.m
      }
    
    points.ecos <- sp.occ %>% st_as_sf() %>%  
      st_intersection(ecoregion, .)
    points.ecos <- unique(points.ecos$ECO_NAME)
    
    area.m <- st_intersection(ecoregion, mcp[[1]]) %>% group_by(ECO_NAME) %>% 
      summarise()
   
   flg.ecos <- 5 <= ((st_area(area.m) * 1E-06) / sum(st_area(area.m) * 1E-06)) %>% 
      as.double() %>% 
      round(5) * 100
    
    ecos.in <- unique(area.m$ECO_NAME) %>% na.omit()
    ecos.in <- unique(c(ecos.in[flg.ecos], points.ecos)) %>% na.omit()
    
    area.m <- ecoregion %>% 
      filter(ECO_NAME  %in% ecos.in) %>% 
      filter(!ECO_NAME  %in% c('Lake', 'Rock and Ice')) %>% 
      st_union()
    
    area.m <- smoothr::fill_holes(area.m, 3e+25)
    
    if(save.area.m.as.shapes) sf::write_sf(area.m, file) else return(area.m)
  }
}


#' get Threshold for given species
#' sp = especie
#'
get.threshold <- function(sp, criteria,
                          th = 'Maximum.training.sensitivity.plus.specificity.Cloglog.threshold',
                          base.dir = 'outputs/final models/'){
  require(dplyr)
  
  stopifnot(criteria  %in% c('minOR', 'minAICc'),
            th  %in% c("Fixed.cumulative.value.1.Cloglog.threshold",
                       "Fixed.cumulative.value.5.Cloglog.threshold",
                       "Fixed.cumulative.value.10.Cloglog.threshold",
                       "Minimum.training.presence.Cloglog.threshold",
                       "X10.percentile.training.presence.Cloglog.threshold",
                       "Equal.training.sensitivity.and.specificity.Cloglog.threshold",
                       "Maximum.training.sensitivity.plus.specificity.Cloglog.threshold",
                       "Balance.training.omission..predicted.area.and.threshold.value.Cloglog.threshold",
                       "Equate.entropy.of.thresholded.and.original.distributions.Cloglog.threshold")
            
  )
  
  tbls <- list.files(paste0(base.dir, sp, '/tables/'), 'contri_permu.csv', full.names = T)
  tbls <- tbls[ grep(criteria, tbls)] %>% read_csv()
  
  
  return(tbls$value[ tbls$variable == th ])
  
}

biaslayer <- function(occs_df,longitude, latitude, raster_mold){
  #function taken and modified from https://luismurao.github.io/ntbox/
  require(terra)
  require(MASS)
  require(sp)
  # occs_df <- x
  
  occs_df <- as.data.frame(occs_df)
  stopifnot(inherits(raster_mold, 'SpatRaster'), ncol(occs_df) > 1, 
            c(longitude, latitude) %in% colnames(occs_df))
  
  rbias <- raster_mold
  
  ll_ras <- rasterize(as.matrix(occs_df[,c(longitude,latitude)]),
                      raster_mold, values=1)
  
  no_na <- which(!is.na(raster_mold[]))
  
  ll_rasIDs <- which(values(ll_ras) == 1)
  occ_T <- sp::coordinates(raster::raster(ll_ras))[ll_rasIDs,]
  
  dens <- MASS::kde2d(occ_T[,1], occ_T[,2],
                      n = c(ncol(ll_ras), nrow(ll_ras)), 
                      lims=ext(raster_mold) %>% as.vector())
  
  biasLayer <- rast(raster::raster(dens))
  
  if(!all(dim(biasLayer)[1:2] == dim(raster_mold)[1:2]))
    rbias <- resample(biasLayer, raster_mold, method='bilinear')
  else{
    rbias[no_na] <- biasLayer[no_na]
  }
  
  names(rbias) <- "BiasLayer"
  
  return(rbias)
}

MyvarImportance <- function (theModel, occSWD, bkgSWD, 
                             responseType = c("link", "exponential", "cloglog", "logistic"), numReplicates = 5) 
{
  varList <- names(theModel$samplemeans)
  importance <- vector("numeric", length(varList))
  names(importance) <- varList
  envData <- rbind(occSWD[, -c(1:2)], bkgSWD[, -c(1:2)])
  fullModelVals <- predict(theModel, envData)
  for (thisVar in varList) {
    correls <- vector("numeric", numReplicates)
    origVals <- envData[, thisVar]
    for (thisRep in 1:numReplicates) {
      permInd <- sample(1:nrow(envData), nrow(envData))
      envData[, thisVar] <- origVals[permInd]
      newModelVals <- predict(theModel, envData)
      correls[thisRep] <- cor(fullModelVals, newModelVals)
    }
    envData[, thisVar] <- origVals
    importance[thisVar] <- mean(correls)
  }
  importance <- 1 - importance
  return(round(100 * importance/sum(importance), 2))
}

update_progress <- function(x, step, completion){
  if(!file.exists('progress-track.csv')) stop('Progress track-file has not been initializae')
  
  stopifnot(any(x$step %in% step), is.logical(completion))
  
  x$completion[ x$step == step ] <- completion
  
  write_csv(x, 'progress-track.csv')
  message(paste0('Progress of step **', step, '** updated to: ', completion))
  return(x %>% distinct())
}

add_step_progress <- function(x, step, completion=FALSE){
  if(any(step %in% x$step)){ 
    message('Step already in progressTrack')
    return(x)
  } else{ 
    stopifnot(is.character(step), is.logical(completion))
    message(paste0('Step: **', step, '** added to track log'))
    x %>% add_row(step=step, completion=completion) %>% return()
  }
}

advance_message <- function(order=NULL, species=NULL, action=NULL, appendLF = FALSE){
  stopifnot(!is.null(order), !is.null(species),
            !is.null(action), is.logical(appendLF))
  
  mess <- str_pad(paste0(order, ': ', species, ' (', action, ')'), width = 120, 'right')
  message('\r', mess, appendLF = appendLF)
  flush.console()
}

get_write_gbif_csv <- function(keys, spp, notingbif, dir){
  stopifnot(length(keys) == length(spp))
  k <- 1
  while(k <= length(keys)){
    if(spp[k] %in% notingbif) {
      k <- k+1
      next
    }
    if(file.exists(paste0(dir,'/', spp[k], '.csv'))) {
      cat(paste0(spp[k], ' already exists \n'))
      k <- k+1
      next
    }
    cat(spp[k], sep = '\n')
    x <- try({occ_data(taxonKey=keys[[k]], hasCoordinate=TRUE, limit=2000)})
    # Set the maximum tries per species
    maxi = 20
    i = 1
    while(class(x) == 'try-error'){
      if(i == 1) cat(paste0('Searching records for: ', spp[k]))
      if(i > 1) cat(paste0('Attempt #', i, '\n\n'))
      Sys.sleep(20); gc()
      x <- try({occ_data(taxonKey=keys[[k]], hasCoordinate=TRUE, limit=2000)})
      i <- i+1
      if(i > maxi) break
    } 
    if(i > maxi) next
    
    if(!is.null(x$data)){
      write_csv(x$data, paste0(dir, '/', spp[k], '.csv'))
    } else if(is.null(x$data)){
      notingbif <- c(notingbif, spp[k])
    } 
    
    k <- k+1
    
  }
  return(notingbif)
}

