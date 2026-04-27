### generate solar radiation ####
library(terra)
library(sf)
library(palinsol)
library(envirem)
library(tidyverse)
library(furrr)
sf_use_s2(F)

aream <- read_sf('outputs/aream/new_area_m.geojson')

rasterTemplate <- rast('inputs/worldclim/2.5m/wc2.1_2.5m_bio_1.tif') %>% crop(aream)

years <- c(71:(2100-1950))

for(y in years) dir.create(paste0('outputs/futureclimate/yearly/', 
                                  as.character(y+1950)), 
                           recursive = T)

library(parallel)
mclapply(years, \(y) {
  # calculate monthly solar radiation, defined for the year 2041
  ETsolradRasters(rasterTemplate = rasterTemplate, year = y, 
                  outputDir = paste0('outputs/futureclimate/yearly/', 
                                     as.character(y+1950)),
                  overwrite = TRUE)
}, mc.cores = 20)

plan(multisession, workers=4)
future_walk(c('2021_2040', '2041_2060', '2061_2080', '2081_2100'), \(case){
  
  dir.create(paste0('outputs/futureclimate/', case, '/solarrad'), 
             recursive = T)
  
  for(i in 1:12){
    # i = 1
    i <- str_pad(i, width = 2, side = 'left', pad = 0)
    c <- str_split(case, '_')[[1]] %>% as.numeric
    
    fls <- paste0('outputs/futureclimate/yearly/', c[1]:c[2], 
                  '/et_solrad_', i, '.tif')
    mr <- rast(fls) %>% mean()
    writeRaster(mr, paste0('outputs/futureclimate/', case, '/solarrad/', 
                           'ET_sr_mon_', i, '_', case, '.tif'), 
                overwrite=T)
  }
})
plan(sequential)
# Make sure all it's OK
r <- rast(list.files('outputs/futureclimate/2021_2040/solarrad/', '.tif', full.names = T))
names(r) <- month.name
plot(r, col=viridis::inferno(10))

unlink('outputs/futureclimate/yearly/', recursive = T)
gc()

ref <- rasterTemplate
# plot(ref)

r <- rast('/mnt/4TB/GIS/Rasters/WorldClim/Future/pr_2_5m/wc2.1_2.5m_prec_CMCC-ESM2_ssp245_2061-2080.tif')
r <- resample(r, ref[[1]]) #%>% crop(r, ref) %>% mask(., ref)

v <- values(r)
v[v < 0] <- NA
values(r) <- v

ref <- r
tifOptions <- c("COMPRESS=DEFLATE", "PREDICTOR=2", "ZLEVEL=6")
# TODO ajustar para nuevos directorios y variables en directorio madre
fls <- list.files('/mnt/4TB/GIS/Rasters/WorldClim/Future/', '.tif$', 
                  full.names = T, recursive = T)
# fls <- fls[ !str_detect(fls, 'bio')]
walk(c('2021_2040', '2041_2060', '2061_2080', '2081_2100'), \(case){
  # case='2021_2040'
  flsin <- fls[str_detect(fls, gsub('_', '-', case))]
  # flsin <- flsin[!str_detect(flsin, 'bio_2_5m')]
  # library(parallel)
  
  mclapply(1:length(flsin), \(i){
    # for (i in 1:length(flsin)) {
    # i=24
    cat(i, ' ')
    r <- rast(flsin[i], ) %>% resample(ref) %>% mask(ref)
    # plot(r)
    bn <- basename(flsin[i]) %>% str_remove('wc2.1_2.5m_') %>% 
      str_remove(paste('_', case)) %>% str_split('_', n = 2, simplify = F) %>% 
      unlist() 
    names(r) <- str_c(bn[1], str_pad(1:nlyr(r),width = 2, pad = '0', side = 'left'), '_', bn[2])
    # if temperature, divide values by 10
    ## recognize temperature files as those containing the terms tasmin or tasmax
    # if (grepl('tn|tx', fls[i])) {
    #   r <- r / 10
    # }
    outdir <- ifelse(
      str_detect(flsin[i], 'pr'), paste0('outputs/futureclimate/', case, '/prec'), 
      ifelse(
        str_detect(flsin[i], 'tn'),paste0('outputs/futureclimate/', case, '/tmin'), 
        ifelse(
          str_detect(flsin[i], 'tx'), paste0('outputs/futureclimate/', case, '/tmax'),
          paste0('outputs/futureclimate/', case, '/bio'))))
    
    dir.create(outdir, recursive = T, showWarnings = F)
    
    outfile <- paste0(outdir, '/', names(r))
    writeRaster(r, filename = outfile, gdal = tifOptions, 
                overwrite = TRUE)
    # }
  }, mc.cores=14)
})

# Generate ENVIREM for future scenarios
library(terra)
library(sf)
library(tidyverse)
library(envirem)
library(parallel)

plan(multisession, workers=4)
future_walk(c('2021_2040', '2041_2060', '2061_2080', '2081_2100'), \(xc){
  # xc = '2061_2080'
  basedir <- paste0('outputs/futureclimate/', xc)
  
  mods <- list.files(paste0(basedir, '/prec/'), pattern = '.tif')
  
  mods <- gsub('prec[0-9]+_|_ssp.*', '', mods) %>% unique()
  rcp <- c('ssp245', 'ssp370','ssp585')
  
  cases <- expand_grid(mods, rcp) %>% reduce(., paste, sep='_')
  
  fls <- list.files(basedir, '.tif$', full.names = T, recursive = T)
  
  solfls <- fls[ str_detect(fls, 'solarrad')]
  solarstack <- rast(solfls)
  names(solarstack) <- paste0('ET_sr_mon_', str_pad(1:12, 2, 'left', 0), '_', xc)
  
  lapply(cases, \(c) {
    # for(c in cases){
    # for(ssp in rcp){
      #   ssp =rcp[1]
      # c = cases[1]
      assignNames(
        tmin = paste0('tmin##_', c, '_', gsub('_','-',xc), '.tif'), 
        tmax =paste0('tmax##_', c, '_', gsub('_','-',xc), '.tif'), 
        # tmean =paste0('CHELSA_tas_mon_', c, '_r1i1p1_g025.nc_##_', xc, '_V1.2'), 
        prec = paste0('prec##_', c,  '_', gsub('_','-',xc), '.tif'),
        solrad = paste0('ET_sr_mon_##_', xc) 
      )
      
      modfls <- fls[ str_detect(fls, c)]
      dir.create(paste0(basedir, '/envirem'), recursive = T, showWarnings = F)
      
      masterstack <- rast(modfls[ !str_detect(modfls, 'bio')])
      
      env <- generateEnvirem(masterstack = masterstack, 
                             solradstack = solarstack, 
                             var = 'all', 
                             tempScale = 1)
      writeRaster(env, 
                  filename = paste0(basedir, '/envirem/',
                                    names(env), '_', c, '_', xc,'.tif'), 
                  overwrite=T)
    # }
  })
})
plan(sequential)

varsIn = c('climaticMoistureIndex',
           'PETColdestQuarter',
           'PETWettestQuarter',
           'bio_17',
           'bio_2',
           'bio_6',
           'PETseasonality')

plan(multisession, workers=4)
future_walk(c('2021_2040', '2041_2060', '2061_2080', '2081_2100'), \(xc){
  
  basedir <- paste0('outputs/futureclimate/', xc)
  
  fls <- list.files(basedir, '.tif$', full.names = T, recursive = T)
  
  rcp <- c('ssp245', 'ssp370','ssp585')
  
  mods <- list.files(paste0(basedir, '/prec'), pattern = '.tif')
  mods <- gsub('prec[0-9]+_|_ssp.*', '', mods) %>% unique()
  
  for(r in rcp){
    # r <- rcp[1]
    fls_rcp <- fls[ str_detect(fls, r )]
    
    fls_rcp_mods <- map(mods, ~fls_rcp[ str_detect(fls_rcp, .x)]) 
    
    vars <- c(  "bio", "envirem")
    fls_rcp_mods_vars <- map(fls_rcp_mods, ~.x[ str_detect(.x, 'bio|envirem')]) 
    names(fls_rcp_mods_vars) <- mods
    
    walk(mods, 
         \(x){
           # x <- mods[[1]]
           dir.create(moddir <- paste0('outputs/futureclimate/',xc,
                                       '/', r, '/', x), 
                      recursive = T, showWarnings = F)
           
           rrs <- rast(fls_rcp_mods_vars[[x]])
        
           var.names <- fls_rcp_mods_vars[[x]] %>% 
             str_remove(., paste0('_', x, '.*')) %>% str_c(., '.tif')
          
           names(rrs) <- names(rrs) %>% str_remove(., '_.*')
           names(rrs) <- names(rrs) %>% str_replace_all('bioc0', 'bio_') %>% str_replace_all('bioc', 'bio_')
           writeRaster(rrs[[varsIn]], filename = paste0(moddir, '/', names(rrs[[varsIn]]), '.asc'), 
                       overwrite=T, NAflag=-9999)
           unlink(fls_rcp_mods_vars[[x]])
         })
  }
})
plan(sequential)

dirs <- dir_ls('outputs/futureclimate/', recurse = T, type='directory')
dirs[ str_detect(dirs, '/b|/p|/t|/so|/en')] %>% 
  map(~unlink(.x, recursive = T))

library(terra)
library(fs)  # For directory and file handling

# Define the main directory containing subdirectories with raster files
main_dir <- "outputs/futureclimate/"

# List all raster files recursively from subdirectories
raster_files <- dir_ls(main_dir, glob = "*.tif", recurse = TRUE)

# Function to rename bands and overwrite the raster file
rename_and_overwrite <- function(file_path) {
  # Open the raster
  raster <- rast(file_path)
  # file.remove(file_path)
  # Update band names by removing file extension
  band_names <- names(raster)
  names(raster) <- gsub("\\.tif$", "", band_names, ignore.case = TRUE)
  # Define a temporary file path
  temp_file <- tempfile(fileext = ".tif")
  
  # Save the raster with updated names to the temporary file
  writeRaster(raster, temp_file, overwrite = TRUE)
  
  # Remove the original file
  file.remove(file_path)
  
  # Rename the temporary file to the original file path
  file.rename(temp_file, file_path)
  
  cat("Updated:", file_path, "\n")
}

# Apply the function to each raster file
lapply(raster_files, rename_and_overwrite)

