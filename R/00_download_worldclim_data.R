library(curl)
library(tidyverse)

case <- paste0('ssp', c('245', '370', '585'))
base <- 'https://geodata.ucdavis.edu/cmip6/2.5m'
mods <- c('CMCC-ESM2', 'EC-Earth3-Veg', 
          'GISS-E2-1-G', 'INM-CM5-0', 
          'ACCESS-CM2', 'MIROC6')
suburl <- 'wc2.1_2.5m' 
vars <- c('tmin', 'tmax', 'prec', 'bioc')
times <- c('2021-2040', '2041-2060', '2061-2080', '2081-2100')

downgrid  <- expand_grid(case, base, mods, suburl, vars, times) %>% 
  mutate(downlink = str_c(base, '/', mods, '/', case, '/', 
                          suburl, '_', vars, '_', mods, '_', 
                          case, '_', times, '.tif')) %>% 
  mutate(basedir = '/mnt/4TB/GIS/Rasters/WorldClim/Future/', 
         outdir = case_when(
           vars == 'tmin' ~ 'tn_2_5m', 
           vars == 'tmax' ~ 'tx_2_5m', 
           vars == 'prec' ~ 'pr_2_5m', 
           vars == 'bioc' ~ 'bio_2_5m', 
         )) %>% 
  mutate(outdir = str_c(basedir, outdir), 
         outfile = str_c(outdir,'/', basename(downlink)))

walk(unique(downgrid$outdir), ~dir.create(.x, recursive = T))

library(furrr)

plan(multisession, workers=15)
future_walk2(downgrid$downlink, downgrid$outfile, \(x,y) curl_download(x, y))
plan(sequential)


fls <- list.files('/mnt/4TB/GIS/Rasters/WorldClim/Future/', 
                  '.tif', recursive = T, full.names = T) 
# unlink(fls[ str_detect(fls, '__')])
stopifnot(length(fls) == nrow(downgrid))
