library(terra)
library(tidyverse)
library(sf)
sf_use_s2(F)

spp <- list.dirs('outputs/models/fitting/', F, F)
fls <- list.files('outputs/models/fitting/', 
                  'tune_args.csv$', 
                  recursive = T, 
                  full.names = T)
res_spp <- map(spp, \(x) {
  map(fls[ str_detect(fls, x)], \(y) {
    read_csv(y, show_col_types = F)}) %>% list_rbind()
}
)
names(res_spp) <- spp
res_spp <- res_spp %>% list_rbind(names_to = 'species')

dds <- list.dirs('outputs/projections', recursive = F)
names(dds) <- spp
mods <- list.dirs('outputs/futureclimate/2061_2080/ssp245/', full.names = F)[-1]

mods <- mods[ -grep('maxent.cache', mods)]

rcp <- c('ssp245', 'ssp585')
times <- c('2061_2080','2081_2100')

aream <- read_sf('outputs/aream/new_area_m.geojson')

for(sp in spp){
  d <- dds[sp]
  ff <- list.files(d, pattern = '.asc$', full.names = T, recursive = T)
  projs <- sapply(mods, \(x) ff[ grep(paste0(x, '.asc'), ff)]) %>% unlist() %>% as.vector()
  
  res_sp <- res_spp %>% filter(species == names(d))
  for(r in 1:nrow(res_sp)){
    # r=1
    pr <- projs[ grep(fcrm <- paste0(res_sp$fc[r], '_', res_sp$rm[r]), projs)]
    for(t in times){
      # t <- times[1]
      prt <- pr[ grep(t, pr)]
      for(s in rcp){
        # s <- rcp[1]
        prts <- prt[ grep(s, prt)] 
        rr <- rast(prts) %>% mean() %>% mask(aream)
        
        on <- paste0(d, '/average_allGCM_', fcrm, '_', t, '_', s, '.tif')
        writeRaster(rr, on, overwrite=T)
      }
    }  
  }
}
