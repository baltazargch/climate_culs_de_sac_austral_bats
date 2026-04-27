library(terra)
library(tidyverse)
library(sf)
# Load individual area M and merge 
modArea <- read_sf('outputs/aream/new_area_m.geojson')

## Load and name WordClim layers (predictors) at 2.5m resolution 
biopreds        <- rast(list.files('inputs/worldclim/2.5m/', '.tif', full.names = T))
names(biopreds) <- gsub('wc2.1_2.5m_', '', names(biopreds))

envirem <- rast(list.files('/mnt/28TB/GIS/Rasters/ENVIREM/2.5 m/', '.tif$', 
                           full.names = T))
names(envirem) <- gsub('current_2-5arcmin_', '', names(envirem))

preds <- c(biopreds%>% crop(modArea), envirem %>% crop(modArea)) %>% mask(modArea)


aream <- modArea

msavi <- rast('/mnt/2TB/GIS/Rasters/Coberturas/Vegetation Indexes/msavi.grd')
msavi <- resample(msavi, envirem)

msavic <- crop(msavi, aream)
enviremc <- crop(envirem, aream)


preds <- preds[[as.character(varsIn$.)]]
plot(preds)
cors <- layerCor(preds, 'pearson', asSample=F, use='pairwise.complete.obs')

personCor <- cors$correlation %>% 
  as.data.frame() %>% 
  rownames_to_column(var='var1') %>% 
  pivot_longer(
    cols=!var1, 
    names_to = 'var2'
  ) %>% 
  filter(var1 != var2)

varsIn
