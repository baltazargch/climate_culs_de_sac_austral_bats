###############################################################
## DATA GATHERING AND CLEANING: OCCURRENCE POINTS AND RANGES ##
###############################################################

#### Load libraries ####
library(rgbif)
library(CoordinateCleaner)
library(redlistr)
library(sf)
library(tidyverse)

#### General parameters ####
sf_use_s2(FALSE)

#### Download taxonomy from GBIF and IUCN ####

ks <- map(c('Lasiurus varius', 
            'Histiotus magellanicus',
            'Myotis chiloensis'), ~name_backbone(name=.x))
ks <- ks %>% map(~.x$usageKey) %>% unlist()

taxonKeys <- c(
  ks
)
for (k in taxonKeys) {
  
  taxGBIF <- occ_download(
    pred_and(
      pred_in("taxonKey", k),
      pred("hasCoordinate", TRUE), 
      pred_in("basisOfRecord", c("Preserved Specimen", "Material Citation"))
    ),
    format = 'SIMPLE_CSV'
  )
  
  if(!dir.exists('inputs/gbif')) dir.create('inputs/gbif', recursive = TRUE)
  
  status = occ_download_meta(taxGBIF)$status
  count = 1
  while(status != 'SUCCEEDED' | count > 20){
    Sys.sleep(120)
    status = occ_download_meta(taxGBIF)$status
    print(count) 
    count = count + 1
  }
  
  x <- occ_download_get(taxGBIF, path = 'inputs/gbif', overwrite =FALSE)
}

x <- as.download(key = "inputs/gbif/0000137-240202131308920")
occs<-occ_download_import(x)

x <- as.download(key='inputs/gbif/0000141-240202131308920') 
occs1<-occ_download_import(x)

x <- as.download(key='inputs/gbif/0000149-240202131308920') 
occs2<-occ_download_import(x)

occs <- rbind(occs, occs1, occs2)

wrld <- read_sf('D:/GIS/Vectores/Mundo/ne_10m_admin_0_countries.shp')

occs_clean <- clean_coordinates(
  occs,seas_ref = wrld,
  tests = c("capitals", "centroids", "equal", "gbif", 
            "institutions", "seas", "zeros")
)

occs_out <- occs[occs_clean$.summary, ] 

colnames(occs_out)<-tolower(colnames(occs_out))

occs_filt <- occs_out %>% 
  filter(coordinateuncertaintyinmeters < 5000)
write_csv(occs_filt, 'inputs/20240202_gbif_records.csv')

##### Filter with IUCN ranges #####
library(sf)
library(tidyverse)
sf_use_s2(FALSE)

db <- read_csv('inputs/20240202_gbif_records.csv')
db$inkownrange <- NA

spat <- st_read(
  'D:/GIS/Vectores/Especies/IUCN/MAMMALS_TERRESTRIAL_ONLY/MAMMALS_TERRESTRIAL_ONLY.shp') 
  
spat <- spat %>% filter(binomial  %in%  c("Lasiurus varius", "Histiotus magellanicus",
            "Myotis chiloensis"))
# plot(spat$geometry)
write_sf(spat, 'inputs/range_maps_IUCN.gpkg')

db.list <- split(db, db$species)

spat.filt <- spat %>% filter(binomial  %in% unique(db$species))
spat.list <- split(spat.filt, spat.filt$binomial)

stopifnot(identical(names(db.list) %>% sort, 
                    names(spat.list) %>% sort))

db.filt <- lapply(
  names(db.list),
  \(x){
    # x <- names(db.list)[1]
    crds <- st_as_sf(db.list[[x]], 
                     coords=c(
                       "decimallongitude",
                       "decimallatitude"
                     )
    )
    st_crs(crds) <- 4326
    pts.in <- st_intersects(crds, spat.list[[x]], sparse=T)
    
    db.list[[x]]$inkownrange <- lengths(pts.in)>0
    
    db.list[[x]]
    
  }
  
)

dir.create('outputs/occs', recursive = T)

write_csv(do.call(rbind, db.filt), 'outputs/occs/20240202_gif_iucn_flagged.csv')
