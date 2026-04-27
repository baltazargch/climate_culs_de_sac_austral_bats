library(terra)
library(sf)
library(tidyterra)
library(tidyverse)
#udf 
suitability_weighted_centroid <- function(suit_map, crs_out = "EPSG:4326") {
  if (!is.lonlat(suit_map)) suit_map <- project(suit_map, crs_out)
  w <- values(suit_map, mat = FALSE)
  ok <- which(!is.na(w) & is.finite(w) & w > 0)
  if (length(ok) == 0) return(data.frame(lon=NA, lat=NA, n_cells=0, weight_sum=0))
  xy <- xyFromCell(suit_map, ok)
  data.frame(
    lon = weighted.mean(xy[,1], w[ok]),
    lat = weighted.mean(xy[,2], w[ok]),
    n_cells = length(ok),
    weight_sum = sum(w[ok], na.rm = TRUE)
  )
}


# Load the data

suitability.files <- list.files('outputs/projections/', '.tif$', full.names = TRUE, recursive = TRUE)
binaries <- list.files('outputs/projections/thresholded_maps/gcm/', '.tif$', full.names = TRUE)

suitability.files <- suitability.files[
  !str_detect(suitability.files, 'nouse') & 
    !str_detect(suitability.files, 'thresholded_maps') & 
    !str_detect(suitability.files, '_binaries.tif$') &
    !str_detect(suitability.files, 'binary_pred_all_thresholds') ]

# analisis grid
data.grid <-
  expand.grid(
    species = list.dirs('outputs/models/fitting', full.names = FALSE, recursive = FALSE), 
    period = c('2021-2040', '2041-2060', '2061-2080', '2081-2100'),
    ssp = c('ssp245', 'ssp370', 'ssp585'), 
    gcm = list.dirs('outputs/futureclimate/2021_2040/ssp245/',recursive = F, 
                    full.names = F)) %>% 
  rbind(tibble(
    species = unique(.$species), 
    period = 'current',
    ssp = 'current', 
    gcm = 'current'
  )) %>% 
  arrange(species, period, ssp) %>% 
  mutate(across(everything(), as.character))

data.grid %>% 
  distinct(gcm) %>% pull(gcm)

# Extract areas x suitability
areas.x.suit <- map(1:nrow(data.grid), \(i){
  # i = 73
  sp.row = data.grid[i,]
  if(sp.row$gcm == 'current') {
    bin.current <- list.files('outputs/projections/thresholded_maps/', '.tif$', full.names = TRUE)
    
    sp.files <- bin.current[ str_detect(bin.current, 'current') & 
                               str_detect(bin.current, gsub(' ', '_', sp.row$species))]
    bin.map <- rast(sp.files)
    
    suit.files <- suitability.files[ str_detect(suitability.files, 'current') & 
                                       str_detect(suitability.files, sp.row$species)]
    suit.map <- rast(suit.files) %>% 
      max() %>% 
      mask(., bin.map > 0, maskvalues = FALSE, updatevalue = NA)
    suitability_weighted_centroid(suit.map) %>% 
      mutate(species = sp.row$species, 
             period = sp.row$period, 
             ssp = sp.row$ssp,
             gcm = sp.row$gcm, .before = 1) 
  } else {
  
  sp.files = binaries[ str_detect(binaries, sp.row$species) & 
                         str_detect(binaries, gsub('-', '_', sp.row$period)) & 
                         str_detect(binaries, sp.row$ssp) & 
                         str_detect(binaries, sp.row$gcm) ]
  
  bin.map <- rast(sp.files) %>% max()
  
  suit.files <- suitability.files[ str_detect(suitability.files, sp.row$species) & 
                                     str_detect(suitability.files, gsub('-', '_', sp.row$period)) & 
                                     str_detect(suitability.files, sp.row$ssp) & 
                                     str_detect(suitability.files, sp.row$gcm) ]
  
  suit.map <- rast(suit.files) %>% max() %>% mask(., bin.map > 0, 
                                                  maskvalues = FALSE, updatevalue = NA)
  
  suitability_weighted_centroid(suit.map) %>% 
    mutate(species = sp.row$species, 
           period = sp.row$period, 
           ssp = sp.row$ssp,
           gcm = sp.row$gcm, .before = 1)
  }
})

centroids <- areas.x.suit %>% 
  list_rbind() %>% 
  mutate(
    species = factor(species),
    period = factor(period, levels = unique(period)),
    ssp = factor(ssp),
    gcm = factor(gcm),
  ) 

centroids %>% write_csv('outputs/centroids_enm_models.csv')

centroids <- read_csv('outputs/centroids_enm_models.csv')

p_lat <- centroids %>%
  mutate(period = factor(period, levels = c("current", "2021-2040", "2041-2060", "2061-2080", "2081-2100"))) %>% 
  ggplot(aes(x = period, y =, group = ssp, colour = ssp)) +
  geom_line() +
  geom_point(size = 2) +
  scale_colour_manual(values=c('black','#2C7FB8', '#BDBDBD', '#CB181D'))+
  facet_wrap(~ species, scales='free_y') +
  labs(
    x = NULL,
    y = "Suitability-weighted centroid latitude (°)"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

p_lat

library(forcats)

km_per_deg_lat <- 111.32

baseline <- centroids %>%
  filter(period == "current", ssp == "current") %>%
  select(species, lat_current = lat)

# 2) join baseline to every row and compute shift (current will be 0 for ssp=="current")
centroids_shift <- centroids %>%
  mutate(
    period = factor(period, levels = c("current", "2021-2040", "2041-2060", "2061-2080", "2081-2100"))
  ) %>%
  left_join(baseline, by = "species") %>%
  mutate(
    shift_km = (lat - lat_current) * km_per_deg_lat
  )

# SSPs you want to plot (exclude the "current" pseudo-ssp)
ssps_to_plot <- centroids_shift %>%
  distinct(ssp) %>%
  filter(ssp != "current") %>%
  pull(ssp)

# build synthetic current rows (one per species x SSP) at 0 km
current_rows <- centroids_shift %>%
  filter(period == "current", ssp == "current") %>%
  select(species, period, lat, lon, lat_current) %>%  # keep whatever cols you have
  crossing(ssp = ssps_to_plot) %>%                               # <-- now ssp doesn't exist yet
  mutate(shift_km = 0)

write_csv(centroids_shift, 'centroid_shift_km.csv')
# final plotting dataset: future SSP rows + synthetic current rows
centroids_shift_plot <- centroids_shift %>%
  filter(ssp != "current") %>%
  bind_rows(current_rows) %>%
  mutate(
    period = factor(period, levels = c("current", "2021-2040", "2041-2060", "2061-2080", "2081-2100")),
    period = fct_drop(period)
  )

p_shift <- 
  centroids_shift_plot %>% 
  mutate(period = as.character(period)) %>% 
  mutate(period = ifelse(period == 'current', 'baseline', period)) %>% 
  mutate(period = factor(period, 
                         levels = c("baseline", "2021-2040", "2041-2060", "2061-2080", "2081-2100"))) %>%
  group_by(ssp, period, species) %>% 
  summarise(shift_km = mean(shift_km)) %>%
  ggplot(aes(x = period, y = shift_km, group = ssp, colour = ssp)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_line(alpha=0.8) +
  geom_point(size = 2, alpha=0.8) +
  scale_colour_manual(values=c('#1f5a8a', '#d66b00', '#b30000'))+
  facet_wrap(~ species, scales = "free_y") +
  labs(x = NULL, y = "Shift in model centroid \nfrom baseline (km)") +
  theme_light()+
  theme(
    
    legend.position = "bottom",
    strip.text =  element_text(face = 'italic'),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

p_shift

ggsave(filename = 'graphical_abstract_shifts_km.svg')

map.arg <- rnaturalearth::ne_countries(
  scale = 10, country = c("Chile","Argentina"), returnclass = "sf"
) %>%
  st_transform(4326)

# optional: keep only future SSPs (recommended)
centroids_map <- centroids %>%
  # filter(ssp != "current") %>%
  mutate(
    period = factor(period, levels = c("current","2021-2040","2041-2060","2061-2080","2081-2100"))
  )

p_map <- ggplot() +
  geom_sf(data = map.arg, fill = "grey95", color = "grey40", linewidth = 0.3) +
  # geom_path(
  #   data = centroids_map,
  #   aes(x = lon, y = lat, group = ssp, colour = ssp),
  #   # arrow = grid::arrow(length = unit(0.10, "inches")),
  #   linewidth = 0.6
  # ) +
  geom_point(
    data = centroids_map,
    aes(x = lon, y = lat, colour = ssp, shape = period),
    size = 2
  ) +
  scale_colour_manual(values = c('black','#2C7FB8', '#BDBDBD', '#CB181D')) +
  coord_sf(expand = FALSE) +
  facet_wrap(~ species) +
  labs(
    x = "Longitude (°)",
    y = "Latitude (°)",
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",legend.byrow = F,
    panel.grid.major = element_line(linewidth = 0.2),
    panel.grid.minor = element_blank()
  )

bb <- centroids_map %>%
  summarise(
    xmin = min(lon, na.rm = TRUE),
    xmax = max(lon, na.rm = TRUE),
    ymin = min(lat, na.rm = TRUE),
    ymax = max(lat, na.rm = TRUE)
  )

p_map_zoom <- p_map +
  coord_sf(
    xlim = c(bb$xmin - 2, bb$xmax + 2),
    ylim = c(bb$ymin - 2, bb$ymax + 2),
    expand = FALSE
  )

p_map_zoom

