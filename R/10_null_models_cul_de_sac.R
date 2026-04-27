#### SETTINGS ####
options(java.parameters = "-Xmx30g")  # Set Java memory limit for MaxEnt.jar

#### LIBRARIES ####
library(terra)       # Raster data handling
library(tidyverse)   # Data wrangling
library(sf)          # Vector data handling
library(ENMeval)       # MaxEnt projections
library(dismo)       # MaxEnt projections
source('R/udf/prePara.R')  # User-defined function for preparing MaxEnt arguments

suitability_weighted_centroid <- function(suit_map, crs_out = "EPSG:4326") {
  if (!is.lonlat(suit_map))
    suit_map <- project(suit_map, crs_out)
  w <- values(suit_map, mat = FALSE)
  ok <- which(!is.na(w) & is.finite(w) & w > 0)
  if (length(ok) == 0)
    return(data.frame(
      lon = NA,
      lat = NA,
      n_cells = 0,
      weight_sum = 0
    ))
  xy <- xyFromCell(suit_map, ok)
  data.frame(
    lon = weighted.mean(xy[, 1], w[ok]),
    lat = weighted.mean(xy[, 2], w[ok]),
    n_cells = length(ok),
    weight_sum = sum(w[ok], na.rm = TRUE)
  )
}

# Build a smoothed sampling-bias raster from observed occurrences
make_bias_raster <- function(occ_df, template_rast, sigma = 2) {
  # occ_df must have lon, lat
  pts <- vect(occ_df[, c("lon", "lat")],
              geom = c("lon", "lat"),
              crs = crs(template_rast))
  
  # count occurrences per raster cell
  r_occ <- rasterize(
    pts,
    template_rast[[1]],
    field = 1,
    fun = "sum",
    background = 0
  )
  
  # gaussian-like kernel
  k <- focalMat(r_occ, d = sigma, type = "Gauss")
  
  # smooth counts to approximate sampling bias
  r_bias <- focal(
    r_occ,
    w = k,
    fun = "sum",
    na.policy = "omit",
    fillvalue = 0
  )
  
  # mask to valid environmental space
  r_bias <- mask(r_bias, template_rast[[1]])
  
  # avoid zero-probability cells inside valid area
  vals <- values(r_bias, mat = FALSE)
  vals[!is.na(vals) & vals <= 0] <- 1e-12
  values(r_bias) <- vals
  
  r_bias
}

# Sample null occurrences preserving sample size and spatial bias
sample_null_occurrences <- function(occ_df, r_bias, template_rast, n = nrow(occ_df)) {
  probs <- values(r_bias, mat = FALSE)
  ok <- which(!is.na(probs) & is.finite(probs) & probs > 0)
  
  samp_cells <- sample(
    x = ok,
    size = n,
    replace = FALSE,
    prob = probs[ok]
  )
  
  xy <- xyFromCell(r_bias, samp_cells)
  
  tibble(lon = xy[, 1], lat = xy[, 2])
}


##### PREPARE OUTPUT DIRECTORIES #####
# List species directories (from model fitting results)
spp <- list.dirs('outputs/models/fitting/',
                 full.names = FALSE,
                 recursive = FALSE)

##### LOAD BEST MODEL SETTINGS #####
# Load previously selected "best" models per species (tune_args.csv files)
fls <- list.files(
  'outputs/models/fitting/',
  'tune_args.csv$',
  recursive = TRUE,
  full.names = TRUE
)

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
mods <- list.dirs(
  'outputs/futureclimate/2061_2080/ssp245/',
  full.names = FALSE,
  recursive = F
)[-1]

# RCP scenarios
rcp <- c('ssp245', 'ssp370', 'ssp585')

# Combination of time x RCP x GCM
cases <- expand_grid(times, rcp, mods) %>% reduce(., paste, sep = '_')

# Load predictors
preds <- rast('outputs/current_predictors.tif')

##### LOAD OCCURRENCE DATA #####
allOccs <- read_csv('outputs/occs/all_occs_flagged.csv')

# Example: split Lasiurus varius by season
lasiurus <- allOccs %>%
  filter(species == 'Lasiurus varius') %>%
  mutate(month = dmy(date) %>% as_date() %>% month()) %>%
  filter(!is.na(month)) %>%
  mutate(month = ifelse(month %in% c(10, 11, 12, 1, 2, 3), 'summer', 'winter')) %>%
  mutate(species = str_c(species, '_', month)) %>%
  split(.$species)

# Split all occurrences by species
occsList <- allOccs %>% split(.$species)
rm(allOccs)

##### COMBINE Lasiurus seasonal splits #####
occsList <- c(occsList)

bg_all <- as.data.frame(preds, xy = T) %>%
  rename(lon = x, lat = y) %>% as_tibble()

##### PARALLELIZATION PREP #####
# Optionally, you could wrap this entire lapply in parallel::mclapply (Linux/Mac only).
##### NULL SETTINGS #####
n_null <- 1000

##### PARALLELIZATION #####
library(furrr)
plan(multisession, workers = 20)

# optional, for reproducible future_map results
furrr_options_seed <- furrr::furrr_options(seed = TRUE)

##### SPECIES LOOP #####
all_results <- map(names(occsList), \(sp) {
  # sp <- names(occsList)[1]
  
  ##### PREPARE DATA FOR THIS SPECIES #####
  sp.dat <- occsList[[sp]]
  r_bias <- make_bias_raster(sp.dat, preds[[1]], sigma = 6)
  
  ftemp <- tempfile(fileext = '.tif')
  writeRaster(r_bias, ftemp)
  ##### LOAD BEST MODEL SETTINGS FOR THIS SPECIES #####
  res_sp <- res_spp %>% filter(species == sp)
  
  ##### NULL REPLICATES #####
  future_map_dfr(1:n_null, \(n) {
    message("Species: ", sp, " | null replicate: ", n)
    # n = 1
    ##### LOAD PREDICTORS #####
    # Load predictors
    preds <- rast('outputs/current_predictors.tif')
    r_bias <- rast(ftemp)
    ##### SAMPLE NULL OCCURRENCES #####
    sp.dat.null <- sample_null_occurrences(
      occ_df = sp.dat,
      r_bias = r_bias,
      template_rast = preds[[1]],
      n = nrow(sp.dat)
    )
    
    ##### SAMPLE BACKGROUND #####
    bgwd <- anti_join(bg_all, sp.dat.null, by = c("lat", "lon")) %>%
      sample_n(10000) %>% mutate(ID = 0) %>% dplyr::select(-lat, -lon)
    
    ##### EXTRACT PREDICTORS #####
    ocwd <- extract(preds, sp.dat.null[, c("lon", "lat")]) %>% na.omit() %>% 
      as_tibble() %>% mutate(ID = 1)
    
    spwd <- rbind(ocwd, bgwd)
    
    ##### FIT AND PROJECT ALL BEST MODELS #####
    map_dfr(1:nrow(res_sp), \(r) {
      # r = 1
      allmods <- list.files(
        "outputs/cmip6_aligned_to_preds/",
        pattern = ".tif$",
        recursive = TRUE,
        full.names = TRUE
      )


      ##### FIT NULL MODEL #####
      mod <- maxent(
        x = spwd[, -1] %>% as.data.frame() %>% remove_rownames(),
        p = spwd$ID %>% as.vector(),
        args = prepPara(
          userfeatures = res_sp$fc[r],
          betamultiplier = res_sp$rm[r],
          doclamp = TRUE,
          randomseed = TRUE,
          outputformat = "cloglog"
        )
      )

      vals_th <- mod@results %>%
        enframe() %>%
        filter(str_detect(name, "Cloglog.threshold")) %>%
        mutate(name = janitor::make_clean_names(name))

      ##### CURRENT CENTROID #####
      cur_mod <- predict(mod,
                         raster::stack(preds),
                         args = c("outputformat=cloglog", "doclamp=true")) %>%
        rast()

      cur_bin <- map(1:nrow(vals_th), ~ cur_mod >= vals_th[[2]][.x]) %>% rast()
      names(cur_bin) <- vals_th$name

      cur_masked <- map(1:nrow(vals_th),
                        ~ mask(cur_mod[[1]], cur_bin[[.x]], maskvalues = FALSE))

      cur_centroids <- map_dfr(cur_masked, suitability_weighted_centroid) %>%
        mutate(threshold = vals_th$name)

      ##### FUTURE PROJECTIONS #####
      map_dfr(times, \(time) {
        # time = times[1]
        map_dfr(rcp, \(rc) {
          # rc = rcp[1]
          modsprj <- allmods[str_detect(allmods, rc)] %>%
            .[str_detect(., gsub("_", "-", time))]

          mods_ord <- list.dirs(
            "outputs/cmip6_aligned_to_preds/",
            recursive = FALSE,
            full.names = FALSE
          )

          names(modsprj) <- mods_ord

          map_dfr(seq_along(modsprj), \(i) {
            # i = 1
            mp <- modsprj[i]
            gcm_name <- names(modsprj)[i]

            fut.preds <- rast(mp)
            fut.preds <- fut.preds[[names(preds)]]
            stopifnot(identical(names(preds), names(fut.preds)))

            fut.mod <- predict(
              mod,
              raster::stack(fut.preds),
              args = c("outputformat=cloglog", "doclamp=true")
            ) %>% rast()

            fut_bin <- map(1:nrow(vals_th), ~ fut.mod >= vals_th[[2]][.x]) %>% rast()
            names(fut_bin) <- vals_th$name

            fut_masked <- map(1:nrow(vals_th),
                              ~ mask(fut.mod[[1]], fut_bin[[.x]], maskvalues = FALSE))

            fut_centroids <- map_dfr(fut_masked, suitability_weighted_centroid) %>%
              mutate(threshold = vals_th$name)

            fut_centroids %>%
              left_join(
                cur_centroids %>%
                  dplyr::select(
                    threshold,
                    lon_current = lon,
                    lat_current = lat
                  ),
                by = "threshold"
              ) %>%
              mutate(
                species = sp,
                fc = res_sp$fc[r],
                rm = res_sp$rm[r],
                rcp = rc,
                time = time,
                model = gcm_name,
                null_n = n,
                shift_km = (lat - lat_current) * 111.32,
                .before = 1
              )
          })
        })
      })
    })
  }, 
  .options = furrr_options_seed)
})
plan(sequential)

saveRDS(all_results, 'outputs/all_centroid_shifts_nulls.rds')
all_results <- readRDS('outputs/all_centroid_shifts_nulls.rds')

list_rbind(all_results) %>% write_csv('outputs/all_centroid_shifts_nulls.csv')

res_all <- list_rbind(all_results)

str(res_all)
res_all$threshold %>% unique()

centroids <- read_csv('outputs/centroids_enm_models.csv')

# 1) clean nulls and keep one latitude per null replicate
null_gcm <- res_all %>%
  rename(ssp = rcp, period = time) %>% 
  filter(!is.na(lat), n_cells > 0, 
         threshold == 'minimum_training_presence_cloglog_threshold') %>%
  dplyr::select(-threshold) %>% 
  group_by(species, model, ssp, period, null_n) %>%
  summarise(
    lat = mean(lat, na.rm = TRUE),
    .groups = "drop"
  ) %>% 
  mutate(period = str_replace(period, '_', '-'))

# 2) observed future centroids
obs_gcm <- centroids %>%   # your real model centroids per GCM
  filter(period!='current') %>% 
  rename(model = gcm) %>% 
  group_by(species, model, ssp, period) %>%
  summarise(
    lat = mean(lat, na.rm = TRUE),
    .groups = "drop"
  )

# 3) optional: p-values to annotate panels
pvals_gcm <- obs_gcm %>%
  rowwise() %>%
  mutate(
    p_value = mean(
      null_gcm$lat[
        null_gcm$species == species &
          null_gcm$model == model &
          null_gcm$ssp == ssp &
          null_gcm$period == period 
      ] <= lat,
      na.rm = TRUE
    )
  ) %>%
  ungroup()

null_summary_gcm <- null_gcm %>%
  group_by(species, model, ssp, period) %>%
  summarise(
    null_mean = mean(lat),
    null_sd   = sd(lat),
    q025 = quantile(lat, 0.025),
    q975 = quantile(lat, 0.975),
    .groups = "drop"
  )

eff_gcm <- obs_gcm %>%
  left_join(null_summary_gcm,
            by = c("species","model","ssp","period")) %>%
  mutate(
    effect = lat - null_mean,
    z_effect = (lat - null_mean) / null_sd
  )

ggplot(eff_gcm, aes(x = period, y = effect, colour = ssp, group = ssp)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point() +
  geom_line() +
  facet_grid(species ~ model, scales = "free_y") +
  theme_bw()

pvals_summary <- pvals_gcm %>%
  mutate(sig = p_value < 0.05) %>%
  group_by(species, ssp, period) %>%
  summarise(
    prop_sig = mean(sig),
    .groups = "drop"
  )


# null distribution per GCM/species/period/ssp
null_gcm %>% 
  filter(ssp == 'ssp370') %>% 
  ggplot(aes(x = lat)) +
  geom_histogram(bins = 35, fill = "grey80", color = "black") +
  geom_vline(
    data = obs_gcm %>% 
      filter(ssp == 'ssp370'),
    aes(xintercept = lat),
    colour = "red",
    linewidth = 0.8
  ) +
  facet_grid(species ~ model + ssp + period, scales = "free_x") +
  theme_bw() +
  labs(
    x = "Future centroid latitude",
    y = "Null frequency"
  )

eff_plot <- eff_gcm %>%
  mutate(
    ymin = q025 - null_mean,
    ymax = q975 - null_mean
  )

ggplot(eff_plot, aes(x = period, y = effect, colour = ssp, group = ssp)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_pointrange(
    aes(ymin = ymin, ymax = ymax),
    position = position_dodge(width = 0.3)
  ) +
  facet_grid(species ~ model, scales = "free_y") +
  theme_bw() +
  labs(
    x = NULL,
    y = "Effect size (observed latitude - mean null latitude)"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

eff_gcm %>% 
  filter(species == 'Myotis chiloensis') %>%
  summary()

ggplot(eff_gcm, aes(x = period, y = z_effect, colour = ssp, group = ssp)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_hline(yintercept = c(-1.96, 1.96), linetype = 3, colour = "grey40") +
  geom_point() +
  geom_line() +
  facet_grid(species ~ model, scales = "free_y") +
  theme_bw() +
  labs(
    x = NULL,
    y = "Standardized effect size ((obs - null mean) / null SD)"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(pvals_gcm, aes(x = period, y = model, fill = p_value)) +
  geom_tile(color = "white") +
  facet_grid(species ~ ssp) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "gray",
    high = "#B2182B",
    midpoint = 0.1,
    limits = c(0, 1)
  ) +
  theme_bw() +
  labs(
    x = NULL,
    y = NULL,
    fill = "p-value"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(pvals_summary, aes(x = period, y = prop_sig, colour = ssp, group = ssp)) +
  geom_line() +
  geom_point(size = 2) +
  facet_wrap(~ species) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_bw() +
  labs(
    x = NULL,
    y = "Proportion of GCMs with p < 0.05"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(null_gcm, aes(x = lat, colour = model)) +
  geom_density() +
  geom_vline(
    data = obs_gcm,
    aes(xintercept = lat, colour = model),
    linewidth = 0.8,
    linetype = 2
  ) +
  facet_grid(species ~ ssp + period, scales = "free_x") +
  theme_bw() +
  labs(
    x = "Future centroid latitude",
    y = "Density"
  )

eff_plot2 <- eff_gcm %>%
  left_join(
    pvals_gcm %>% dplyr::select(species, model, ssp, period, p_value),
    by = c("species", "model", "ssp", "period")
  ) %>%
  mutate(
    sig_lab = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ ""
    ),
    ymin = q025 - null_mean,
    ymax = q975 - null_mean
  )

ggplot(eff_plot2, aes(x = period, y = effect, colour = ssp, group = ssp)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_pointrange(
    aes(ymin = ymin, ymax = ymax),
    position = position_dodge(width = 0.3)
  ) +
  geom_text(
    aes(label = sig_lab),size=15,
    position = position_dodge(width = 0.3),
    vjust = -0.8,
    show.legend = FALSE
  ) +
  facet_grid(species ~ model, scales = "free_y") +
  theme_bw() +
  labs(
    x = NULL,
    y = "Effect size (observed - null mean latitude)"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
