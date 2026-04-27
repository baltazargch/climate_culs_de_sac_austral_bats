library(tidyverse)
library(tidysdm)
library(terra)
# UDF calculate eval summary metrics
calc_eval_metrics <- function(list.csv){
  reduce(list.csv, \(x, y) left_join(x, y, by = "name")) %>%
    rowwise() %>%
    mutate(
      mean.val   = mean(c_across(contains('value')), na.rm = TRUE),
      sd.val     = sd(c_across(contains('value')), na.rm = TRUE),
      min.val    = min(c_across(contains('value')), na.rm = TRUE),
      max.val    = max(c_across(contains('value')), na.rm = TRUE),
      n.val      = sum(!is.na(c_across(contains('value')))),
      se.val     = sd.val / sqrt(n.val),
      cv.val     = sd.val / mean.val,                 # coefficient of variation
      range.val  = max.val - min.val,
      median.val = median(c_across(contains('value')), na.rm = TRUE),
      iqr.val    = IQR(c_across(contains('value')), na.rm = TRUE),
      q05.val    = quantile(c_across(contains('value')), 0.05, na.rm = TRUE, names = FALSE),
      q95.val    = quantile(c_across(contains('value')), 0.95, na.rm = TRUE, names = FALSE)
    ) %>%
    ungroup() %>% 
    select(-starts_with('value')) %>% return()
}


dir.mods <-  list.files('outputs/projections/', pattern = 'Results.csv$', full.names = T, recursive = T)

dir.mods <- dir.mods[!str_detect(dir.mods, 'nouse')]

grid <- expand.grid(
  species = list.dirs('outputs/models/fitting/',recursive = F, full.names = F), 
  period = c('current', list.dirs('outputs/cmip6_aligned_to_preds/ACCESS-CM2/ssp245/', F, F)), 
  ssp = str_c('ssp', c('245', '370', '585'))
) %>% as_tibble %>% arrange(species, period, ssp) %>% 
  mutate(across(everything(), as.character)) %>% 
  mutate(ssp = ifelse(period == 'current', 'current', ssp)) %>% 
  distinct(species, period, ssp, .keep_all = T) %>% 
  arrange(species, ssp, period)

curr.grid <- grid %>% 
  filter(period =='current')

current.metrics <-  map(1:nrow(curr.grid), \(i){
  # i=1
  sp.grid = curr.grid[i,]
  
  sp.files = dir.mods[ str_detect(dir.mods, sp.grid$species)]
  sp.files = sp.files[ ! str_detect(sp.files, 'ssp') ]
  sp.files = sp.files[ ! str_detect(sp.files, 'nulls') ]
  
  csv <- map(sp.files, \(x){
    # x = sp.files[3]
    read_csv(x, show_col_types = F) %>% select(-Species) %>% pivot_longer(cols = everything()) %>%
      filter(str_detect(name, 'AUC|contribution|importance') | name == 'Minimum training presence Cloglog threshold')
  })
  
  eval.sp <- calc_eval_metrics(csv) %>% 
    mutate(species = sp.grid$species, 
           period = sp.grid$period, 
           ssp = sp.grid$ssp, .before = 1)
  
  csv <- map(sp.files, \(x){
    # x = sp.files[1]
    read_csv(x, show_col_types = F) %>% select(-Species) %>% pivot_longer(cols = everything()) %>%
      filter(name == 'Minimum training presence Cloglog threshold')
  }) 
  
  fcrm <- str_extract(sp.files, '/LQ.*./') %>% 
    str_remove('/current/') %>% 
    str_remove('/') 
  
  list(eval.sp = eval.sp, 
       th = csv %>% list_rbind %>%
         mutate(species = sp.grid$species, 
                period = sp.grid$period, 
                ssp = sp.grid$ssp, 
                fcrm = fcrm, .before = 1))
})

current.metrics %>% map(~.x$eval.sp) %>% list_rbind() %>% 
  select(-ssp) %>% 
  write_csv('current_eval_full_model_metrics.csv')

th <- current.metrics %>%  map(~.x$th) %>% list_rbind()

future.grid <- anti_join(grid, curr.grid, by = 'period')

fut.files <-  list.files('outputs/projections/', pattern = '.tif$', full.names = T, recursive = T)
fut.files <- fut.files[!str_detect(fut.files, 'nouse|current')]

curr.files <- list.files("outputs/projections/", "binary_pred_all_thresholds.tif", full.names = T, recursive = T)
dir.create('outputs/projections/thresholded_maps/', recursive = T, showWarnings = F)

curr.areas <- map(curr.grid$species, \(sp){
  # sp = curr.grid$species[1]
  curr.map <- curr.files[str_detect(curr.files, sp)]
  curr.map <- curr.map[ !str_detect(curr.map, 'nouse')]
  
  curr.map <- map(curr.map, ~rast(.x, lyrs = 'minimum_training_presence_cloglog_threshold'))
  
  rout <- mean(rast(curr.map))
  writeRaster(rout, filename = paste0('outputs/projections/thresholded_maps/', gsub(' ', '_', sp), '_current.tif'), overwrite = T)
  
  curr.vars <- map_dfr(curr.map, ~expanse(.x, unit="km", byValue=TRUE, wide=TRUE) %>% 
                         rename(model = layer) %>% 
                         mutate(model = 'current') %>% 
                         select(-`0`))
  
  tibble(
    species = sp, 
    period ='current', 
    ssp='current', 
    mean.val   = mean(curr.vars[,2], na.rm = TRUE),
    sd.val     = sd(curr.vars[,2], na.rm = TRUE),
    min.val    = min(curr.vars[,2], na.rm = TRUE),
    max.val    = max(curr.vars[,2], na.rm = TRUE),
    n.val      = sum(!is.na(curr.vars[,2])),
    se.val     = sd.val / sqrt(n.val),
    cv.val     = sd.val / mean.val,                 # coefficient of variation
    range.val  = max.val - min.val,
    median.val = median(curr.vars[,2], na.rm = TRUE),
    iqr.val    = IQR(curr.vars[,2], na.rm = TRUE),
    q05.val    = quantile(curr.vars[,2], 0.05, na.rm = TRUE, names = FALSE),
    q95.val    = quantile(curr.vars[,2], 0.95, na.rm = TRUE, names = FALSE)
  ) %>% return()
  
}) %>% list_rbind()

curr.areas              

future.areas <- map(1:nrow(future.grid), \(i){
  # i=1
  sp.grid = future.grid[i,]
  sp.files <- fut.files[ str_detect(fut.files, sp.grid$species)] %>% 
    .[str_detect(., gsub('-', '_', sp.grid$period))] %>%  .[str_detect(., sp.grid$ssp)] %>% 
    .[str_detect(., 'binaries')]
  
  names(sp.files) <- sp.files %>% basename() %>% str_remove('_.*')
  sp.list <- split(sp.files, names(sp.files))
  
  sp.th <- th %>% filter(species == sp.grid$species) %>% split(.$fcrm)
  
  mod.maps <- map(sp.list, \(x){
    # x = sp.list[[1]]
    fc_rm <- stringr::str_match(
      basename(x),
      "^[^_]+_([A-Z]+_[0-9.]+)_"
    )[,2]
    
    map2(x, sp.th[ fc_rm ], \(f, t){
      # f = x[1]
      # t = sp.th[[1]]
      
      r = rast(f, lyrs = 'minimum_training_presence_cloglog_threshold') 
      bin_r = r >= t$value
      names(bin_r) <- t$fcrm
      bin_r
    })
  })
  
  mod.maps <- map(mod.maps, \(x){
    names(x) <- map_chr(x, varnames)
    x
  })
  dir.create('outputs/projections/thresholded_maps/gcm/')
  mod.maps.flat <- purrr::flatten(mod.maps)
  
  names(mod.maps.flat) <- paste0(sp.grid$species, '_', names(mod.maps.flat))
  writeRaster(rast(mod.maps.flat), 
              filename = paste0('outputs/projections/thresholded_maps/gcm/',
                                names(mod.maps.flat), '.tif'), 
              overwrite = T)
  
  # rout <- mean(rast(mod.maps))
  # names(rout) <- paste0(gsub(' ', '_', sp.grid$species), '_', sp.grid$period, '_', sp.grid$ssp)
  # 
  # writeRaster(rout, filename = paste0('outputs/projections/thresholded_maps/', names(rout), '.tif'), overwrite = T)
  # 
  # fut.vars <- imap_dfr(mod.maps, \(x, n){
  #   expanse(x, unit="km", byValue=TRUE, wide=TRUE) %>% 
  #     rename(name = layer) %>%
  #     mutate(name = n) %>%
  #     select(!all_of('0'))
  # })
  # 
  # tibble(
  #   species = sp.grid$species, 
  #   period =sp.grid$period, 
  #   ssp=sp.grid$ssp, 
  #   mean.val   = mean(fut.vars[,2], na.rm = TRUE),
  #   sd.val     = sd(fut.vars[,2], na.rm = TRUE),
  #   min.val    = min(fut.vars[,2], na.rm = TRUE),
  #   max.val    = max(fut.vars[,2], na.rm = TRUE),
  #   n.val      = sum(!is.na(fut.vars[,2])),
  #   se.val     = sd.val / sqrt(n.val),
  #   cv.val     = sd.val / mean.val,                 # coefficient of variation
  #   range.val  = max.val - min.val,
  #   median.val = median(fut.vars[,2], na.rm = TRUE),
  #   iqr.val    = IQR(fut.vars[,2], na.rm = TRUE),
  #   q05.val    = quantile(fut.vars[,2], 0.05, na.rm = TRUE, names = FALSE),
  #   q95.val    = quantile(fut.vars[,2], 0.95, na.rm = TRUE, names = FALSE)
  # ) %>% return()

}) 

future.areas %>% list_rbind() %>% rbind(., curr.areas) %>% 
  arrange(species, ssp, period) %>% 
  write_csv('outputs/projections/thresholded_maps/all_areas_summary.csv')



