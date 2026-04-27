library(terra)
library(geodata)
library(envirem)
library(sf)
library(tidyverse)
# =========================
# INPUT: your template raster
# =========================
aream <- read_sf('outputs/aream/new_area_m.geojson')
preds <- rast("inputs/worldclim/2.5m/wc2.1_2.5m_bio_2.tif")  %>% # you already have this in memory
  crop(aream, mask=T)
stopifnot(inherits(preds, "SpatRaster"))

# Use first layer as template grid + mask footprint
preds_template <- preds[[1]]

# A binary mask: 1 where preds exists, NA where outside
preds_mask <- ifel(is.na(preds_template), NA, 1)

# =========================
# USER SETTINGS
# =========================
res_mins <- 2.5

gcms <- c('CMCC-ESM2', 'EC-Earth3-Veg',
          'GISS-E2-1-G', 'INM-CM5-0',
          'ACCESS-CM2', 'MIROC6')

periods <- c('2021-2040', '2041-2060', '2061-2080', '2081-2100')
ssps <- c('245', '370', '585')

base_dir <- "outputs/cmip6_aligned_to_preds"
dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)

# Variables you requested
env_vars <- c("climaticMoistureIndex", "PETColdestQuarter", "PETWettestQuarter", "PETseasonality")
bio_keep <- c("bio2", "bio6", "bio17")  # WorldClim names; will rename to bio_*

# =========================
# HELPERS
# =========================

period_to_envirem_year <- function(period) {
  yrs <- as.integer(strsplit(period, "-", fixed = TRUE)[[1]])
  mid <- round(mean(yrs))
  mid - 1950
}

# Align any raster stack to preds: project -> resample -> mask
align_to_preds <- function(r, template, mask_rast, method = "bilinear") {
  # project to template CRS
  r2 <- project(r, template, method = method)
  # force exact grid match (origin/resolution/extent)
  r3 <- resample(r2, template, method = method)
  # mask to preds footprint (NA outside)
  mask(r3, mask_rast)
}

prepare_master_for_envirem <- function(tmin, tmax, prec) {
  tmean <- (tmin + tmax) / 2
  
  names(tmin)  <- paste0("tmin_",  1:12)
  names(tmax)  <- paste0("tmax_",  1:12)
  names(tmean) <- paste0("tmean_", 1:12)
  names(prec)  <- paste0("prec_",  1:12)
  
  master <- c(tmin, tmax, tmean, prec)
  
  assignNames(
    tmin   = "tmin_##",
    tmax   = "tmax_##",
    tmean  = "tmean_##",
    precip = "prec_##"
  )
  verifyRasterNames(master)
  master
}

# library(furrr)

# plan(multisession, workers=8)
# =========================
# MAIN LOOP
# =========================
walk(gcms, \(gcm) {
  # for (gcm in gcms) {
  for (ssp in ssps) {
    for (period in periods) {
      # gcm =gcms[1]; ssp=ssps[1]; period=periods[1] # for testing
      
      cat("\n----------------------------------------\n")
      cat("GCM:", gcm, " SSP:", ssp, " Period:", period, "\n")
      cat("----------------------------------------\n")
      
      combo_dir <- file.path(base_dir, gcm, paste0("ssp", ssp), period)
      dir.create(combo_dir, recursive = TRUE, showWarnings = FALSE)
      
      # ---- 1) BIOCLIM (download -> subset -> align)
      bioc <- geodata::cmip6_world(
        model = gcm, ssp = ssp, time = period,
        var = "bioc", res = res_mins,
        path = combo_dir
      )
      names(bioc) <- paste0('bio', 1:19)
      bio_sel <- bioc[[bio_keep]]
      names(bio_sel) <- c("bio_2", "bio_6", "bio_17")
      
      # Align bio to preds (bilinear is fine for continuous)
      bio_sel_a <- align_to_preds(bio_sel, preds_template, preds_mask, method = "bilinear")
      
      # ---- 2) Monthly climate for ENVIREM (download -> align)
      tmin <- geodata::cmip6_world(gcm, ssp, period, var="tmin", res=res_mins, path=combo_dir)
      tmax <- geodata::cmip6_world(gcm, ssp, period, var="tmax", res=res_mins, path=combo_dir)
      prec <- geodata::cmip6_world(gcm, ssp, period, var="prec", res=res_mins, path=combo_dir)
      
      tmin_a <- align_to_preds(tmin, preds_template, preds_mask, method = "bilinear")
      tmax_a <- align_to_preds(tmax, preds_template, preds_mask, method = "bilinear")
      prec_a <- align_to_preds(prec, preds_template, preds_mask, method = "bilinear")
      
      # ---- 3) ENVIREM generation on the aligned monthly stacks
      master <- prepare_master_for_envirem(tmin_a, tmax_a, prec_a)
      
      env_year <- period_to_envirem_year(period)
      solrad <- ETsolradRasters(rasterTemplate = master[[1]], year = env_year, outputDir = NULL)
      
      env_stack <- generateEnvirem(
        masterstack = master,
        solradstack = solrad,
        var = env_vars,
        tempScale = 1,
        precipScale = 1
      )
      
      # ENVIREM output should already be aligned (since master/template are),
      # but we’ll enforce mask again just in case:
      env_stack_a <- mask(env_stack, preds_mask)
      
      # ---- 4) Combine + write
      out <- c(env_stack_a, bio_sel_a)
      
      out_file <- file.path(
        combo_dir,
        paste0("aligned_envirem_and_bio_preds.tif")
      )
      writeRaster(out, out_file, overwrite = TRUE)
      cat("Wrote:", out_file, "\n")
      
      rm(bioc, bio_sel, bio_sel_a, tmin, tmax, prec, tmin_a, tmax_a, prec_a,
         master, solrad, env_stack, env_stack_a, out)
      gc()
      unlink(paste0(combo_dir, '/climate'), recursive = TRUE) # clean up raw downloads
    }
  }
})
# plan(sequential) # back to sequential for any remaining code
cat("\nDONE.\n")
