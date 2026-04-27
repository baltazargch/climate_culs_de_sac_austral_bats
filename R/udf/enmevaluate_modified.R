ENMevaluate_modified <- function (occs, envs = NULL, bg = NULL, tune.args = NULL, partitions = NULL, 
          algorithm = NULL, partition.settings = NULL, other.settings = NULL, 
          categoricals = NULL, doClamp = TRUE, clamp.directions = NULL, 
          user.enm = NULL, user.grp = NULL, occs.testing = NULL, taxon.name = NULL, 
          n.bg = 10000, overlap = FALSE, overlapStat = c("D", "I"), 
          user.val.grps = NULL, user.eval = NULL, rmm = NULL, parallel = FALSE, 
          numCores = NULL, parallelType = "doSNOW", updateProgress = FALSE, 
          quiet = FALSE, occ = NULL, env = NULL, bg.coords = NULL, 
          RMvalues = NULL, fc = NULL, occ.grp = NULL, bg.grp = NULL, 
          method = NULL, bin.output = NULL, rasterPreds = NULL, clamp = NULL, 
          progbar = NULL) 
{
  all.legacy <- list(occ, env, bg.coords, RMvalues, fc, occ.grp, 
                     bg.grp, method, bin.output, rasterPreds)
  if (sum(sapply(all.legacy, function(x) !is.null(x))) > 0) {
    if (quiet != TRUE) 
      message("* Running ENMeval v2.0.4 with legacy arguments. These will be phased out in the next version.")
  }
  if (!is.null(occ)) 
    occs <- occ
  if (!is.null(env)) 
    envs <- env
  if (!is.null(bg.coords)) 
    bg <- bg.coords
  if (!is.null(method)) 
    partitions <- method
  if (!is.null(clamp)) 
    doClamp <- clamp
  if (!is.null(rasterPreds)) {
    stop("The argument \"rasterPreds\" was deprecated. Please leave this argument at its default NULL. If you want to avoid generating model prediction rasters, include predictor variable values in the occurrence and background data frames (SWD format). See Details in ?ENMevaluate for more information.")
  }
  if (!is.null(RMvalues)) 
    tune.args$rm <- RMvalues
  if (!is.null(fc)) 
    tune.args$fc <- fc
  if (!is.null(occ.grp) & !is.null(bg.grp)) {
    user.grp <- list(occs.grp = occ.grp, bg.grp = bg.grp)
    if (quiet != TRUE) 
      message("Warning: These arguments were deprecated and replaced with the argument \"user.grp\".")
  }
  if ((!is.null(occ.grp) & is.null(bg.grp)) | (is.null(occ.grp) & 
                                               !is.null(bg.grp))) {
    stop("For user partitions, please input both arguments \"occ.grp\" and \"bg.grp\". Warning: These are legacy arguments that were replaced with the argument \"user.grp\".")
  }
  if (is.null(algorithm) & is.null(user.enm)) {
    stop("* Please select a model name (argument \"algorithm\") or specify a user model (argument \"user.enm\").")
  }
  start.time <- proc.time()
  if ((requireNamespace("ecospat", quietly = TRUE))) {
    ecospat.use <- TRUE
  }
  else {
    message("Package ecospat is not installed, so Continuous Boyce Index (CBI) cannot be calculated.")
    ecospat.use <- FALSE
  }
  occs <- as.data.frame(occs)
  if (!is.null(bg)) 
    bg <- as.data.frame(bg)
  if (is.null(partition.settings)) 
    partition.settings <- list(orientation = "lat_lon", 
                               aggregation.factor = 2, kfolds = 5)
  other.settings <- c(other.settings, list(abs.auc.diff = TRUE, 
                                           pred.type = "cloglog", validation.bg = "full"))
  other.settings <- c(other.settings, ecospat.use = ecospat.use)
  if (inherits(occs[, 1], "character") | inherits(bg[, 1], 
                                                  "character")) 
    stop("* If first column of input occurrence or background data is the taxon name, remove it and instead include the 'taxon.name' argument. The first two columns must be the longitude and latitude of the occurrence/background localities.")
  if (is.null(taxon.name)) {
    if (quiet != TRUE) 
      message(paste0("*** Running initial checks... ***\n"))
  }
  else {
    if (quiet != TRUE) 
      message(paste0("*** Running initial checks for ", 
                     taxon.name, " ... ***\n"))
  }
  all.partitions <- c("jackknife", "randomkfold", "block", 
                      "checkerboard1", "checkerboard2", "user", "testing", 
                      "none")
  if (!(partitions %in% all.partitions)) {
    stop("Please enter an accepted partition method.")
  }
  if (partitions == "testing" & is.null(occs.testing)) {
    stop("If doing testing evaluations, please provide testing data (occs.testing).")
  }
  if ((partitions == "checkerboard1" | partitions == "checkerboard2") & 
      is.null(envs)) {
    stop("For checkerboard partitioning, predictor variable rasters \"envs\" are required.")
  }
  if (partitions == "randomkfold") {
    if (is.null(partition.settings$kfolds)) {
      stop("For random k-fold partitioning, a numeric, non-zero value of \"kfolds\" is required.")
    }
    else {
      if (partition.settings$kfolds == 0) {
        stop("For random k-fold partitioning, a numeric, non-zero value of \"kfolds\" is required.")
      }
    }
  }
  if (partitions == "testing") {
    if (is.null(occs.testing)) {
      stop("If performing fully withheld testing, enter occs.testing dataset and assign partitions to \"testing\".")
    }
    if (nrow(occs.testing) == 0) {
      stop("If performing fully withheld testing, enter occs.testing dataset and assign partitions to \"testing\".")
    }
  }
  if (is.null(tune.args) & overlap == TRUE) {
    if (quiet != TRUE) 
      message("* As no tuning arguments were specified, turning off niche overlap.")
    overlap <- FALSE
  }
  if (all(names(occs) != names(bg)) & !is.null(bg)) {
    stop("Datasets \"occs\" and \"bg\" have different column names. Please make them identical and try again.")
  }
  if (doClamp == FALSE & !is.null(clamp.directions)) {
    stop("If specifying clamp directions, please make doClamp = TRUE.")
  }
  tune.args.num <- which((sapply(tune.args, class) %in% c("numeric", 
                                                          "integer")) & sapply(tune.args, length) > 1)
  for (i in tune.args.num) {
    tune.args[[i]] <- sort(tune.args[[i]])
  }
  if (!(partitions %in% c("block", "checkerboard1", "checkerboard2", 
                          "user")) & other.settings$validation.bg == "partition") {
    stop("If using non-spatial partitions, please set validation.bg to \"full\". The \"partition\" option only makes sense when partitions represent different regions of the study extent. See ?ENMevaluate for details.")
  }
  else if (partitions == "user" & other.settings$validation.bg == 
           "partition") {
    message("* Please make sure that the user-specified partitions are spatial, else validation.bg should be set to \"full\". The \"partition\" option only makes sense when partitions represent different regions of the study extent. See ?ENMevaluate for details.")
  }
  if (is.null(user.enm)) {
    enm <- ENMeval:::lookup.enm(algorithm)
  }
  else {
    enm <- user.enm
  }
  if (!is.null(user.grp)) {
    user.grp.uniq <- unique(c(user.grp$occs.grp, user.grp$bg.grp))
  }
  enm@errors(occs, envs, bg, tune.args, partitions, algorithm, 
             partition.settings, other.settings, categoricals, doClamp, 
             clamp.directions)
  if (!is.null(envs)) {
     envs <- raster::stack(envs)
  #   envs.z <- raster::values(envs)
  #   envs.naMismatch <- sum(apply(envs.z, 1, function(x) !all(is.na(x)) & 
  #                                  !all(!is.na(x))))
  #   if (envs.naMismatch > 0) {
  #     if (quiet != TRUE) 
  #       message(paste0("* Found ", envs.naMismatch, 
  #                      " raster cells that were NA for one or more, but not all, predictor variables. Converting these cells to NA for all predictor variables."))
  #     envs.names <- names(envs)
  #     envs <- raster::stack(raster::calc(envs, fun = function(x) if (sum(is.na(x)) > 
  #                                                                    0) 
  #       x * NA
  #       else x))
  #     names(envs) <- envs.names
  #   }
    if (is.null(bg)) {
      if (quiet != TRUE)
        message(paste0("* Randomly sampling ", n.bg,
                       " background points ..."))
      bg <- as.data.frame(dismo::randomPoints(envs, n = n.bg))
      names(bg) <- names(occs)
    }
    occs.cellNo <- raster::extract(envs, occs, cellnumbers = TRUE)
    occs.dups <- duplicated(occs.cellNo[, 1])
    if (sum(occs.dups) > 0)
      if (quiet != TRUE)
        message(paste0("* Removed ", sum(occs.dups),
                       " occurrence localities that shared the same grid cell."))
    occs <- occs[!occs.dups, ]
    if (!is.null(user.grp))
      user.grp$occs.grp <- user.grp$occs.grp[!occs.dups]
    occs.z <- raster::extract(envs, occs)
    bg.z <- raster::extract(envs, bg)
    occs <- cbind(occs, occs.z)
    bg <- cbind(bg, bg.z)
  }
  else {
    
    if (is.null(bg)) {
      stop("* If inputting variable values without rasters, please make sure to input background coordinates with values as well as occurrences.")
    }
    if (quiet != TRUE) 
      message("* Variable values were input along with coordinates and not as raster data, so no raster predictions can be generated and AICc is calculated with background data for Maxent models.")
    if (ncol(occs) < 3 | ncol(bg) < 3) 
      stop("* If inputting variable values without rasters, please make sure these values are included in the occs and bg tables proceeding the coordinates.")
  }
  occs.z.na <- which(rowSums(is.na(occs)) > 0)
  if (length(occs.z.na) > 0) {
    if (quiet != TRUE) 
      message(paste0("* Removed ", length(occs.z.na), 
                     " occurrence points with NA predictor variable values."))
    occs <- occs[-occs.z.na, ]
    if (!is.null(user.grp)) 
      user.grp$occs.grp <- user.grp$occs.grp[-occs.z.na]
  }
  bg.z.na <- which(rowSums(is.na(bg)) > 0)
  if (length(bg.z.na) > 0) {
    if (quiet != TRUE) 
      message(paste0("* Removed ", length(bg.z.na), " background points with NA predictor variable values."))
    bg <- bg[-bg.z.na, ]
    if (!is.null(user.grp)) 
      user.grp$bg.grp <- user.grp$bg.grp[-bg.z.na]
  }
  d <- rbind(occs, bg)
  d$pb <- c(rep(1, nrow(occs)), rep(0, nrow(bg)))
  if (!is.null(user.grp)) {
    d[d$pb == 1, "grp"] <- as.numeric(as.character(user.grp$occs.grp))
    d[d$pb == 0, "grp"] <- as.numeric(as.character(user.grp$bg.grp))
    if (!all(user.grp.uniq %in% d$grp)) 
      stop("Removal of cell duplicates caused one or more user partition groups to be missing. Please make sure all partition groups are represented by at least one non-duplicate occurrence record.")
    d$grp <- factor(d$grp)
  }
  if (partitions == "testing") {
    if (!is.null(envs)) {
      occs.testing.z <- as.data.frame(raster::extract(envs, 
                                                      occs.testing))
      occs.testing.z <- cbind(occs.testing, occs.testing.z)
    }
    else {
      occs.testing.z <- occs.testing
    }
  }
  else {
    occs.testing.z <- NULL
  }
  if (!is.null(envs)) {
    categoricals <- unique(c(categoricals, names(envs)[which(raster::is.factor(envs))]))
  }
  else {
    categoricals <- unique(c(categoricals, names(occs)[which(sapply(occs, 
                                                                    is.factor))]))
  }
  if (length(categoricals) == 0) 
    categoricals <- NULL
  if (!is.null(categoricals)) {
    for (i in 1:length(categoricals)) {
      if (quiet != TRUE) 
        message(paste0("* Assigning variable ", categoricals[i], 
                       " to categorical ..."))
      d[, categoricals[i]] <- as.factor(d[, categoricals[i]])
      if (!is.null(user.val.grps)) 
        user.val.grps[, categoricals[i]] <- factor(user.val.grps[, 
                                                                 categoricals[i]], levels = levels(d[, categoricals[i]]))
      if (!is.null(occs.testing.z)) 
        occs.testing.z[, categoricals[i]] <- factor(occs.testing.z[, 
                                                                   categoricals[i]], levels = levels(d[, categoricals[i]]))
    }
  }
  other.settings$categoricals <- categoricals
  if (doClamp == TRUE) {
    if (!is.null(envs)) {
      if (is.null(clamp.directions)) {
        clamp.directions$left <- names(envs)
        clamp.directions$right <- names(envs)
      }
      other.settings$clamp.directions <- clamp.directions
      envs <- clamp.vars(orig.vals = envs, ref.vals = rbind(occs.z, 
                                                            bg.z), left = clamp.directions$left, right = clamp.directions$right, 
                         categoricals = categoricals)
      if (quiet != TRUE) 
        message("* Clamping predictor variable rasters...")
    }
    else {
      if (is.null(clamp.directions)) {
        clamp.directions$left <- names(d[, 3:(ncol(d) - 
                                                1)])
        clamp.directions$right <- names(d[, 3:(ncol(d) - 
                                                 1)])
      }
    }
  }
  other.settings$doClamp <- FALSE
  d.occs <- d %>% dplyr::filter(pb == 1) %>% dplyr::select(1:2)
  d.bg <- d %>% dplyr::filter(pb == 0) %>% dplyr::select(1:2)
  grps <- switch(partitions, jackknife = get.jackknife(d.occs, 
                                                       d.bg), randomkfold = get.randomkfold(d.occs, d.bg, partition.settings$kfolds), 
                 block = get.block(d.occs, d.bg, partition.settings$orientation), 
                 checkerboard1 = get.checkerboard1(d.occs, envs, d.bg, 
                                                   partition.settings$aggregation.factor), checkerboard2 = get.checkerboard2(d.occs, 
                                                                                                                             envs, d.bg, partition.settings$aggregation.factor), 
                 user = NULL, testing = list(occs.grp = rep(0, nrow(d.occs)), 
                                             bg.grp = rep(0, nrow(d.bg))), none = list(occs.grp = rep(0, 
                                                                                                      nrow(d.occs)), bg.grp = rep(0, nrow(d.bg))))
  parts.message <- switch(partitions, jackknife = "* Model evaluations with k-1 jackknife (leave-one-out) cross validation...", 
                          randomkfold = paste0("* Model evaluations with random ", 
                                               partition.settings$kfolds, "-fold cross validation..."), 
                          block = paste0("* Model evaluations with spatial block (4-fold) cross validation and ", 
                                         partition.settings$orientation, " orientation..."), 
                          checkerboard1 = "* Model evaluations with checkerboard (2-fold) cross validation...", 
                          checkerboard2 = "* Model evaluations with hierarchical checkerboard (4-fold) cross validation...", 
                          user = paste0("* Model evaluations with user-defined ", 
                                        length(unique(user.grp$occs.grp)), "-fold cross validation..."), 
                          testing = "* Model evaluations with testing data...", 
                          none = "* Skipping model evaluations (only calculating full model statistics)...")
  if (quiet != TRUE) 
    message(parts.message)
  if (partitions == "jackknife") 
    other.settings$cbi.cv <- FALSE
  else other.settings$cbi.cv <- TRUE
  if (partitions == "user") {
    user.nk <- length(unique(user.grp$occs.grp))
    partition.settings$kfolds <- user.nk
    if (user.nk == nrow(d[d$pb == 1, ])) 
      other.settings$cbi.cv <- FALSE
    else other.settings$cbi.cv <- TRUE
  }
  if (!is.null(grps)) 
    d$grp <- factor(c(grps$occs.grp, grps$bg.grp))
  if (is.null(taxon.name)) {
    if (quiet != TRUE) 
      message(paste("\n*** Running ENMeval v2.0.4 with", 
                    enm@msgs(tune.args, other.settings), "***\n"))
  }
  else {
    if (quiet != TRUE) 
      message(paste("\n*** Running ENMeval v2.0.4 for", 
                    taxon.name, "with", enm@msgs(tune.args, other.settings), 
                    "***\n"))
  }
  tune.tbl <- expand.grid(tune.args, stringsAsFactors = FALSE) %>% 
    tibble::as_tibble()
  if (nrow(tune.tbl) == 0) 
    tune.tbl <- NULL
  if (parallel) {
    results <- ENMeval:::tune.parallel(d, envs, enm, partitions, tune.tbl, 
                             doClamp, other.settings, partition.settings, user.val.grps, 
                             occs.testing.z, numCores, parallelType, user.eval, 
                             algorithm, quiet)
  }
  else {
    results <- ENMeval:::tune.regular(d, envs, enm, partitions, tune.tbl, 
                            doClamp, other.settings, partition.settings, user.val.grps, 
                            occs.testing.z, updateProgress, user.eval, algorithm, 
                            quiet)
  }
  train.stats.all <- dplyr::bind_rows(lapply(results, function(x) x$train.stats))
  val.stats.all <- dplyr::bind_rows(lapply(results, function(x) x$cv.stats))
  if (is.null(tune.tbl)) {
    tune.names <- enm@name
  }
  else {
    tune.tbl <- dplyr::mutate_all(tune.tbl, as.factor)
    tune.names <- train.stats.all$tune.args
    tune.tbl$tune.args <- factor(tune.names, levels = tune.names)
  }
  mod.full.all <- lapply(results, function(x) x$mod.full)
  names(mod.full.all) <- tune.names
  if (!is.null(envs)) {
    mod.full.pred.all <- raster::stack(sapply(results, function(x) x$mod.full.pred))
    names(mod.full.pred.all) <- tune.names
  }
  else {
    mod.full.pred.all <- raster::stack()
  }
  if (partitions %in% c("testing", "none")) {
    nk <- 0
  }
  else {
    nk <- length(unique(d[d$pb == 1, "grp"]))
  }
  if (nk > 0) {
    nset <- ifelse(!is.null(tune.tbl), ncol(tune.tbl), 0)
    if (partitions == "jackknife") {
      sum.list <- list(avg = mean, sd = ~sqrt(corrected.var(., 
                                                            nk)))
    }
    else {
      sum.list <- list(avg = mean, sd = sd)
    }
    if (nk == 1 | partitions == "testing") 
      sum.list <- list(function(x) {
        x
      })
    if (!is.null(tune.tbl)) 
      val.stats.all$tune.args <- factor(val.stats.all$tune.args, 
                                        levels = tune.names)
    cv.stats.sum <- val.stats.all %>% dplyr::group_by(tune.args) %>% 
      dplyr::select(-fold) %>% dplyr::summarize_all(sum.list) %>% 
      dplyr::ungroup()
    names(cv.stats.sum) <- gsub("(.*)_(.*)", "\\1.\\2", 
                                names(cv.stats.sum))
    cv.stats.sum <- cv.stats.sum[, order(colnames(cv.stats.sum))]
    if (!is.null(tune.tbl)) {
      train.stats.all$tune.args <- factor(train.stats.all$tune.args, 
                                          levels = tune.names)
      eval.stats <- tune.tbl %>% dplyr::left_join(train.stats.all, 
                                                  by = "tune.args") %>% dplyr::left_join(cv.stats.sum, 
                                                                                         by = "tune.args")
    }
    else {
      train.stats.all$tune.args <- NULL
      cv.stats.sum$tune.args <- NULL
      eval.stats <- dplyr::bind_cols(train.stats.all, 
                                     cv.stats.sum)
    }
  }
  else {
    train.stats.all$tune.args <- factor(train.stats.all$tune.args, 
                                        levels = tune.names)
    eval.stats <- dplyr::left_join(tune.tbl, train.stats.all, 
                                   by = "tune.args")
    if (nrow(val.stats.all) > 0) 
      eval.stats <- dplyr::left_join(eval.stats, val.stats.all, 
                                     by = "tune.args")
    if ("fold" %in% names(eval.stats)) 
      eval.stats <- eval.stats %>% dplyr::select(-fold)
  }
  ncoefs <- sapply(mod.full.all, enm@ncoefs)
  if ((enm@name == "maxnet" | enm@name == "maxent.jar")) {
    pred.type.raw <- switch(enm@name, maxnet = "exponential", 
                            maxent.jar = "raw")
    aic.settings <- other.settings
    aic.settings$pred.type <- pred.type.raw
    if (!is.null(envs)) {
      pred.all.raw <- raster::stack(lapply(mod.full.all, 
                                           enm@predict, envs, aic.settings))
    }
    else {
      pred.all.raw <- NULL
    }
    occs.pred.raw <- dplyr::bind_rows(lapply(mod.full.all, 
                                             enm@predict, occs[, -c(1, 2)], aic.settings))
    aic <- aic.maxent(occs.pred.raw, ncoefs, pred.all.raw)
    eval.stats <- dplyr::bind_cols(eval.stats, aic)
  }
  eval.stats$ncoef <- ncoefs
  if (is.null(taxon.name)) 
    taxon.name <- ""
  if (is.null(tune.tbl)) 
    tune.tbl <- data.frame()
  if (is.null(occs.testing.z)) 
    occs.testing.z <- data.frame()
  if (partitions != "block") 
    partition.settings$orientation <- NULL
  if (partitions != "checkerboard1" | partitions != "checkerboard2") 
    partition.settings$aggregation.factor <- NULL
  if (partitions != "randomkfold") 
    partition.settings$kfolds <- NULL
  if (is.null(partition.settings) | length(partition.settings) == 
      0) 
    partition.settings <- list()
  if (is.null(clamp.directions)) 
    clamp.directions <- list()
  variable.importance.all <- lapply(mod.full.all, enm@variable.importance)
  other.settings$doClamp <- NULL
  e <- ENMevaluation(algorithm = enm@name, tune.settings = as.data.frame(tune.tbl), 
                     results = as.data.frame(eval.stats), results.partitions = val.stats.all, 
                     predictions = mod.full.pred.all, models = mod.full.all, 
                     variable.importance = variable.importance.all, partition.method = partitions, 
                     partition.settings = partition.settings, other.settings = other.settings, 
                     doClamp = doClamp, clamp.directions = clamp.directions, 
                     taxon.name = as.character(taxon.name), occs = d[d$pb == 
                                                                       1, 1:(ncol(d) - 2)], occs.testing = occs.testing.z, 
                     occs.grp = factor(d[d$pb == 1, "grp"]), bg = d[d$pb == 
                                                                      0, 1:(ncol(d) - 2)], bg.grp = factor(d[d$pb == 0, 
                                                                                                             "grp"]), rmm = list())
  e@rmm <- buildRMM(e, envs, rmm)
  if (overlap == TRUE) {
    nr <- raster::nlayers(e@predictions)
    if (nr == 0) {
      if (quiet != TRUE) 
        message("Warning: calculate niche overlap without model prediction rasters.")
    }
    else if (nr == 1) {
      if (quiet != TRUE) 
        message("Warning: only 1 model prediction raster found. Need at least 2 rasters to calculate niche overlap. Increase number of tuning arguments and run again.")
    }
    else {
      for (ovStat in overlapStat) {
        if (quiet != TRUE) 
          message(paste0("Calculating niche overlap for statistic ", 
                         ovStat, "..."))
        predictions.noNegs <- raster::calc(e@predictions, 
                                           function(x) {
                                             x[x < 0] <- 0
                                             x
                                           })
        overlap.mat <- calc.niche.overlap(predictions.noNegs, 
                                          ovStat, quiet)
        e@overlap[[ovStat]] <- overlap.mat
      }
    }
  }
  timed <- proc.time() - start.time
  t.min <- floor(timed[3]/60)
  t.sec <- timed[3] - (t.min * 60)
  if (quiet != TRUE) 
    message(paste("ENMevaluate completed in", t.min, "minutes", 
                  round(t.sec, 1), "seconds."))
  return(e)
}
