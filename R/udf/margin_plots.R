my_margin_plot <- function (model, env, layer, bg, points,
                            standardize = TRUE, verbose = FALSE) 
{
  test=F
  if(test){
    env = spM1
    layer = 'bio_2'
    points = occs[,c('lon', 'lat')]
    standardize = TRUE
    verbose = FALSE
    bg = bg
    model = enm@models[[ best$tune.args[1] ]]
  }
  if (!layer %in% names(env)) {
    stop(paste("Couldn't find layer named", layer, "in environmental rasters!"))
  }
  
  minmax <- terra::minmax(env)
  if (any(is.na(minmax))) {
    env <- terra::setMinMax(env)
    message("\n\nSetting min and max for environment layers...\n\n")
  }
  names <- layer
  minmax <- terra::minmax(env[[layer]])
  plot.df <- seq(minmax[1, ], minmax[2, ], length = 100)
  
  for (i in names(env)) {
    if (i != layer) {
      layer.values <- terra::extract(env[[i]], points, 
                                     ID = FALSE)
      plot.df <- cbind(plot.df, rep(mean(unlist(layer.values), 
                                         na.rm = TRUE), 100))
      names <- c(names, i)
    }
  }
  
  if (standardize == TRUE) {
    plot.df <- data.frame(plot.df)
  }
  colnames(plot.df) <- names
  minmax <- terra::minmax(env[[layer]])
  
  pres.env <- unlist(terra::extract(env[[layer]], points, ID = FALSE))
  pres.dens <- density(pres.env, from = minmax[1, ], to = minmax[2, 
  ], n = 100, na.rm = TRUE)$y
  pres.dens <- pres.dens/max(pres.dens)
  
  bg.env <- unlist(terra::extract(env[[layer]], bg, ID = FALSE))
  bg.dens <- density(bg.env, from = minmax[1, ], to = minmax[2, 
  ], n = 100, na.rm = TRUE)$y
  bg.dens <- bg.dens/max(bg.dens)
  
  pred <- predict(model, x = plot.df, type = "response") 
  
  pred <- pred/max(pred)
  
  plot.df.long <- data.frame(
    layer = c(plot.df[, layer], 
              plot.df[, layer], 
              plot.df[, layer]), 
    value = c(pred, 
              pres.dens, 
              bg.dens), 
    source = c(rep("Suitability", 100), 
               rep("Presence", 100), 
               rep("Background", 100)))
  
  response.plot <- ggplot(data = plot.df.long, aes(x = layer, y = value)) + 
    geom_line(aes(colour = source, linetype = source)) + 
    xlab(layer) + ylab("Value") + theme_bw() + 
    scale_color_manual(values = c("green4", "red", "blue")) + 
    scale_linetype_manual(values = c("dashed", "twodash", "solid")) + 
    theme(plot.title = element_text(hjust = 0.5)) + 
    theme(plot.title = element_text(hjust = 0.5)) + 
    theme(legend.title = element_blank())
  
  return(response.plot)
}
