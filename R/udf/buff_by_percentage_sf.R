calculate_buffer_value <- function(map, f=NULL, d=NULL, 
                                   min.train=4, p=1.1,
                                   max.train=100,
                                   calc.new.area=FALSE, 
                                   verbose=FALSE, 
                                   plot=FALSE, 
                                   transform = TRUE){
  require(sf)
  require(tidyverse)
  # listmaps <- read_sf('D:/GIS/Vectores/Paises/Colombia/AP/runap2/runap2Polygon.shp')
  # listmaps <- split(listmaps, listmaps$nombre)
  # 
  # map = listmaps[[1]]; min.train=5; max.train=50; p = 1.5;
  # 
  # f=NULL; d=NULL; calc.new.area = TRUE; verbose=TRUE; transform=TRUE;plot=TRUE
  if(verbose) calc.new.area <- TRUE
  stopifnot(inherits(map, 'sf'), nrow(map)>0)
  
  if(transform) map <- st_transform(map, "ESRI:54009")
  
  a <- ini <- st_area(map) %>% sum() %>% as.numeric()
  tar <- (a * p) %>% sum() %>% as.numeric()
  
  if(is.null(f)) f = sqrt((tar-a) * 0.001)
  if(is.null(d)) d = f/2
  
  buff <- 0
  area <- a %>% as.numeric()
  
  n = 1
  
  if(st_crs(map) == st_crs('ESRI:54009')){
    while(n < min.train){ #tar < max(area)
      a <- st_area(st_buffer(map, f)) %>% sum()
      # print👎
      n = n+1
      buff[n] <- f
      area[n] <- a %>% as.numeric()
      f <- f+d
    }
  } else {
    while((!all(n > min.train, tar < max(area))) & n < max.train){
      a <- st_area(st_buffer(map, f)) %>% sum()
      # print👎
      n = n+1
      buff[n] <- f
      area[n] <- a %>% as.numeric()
      f <- f+d
    }
  }
  
  if(plot) {
    plot(buff, area, ylim = c(ini, max(tar+10000, max(area))))
    lines(buff, area, ylim = c(ini, max(tar+10000, max(area))), type='l')
    abline(h=tar, col='red')
    text(x=median(buff), y=tar-(tar*0.005), labels = 'Target area')
  }
  
  fit <- lm(buff ~ area)
  nbuff <- coef(fit)[1] + tar * coef(fit)[2]
  
  if(calc.new.area) newarea <- st_area(st_buffer(map, nbuff)) %>% sum()
  
  if(verbose){
    cat(
      paste( 'Area inicial:', ini, '\n', 
             'Area objetivo:', tar, '\n', 
             'Area lograda:', newarea, '\n', 
             'Diferencia objetivo lograda(%):', 
             round(100 - as.numeric(abs((newarea/tar)*100)),2), '%\n\n'
      )
    )
  }
  return(nbuff %>% as.numeric())
}


buffer_by_percentage <- function(x, p=1.1, ...){
  require(sf)
  if(st_crs(x) == st_crs('ESRI:54009')) transform = FALSE else transform = TRUE
  
  buff_value54009 <- calculate_buffer_value(x, p=p, transform = transform, ...)
  
  if(st_crs(x) != st_crs('ESRI:54009')) x <- st_transform(x, "ESRI:54009")
  
  x <- st_buffer(x, buff_value54009)
  # x <- st_transform(x, 4326)
  return(x)
}