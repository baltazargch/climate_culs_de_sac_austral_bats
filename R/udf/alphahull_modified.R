ahull_modified <- function (x, fraction = 0.95, partCount = 1, buff = 10000, initialAlpha = 3, 
                            coordHeaders = c("Longitude", "Latitude"), clipToCoast = "terrestrial", 
                            alphaIncrement = 1, verbose = FALSE, alphaCap = 400) 
{
  require(rangeBuilder)
  # x = sp.occ@coords
  # clipToCoast = "terrestrial"
  # verbose=T
  # alphaCap = 400
  # alphaIncrement = 4
  # fraction = 1
  # buff = 10000
  # initialAlpha = 3
  # partCount = 1
  
  if (clipToCoast == FALSE) clipToCoast <- "no"
  clipToCoast <- match.arg(clipToCoast, c("no", "terrestrial", "aquatic"))
  
  if (ncol(x) == 2) coordHeaders <- c(1, 2)
  x <- x[!duplicated(x[, coordHeaders]), coordHeaders]
  x <- x[stats::complete.cases(x), ]
  
  if (nrow(x) < 3) {
    stop("This function requires a minimum of 3 unique coordinates (after removal of duplicates).")
  }
  
  while ((x[1, coordHeaders[1]] == x[2, coordHeaders[1]] & 
          x[2, coordHeaders[1]] == x[3, coordHeaders[1]]) | 
         (x[1, 2] == x[2, coordHeaders[2]] & 
          x[2, coordHeaders[2]] == x[3, coordHeaders[2]])) {
    x <- x[sample(1:nrow(x), size = nrow(x)), ]
  }
  x <- sf::st_as_sf(as.data.frame(x), coords = 1:2, crs = 4326)
  
  if (nrow(x) < 3) stop("This function requires a minimum of 3 unique coordinates.")
  
  if( nrow(x) == 3) {
    message("Only 3 non-duplicated coordinates. Returning a MCH")
    hull <- sf::st_convex_hull(st_union(x))
    hull <- sf::st_transform(hull, crs = "+proj=eqearth")
    hull <- sf::st_buffer(hull, dist = buff)
    hull <- sf::st_transform(hull, crs = 4326)
    alphaVal = "MCH"
    return(list(hull, alpha = paste("alpha", alphaVal, sep = "")))
    break
  }
  
  alpha <- initialAlpha
  problem <- FALSE
  
  if (verbose) cat("\talpha:", alpha, "\n")
  
  hull <- try(alphahull::ahull(sf::st_coordinates(x), alpha = alpha), silent = TRUE)

  if(inherits(hull, "try-error") & any(grepl("duplicate points", hull)) |
     inherits(hull, "try-error") & hull[1]== "Error in shull.deltri(x1, y1) : error counting arcs!\n") {
    ptDist <- sf::st_distance(x, x)
    diag(ptDist) <- NA
    units(ptDist) <- NULL
  }
  dropPt <- c()
  while(inherits(hull, "try-error") & any(grepl("duplicate points", hull)) | 
        inherits(hull, "try-error") & hull[1]== "Error in shull.deltri(x1, y1) : error counting arcs!\n") {
    closest <- which(ptDist == min(ptDist, na.rm = TRUE), arr.ind = TRUE)
    
    dropPt <- c(dropPt, closest[1,1] %>% unname())
    hull <- try(alphahull::ahull(sf::st_coordinates(x)[-dropPt], alpha = alpha), silent = TRUE)
    print(dropPt)
    if (inherits(hull, "ahull")) {
      break
    } else {
      ptDist[closest[1,1],closest[1,2]] <- NA
    }
  }
  
  if (length(dropPt) > 0) x <- x[-dropPt, ]
  
  hull <- try(alphahull::ahull(sf::st_coordinates(x), alpha = alpha), silent = TRUE)
  
  while (inherits(hull, "try-error")) {
    if (verbose) cat("\talpha:", alpha, "\n")
    alpha <- alpha + alphaIncrement
    hull <- try(alphahull::ahull(sf::st_coordinates(x), alpha = alpha), silent = TRUE)
    
    if (alpha > alphaCap) {
      problem <- TRUE
      break
    }
  }
  
  while (any(hull$arcs[,c(3:6)] == 0)) {
    if (verbose) cat("\talpha:", alpha, "\n")
    alpha <- alpha + alphaIncrement
    hull <- try(alphahull::ahull(sf::st_coordinates(x), alpha = alpha), silent = TRUE)
    
    if (alpha > alphaCap) {
      problem <- TRUE
      break
    }
  }
  
  
  suppressMessages({
    if (!problem) {
      hull <- ahull2poly(hull)
      st_crs(hull) <- 4326
      
      while (is.null(hull) | 
             inherits(hull, "try-error") | 
             !all(sf::st_is_valid(hull))) {
        
        alpha <- alpha + alphaIncrement
        
        if (verbose) cat("\talpha:", alpha, "\n")
        
        hull <- try(ahull2poly(alphahull::ahull(sf::st_coordinates(x), 
                                                alpha = alpha)), silent = TRUE)
        st_crs(hull) <- 4326
      }
      
      pointWithin <- sf::st_intersects(x, hull)
      alphaVal <- alpha
      buffered <- FALSE
      buff <- units::set_units(buff, "m")
      
      while (any(
        length(hull) > partCount, 
        (sum(lengths(pointWithin))/nrow(x)) < fraction, 
        !all(sf::st_is_valid(hull)), 
        inherits(hull, "try-error"),       
        alpha <= alphaCap
        )
        ) {
        
        if(alpha > 150) alphaIncrement <- 2
        alpha <- alpha + alphaIncrement
        if (verbose) cat("\talpha:", alpha, "\n") 
        
        hull <- try(alphahull::ahull(sf::st_coordinates(x), alpha = alpha), silent = TRUE)
        
        hull <- ahull2poly(hull)
        st_crs(hull) <- 4326
        
        pointWithin <- sf::st_intersects(x, hull)
        if((sum(lengths(pointWithin))/nrow(x)) < fraction & length(hull) > partCount) break
        if(alpha > alphaCap) break
      }
      
      if (!inherits(hull, "try-error") & alpha <= alphaCap) {
        hull <- sf::st_transform(hull, crs = "+proj=eqearth")
        hull <- sf::st_buffer(hull, dist = buff)
        hull <- sf::st_transform(hull, crs = 4326)
      }
      
      alphaVal <- alpha
      if (alpha >= alphaCap) {
        hull <- sf::st_convex_hull(st_union(x))
        hull <- sf::st_transform(hull, crs = "+proj=eqearth")
        hull <- sf::st_buffer(hull, dist = buff)
        hull <- sf::st_transform(hull, crs = 4326)
        buffered <- TRUE
        alphaVal = "MCH"
      }
    } else {
      hull <- sf::st_convex_hull(st_union(x))
      hull <- sf::st_transform(hull, crs = "+proj=eqearth")
      hull <- sf::st_buffer(hull, dist = buff)
      hull <- sf::st_transform(hull, crs = 4326)
      buffered <- TRUE
      alphaVal = "MCH"
    }
  })
  
  if (!buffered) {
    hull <- sf::st_transform(hull, crs = "+proj=eqearth")
    hull <- sf::st_buffer(hull, dist = buff)
    hull <- sf::st_transform(hull, crs = 4326)
  }
  if (clipToCoast != "no") {
    world <- rangeBuilder:::loadWorldMap()
    if (clipToCoast == "terrestrial") {
      hull <- sf::st_intersection(hull, world)
    }
    else {
      hull <- sf::st_difference(hull, world)
    }
  }
  
  return(list(hull, alpha = paste("alpha", alphaVal, sep = "")))
}
