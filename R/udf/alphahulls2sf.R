library(igraph)
library(alphahull)

#taken and modified from: https://babichmorrowc.github.io/post/2019-03-18-alpha-hull/

arc2line <- function(center, r, vector, theta, npoints = 100) {
  # Get the angles at the extremes of the arcs
  angles <- anglesArc(vector, theta)
  # Generate sequence of angles along the arc to determine the points
  seqang <- seq(angles[1], angles[2], length = npoints)
  # Generate x coordinates for points along the arc
  x <- center[1] + r * cos(seqang)
  # Generate y coordinates for points along the arc
  y <- center[2] + r * sin(seqang)
  coords.xy <- cbind(x,y)
  line <- Line(coords = coords.xy)
  return(line)
}

ahull2lines <- function(x){
 
  arclist <- x$arcs
  # arclist <- arclist[!apply(arclist[,c(3:6)], 1, \(x) all(x == 0) ),]
  lines <- list()
  for (i in 1:nrow(arclist)) {
    # Extract the attributes of arc i
    center_i <- arclist[i, 1:2]
    radius_i <- arclist[i, 3]
    vector_i <- arclist[i, 4:5]
    theta_i <- arclist[i, 6]
    # Convert arc i into a Line object
    line_i <- arc2line(center = center_i, r = radius_i, vector = vector_i, theta = theta_i)
    list_length <- length(lines)
    if(list_length > 0){
      # If a line has already been added to the list of lines
      # Define last_line_coords as the coordinates of the last line added to the list before the ith line
      last_line_coords <- lines[[list_length]]@coords
    }
    if(i == 1){
      # Add the first line to the list of lines
      lines[[i]] <- line_i
    } else if(isTRUE(all.equal(line_i@coords[1,], last_line_coords[nrow(last_line_coords),]))){
      # If the first coordinate in the ith line is equal to the last coordinate in the previous line
      # then those lines should be connected
      # Row bind the coordinates for the ith line to the coordinates of the previous line in the list
      lines[[list_length]]@coords <- rbind(last_line_coords, line_i@coords[2:nrow(line_i@coords),])
    } else {
      # If the first coordinate in the ith line does not match the last coordinate in the previous line
      # then the ith line represents a new line
      # Add the ith line to the list as a new element
      lines[[length(lines) + 1]] <- line_i
    }
  }
  # Convert the list of lines to a Line object
  lines <- Lines(lines, ID = 'l')
  # Convert the Line object to a SpatialLines object
  sp_lines <- SpatialLines(list(lines))
  return(sp_lines)
}

spLines2poly <- function(x){
  
  # x <-  ahull2lines(hull)
  # Extract the lines slot
  lines_slot <- x
  
  # Create SpatialPolygons
  sp_polys <- try(st_as_sf(lines_slot))
  
  if(!inherits(sp_polys, 'try-error')) {
    
    sp_polys <-  st_make_valid(sp_polys) %>% st_cast(, 'MULTLINESTRING', do_split = TRUE)
    sp_polys <- sp_polys[ !st_is(sp_polys, 'POINT'), ] %>% st_union()
    st_crs(sp_polys) <- 4326
    
    crds <- st_coordinates(sp_polys)
    crds <- crds[ !duplicated(crds[,1:2]),1:2]
    lon <- crds[,1]
    lat <- crds[,2]
    
    pol = st_polygon(
      list(
        cbind(c(lon, lon[1]), c(lat, lat[1]))
      )
    )
    polc = st_sfc(pol, crs=4326)
    
    return(polc)
  } else {
    stop()
  # 
  # # Create a list of booleans indicating whether a given Line represents a polygon
  # poly_bool <- sapply(lines_slot@Lines, \(x){
  #   # x <- lines_slot@Lines[[2]]
  #   coords <- x@coords
  #   # Check if the first coordinate in the line is the same as the last
  #   dist <- st_distance(st_sfc(st_point(coords[1,])), st_sfc(st_point(coords[nrow(coords),]))) %>% round(6)
  #   out <- dist <= 0
  # 
  #   if(all(coords[-1,1] == coords[1,1])) return(FALSE) else out
  # 
  #   })
  # # Pull out the lines that form polygons
  # poly_lines <- sp_lines
  # poly_lines@lines[[1]]@Lines <- poly_lines@lines[[1]]@Lines[ poly_bool ]
  # 
  # # Create SpatialPolygons
  # sp_polys <- SpatialPolygons(
  #   list(
  #     Polygons(
  #       lapply(poly_lines@lines[[1]]@Lines, \(x) {Polygon(x@coords)}), 
  #       ID = "1")
  #     )
  #   )
  # sp_polys <- st_as_sfc(sp_polys) %>% st_cast(., 'POLYGON') 
  # 
  # return(sp_polys)
  }
}

ahull2poly <- function(x){
  # Convert the alpha hull to SpatialLines
  # x = try(alphahull::ahull(sf::st_coordinates(x), alpha = alpha), silent = TRUE)
  # x=hull
  hull2SpatialLines <- ahull2lines(x)
  # Convert SpatialLines to SpatialPolygon
  SpatialLines2SpatialPolygon <- spLines2poly(hull2SpatialLines) %>%  st_make_valid() %>% st_cast('POLYGON')
  zeroarea <- round(st_area(SpatialLines2SpatialPolygon) %>% as.numeric, 2) == 0
  
  if(any( zeroarea) )SpatialLines2SpatialPolygon <- SpatialLines2SpatialPolygon[ !zeroarea]
  
  return(SpatialLines2SpatialPolygon)
}
