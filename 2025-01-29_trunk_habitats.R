##################
## Calculate stem volume
## from presegmented point clouds
## by Julian Frey
## 2025-01-29
##################

# requires RTools https://cran.r-project.org/bin/windows/Rtools/

## necessary packages 
# concaveman
# lidR
# conicfit
# sf
# Rcpp

# ToDos
# - mean coordinate of evry tree
# - mean surface variation l3 / (l1 + l2 + l3)
# - mean planarity (l2 - l3) / l1
# - mean stem volume
# - mean surface area
# - mean concave perimeter vs. convex perimeter
# - write a loop to process multiple trees

library(lidR)

# read the cpp script
Rcpp::sourceCpp("./eigen_decomposition.cpp")

# function to calculate the area of a convex hull
convhull_area <- function(xy, method = c("convex", "concave")){
  method <- method[1]
  xy <- xy |> as.data.frame()
  if(nrow(xy) < 3){
    return(NA)
  }
  if(method == "concave"){
    ch <- concaveman::concaveman(as.matrix(xy), concavity = 0.1)
    return(abs(0.5 * sum(ch[,1] * c(tail(ch[,2], -1), head(ch[,2], 1)) - c(tail(ch[,1], -1), head(ch[,1], 1)) * ch[,2])))
  } else if(method == "convex"){
    ch <- chull(xy)
    return(abs(0.5 * sum(xy[ch,1] * c(tail(xy[ch,2], -1), head(xy[ch,2], 1)) - c(tail(xy[ch,1], -1), head(xy[ch,1], 1)) * xy[ch,2])))
  } else {
    stop("method must be either 'convex' or 'concave'")
  }
  
}

# na to zero
na2zero <- function(x){
  x[is.na(x)] <- 0
  return(x)
}

# function to calculate stem volume
stem_volume <- function(las, slice_height = 0.1, method = c("convex", "concave")) {
  z_range <- range(las$Z)
  z_slices <- seq(z_range[1], z_range[2], slice_height)
  stem_volume <- 0
  for(i in 1:(length(z_slices) - 1)){
    slice <- filter_poi(las, Z >= z_slices[i] & Z < z_slices[i + 1])@data
    if(nrow(slice) > 0){
      stem_volume <- stem_volume + na2zero(convhull_area(slice[,c("X","Y")], method = method[1]) * slice_height)
    }
  }
  return(stem_volume)
}

# function to compute convex/concave perimeter
convex_perimeter <- function(slice, method = c("convex", "concave"), dTolerance = 0.05, concavity = 1){
  method <- method[1]
  xy <- slice@data[,c("X","Y")] |> as.data.frame()
  if(nrow(xy) < 3){
    return(NA)
  }
  if(method == "concave"){
    peri <- xy |> sf::st_as_sf(coords = c("X","Y")) |> 
      concaveman::concaveman(concavity = concavity) |>
      sf::st_cast("MULTILINESTRING") |> 
      sf::st_simplify(dTolerance = dTolerance) |> 
      sf::st_length() |> sum()
    return(peri)
  } else if(method == "convex"){
    ch <- chull(xy)
    ch <- c(ch, ch[1])
    return(sum(sqrt(diff(xy[ch,1])^2 + diff(xy[ch,2])^2)))
  } else {
    stop("method must be either 'convex' or 'concave'")
  }
}

# function to compute the concave/convex stem surface area
stem_surface_area <- function(las, slice_height = 0.1, method = c("convex", "concave"), dTolerance = 0.01, concavity = 2){
  z_range <- range(las$Z)
  z_slices <- seq(z_range[1], z_range[2], slice_height)
  stem_surface_area <- 0
  for(i in 1:(length(z_slices) - 1)){
    slice <- filter_poi(las, Z >= z_slices[i] & Z < z_slices[i + 1])
    if(nrow(slice) > 0){
      stem_surface_area <- stem_surface_area + na2zero(convex_perimeter(slice, method = method[1], dTolerance = dTolerance, concavity = concavity) * slice_height)
    }
  }
  return(stem_surface_area)
}

# function to compute the mean coordinate
mean_coordinate <- function(las){
  return(c(mean(las$X), mean(las$Y), mean(las$Z)))
}


# function to compute the eigenvalues
trunk_features <- function(las, k = 10L, n_cores = 1) {
  # check if inputs of the right type
  if (!lidR::is(las,"LAS")) {
    stop('las has to be a LAS object.')
  }
  if(!(as.integer(k) == k & length(k) == 1 & k > 0 )) {
    stop('k has to be one positive integer.')
  }
  # necessary for raster_geometry
  # returns geometric features based on eigenvalues
  eigen <- eigen_decomposition(las, k, n_cores) # k neighbours, n cores
  las <- las |>
    add_lasattribute(eigen[,3] / (eigen[,1] + eigen[,2] + eigen[, 3]), 'Curvature', 'curvature') |>
    add_lasattribute((eigen[,2] - eigen[,3]) / eigen[,1], 'Planarity', 'planarity') |> 
    add_lasattribute(1 - abs(eigen[,4]) ,'Verticality','verticality')
  return(las)
}

# test

# # Load presegmented point cloud
# las = lidR::readLAS("X:/Weidbuchen/Edited/pointclouds/2024-12-13weidbuche9_greinwald.las", select = "xyz")
# 
# las |> mean_coordinate()
# 
# las |> stem_volume(method = "convex")
# las |> stem_volume(method = "concave")
# 
# las |> stem_surface_area(method = "convex")
# las |> stem_surface_area(method = "concave")
# 
# las |> trunk_features() |> plot(color = "Curvature")
# 
# test_cylinder <- as.data.frame(cbind(conicfit::calculateCircle(0,0,1, steps = 10000), runif(10000)))
# names(test_cylinder) <- c("X","Y","Z")
# test_cylinder <- LAS(test_cylinder) 
# 
# # plot(test_cylinder)
# # 
# test_cylinder |> stem_volume(method = "concave")
# test_cylinder |> stem_volume(method = "convex")
# 
# test_cylinder |> convex_perimeter(method = "convex")
# test_cylinder |> convex_perimeter(method = "concave")
# 
# test_cylinder |> stem_surface_area(method = "convex")
# test_cylinder |> stem_surface_area(method = "concave")

# final loop
files <- list.files("X:/Weidbuchen/Edited/pointclouds/", pattern = ".las", full.names = TRUE)
results <- data.frame(file = character(), X = numeric(), Y = numeric(), Z = numeric(), volume_convex = numeric(), volume_concave = numeric(), surface_area_convex = numeric(), surface_area_concave = numeric(), mean_curvature = numeric(), mean_planarity = numeric(), mean_verticality = numeric(), stringsAsFactors = FALSE)

for( f in files){
  las <- lidR::readLAS(f, select = "xyz")
  las <- trunk_features(las)
  results <- rbind(results, data.frame(file = f, mean_coordinate(las), stem_volume(las, method = "convex"), stem_volume(las, method = "concave"), stem_surface_area(las, method = "convex"), stem_surface_area(las, method = "concave"), mean(las$Curvature), mean(las$Planarity), mean(las$Verticality)))
}

