#' Predict values from the mark distribution
#'
#' @param sim_realization a realization of simulated values
#' @param raster_list a list of rasters
#' @param size_model a predictive model
#' @param bounds a vector of domain bounds (2 for x, 2 for y)
#' @param correction type of correction to apply ("none" or "toroidal")
#'
#' @return a vector of predicted mark values
#' @export
#'
predict_marks <- function(sim_realization, raster_list, size_model, bounds, correction = "none"){

  s <- sim_realization[,c("x", "y")]
  raster_trans <- scale_rasters(raster_list)
  X <- extract_covars(x = s, raster_list = raster_trans)
  X$x <- s[,1]
  X$y <- s[,2]
  X$time <- sim_realization[,1]

  X$near.nbr.dist <- NA
  X$near.nbr.num <- NA
  X$avg.nbr.dist.15 <- NA
  X$near.nbr.time <- NA
  X$near.nbr.time.all <- NA
  X$near.nbr.time.dist.ratio <- NA
  if((correction == "none")) {
    distance.matrix <- base::as.matrix(stats::dist(s, method = "euclidean"))
  }else if(correction == "toroidal") {
    distance.matrix <- toroidal_dist_matrix_optimized(s, bounds[2] - bounds[1], bounds[4] - bounds[3])
  }

  for(i in 1:base::nrow(X)){
    close.points.15 <- base::unique(base::which(distance.matrix[i,] < 15 & distance.matrix[i,] != 0))
    close.times.15 <- X$time[close.points.15]
    X$near.nbr.dist[i] <- base::min(distance.matrix[i,][-i])
    X$near.nbr.num[i] <- base::length(close.points.15)
    X$avg.nbr.dist.15[i] <- base::mean(distance.matrix[i,][close.points.15])
    if(base::length(close.points.15) == 0){
      X$avg.nbr.dist.15[i] <- base::min(distance.matrix[i,][-i])
    }
    X$near.nbr.time[i] <- X$time[base::unique(base::which(distance.matrix[i,] == X$near.nbr.dist[i]))]
    X$near.nbr.time.all[i] <- mean(close.times.15)
    if(base::length(close.points.15) == 0){
      X$near.nbr.time.all[i] <- X$time[base::unique(base::which(distance.matrix[i,] == X$near.nbr.dist[i]))]
    }
    X$near.nbr.time.dist.ratio[i] <- X$near.nbr.time[i]/X$near.nbr.dist[i]
  }

  return(stats::predict(size_model, X))
}
