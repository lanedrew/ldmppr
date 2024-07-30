#' Predict values from the mark distribution
#'
#' @param sim_realization a realization of simulated values
#' @param raster_list list of rasters
#' @param size_model a predictive model
#'
#' @return a vector of predictions
#' @export
#'
predict_marks <- function(sim_realization, raster_list, size_model){

  s <- sim_realization[,c("x", "y")]
  raster_trans <- scale_rasters(raster_list)
  X <- extract_covars(x = s, raster_list = raster_trans)
  X$x <- s[,1]
  X$y <- s[,2]
  X$time <- sim_realization[,1]

  X$near.nbr.dist <- NA
  X$near.nbr.num <- NA
  X$avg.nbr.dist.15 <- NA
  X$near.nbr.size <- NA
  X$near.nbr.size.all <- NA
  X$near.nbr.size.dist.ratio <- NA
  distance.matrix <- base::as.matrix(stats::dist(s, method = "euclidean"))

  for(i in 1:base::nrow(X)){
    close.points.15 <- base::unique(base::which(distance.matrix[i,] < 15 & distance.matrix[i,] != 0))
    close.times.15 <- X$time[close.points.15]
    X$near.nbr.dist[i] <- min(distance.matrix[i,][-i])
    X$near.nbr.num[i] <- length(close.points.15)
    X$avg.nbr.dist.15[i] <- mean(distance.matrix[i,][close.points.15])
    if(length(close.points.15) == 0){
      X$avg.nbr.dist.15[i] <- min(distance.matrix[i,][-i])
    }
    X$near.nbr.time[i] <- X$time[base::unique(base::which(distance.matrix[i,] == X$near.nbr.dist[i]))]
    X$near.nbr.time.all[i] <- mean(close.times.15)
    if(length(close.points.15) == 0){
      X$near.nbr.time.all[i] <- X$time[base::unique(base::which(distance.matrix[i,] == X$near.nbr.dist[i]))]
    }
    X$near.nbr.size.dist.ratio[i] <- X$near.nbr.time[i]/X$near.nbr.dist[i]
  }

  return(stats::predict(size_model, X))
}
