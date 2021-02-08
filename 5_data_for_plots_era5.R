
# A refaire et à connecter avec le script précédent
args = commandArgs(trailingOnly=TRUE)
stopifnot(length(args) == 1)
var = args[1]
print(var)


var_future <- readRDS(
  file = paste0("output/", var, "_var_futur_era5.rds")
)

mmcombi <- readRDS(file = paste0("output/", var, "_mmcombi_era5.rds"))


bias_map <- 
  lapply(mmcombi, function(mmc){
    mmc$Bias
  })

mae_bias_map <- 
  lapply(mmcombi, function(mmc){
    abs(mmc$Bias)
  })

mae_bias <-
  sapply(mmcombi, function(mmc){
    mean(abs(mmc$Bias))
  })


gradient_mse <- function(ref, var){
  mse <- matrix(0, nrow(ref), ncol(ref))

  ileft <- 1:(nrow(ref)-1)
  iright <- 2:nrow(ref)
  itop <- 1:(ncol(ref)-1)
  ibottom <- 2:ncol(ref)
  
  mse[ileft, ] <- mse[ileft,] + (var[ileft, ] - var[iright, ] - ref[ileft, ] + ref[iright, ])^2
  mse[iright, ] <- mse[iright,] + (var[iright, ] - var[ileft, ] - ref[iright, ] + ref[ileft, ])^2
  mse[, itop] <- mse[, itop] + (var[, itop] - var[, ibottom] - ref[, itop] + ref[, ibottom ])^2
  mse[, ibottom] <- mse[, ibottom] + (var[, ibottom] - var[, itop] - ref[, ibottom] + ref[, itop])^2
  
  return(mse)
}


mse_gradient_map <- 
  lapply(mmcombi, function(mmc){
    gradient_mse(var_future[,, 1], mmc$Var)
  })

mse_gradient <- 
  sapply(mse_gradient_map, function(mmc){
    mean(mmc)
  })




saveRDS(mae_bias,  paste0("output/", var, "_mae_bias_era5.rds"))
saveRDS(mae_bias_map,  paste0("output/", var, "_mae_bias_map_era5.rds"))
saveRDS(bias_map,  paste0("output/", var, "_bias_map_era5.rds"))
saveRDS(mse_gradient,  paste0("output/", var, "_mse_gradient_era5.rds"))
saveRDS(mse_gradient_map,  paste0("output/", var, "_mse_gradient_map_era5.rds"))
