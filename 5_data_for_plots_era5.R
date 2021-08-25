source("functions/gradients.R")

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


mae_gradient_map <- 
  lapply(mmcombi, function(mmc){
    gradient_mae(var_future[,, 1], mmc$Var)
  })

mae_gradient <- 
  sapply(mae_gradient_map, function(mmc){
    mean(mmc)
  })




saveRDS(mae_bias,  paste0("output/", var, "_mae_bias_era5.rds"))
saveRDS(mae_bias_map,  paste0("output/", var, "_mae_bias_map_era5.rds"))
saveRDS(bias_map,  paste0("output/", var, "_bias_map_era5.rds"))
saveRDS(mae_gradient,  paste0("output/", var, "_mae_gradient_era5.rds"))
saveRDS(mae_gradient_map,  paste0("output/", var, "_mae_gradient_map_era5.rds"))
