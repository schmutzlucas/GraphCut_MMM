mmm <- function(ref, var) {

  width = ncol(ref) ### Number of longitudes
  height = nrow(ref) ### Number of latitudes
  nlabs = dim(var)[3] ### Number of longitudes
  
  # Preparing the data to perform the MMM
  var_mmm <- apply(var, 1:2, mean)
  bias_mmm <- var_mmm - ref
  bias_var_mmm <- list("Var" = var_mmm, "Bias" = bias_mmm)
  
  return(bias_var_mmm)
}
