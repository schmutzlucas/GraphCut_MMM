# 


bias_var <- function(ref_future, var_future, labeling){
  width <- ncol(var_future)   ### Number of longitudes
  height <- nrow(var_future)
  nlabs <- dim(var_future)[3] ### Number of longitudes
  
  # allocating the memory 
  var_gc <- bias_gc <- matrix(0,nrow = height,ncol = width)
  
  for(l in 0:(nlabs-1)){
    islabel <- which(labeling == l)
    var_gc[islabel] <- var_future[,,(l+1)][islabel]
  }
  bias_gc <- if(is.null(ref_future)) NULL else var_gc -  ref_future
  bias_var <- list("Var" = var_gc, "Bias" = bias_gc)
  return(bias_var)
}
