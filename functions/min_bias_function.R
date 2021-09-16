# This functions takes a list of models and maps for each pixel the model with 
# the least difference compared to a model of reference
# returns the label attribution of the method "naive pixel perfect"
min_bias <- function(ref, var){
  
  width <- ncol(ref) ### Number of longitudes
  height <- nrow(ref) ### Number of latitudes
  nlabs <- dim(var)[3] # Number of labels used in the GC

  bias <- array(0, c(height, width, nlabs)) 
  
  for(i in 1:nlabs){
    bias[,,i] <- var[,,i] - ref
  }
  
  label_attribution <- matrix(0, height, width)
  for(x in 1:height){
    for(y in 1:width){
      label_attribution[x,y] = which.min(abs(bias[x,y,]))-1
    }
  }
  
  return(label_attribution)
}
