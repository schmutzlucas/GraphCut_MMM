list.of.packages <- c("RcppXPtrUtils","devtools")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
# if(length(new.packages)) install.packages(new.packages,repos = "http://cran.us.r-project.org")
lapply(list.of.packages, library, character.only = TRUE)
# install_github("thaos/gcoWrapR")
library(gcoWrapR)


graphcut <- function(
  ref.datacost,
  models.datacost,
  models.smoothcost,
  weight.data,
  weight.smooth
){
  
  nlabs = dim(models.datacost)[3] # Number of labels used in the GC
  width = ncol(ref.datacost) # Number of latitudes
  height = nrow(ref.datacost) # Number of longitudes
  
    # Instanciation of the GraphCut environment
  
  gco <- new(GCoptimizationGridGraph, width, height, nlabs)
  
  # Preparing the DataCost and SmoothCost functions of the GraphCut in C++
  
  cat("Creating DataCost function...  ")
  ptrDataCost <- cppXPtr(
    code = 'float dataFn(int p, int l, Rcpp::List extraData)
{
  int numPix = extraData["numPix"];
  float weight = extraData["weight"];
  NumericVector data = extraData["data"];
  return(weight * std::abs(data[p + numPix * l]) );
}',
  includes = c("#include <math.h>", "#include <Rcpp.h>"),
  rebuild = FALSE, showOutput = FALSE, verbose = FALSE
  )
  
  cat("Creating SmoothCost function...  ")
  ptrSmoothCost <- cppXPtr(
    code = 'float smoothFn(int p1, int p2, int l1, int l2, Rcpp::List extraData)
{
  int numPix = extraData["numPix"];
  float weight = extraData["weight"];
  NumericVector data = extraData["data"];
  float cost = std::abs(data[p1 + numPix * l1]-data[p1 + numPix * l2])
  + std::abs(data[p2 + numPix * l1] - data[p2 + numPix * l2]) ;
  return(weight * cost);
}',
  includes = c("#include <math.h>", "#include <Rcpp.h>"),
  rebuild = FALSE, showOutput = FALSE, verbose = FALSE
  )
  
  # Preparing the data to perform GraphCut
  bias <- array(0, c(height, width, nlabs)) 
  for(i in 1:nlabs){
    bias[,,i] <- models.datacost[,, i] - ref.datacost
  }
  
  # Permuting longitude and latitude since the indexing isn't the same in R and in C++
  bias_cpp <- c(aperm(bias, c(2, 1, 3)))
  smooth_cpp <- c(aperm(models.smoothcost, c(2, 1, 3)))
  
  gco$setDataCost(ptrDataCost, list(numPix = width * height, data = bias_cpp, weight = weight.data))
  gco$setSmoothCost(ptrSmoothCost, list(numPix = width * height, data = smooth_cpp, weight = weight.smooth))
  
  # Creating the initialization matrix based on the best model from the previous list
  mae_list <- numeric(nlabs)   
  for(i in 1:length(mae_list)){
    mae_list[[i]] <- mean(abs(bias[,,i]))
  }
  best_label <- which.min(mae_list)-1 # in C++ label indices start at 0
  for(z in 0:(length(ref.datacost)-1)){
    gco$setLabel(z, best_label) 
  }
  
  # Optimizing the MRF energy with alpha-beta swap (-1 refers to the optimization until convergence)
  cat("Starting GraphCut optimization...  ")
  begin <- Sys.time()
  gco$swap(-1)
  time_spent <- Sys.time()-begin
  cat("GraphCut optimization done :  ")
  print(time_spent)
  
  
  data_cost <- gco$giveDataEnergy()
  smooth_cost <- gco$giveSmoothEnergy()
  data_smooth_list <- list("Data cost" = data_cost, "Smooth cost" = smooth_cost)
  
  label_attribution <- matrix(0,nrow = height,ncol = width)
  for(j in 1:height){
    for(i in 1:width){
      label_attribution[j,i] <- gco$whatLabel((i - 1) + width * (j - 1)) ### Permuting from the C++ indexing to the R indexing
    }
  }
  
  gc_result <- vector("list",length=2)
  gc_result <- list("Label attribution" = label_attribution, "Data and smooth cost" = data_smooth_list)
  
  return(gc_result)
}
