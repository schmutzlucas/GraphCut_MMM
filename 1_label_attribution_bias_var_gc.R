# Loading libraries
# netCDF is the file format the climate models use
library(ncdf4)
library(ncdf4.helpers)

# Loading the functions written by Mathieu and Soulivanh see /functions
source("functions/gc_function.R")
source("functions/bias_var_function.R")
source("functions/mmm_function.R")
source("functions/mmmo_function.R")
source("functions/min_bias_function.R")

# Variable utilisée : "tas" ou "pr"
# args = commandArgs(trailingOnly=TRUE)
# stopifnot(length(args) == 1)
# var = args[1]
# print(var)

# Alternative method 
var = "tas"
# var = "pr" 
stopifnot(length(args) == 1)
print(var)


# Defining time periods 
year_present <- 1979:2008
year_future <- 2071:2100


# Names of models used in the ensemble
# Change the imodels numbers to change model
imodels <- c(
  1, 2, 3, 4, 6,
  8, 10, 12, 13, 14, 
  17, 20, 22, 24, 26,
  30, 31, 33, 34, 37
)
# Store the corresponding model names
model_names <- readRDS("input/models_list.rds")[imodels, 2]


# Read and store the models and metadata
for(i in 1:length(model_names)){
  # print model i name
  print(model_names[i])
  # Load the model with the ncd4 libarby function
  nc <- nc_open(paste0("netcdf/", var, "_", model_names[i], ".nc"))
  # For the first model, load latitude, longitude and create 0-filled arrays 
  #   of dim [longitude, latitude, nb of models]
  if(i == 1) {
    lat <- ncvar_get(nc, "lat")
    lon <- ncvar_get(nc, "lon")
    var_present <- var_future <- array(0,c(length(lon), length(lat), length(model_names)))
  }
  # vector of the time series in the models (1850-2100)
  yyyy <- substr(as.character(nc.get.time.series(nc)), 1, 4)
  # vector of the years defined as present
  iyyyy <- which(yyyy %in% year_present)
  
  # Average of [var] for the present time period 
  ## in apply() 1:2 == margins how?
  var_present[,,i] <- apply(
    ncvar_get(nc, var, start = c(1, 1, min(iyyyy)), count = c(-1, -1, length(iyyyy))),
    1:2, 
    mean
  )
  
  # doing the same for the future time period
  iyyyy <- which(yyyy %in% year_future)
  var_future[,,i] <- apply(
    ncvar_get(nc, var, start = c(1, 1, min(iyyyy)), count = c(-1, -1, length(iyyyy))),
    1:2, 
    mean
  )
  
  # freeing memory of [nc] variable
  nc_close(nc)
}

# naming the dimensions of the data arrays 
dimnames(var_future) <- dimnames(var_present) <- list(lon = lon, lat = lat, model = model_names)


# List of label attribution matrices obtained with GraphCut
# memory allocation for mmcombi
mmcombi <- list()

# loop going through the list of models
# for the min_bias method, each model is successively used as reference
# results are store in [mmcombi->min_bias]

for(m in 1:length(model_names)){
  # allocation for the mth element
  mmcombi[[m]] <- list()

  ##### Min bias #####
  # call of the min_bias function --> see /functions
  # the minus index means we drop the index data
  minbias_result <- min_bias(
    ref = var_present[,, m],
    var = var_present[,, -m]
  )
  
  # This function returns two maps, the results for the variable and the bias 
  # compared with the reference
  minbias_bias_var <- bias_var(
    ref_future = var_future[,, m],
    var_future = var_future[,, -m],
    labeling = minbias_result
  )
  
  # storing the results for the min_bias method for all models and label attribution
  ## 
  mmcombi[[m]][["min_bias"]] <- c(minbias_bias_var, list("Label attribution" = minbias_result))
  
  
  
  #### GraphCut ####
  gc_result <- graphcut(
    ref.datacost = var_present[,, m],
    models.datacost = var_present[,, -m],
    models.smoothcost = var_present[,, -m],
    weight.data = 1,
    weight.smooth = 1 
  )
  gc_bias_var <- bias_var(
    ref_future = var_future[,, m],
    var_future = var_future[,, -m],
    labeling = gc_result$`Label attribution`
  )
  mmcombi[[m]][["gc_present"]] <- c(gc_bias_var, gc_result) 

  gc_result <- graphcut(
    ref.datacost = var_present[,, m],
    models.datacost = var_present[,, -m],
    models.smoothcost = var_future[,, -m],
    weight.data = 1,
    weight.smooth = 1
  )
  gc_bias_var <- bias_var(
    ref_future = var_future[,, m],
    var_future = var_future[,, -m],
    labeling = gc_result$`Label attribution`
  )
  mmcombi[[m]][["gc_hybrid"]] <- c(gc_bias_var, gc_result)

  gc_result <- graphcut(
    ref.datacost = var_future[,, m],
    models.datacost = var_future[,, -m],
    models.smoothcost = var_future[,, -m],
    weight.data = 1,
    weight.smooth = 1
  )
  gc_bias_var <- bias_var(
    ref_future = var_future[,, m],
    var_future = var_future[,, -m],
    labeling = gc_result$`Label attribution`
  )
  mmcombi[[m]][["gc_future"]] <- c(gc_bias_var, gc_result)


  ##### MMM #####
  mmcombi[[m]][["mmm"]] <- mmm(
    ref = var_future[,, m],
    var = var_future[,, -m]
  )

  ##### MMM optimized #####
  mmcombi[[m]][["om_present"]] <- mmmo(
    ref_present = var_present[,, m],
    ref_future = var_future[,, m],
    var_present = var_present[,, -m],
    var_future = var_future[,, -m]
  )
  mmcombi[[m]][["om_future"]] <- mmmo(
    ref_present = var_future[,, m],
    ref_future = var_future[,, m],
    var_present = var_future[,, -m],
    var_future = var_future[,, -m]
  )

}

names(mmcombi) <- model_names
saveRDS(
  mmcombi,
  file = paste0("output/", var, "_mmcombi_perfectmodel.rds")
)
saveRDS(
  var_present,
  file = paste0("output/", var, "_var_present_perfectmodel.rds")
)
saveRDS(
  var_future,
  file = paste0("output/", var, "_var_futur_perfectmodel.rds")
)

