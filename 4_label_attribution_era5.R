library(ncdf4)
library(ncdf4.helpers)

source("functions/gc_function.R")
source("functions/bias_var_function.R")
source("functions/mmm_function.R")
source("functions/mmmo_function.R")
source("functions/min_bias_function.R")

# Variable utilisée : "tas" ou "pr"
args = commandArgs(trailingOnly=TRUE)
stopifnot(length(args) == 1)
var = args[1]
print(var)

year_present <- 1979:1998
year_future <- 1999:2019

# Names of models used in the ensemble
imodels <- c(
  1, 2, 3, 4, 7,
  8, 10, 12, 13, 14, 
  17, 20, 22, 24, 26,
  30, 31, 33, 34, 37
)

model_names <- readRDS("input/models_list.rds")[imodels, 2]
model_names <- c("era5", model_names)

for(i in 1:length(model_names)){
  nc <- nc_open(paste0("netcdf/", var, "_", model_names[i], ".nc"))
  if(i == 1) {
    lat <- ncvar_get(nc, "lat")
    lon <- ncvar_get(nc, "lon")
    var_present <- var_future <- array(0,c(length(lon), length(lat), length(model_names)))
  }
  yyyy <- substr(as.character(nc.get.time.series(nc)), 1, 4)
  iyyyy <- which(yyyy %in% year_present)
  var_present[,,i] <- apply(
    ncvar_get(nc, var, start = c(1, 1, min(iyyyy)), count = c(-1, -1, length(iyyyy))),
    1:2, 
    mean
  )
  iyyyy <- which(yyyy %in% year_future)
  var_future[,,i] <- apply(
    ncvar_get(nc, var, start = c(1, 1, min(iyyyy)), count = c(-1, -1, length(iyyyy))),
    1:2, 
    mean
  )
  nc_close(nc)
}

dimnames(var_future) <- dimnames(var_present) <- list(lon = lon, lat = lat, model = model_names)

# List of label attribution matrices obtained with GraphCut
mmcombi <- list()

##### Min bias #####
minbias_result <- min_bias(
  ref = var_present[,, 1],
  var = var_present[,, -1]
)
minbias_bias_var <- bias_var(
  ref_future = var_future[,, 1],
  var_future = var_future[,, -1],
  labeling = minbias_result
  )
mmcombi[["min_bias"]] <- c(minbias_bias_var, list("Label attribution" = minbias_result))

#### GraphCut ####
gc_result <- graphcut(
  ref.datacost = var_present[,, 1],
  models.datacost = var_present[,, -1],
  models.smoothcost = var_present[,, -1],
  weight.data = 1,
  weight.smooth = 1
)
gc_bias_var <- bias_var(
  ref_future = var_future[,, 1],
  var_future = var_future[,, -1],
  labeling = gc_result$`Label attribution`
)
mmcombi[["gc_present"]] <- c(gc_bias_var, gc_result) 

gc_result <- graphcut(
  ref.datacost = var_present[,, 1],
  models.datacost = var_present[,, -1],
  models.smoothcost = var_future[,, -1],
  weight.data = 1,
  weight.smooth = 1
)
gc_bias_var <- bias_var(
  ref_future = var_future[,, 1],
  var_future = var_future[,, -1],
  labeling = gc_result$`Label attribution`
)
mmcombi[["gc_hybrid"]] <- c(gc_bias_var, gc_result)

gc_result <- graphcut(
  ref.datacost = var_future[,, 1],
  models.datacost = var_future[,, -1],
  models.smoothcost = var_future[,, -1],
  weight.data = 1,
  weight.smooth = 1
)
gc_bias_var <- bias_var(
  ref_future = var_future[,, 1],
  var_future = var_future[,, -1],
  labeling = gc_result$`Label attribution`
)
mmcombi[["gc_future"]] <- c(gc_bias_var, gc_result)

##### MMM #####
mmcombi[["mmm"]] <- mmm(
  ref = var_future[,, 1],
  var = var_future[,, -1]
)

##### MMM optimized #####
mmcombi[["om_present"]] <- mmmo(
  ref_present = var_present[,, 1],
  ref_future = var_future[,, 1],
  var_present = var_present[,, -1],
  var_future = var_future[,, -1]
)
mmcombi[["om_future"]] <- mmmo(
  ref_present = var_future[,, 1],
  ref_future = var_future[,, 1],
  var_present = var_future[,, -1],
  var_future = var_future[,, -1]
)


saveRDS(
mmcombi,
file = paste0("output/", var, "_mmcombi_era5.rds")
)
saveRDS(
var_present,
file = paste0("output/", var, "_var_present_era5.rds")
)
saveRDS(
var_future,
file = paste0("output/", var, "_var_futur_era5.rds")
)

