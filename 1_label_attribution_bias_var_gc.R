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


year_present <- 1979:2008
year_future <- 2071:2100

# Names of models used in the ensemble
imodels <- c(
  1, 2, 3, 4, 7,
  8, 10, 12, 13, 14, 
  17, 20, 22, 24, 26,
  30, 31, 33, 34, 37
)

model_names <- readRDS("input/models_list.rds")[imodels, 2]

for(i in 1:length(model_names)){
  print(model_names[i])
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
for(m in 1:length(model_names)){
  mmcombi[[m]] <- list()
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

  ##### Min bias #####
  minbias_result <- min_bias(
    ref = var_present[,, m],
    var = var_present[,, -m]
  )
  minbias_bias_var <- bias_var(
    ref_future = var_future[,, m],
    var_future = var_future[,, -m],
    labeling = minbias_result
    )
  mmcombi[[m]][["min_bias"]] <- c(minbias_bias_var, list("Label attribution" = minbias_result))

  ##### MMM #####
  mmcombi[[m]][["mmm"]] <- mmm(
    ref = var_future[,, m],
    var = var_future[,, -m]
  )

  ##### MMM optimized #####
  mmcombi[[m]][["mmmo"]] <- mmmo(
    ref_present = var_present[,, m],
    ref_future = var_future[,, m],
    var_present = var_present[,, -m],
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

