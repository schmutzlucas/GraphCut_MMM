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

year_future <- 2071:2100

# Names of models used in the ensemble
imodels <- c(
  1, 2, 3, 4, 7,
  8, 10, 12, 13, 14, 
  17, 20, 22, 24, 26,
  30, 31, 33, 34, 37
)

model_names <- readRDS("input/models_list.rds")[imodels, 2]

var_present <- readRDS(paste0("output/", var, "_var_present_era5.rds"))
ref_present <- var_present[,, 1]
var_present <- var_present[,, -1]
var_future <- array(NA, dim(var_present))
dimnames(var_future) <- dimnames(var_present)

for(i in 1:length(model_names)){
  nc <- nc_open(paste0("netcdf/", var, "_", model_names[i], ".nc"))
  if(i == 1) {
    lat <- ncvar_get(nc, "lat")
    lon <- ncvar_get(nc, "lon")
  }
  yyyy <- substr(as.character(nc.get.time.series(nc)), 1, 4)
  iyyyy <- which(yyyy %in% year_future)
  var_future[,,i] <- apply(
    ncvar_get(nc, var, start = c(1, 1, min(iyyyy)), count = c(-1, -1, length(iyyyy))),
    1:2, 
    mean
  )
  nc_close(nc)
}


# List of label attribution matrices obtained with GraphCut
mmcombi <- list()
#### GraphCut ####
gc_result <- graphcut(
  ref.datacost = ref_present,
  models.datacost = var_present,
  models.smoothcost = var_present,
  weight.data = 1,
  weight.smooth = 1
)
gc_bias_var <- bias_var(
  ref_future = NULL,
  var_future = var_future,
  labeling = gc_result$`Label attribution`
)
mmcombi[["gc_present"]] <- c(gc_bias_var, gc_result) 

gc_result <- graphcut(
  ref.datacost = ref_present,
  models.datacost = var_present,
  models.smoothcost = var_future,
  weight.data = 1,
  weight.smooth = 1
)
gc_bias_var <- bias_var(
  ref_future = NULL,
  var_future = var_future,
  labeling = gc_result$`Label attribution`
)
mmcombi[["gc_hybrid"]] <- c(gc_bias_var, gc_result)

##### Min bias #####
minbias_result <- min_bias(
  ref = ref_present,
  var = var_present
)
minbias_bias_var <- bias_var(
  ref_future = NULL,
  var_future = var_future,
  labeling = minbias_result
  )
mmcombi[["min_bias"]] <- c(minbias_bias_var, list("Label attribution" = minbias_result))

##### MMM #####
mmcombi[["mmm"]] <- mmm(
  ref = NULL,
  var = var_future
)

##### MMM optimized #####
mmcombi[["mmmo"]] <- mmmo(
  ref_present = ref_present,
  ref_future = NULL,
  var_present = var_present,
  var_future = var_future
)


saveRDS(
mmcombi,
file = paste0("output/", var, "_mmcombi_era5_proj.rds")
)

combi_names <- names(mmcombi)
var_map <- 
  lapply(mmcombi, function(mmc){
    mmc$Var
  })

##### Maps #####
colors <- c('#ffffcc','#ffeda0','#fed976','#feb24c','#fd8d3c','#fc4e2a','#e31a1c','#bd0026','#800026')
if(var == "tas"){
  colors <- rev(c('#67001f','#b2182b','#d6604d','#f4a582','#fddbc7','#f7f7f7','#d1e5f0','#92c5de','#4393c3','#2166ac','#053061'))
} else {
  colors <- c('#fff7fb','#ece7f2','#d0d1e6','#a6bddb','#74a9cf','#3690c0','#0570b0','#045a8d','#023858')
}
ncolors <- length(colors)

# Var projections

pdf(paste0("figures/", var, "_var_map_era5_proj.pdf"))
zlim <- range(unlist(var_map))
zmin <- zlim[1]
zmax <- zlim[2]
if(var == "tas"){
  scale <- function(x) (x - zmin) / (zmax - zmin)
  unscale <- function(x) x * (zmax - zmin) + zmin
} else {
  scale <- function(x) log(1 + x) / log(1 + zmax)
  unscale <- function(x) exp(x * log(1 + zmax)) - 1 
}
bias <- lapply(var_map, function(x) atan(x)/(pi/2)/zmax)
zlim <- c(0, 1)
breaks <- seq(0, 1, length.out = ncolors + 1)

par(mfrow = c(3, 2))
for(icombi in seq_along(var_map)){
  var <- scale(var_map[[icombi]])
  fields::image.plot(lon, lat, var,
             xlab = '',
             ylab = '',
             main = combi_names[icombi],
             zlim = zlim,
             col=colors,
             breaks = breaks,
             axis.args = list(at=breaks, labels = sprintf("%.2f", unscale(breaks)), cex.axis = 0.5)
  )
  maps::map("world2",add=T)
}
dev.off()
