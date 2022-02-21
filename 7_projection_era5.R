library(ncdf4)
library(ncdf4.helpers)

source("functions/gc_function.R")
source("functions/bias_var_function.R")
source("functions/mmm_function.R")
source("functions/mmmo_function.R")
source("functions/min_bias_function.R")
source("functions/multimap.R")

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


##### MMM #####
mmcombi[["mmm"]] <- mmm(
  ref = NULL,
  var = var_future
)

##### MMM optimized #####
mmcombi[["om_present"]] <- mmmo(
  ref_present = ref_present,
  ref_future = NULL,
  var_present = var_present,
  var_future = var_future
)

saveRDS(
mmcombi,
file = paste0("output/", var, "_mmcombi_era5_proj.rds")
)

mmcombi <- mmcombi[c("min_bias", "mmm", "gc_hybrid", "om_present")]
combi_names <- names(mmcombi)
var_map_7100 <- 
  lapply(mmcombi, function(mmc){
    mmc$Var
  })

mmcombi_0019 <- readRDS(paste0("output/", var, "_mmcombi_era5.rds"))
var_map_0019 <- 
  lapply(mmcombi_0019[combi_names], function(mmc){
    mmc$Var
  })

var_map_diff <- mapply(
  FUN = function(x, y) x-y,
  x = var_map_7100,
  y = var_map_0019,
  SIMPLIFY = FALSE
)

##### Maps #####
colors <- c('#ffffcc','#ffeda0','#fed976','#feb24c','#fd8d3c','#fc4e2a','#e31a1c','#bd0026','#800026')
if(var == "tas"){
  colors <- rev(c('#67001f','#b2182b','#d6604d','#f4a582','#fddbc7','#f7f7f7','#d1e5f0','#92c5de','#4393c3','#2166ac','#053061'))
} else {
  colors <- c('#fff7fb','#ece7f2','#d0d1e6','#a6bddb','#74a9cf','#3690c0','#0570b0','#045a8d','#023858')
}
ncolors <- length(colors)

# Var projections

pdf(paste0("figures/", var, "_var_map_era5_proj.pdf"), height = 2/3 * 7 + 1)
zlim <- range(unlist(var_map_7100))
zmin <- zlim[1]
zmax <- zlim[2]
if(var == "tas"){
  scale <- function(x) (x - zmin) / (zmax - zmin)
  unscale <- function(x) x * (zmax - zmin) + zmin
} else {
  scale <- function(x) log(1 + x) / log(1 + zmax)
  unscale <- function(x) exp(x * log(1 + zmax)) - 1 
}
maintitle <- sprintf(
  "(%s) %s",
  letters[seq.int(length(combi_names))],
  combi_names
)
subtitle <- rep("", length(combi_names))
lvalues <- lapply(var_map_7100, scale)
zlim <- c(0, 1)
breaks <- seq(zlim[1], zlim[2], length.out = length(colors) + 1)
labels <- sprintf("%.2f", unscale(breaks))
multimap(
  lon, lat, lvalues,
  maintitle, subtitle,
  zlim, colors,
  breaks, labels,
  outersplit = rbind(c(0, 1, 3/17, 1), c(0, 1, 0, 3/17)),
  innersplit = c(2,2),
  mai = c(0.6732, 0.5412, 0.5412, 0.2772),
  legend_par = list(mar=c(8, 1, 1, 1), plt = c(0.1, 0.9, 0.9, 1))
)
dev.off()

##### Diff Maps #####
pal_redwhiteblue <- c('#67001f','#b2182b','#d6604d','#f4a582','#fddbc7','#f7f7f7','#d1e5f0','#92c5de','#4393c3','#2166ac','#053061')
colors <- rev(pal_redwhiteblue)
ncolors <- length(colors)

# Var evolution
pdf(paste0("figures/", var, "_evo_map_era5_proj.pdf"), height = 2/3 * 7 + 1)
zmax <- max(abs(unlist(var_map_diff)))
scale <- function(x) x / zmax
unscale <- function(x) x * zmax
subtitle <- rep("", length(combi_names))
lvalues <- lapply(var_map_diff, scale)
zlim <- c(-1, 1)
breaks <- seq(zlim[1], zlim[2], length.out = length(colors) + 1)
labels <- sprintf("%.2f", unscale(breaks))
multimap(
  lon, lat, lvalues,
  maintitle, subtitle,
  zlim, colors,
  breaks, labels,
  outersplit = rbind(c(0, 1, 3/17, 1), c(0, 1, 0, 3/17)),
  innersplit = c(2,2),
  mai = c(0.6732, 0.5412, 0.5412, 0.2772),
  legend_par = list(mar=c(8, 1, 1, 1), plt = c(0.1, 0.9, 0.9, 1))
)
dev.off()
