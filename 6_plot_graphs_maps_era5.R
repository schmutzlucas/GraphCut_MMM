library(fields)
library(ncdf4)
library(maps)

args = commandArgs(trailingOnly=TRUE)
stopifnot(length(args) == 1)
var = args[1]
print(var)

imodels <- c(
  1, 2, 3, 4, 7,
  8, 10, 12, 13, 14, 
  17, 20, 22, 24, 26,
  30, 31, 33, 34, 37
)
model_names <- readRDS("input/models_list.rds")[imodels, 2]

al_trans <- function(x) atan(log(fartorr(x)))/(pi/2)
al_inv <- function(x) rrtofar(exp(tan((pi/2)*x)))
rrtofar <- function(x) 1 - 1/x
fartorr <- function(x) 1 /(1 - x) 

mmcombi <- readRDS(paste0("output/", var, "_mmcombi_era5.rds"))
mae_bias <- readRDS(paste0("output/", var, "_mae_bias_era5.rds"))
mae_bias_map <- readRDS(paste0("output/", var, "_mae_bias_map_era5.rds"))
bias_map <- readRDS(paste0("output/", var, "_bias_map_era5.rds"))
mse_gradient <- readRDS(paste0("output/", var, "_mse_gradient_era5.rds"))
mse_gradient_map  <- readRDS(paste0("output/", var, "_mse_gradient_map_era5.rds"))


nc <- nc_open("netcdf/tas_era5.nc")
lon <- ncvar_get(nc, "lon")
lat <- ncvar_get(nc, "lat")
nc_close(nc)

combi_names <- names(mae_bias)
icombi_ref <- which(combi_names == "mmm")
col_combi <- c("red", "darkolivegreen1", "blue", "purple", "grey", "black")

##### MAE Bias #####

op <- par(no.readonly = T)

pdf(paste0("figures/", var, "_bias_era5.pdf"))

barplot(mae_bias,  
     names.arg = combi_names,
     main="",
     ylab="MAE_b",
     col=col_combi
)

dev.off()
par(op)

# MSE Gradient
pdf(paste0("figures/", var, "_gradient_era5.pdf"))

barplot(mse_gradient,  
     names.arg = combi_names,
     main="",
     ylab="MSE_g",
     col=col_combi
)

dev.off()
par(op)




##### Maps #####
##### Maps #####
# Labels
colors <- rainbow(length(model_names)) 
ncolors <- length(colors)
breaks <- seq.int(0, ncolors)
at <- breaks[-1] - 0.5
igc <- which(names(mmcombi) %in% c("gc_present", "gc_hybrid", "gc_future", "min_bias"))
pdf(paste0("figures/", var, "_label_map_era5.pdf"))
set.panel() # reset plotting device
iouter <- split.screen(rbind(c(0, 1, 0.2, 1), c(0, 1, 0, 0.25)))
iinner <- split.screen(c(2,2), screen= iouter[1])
for(icombi in seq_along(igc)){
  screen(iinner[icombi])
  image(
    lon, lat, mmcombi[[igc[icombi]]][["Label attribution"]],
    xlab = '',
    ylab = '',
    main = names(mmcombi)[igc[icombi]],
    col= colors
  )
  map("world2",add=T)
}
screen(iouter[2])
oldpar <- par(mar=c(8, 1, 1, 1), plt = c(0.1, 0.9, 0.9, 1))
image(
  x = at,
  y = 1:2, 
  z = matrix(at, ncol = 1),
  col = colors,
  breaks = breaks,
  xaxt = "n", yaxt = "n", xlab = "", ylab = ""
)
axis(1, at = at, labels = model_names, las = 3)
mtext("labels", side = 3)
close.screen(all=TRUE)
par(oldpar)
dev.off()

# Label frequency
freq <- matrix(0, nrow = length(model_names), ncol = length(igc))
dimnames(freq) <- list("model" = model_names, "icombi" = rownames(mae_bias)[igc])
pdf(paste0("figures/", var, "_label_frenquency_era5.pdf"))
par(mfrow = c(2, 2), mar = c(8, 3, 1, 1))
for(icombi in seq_along(igc)){
  labels <- mmcombi[[igc[icombi]]][["Label attribution"]] + 1
  count <- table(model_names[labels])
  freq[names(count), igc[icombi]] <- count / sum(count)
  barplot(freq[, igc[icombi]], main = names(mae_bias)[igc[icombi]], xlab="", ylim= c(0, max(freq)), las=2)
}
dev.off()

colors <- c('#ffffcc','#ffeda0','#fed976','#feb24c','#fd8d3c','#fc4e2a','#e31a1c','#bd0026','#800026')
colors <- c('#67001f','#b2182b','#d6604d','#f4a582','#fddbc7','#f7f7f7','#d1e5f0','#92c5de','#4393c3','#2166ac','#053061')
ncolors <- length(colors)
# MAE bias

pdf(paste0("figures/", var, "_bias_map_era5.pdf"))
zmax <- max(abs(atan(unlist(bias_map))/(pi/2)))
zlim <- c(-1, 1)
breaks <- seq(-1, 1, length.out = ncolors + 1)
par(mfrow = c(3, 2))
for(icombi in seq_along(bias_map)){
  bias <- atan(bias_map[[icombi]]) / (pi/2) / zmax
  image.plot(lon, lat, bias,
             xlab = '',
             ylab = '',
             main = combi_names[icombi],
             sub = sprintf("MAE_b = %.2f ", mae_bias[icombi]),
             zlim = zlim,
             col=colors,
             breaks = breaks,
             axis.args = list(at=breaks, labels = sprintf("%.2f", tan(breaks * zmax * (pi/2))), cex.axis = 0.5)
  )
  map("world2",add=T)
}
dev.off()

colors <- c('#fff5eb','#fee6ce','#fdd0a2','#fdae6b','#fd8d3c','#f16913','#d94801','#a63603','#7f2704')
ncolors <- length(colors)
pdf(paste0("figures/", var, "_mae_bias_map_era5.pdf"))
# zlim <- max(abs(range(unlist(mae_bias_map[[imodel]])))) * c(-1, 1)
zmax <- max(abs(atan(unlist(mae_bias_map))/(pi/2)))
zlim <- c(0, 1)
breaks <- seq(0, 1, length.out = ncolors + 1)
par(mfrow = c(3, 2))
for(icombi in seq_along(mae_bias_map)){
  mae <- atan(mae_bias_map[[icombi]]) / (pi/2) / zmax
  image.plot(lon, lat, mae,
             xlab = '',
             ylab = '',
             main = combi_names[icombi],
             sub = sprintf("MAE_b = %.2f ", mae_bias[icombi]),
             zlim = zlim,
             col=colors,
             breaks = breaks,
             axis.args = list(at=breaks, labels = sprintf("%.2f", tan(breaks * zmax * (pi/2))), cex.axis = 0.5)
  )
  map("world2",add=T)
}
dev.off()

# MAE bias Skill Score
colors <- c('#9e0142','#d53e4f','#f46d43','#fdae61','#fee08b','#ffffbf','#e6f598','#abdda4','#66c2a5','#3288bd','#5e4fa2')
ncolors <- length(colors)
pdf(paste0("figures/", var, "_mae_bias_skillscore_map_era5.pdf"))
zlim <- c(-1, 1)
breaks <- seq(-1, 1, length.out = ncolors + 1)
par(mfrow = c(3, 2))
for(icombi in seq_along(mae_bias_map)){
  skill <-  al_trans(1 - mae_bias_map[[icombi]] / mae_bias_map[[icombi_ref]])
  image.plot(lon, lat, skill,
             xlab = '',
             ylab = '',
             main = combi_names[icombi],
             sub = sprintf("MAE_b = %.2f ", mae_bias[icombi]),
             zlim = zlim,
             col = colors,
             breaks = breaks,
             axis.args = list(at=breaks, labels = sprintf("%.2f", al_inv(breaks)), cex.axis = 0.5)
  )
  map("world2",add=T)
}
dev.off()

# MSE gradient
colors <- rev(viridis::viridis(11))
colors <- c('#fcfbfd','#efedf5','#dadaeb','#bcbddc','#9e9ac8','#807dba','#6a51a3','#54278f','#3f007d')
ncolors <- length(colors)
pdf(paste0("figures/", var, "_mse_gradient_map_era5.pdf"))
zmax <- max(abs(atan(unlist(mse_gradient_map))/(pi/2)))
zlim <- c(0, 1)
breaks <- seq(0, 1, length.out = ncolors + 1)
par(mfrow = c(3, 2))
for(icombi in seq_along(mse_gradient_map)){
  mse <- atan(mse_gradient_map[[icombi]]) / (pi/2) / zmax
  image.plot(lon, lat, mse,
             xlab = '',
             ylab = '',
             main = combi_names[icombi],
             sub = sprintf("MAE_g = %.2f ", mse_gradient[icombi]),
             zlim = zlim,
             breaks = breaks,
             col = colors,
             axis.args = list(at=breaks, labels = sprintf("%.2f", tan(breaks * zmax * (pi/2))), cex.axis = 0.5)
  )
  map("world2",add=T)
}
dev.off()

# MSE gradient Skill Score
colors <- c('#9e0142','#d53e4f','#f46d43','#fdae61','#fee08b','#ffffbf','#e6f598','#abdda4','#66c2a5','#3288bd','#5e4fa2')
ncolors <- length(colors)
pdf(paste0("figures/", var, "_mse_gradient_skillscore_map_era5.pdf"))
zlim <- c(-1, 1)
breaks <- seq(-1, 1, length.out = ncolors + 1)
par(mfrow = c(3, 2))
for(icombi in seq_along(mse_gradient_map)){
  skill <-  al_trans(1 - mse_gradient_map[[icombi]] / mse_gradient_map[[icombi_ref]])
  image.plot(lon, lat, skill,
             xlab = '',
             ylab = '',
             main = combi_names[icombi],
             sub = sprintf("MAE_b = %.2f ", mse_gradient[icombi]),
             zlim = zlim,
             col = colors,
             breaks = breaks,
             axis.args = list(at=breaks, labels = sprintf("%.2f", al_inv(breaks)), cex.axis = 0.5)
  )
  map("world2",add=T)
}
dev.off()
