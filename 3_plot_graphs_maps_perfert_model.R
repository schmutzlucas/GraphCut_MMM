library(fields)
library(ncdf4)
library(maps)


args = commandArgs(trailingOnly=TRUE)
stopifnot(length(args) == 1)
var = args[1]
print(var)
units <- switch(var, "pr" = "mm/s", "tas" = "°C")
al_trans <- function(x) atan(log(fartorr(x)))/(pi/2)
al_inv <- function(x) rrtofar(exp(tan((pi/2)*x)))
rrtofar <- function(x) 1 - 1/x
fartorr <- function(x) 1 /(1 - x) 



mmcombi <- readRDS(paste0("output/", var, "_mmcombi_perfectmodel.rds"))
mae_bias <- readRDS(paste0("output/", var, "_mae_bias_perfectmodel.rds"))
mae_bias_map <- readRDS(paste0("output/", var, "_mae_bias_map_perfectmodel.rds"))
bias_map <- readRDS(paste0("output/", var, "_bias_map_perfectmodel.rds"))
mse_gradient <- readRDS(paste0("output/", var, "_mse_gradient_perfectmodel.rds"))
mse_gradient_map  <- readRDS(paste0("output/", var, "_mse_gradient_map_perfectmodel.rds"))

model_names <- colnames(mae_bias)
icombi_ref <- which(names(mmcombi[[1]]) == "mmm")
col_combi <- c("red", "green", "blue", "purple", "grey48", "black")

nc <- nc_open("netcdf/tas_era5.nc")
lon <- ncvar_get(nc, "lon")
lat <- ncvar_get(nc, "lat")
nc_close(nc)


##### MAE Bias #####


opt <- options("scipen" = 3)
pdf(paste0("figures/", var, "_bias_perfectmodel.pdf"))
par(mar = c(6.1, 6.1, 4.1, 2.1), mgp = c(5, 1, 0))
plot(NA,  
     xlim=range(seq_along(model_names)),
     main="",
     type="n",
     ylim=range(mae_bias),
     xlab ="",
     ylab="MAE_b",
     col="red",
     xaxt="n",
     yaxt='n')

axis(1,
     at=seq_along(model_names),
     labels = model_names,
     las=3,
     cex.axis=0.75,
     srt=45)
axis(2,
     las=1)


for(icombi in seq.int(nrow(mae_bias))){
  points(mae_bias[icombi, ], col = col_combi[icombi])
  lines(mae_bias[icombi, ], col = col_combi[icombi], lty = 2)
}

legend("topright", legend = rownames(mae_bias),
       col = col_combi, pch= 1, cex=0.8, bg = "white"
)

dev.off()

# MSE Gradient
pdf(paste0("figures/", var, "_gradient_perfectmodel.pdf"))
par(mar = c(6.1, 6.1, 4.1, 2.1), mgp = c(5, 1, 0))

plot(NA,  
     xlim=range(seq_along(model_names)),
     main="",
     type="n",
     ylim=range(mse_gradient),
     xlab ="",
     ylab="MSE_g",
     col="red",
     xaxt="n",
     yaxt='n')

axis(1,
     at=seq_along(model_names),
     labels = model_names,
     las=3,
     cex.axis=0.75,
     srt=45)
axis(2, las=1)

for(icombi in seq.int(nrow(mse_gradient))){
  points(mse_gradient[icombi, ], col = col_combi[icombi])
  lines(mse_gradient[icombi, ], col = col_combi[icombi], lty = 2)
}

legend("topright", legend = rownames(mse_gradient),
       col = col_combi, pch= 1, cex=0.8, bg = "white")

dev.off()

options(opt)

##### Maps #####
# Labels
colors <- rainbow(length(mmcombi)) 
ncolors <- length(colors)
breaks <- seq.int(0, ncolors)
at <- breaks[-1] - 0.5
igc <- which(names(mmcombi[[1]]) %in% c("gc_present", "gc_hybrid", "gc_future", "min_bias"))
for(imodel in seq_along(bias_map)){
  pdf(paste0("figures/", var, "_label_map_", model_names[imodel], "_perfectmodel.pdf"))
  set.panel() # reset plotting device
  iouter <- split.screen(rbind(c(0, 1, 0.2, 1), c(0, 1, 0, 0.25)))
  iinner <- split.screen(c(2,2), screen= iouter[1])
  for(icombi in seq_along(igc)){
    screen(iinner[icombi])
    image(
      lon, lat, mmcombi[[imodel]][[igc[icombi]]][["Label attribution"]],
      xlab = '',
      ylab = '',
      main = rownames(mae_bias)[igc[icombi]],
      col=colors[-imodel]
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
  axis(1, at = at, labels = names(mmcombi), las = 3)
  mtext("labels", side = 3)
  close.screen(all=TRUE)
  par(oldpar)
  dev.off()
}

# Label frequency
freq <- array(0, dim = c(length(model_names) - 1, length(igc), length(model_names) - 1))
dimnames(freq) <- list("model" = head(model_names, -1), "icombi" = rownames(mae_bias)[igc], "reference" = head(model_names, -1))
for(imodel in seq.int(length(bias_map))){
  pdf(paste0("figures/", var, "_label_frenquency_", model_names[imodel], "_perfectmodel.pdf"))
  par(mfrow = c(2, 2), mar = c(8, 3, 1, 1))
  for(icombi in seq_along(igc)){
    labels <- mmcombi[[imodel]][[igc[icombi]]][["Label attribution"]] + 1
    labels[labels >= imodel] <- labels[labels >= imodel] + 1
    count <- table(model_names[labels])
    freq[names(count), igc[icombi], imodel] <- count / sum(count)
    barplot(freq[-imodel, igc[icombi], imodel], main = rownames(mae_bias)[igc[icombi]], xlab="", ylim= c(0, max(freq[, , imodel])), las=2)
  }
  dev.off()
}
pdf(paste0("figures/", var, "_label_frenquency_summary_perfectmodel.pdf"))
freq_summary <- apply(freq, 1:2, mean)
par(mfrow = c(2, 2), mar = c(8, 3, 1, 1))
for(icombi in seq_along(igc)){
    barplot(freq_summary[, igc[icombi]],  main=rownames(mae_bias)[igc[icombi]], xlab="", ylim=c(0, max(freq_summary)), las=2)
}
dev.off()

# MAE bias
colors <- c('#ffffcc','#ffeda0','#fed976','#feb24c','#fd8d3c','#fc4e2a','#e31a1c','#bd0026','#800026')
colors <- c('#67001f','#b2182b','#d6604d','#f4a582','#fddbc7','#f7f7f7','#d1e5f0','#92c5de','#4393c3','#2166ac','#053061')
ncolors <- length(colors)

for(imodel in seq_along(bias_map)){
  pdf(paste0("figures/", var, "_bias_map_", model_names[imodel], "_perfectmodel.pdf"))
  zmax <- max(abs(atan(unlist(bias_map[[imodel]]))/(pi/2)))
  zlim <- c(-1, 1)
  breaks <- seq(-1, 1, length.out = ncolors + 1)
  par(mfrow = c(3, 2))
  for(icombi in seq_along(bias_map[[imodel]])){
    bias <- atan(bias_map[[imodel]][[icombi]]) / (pi/2) / zmax
    image.plot(lon, lat, bias,
               xlab = '',
               ylab = '',
               main = rownames(mae_bias)[icombi],
               sub = sprintf("MAE_b = %.2f ", mae_bias[icombi, imodel]),
               zlim = zlim,
               col=colors,
               breaks = breaks,
               axis.args = list(at=breaks, labels = sprintf("%.2f", tan(breaks * zmax * (pi/2))), cex.axis = 0.5)
    )
    map("world2",add=T)
  }
  dev.off()
}

colors <- rev(viridis::viridis(11))
colors <- c('#fff5eb','#fee6ce','#fdd0a2','#fdae6b','#fd8d3c','#f16913','#d94801','#a63603','#7f2704')
ncolors <- length(colors)
for(imodel in seq_along(mae_bias_map)){
  pdf(paste0("figures/", var, "_mae_bias_map_", model_names[imodel], "_perfectmodel.pdf"))
  zmax <- max(abs(atan(unlist(mae_bias_map[[imodel]]))/(pi/2)))
  zlim <- c(0, 1)
  breaks <- seq(0, 1, length.out = ncolors + 1)
  par(mfrow = c(3, 2))
  for(icombi in seq_along(mae_bias_map[[imodel]])){
    mae <- atan(mae_bias_map[[imodel]][[icombi]]) / (pi/2) / zmax
    image.plot(lon, lat, mae,
               xlab = '',
               ylab = '',
               main = rownames(mae_bias)[icombi],
               sub = sprintf("MAE_b = %.2f ", mae_bias[icombi, imodel]),
               zlim = zlim,
               col=colors,
               breaks = breaks,
               axis.args = list(at=breaks, labels = sprintf("%.2f", tan(breaks * zmax * (pi/2))), cex.axis = 0.5)
    )
    map("world2",add=T)
  }
  dev.off()
}

# MAE bias Skill Score
colors <- rev(viridis::viridis(11))
colors <- c('#9e0142','#d53e4f','#f46d43','#fdae61','#fee08b','#ffffbf','#e6f598','#abdda4','#66c2a5','#3288bd','#5e4fa2')
ncolors <- length(colors)
for(imodel in seq_along(mae_bias_map)){
  pdf(paste0("figures/", var, "_mae_bias_skillscore_map_", model_names[imodel], "_perfectmodel.pdf"))
  # zlim <- max(abs(range(unlist(mae_bias_map[[imodel]])))) * c(-1, 1)
  zlim <- c(-1, 1)
  breaks <- seq(-1, 1, length.out = ncolors + 1)
  par(mfrow = c(3, 2))
  for(icombi in seq_along(mae_bias_map[[imodel]])){
    skill <-  al_trans(1 - mae_bias_map[[imodel]][[icombi]] / mae_bias_map[[imodel]][[icombi_ref]])
    image.plot(lon, lat, skill,
               xlab = '',
               ylab = '',
               main = rownames(mae_bias)[icombi],
               sub = sprintf("MAE_b = %.2f ", mae_bias[icombi, imodel]),
               zlim = zlim,
               col = colors,
               breaks = breaks,
               axis.args = list(at=breaks, labels = sprintf("%.2f", al_inv(breaks)), cex.axis = 0.5)
    )
    map("world2",add=T)
  }
  dev.off()
}

# MSE gradient
colors <- rev(viridis::viridis(11))
colors <- c('#fcfbfd','#efedf5','#dadaeb','#bcbddc','#9e9ac8','#807dba','#6a51a3','#54278f','#3f007d')
ncolors <- length(colors)
for(imodel in seq_along(mse_gradient_map)){
  pdf(paste0("figures/", var, "_mse_gradient_map_", model_names[imodel], "_perfectmodel.pdf"))
  zmax <- max(abs(atan(unlist(mse_gradient_map[[imodel]]))/(pi/2)))
  zlim <- c(0, 1)
  breaks <- seq(0, 1, length.out = ncolors + 1)
  par(mfrow = c(3, 2))
  for(icombi in seq_along(mse_gradient_map[[imodel]])){
    mse <- atan(mse_gradient_map[[imodel]][[icombi]]) / (pi/2) / zmax
    image.plot(lon, lat, mse,
               xlab = '',
               ylab = '',
               main = rownames(mse_gradient)[icombi],
               sub = sprintf("MSE_g = %.2f ", mse_gradient[icombi, imodel]),
               zlim = zlim,
               breaks = breaks,
               col = colors,
               axis.args = list(at=breaks, labels = sprintf("%.2f", tan(breaks * zmax * (pi/2))), cex.axis = 0.5)
    )
    map("world2",add=T)
  }
  dev.off()
}

# MSE gradient Skill Score
colors <- rev(viridis::viridis(11))
colors <- c('#9e0142','#d53e4f','#f46d43','#fdae61','#fee08b','#ffffbf','#e6f598','#abdda4','#66c2a5','#3288bd','#5e4fa2')
ncolors <- length(colors)
for(imodel in seq_along(mse_gradient_map)){
  pdf(paste0("figures/", var, "_mse_gradient_skillscore_map_", model_names[imodel], "_perfectmodel.pdf"))
  zlim <- c(-1, 1)
  breaks <- seq(-1, 1, length.out = ncolors + 1)
  par(mfrow = c(3, 2))
  for(icombi in seq_along(mse_gradient_map[[imodel]])){
    skill <-  al_trans(1 - mse_gradient_map[[imodel]][[icombi]] / mse_gradient_map[[imodel]][[icombi_ref]])
    image.plot(lon, lat, skill,
               xlab = '',
               ylab = '',
               main = rownames(mse_gradient)[icombi],
               sub = sprintf("MSE_g = %.2f ", mse_gradient[icombi, imodel]),
               zlim = zlim,
               col = colors,
               breaks = breaks,
               axis.args = list(at=breaks, labels = sprintf("%.2f", al_inv(breaks)), cex.axis = 0.5)
    )
    map("world2",add=T)
  }
  dev.off()
}
