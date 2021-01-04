library(fields)
library(ncdf4)
library(maps)

var <- "tas"

mae_bias <- readRDS(paste0("output/", var, "_mae_bias.rds"))
mae_bias_map <- readRDS(paste0("output/", var, "_mae_bias_map.rds"))
mse_gradient <- readRDS(paste0("output/", var, "_mse_gradient.rds"))
mse_gradient_map  <- readRDS(paste0("output/", var, "_mse_gradient_map.rds"))

model_names <- colnames(mae_bias)

nc <- nc_open("input/tas_30Ayr_ERA5_197901_200812.nc")
lon <- ncvar_get(nc, "lon")
lat <- ncvar_get(nc, "lat")
nc_close(nc)


##### MAE Bias #####

op <- par(no.readonly = T)

pdf(paste0("output/", var, "_bias.pdf"))

plot(NA,  
     xlim=range(seq_along(model_names)),
     main="",
     type="n",
     ylim=range(mae_bias),
     xlab ="",
     ylab="MAE_b (°C)",
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

col_combi <- c("red", "darkolivegreen1", "blue", "purple", "grey", "black")

for(icombi in seq.int(nrow(mae_bias))){
  points(mae_bias[icombi, ], col = col_combi[icombi])
  lines(mae_bias[icombi, ], col = col_combi[icombi], lty = 2)
}

legend("topright", legend = rownames(mae_bias),
       col = col_combi, pch= 1, cex=0.8)

dev.off()
par(op)

# MSE Gradient
pdf(paste0("output/", var, "_gradient.pdf"))

plot(NA,  
     xlim=range(seq_along(model_names)),
     main="",
     type="n",
     ylim=range(mse_gradient),
     xlab ="",
     ylab="MAE_g (°C)",
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
       col = col_combi, pch= 1, cex=0.8)

dev.off()
par(op)




##### Maps #####
colors <- c('#ffffcc','#ffeda0','#fed976','#feb24c','#fd8d3c','#fc4e2a','#e31a1c','#bd0026','#800026')
colors <- rev(viridis::viridis(11))
ncolors <- length(colors)
# MAE bias

for(imodel in seq_along(mae_bias_map)){
  pdf(paste0("output/", var, "_bias_map_", model_names[imodel], ".pdf"))
  # zlim <- max(abs(range(unlist(mae_bias_map[[imodel]])))) * c(-1, 1)
  zlim <- c(0, max(unlist(mae_bias_map[[imodel]])))
  par(mfrow = c(3, 2))
  for(icombi in seq_along(mae_bias_map[[imodel]])){
    image.plot(lon, lat, mae_bias_map[[imodel]][[icombi]],
               xlab = '',
               ylab = '',
               main = rownames(mae_bias)[icombi],
               sub = sprintf("MAE_b = %.2f ", mae_bias[icombi, imodel]),
               zlim = zlim,
               col=colors
    )
               # axis.args=list( at=1:length(colorBreaks), labels=paste0("> ",colorBreaks)))
               # lab.breaks = names(colorBreaks))
    #mtext(, side = 4)
    map("world2",add=T)
  }
  dev.off()
}

# MSE gradient
for(imodel in seq_along(mse_gradient_map)){
  pdf(paste0("output/", var, "_gradient_map_", model_names[imodel], ".pdf"))
  zlim <- c(0, max(unlist(mse_gradient_map[[imodel]])))
  print(zlim)
  breaks <- c(0, zlim[2]/2^((ncolors - 1):0))
  par(mfrow = c(3, 2))
  for(icombi in seq_along(mse_gradient_map[[imodel]])){
    z <- mse_gradient_map[[imodel]][[icombi]]
    zdiscretized <- matrix(as.integer(cut(z, breaks = breaks)), nrow = nrow(z), ncol = ncol(z))
    image.plot(lon, lat, zdiscretized,
               xlab = '',
               ylab = '',
               main = rownames(mse_gradient)[icombi],
               sub = sprintf("MAE_g = %.2f ", mse_gradient[icombi, imodel]),
               zlim = c(0, ncolors),
               breaks = 0:ncolors,
               col = colors,
               axis.args = list(at=seq_along(breaks) - 1, labels = sprintf("%.3f", breaks))
    )
               # lab.breaks = names(colorBreaks))
    #mtext(, side = 4)
    map("world2",add=T)
  }
  dev.off()
}
