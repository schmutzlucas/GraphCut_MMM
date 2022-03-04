library(fields)
library(ncdf4)
library(maps)
library(magrittr)

source("functions/multimap.R")


args = commandArgs(trailingOnly=TRUE)
stopifnot(length(args) == 1)
var = args[1]
print(var)
units <- switch(var, "pr" = "mm/s", "tas" = "°C")
al_trans <- function(x) atan(log(fartorr(x)))/(pi/2)
al_inv <- function(x) rrtofar(exp(tan((pi/2)*x)))
rrtofar <- function(x) 1 - 1/x
fartorr <- function(x) 1 /(1 - x) 

atan_trans <- function(x, zmax = 1) atan(x)/ (pi/2) / zmax
atan_inv <- function(x, zmax = 1) tan((pi/2) * x * zmax)

pal_redwhiteblue <- c('#67001f','#b2182b','#d6604d','#f4a582','#fddbc7','#f7f7f7','#d1e5f0','#92c5de','#4393c3','#2166ac','#053061')
pal_violetwhitegreen <- c('#8e0152','#c51b7d','#de77ae','#f1b6da','#fde0ef','#f7f7f7','#e6f5d0','#b8e186','#7fbc41','#4d9221','#276419')
pal_whitered <- c('#fff5eb','#fee6ce','#fdd0a2','#fdae6b','#fd8d3c','#f16913','#d94801','#a63603','#7f2704')
pal_whiteredblack <- c('#fff5f0','#fee0d2','#fcbba1','#fc9272','#fb6a4a','#ef3b2c','#cb181d','#a50f15','#000000')
pal_whitepurple <- c('#fcfbfd','#efedf5','#dadaeb','#bcbddc','#9e9ac8','#807dba','#6a51a3','#54278f','#3f007d')

mmcombi <- readRDS(paste0("output/", var, "_mmcombi_perfectmodel.rds"))
mae_bias <- readRDS(paste0("output/", var, "_mae_bias_perfectmodel.rds"))
mae_bias_map <- readRDS(paste0("output/", var, "_mae_bias_map_perfectmodel.rds"))
bias_map <- readRDS(paste0("output/", var, "_bias_map_perfectmodel.rds"))
mae_gradient <- readRDS(paste0("output/", var, "_mae_gradient_perfectmodel.rds"))
mae_gradient_map  <- readRDS(paste0("output/", var, "_mae_gradient_map_perfectmodel.rds"))

model_names <- colnames(mae_bias)
icombi_ref <- which(names(mmcombi[[1]]) == "mmm")
col_combi <- c("red", "green", "blue", "purple", "grey48", "orange", "yellow")
pch_combi <- c(21, 22, 23, 8, 24, 25, 25)

nc <- nc_open("netcdf/tas_era5.nc")
lon <- ncvar_get(nc, "lon")
lat <- ncvar_get(nc, "lat")
nc_close(nc)


##### MAE Bias #####


opt <- options("scipen" = 3)
<<<<<<< HEAD
pdf(paste0("figures/", var, "_bias_perfectmodel.pdf"), width = 9)
par(mar = c(6.1, 6.1, 4.1, 8.1), mgp = c(5, 1, 0))
=======
pdf(paste0("figures/", var, "_bias_perfectmodel.pdf"))
par(mar = c(6.1, 6.1, 4.1, 2.1), mgp = c(5, 1, 0))
>>>>>>> dfee4809dbf5fef46729f8272b5ecbd0c5d1fb53
xlim <- c(1, length(model_names) - 1 + nrow(mae_bias)/2 + .5)
plot(
  NA,  
  xlim = xlim,
  main = "",
  type = "n",
  ylim = range(mae_bias),
  xlab = "",
  ylab = "MAE_b",
  col = "red",
  xaxt = "n",
  yaxt = 'n'
)
lapply(
  seq.int(length(model_names) - 1),
  function(v) abline(v = v, col = "lightgrey", lwd = 0.4)
)

axis(
  1,
  at = seq_along(head(model_names, -1)),
  labels = head(model_names, -1),
  las = 3,
  cex.axis = 0.75,
  srt = 45
)
axis(2, las = 1)

imodels <- (model_names != "summary")

for (icombi in seq.int(nrow(mae_bias))) {
  # lines(mae_bias[icombi, imodels], col = col_combi[icombi], lty = 1)
  points(
    mae_bias[icombi, imodels],
    col = col_combi[icombi],
    bg = col_combi[icombi],
    pch = pch_combi[icombi]
  )
  boxplot(
    mae_bias[icombi, imodels],
    at = length(model_names) + (icombi - 1) / 2,
    col = col_combi[icombi],
    ouline = FALSE,
    add = TRUE,
    axes = FALSE,
    boxwex = 0.4,
    range = 0
  )
}

axis(
  1,
  at = c(length(model_names), length(model_names) + nrow(mae_bias)/2),
  labels = rep("", 2),
  line = 0.5,
  las = 3,
  cex.axis = 0.75,
  srt = 45
)
mtext(
  "summary",
  side = 1,
  line = 1.5, 
  adj = 0,
  at = length(model_names)
)

<<<<<<< HEAD
legend(xlim[2] + 2, max(range(mae_bias)),  xpd = TRUE, 
       legend = rownames(mae_gradient), bty = "n",
       col = col_combi, pt.bg = col_combi, pch = pch_combi, cex=0.8, bg = "white") 
dev.off()

# MAE Gradient
pdf(paste0("figures/", var, "_gradient_perfectmodel.pdf"), width = 9)
par(mar = c(6.1, 6.1, 4.1, 8.1), mgp = c(5, 1, 0))
=======
legend("topright", legend = rownames(mae_gradient),
       col = col_combi, pt.bg = col_combi, pch = pch_combi, cex=0.8, bg = "white")

dev.off()

# MAE Gradient
pdf(paste0("figures/", var, "_gradient_perfectmodel.pdf"))
par(mar = c(6.1, 6.1, 4.1, 2.1), mgp = c(5, 1, 0))
>>>>>>> dfee4809dbf5fef46729f8272b5ecbd0c5d1fb53

plot(NA,  
     xlim=xlim,
     main="",
     type="n",
     ylim=range(mae_gradient),
     xlab ="",
     ylab="MAE_g",
     col="red",
     xaxt="n",
     yaxt='n')
lapply(
  seq.int(length(model_names) - 1),
  function(v) abline(v = v, col = "lightgrey", lwd = 0.4)
)

axis(
  1,
  at = seq_along(head(model_names, -1)),
  labels = head(model_names, -1),
  las = 3,
  cex.axis = 0.75,
  srt = 45
)
axis(2, las = 1)

for(icombi in seq.int(nrow(mae_gradient))){
  points(
    mae_gradient[icombi, imodels],
    col = col_combi[icombi],
    bg = col_combi[icombi],
    pch = pch_combi[icombi]
  )
  boxplot(
    mae_gradient[icombi, imodels],
    at = length(model_names) + (icombi - 1) / 2,
    col = col_combi[icombi],
    ouline = FALSE,
    add = TRUE,
    axes = FALSE,
    boxwex = 0.4,
    range = 0
  )
}

axis(
  1,
  at = c(length(model_names), length(model_names) + nrow(mae_bias)/2),
  labels = rep("", 2),
  line = 0.5,
  las = 3,
  cex.axis = 0.75,
  srt = 45
)
mtext(
  "summary",
  side = 1,
  line = 1.5, 
  adj = 0,
  at = length(model_names)
)

legend(xlim[2] + 2, max(range(mae_gradient)),  xpd = TRUE, 
       legend = rownames(mae_gradient), bty = "n",
       col = col_combi, pt.bg = col_combi, pch = pch_combi, cex=0.8, bg = "white")

dev.off()

options(opt)

##### Maps #####
# Labels
par(mfrow = c(3, 2))
mai = par("mai")
dev.off()

igc <- sapply(c("min_bias", "gc_present", "gc_hybrid", "gc_future"), function(x) which(x == names(mmcombi[[1]])))
imaps <- sapply(c("min_bias", "mmm", "gc_hybrid", "om_present", "gc_future", "om_future"), function(x) which(x == names(mmcombi[[1]])))
for(imodel in seq_along(bias_map)){
pdf(paste0("figures/", var, "_label_map_", model_names[imodel], "_perfectmodel.pdf"))
  mae <- mae_bias[igc, imodel]
  maintitle <- sprintf(
    "(%s) %s",
    letters[seq.int(length(igc))],
    names(mae)
  )
  subtitle <- rep("", length(igc))
  lvalues <- lapply(
    mmcombi[[imodel]][igc],
    function(x){
      x <- x[["Label attribution"]] + 1
      x[x >= imodel] <- x[x >= imodel] + 1
      return(x)
    }
  )
  colors <- rainbow(length(mmcombi)) 
  ncolors <- length(colors)
  breaks <- seq.int(ncolors + 1) - 0.5
  zlim <- range(breaks)
  labels <- names(mmcombi)
  multimap(
    lon, lat, lvalues,
    maintitle, subtitle,
    zlim, colors,
    breaks, labels,
    outersplit = rbind(c(0, 1, 0.33, 1), c(0, 1, 0, 0.33)), 
    innersplit = c(2,2),
    mai = c(0.6732, 0.5412, 0.5412, 0.2772),
    legend_par = list(mar=c(8, 1, 1, 1), plt = c(0.1, 0.9, 0.9, 1))
  )
  dev.off()
}

# Label frequency
freq <- array(0, dim = c(length(model_names) - 1, length(igc), length(model_names) - 1))
dimnames(freq) <- list("model" = head(model_names, -1), "icombi" = rownames(mae_bias)[igc], "reference" = head(model_names, -1))
for(imodel in seq.int(length(bias_map))){
  pdf(paste0("figures/", var, "_label_frenquency_", model_names[imodel], "_perfectmodel.pdf"))
  par(mfrow = c(2, 2), mar = c(8, 3, 1, 1))
  count <- list()
  for(icombi in seq_along(igc)){
    labels <- mmcombi[[imodel]][[igc[icombi]]][["Label attribution"]] + 1
    labels[labels >= imodel] <- labels[labels >= imodel] + 1
    count <- table(model_names[labels])
    freq[names(count), igc[icombi], imodel] <- count / sum(count)
  }
  for(icombi in seq_along(igc)){
    barplot(
      freq[-imodel, igc[icombi], imodel],
      main = paste0("(", letters[icombi], ") ", rownames(mae_bias)[igc[icombi]]),
      xlab="", ylim= c(0, max(freq[, , imodel])),
      las=2
    )
  }
  dev.off()
}

pdf(paste0("figures/", var, "_label_frenquency_summary_perfectmodel.pdf"))
freq_summary <- apply(freq, 1:2, mean)
par(mfrow = c(2, 2), mar = c(8, 3, 1, 1))
for(icombi in seq_along(igc)){
    barplot(
      freq_summary[, igc[icombi]],
      main = paste0("(", letters[icombi], ") ", rownames(mae_bias)[igc[icombi]]),
      xlab="", ylim=c(0, max(freq_summary)),
      las=2
    )
}
dev.off()

# Bias
colors <- if(var == "tas") rev(pal_redwhiteblue) else pal_redwhiteblue

for(imodel in seq_along(bias_map)){
  pdf(paste0("figures/", var, "_bias_map_", model_names[imodel], "_perfectmodel.pdf"), height = 8)
  mae <- mae_bias[imaps, imodel]
  maintitle <- sprintf(
    "(%s) %s",
    letters[seq.int(length(imaps))],
    names(mae)
  )
  subtitle <- sprintf("MAE_b = %.2f ", mae)
  lvalues <- lapply(bias_map[[imodel]][imaps], atan_trans)
  zmax <- unlist(lvalues) %>% abs %>% max
  lvalues <- lapply(lvalues, "/", y = zmax)
  zlim <- c(-1, 1)
  breaks <- seq(zlim[1], zlim[2], length.out = length(colors) + 1)
  labels <- sprintf("%.2f", atan_inv(breaks, zmax))
  multimap(
    lon, lat, lvalues,
    maintitle, subtitle,
    zlim, colors,
    breaks, labels
  )
  dev.off()
}

# MAE bias
colors <- pal_whitered
ncolors <- length(colors)
for(imodel in seq_along(mae_bias_map)){
  pdf(paste0("figures/", var, "_mae_bias_map_", model_names[imodel], "_perfectmodel.pdf"), height = 8)
  mae <- mae_bias[imaps, imodel]
  maintitle <- sprintf(
    "(%s) %s",
    letters[seq.int(length(imaps))],
    names(mae)
  )
  subtitle <- sprintf("MAE_b = %.2f ", mae)
  lvalues <- lapply(mae_bias_map[[imodel]][imaps], atan_trans)
  zmax <- unlist(lvalues) %>% abs %>% max
  lvalues <- lapply(lvalues, "/", y = zmax)
  zlim <- c(0, 1)
  breaks <- seq(zlim[1], zlim[2], length.out = length(colors) + 1)
  labels <- sprintf("%.2f", atan_inv(breaks, zmax))
  multimap(
    lon, lat, lvalues,
    maintitle, subtitle,
    zlim, colors,
    breaks, labels
  )
  dev.off()
}

# MAE bias Skill Score
colors <- rev(pal_violetwhitegreen)
ncolors <- length(colors)
for(imodel in seq_along(mae_bias_map)){
  pdf(paste0("figures/", var, "_mae_bias_diff_map_", model_names[imodel], "_perfectmodel.pdf"), height = 8)
  mae <- mae_bias[imaps, imodel]
  maintitle <- sprintf(
    "(%s) %s",
    letters[seq.int(length(imaps))],
    names(mae)
  )
  subtitle <- sprintf("MAE_b = %.2f ", mae)
  lvalues <- lapply(
    mae_bias_map[[imodel]][imaps], 
    function(x) atan_trans(x - mae_bias_map[[imodel]][[icombi_ref]])
  )
  zmax <- unlist(lvalues) %>% abs %>% max
  lvalues <- lapply(lvalues, "/", y = zmax)
  zlim <- c(-1, 1)
  breaks <- seq(zlim[1], zlim[2], length.out = length(colors) + 1)
  labels <- sprintf("%.2f", atan_inv(breaks, zmax))
  multimap(
    lon, lat, lvalues,
    maintitle, subtitle,
    zlim, colors,
    breaks, labels
  )
  dev.off()
}

# MAE gradient
colors <- pal_whiteredblack
ncolors <- length(colors)
for(imodel in seq_along(mae_gradient_map)){
  pdf(paste0("figures/", var, "_mae_gradient_map_", model_names[imodel], "_perfectmodel.pdf"), height = 8)
  mae <- mae_gradient[imaps, imodel]
  maintitle <- sprintf(
    "(%s) %s",
    letters[seq.int(length(imaps))],
    names(mae)
  )
  subtitle <- sprintf("MAE_g = %.2f ", mae)
  lvalues <- lapply(mae_gradient_map[[imodel]][imaps], atan_trans)
  zmax <- unlist(lvalues) %>% abs %>% max
  lvalues <- lapply(lvalues, "/", y = zmax)
  zlim <- c(0, 1)
  breaks <- seq(zlim[1], zlim[2], length.out = length(colors) + 1)
  labels <- sprintf("%.2f", atan_inv(breaks, zmax))
  multimap(
    lon, lat, lvalues,
    maintitle, subtitle,
    zlim, colors,
    breaks, labels
  )
  dev.off()
}

# MAE gradient Skill Score
colors <- rev(pal_violetwhitegreen)
ncolors <- length(colors)
for(imodel in seq_along(mae_gradient_map)){
  pdf(paste0("figures/", var, "_mae_gradient_diff_map_", model_names[imodel], "_perfectmodel.pdf"), height = 8)
  mae <- mae_gradient[imaps, imodel]
  maintitle <- sprintf(
    "(%s) %s",
    letters[seq.int(length(imaps))],
    names(mae)
  )
  subtitle <- sprintf("MAE_g= %.2f ", mae)
  lvalues <- lapply(
    mae_gradient_map[[imodel]][imaps], 
    function(x) atan_trans(x - mae_gradient_map[[imodel]][[icombi_ref]])
  )
  zmax <- unlist(lvalues) %>% abs %>% max
  lvalues <- lapply(lvalues, "/", y = zmax)
  zlim <- c(-1, 1)
  breaks <- seq(zlim[1], zlim[2], length.out = length(colors) + 1)
  labels <- sprintf("%.2f", atan_inv(breaks, zmax))
  multimap(
    lon, lat, lvalues,
    maintitle, subtitle,
    zlim, colors,
    breaks, labels
  )
  dev.off()
}
