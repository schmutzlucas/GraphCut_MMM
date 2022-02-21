library(fields)
library(ncdf4)
library(maps)
library(magrittr)

source("functions/multimap.R")

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

atan_trans <- function(x, zmax = 1) atan(x)/ (pi/2) / zmax
atan_inv <- function(x, zmax = 1) tan((pi/2) * x * zmax)

pal_redwhiteblue <- c('#67001f','#b2182b','#d6604d','#f4a582','#fddbc7','#f7f7f7','#d1e5f0','#92c5de','#4393c3','#2166ac','#053061')
pal_violetwhitegreen <- c('#8e0152','#c51b7d','#de77ae','#f1b6da','#fde0ef','#f7f7f7','#e6f5d0','#b8e186','#7fbc41','#4d9221','#276419')
pal_whitered <- c('#fff5eb','#fee6ce','#fdd0a2','#fdae6b','#fd8d3c','#f16913','#d94801','#a63603','#7f2704')
pal_whiteredblack <- c('#fff5f0','#fee0d2','#fcbba1','#fc9272','#fb6a4a','#ef3b2c','#cb181d','#a50f15','#000000')
pal_whitepurple <- c('#fcfbfd','#efedf5','#dadaeb','#bcbddc','#9e9ac8','#807dba','#6a51a3','#54278f','#3f007d')

mmcombi <- readRDS(paste0("output/", var, "_mmcombi_era5.rds"))
mae_bias <- readRDS(paste0("output/", var, "_mae_bias_era5.rds"))
mae_bias_map <- readRDS(paste0("output/", var, "_mae_bias_map_era5.rds"))
bias_map <- readRDS(paste0("output/", var, "_bias_map_era5.rds"))
mae_gradient <- readRDS(paste0("output/", var, "_mae_gradient_era5.rds"))
mae_gradient_map  <- readRDS(paste0("output/", var, "_mae_gradient_map_era5.rds"))


nc <- nc_open("netcdf/tas_era5.nc")
lon <- ncvar_get(nc, "lon")
lat <- ncvar_get(nc, "lat")
nc_close(nc)

combi_names <- names(mae_bias)
icombi_ref <- which(combi_names == "mmm")
col_combi <- c("red", "green", "blue", "purple", "grey48", "orange", "yellow")

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

# MAE Gradient
pdf(paste0("figures/", var, "_gradient_era5.pdf"))

barplot(mae_gradient,  
     names.arg = combi_names,
     main="",
     ylab="MAE_g",
     col=col_combi
)

dev.off()




##### Maps #####
# Labels
igc <- sapply(c("min_bias", "gc_present", "gc_hybrid", "gc_future"), function(x) which(x == names(mmcombi)))
imaps <- sapply(c("min_bias", "mmm", "gc_hybrid", "om_present", "gc_future", "om_future"), function(x) which(x == names(mmcombi)))
pdf(paste0("figures/", var, "_label_map_era5.pdf"))
mae <- mae_bias[igc]
maintitle <- sprintf(
  "(%s) %s",
    letters[seq.int(length(igc))],
    names(mae)
)
subtitle <- rep("", length(igc))
lvalues <- lapply(
  mmcombi[igc],
  function(x){
    x <- x[["Label attribution"]] + 1
    return(x)
  }
)
colors <- rainbow(length(model_names)) 
ncolors <- length(colors)
breaks <- seq.int(ncolors + 1) - 0.5
zlim <- range(breaks)
labels <- model_names
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


# Label frequency
freq <- matrix(0, nrow = length(model_names), ncol = length(igc))
dimnames(freq) <- list("model" = model_names, "icombi" = rownames(mae_bias)[igc])
pdf(paste0("figures/", var, "_label_frenquency_era5.pdf"))
par(mfrow = c(2, 2), mar = c(8, 3, 1, 1))
for(icombi in seq_along(igc)){
  labels <- mmcombi[[igc[icombi]]][["Label attribution"]] + 1
  count <- table(model_names[labels])
  freq[names(count), igc[icombi]] <- count / sum(count)
}
for(icombi in seq_along(igc)){
  barplot(
    freq[, igc[icombi]],
    main = paste0("(", letters[icombi], ") ", names(mae_bias)[igc[icombi]]),
    xlab="", ylim= c(0, max(freq)),
    las=2
  )
}
dev.off()

# MAE bias

pdf(paste0("figures/", var, "_bias_map_era5.pdf"), height = 8)
colors <- if(var == "tas") rev(pal_redwhiteblue) else pal_redwhiteblue
ncolors <- length(colors)
mae <- mae_bias[imaps]
maintitle <- sprintf(
  "(%s) %s",
  letters[seq.int(length(imaps))],
  names(mae)
)
subtitle <- sprintf("MAE_b = %.2f ", mae)
lvalues <- lapply(bias_map[imaps], atan_trans)
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

pdf(paste0("figures/", var, "_bias_hist_era5.pdf"), height = 8)
xrange <- range(sapply(bias_map[imaps], range))
breaks <- seq(xrange[1], xrange[2], length.out = 41)
lhist <- lapply(bias_map[imaps], hist, breaks = breaks, plot = FALSE)
ymax <- max(sapply(lhist, function(h) max(h$count))) * 1.1
par(mfrow = c(3, 2))
for (i in seq_along(imaps)) {
    hist(
        bias_map[[imaps[i]]],
        breaks = breaks, ylim = c(0, ymax), col = "grey",
        xlab = "", main = maintitle[i], sub = subtitle[i]
    )
}
dev.off()


pdf(paste0("figures/", var, "_mae_bias_map_era5.pdf"), height = 8)
colors <- pal_whitered
ncolors <- length(colors)
lvalues <- lapply(mae_bias_map[imaps], atan_trans)
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

# MAE bias Skill Score
pdf(paste0("figures/", var, "_mae_bias_diff_map_era5.pdf"), height = 8)
lvalues <- lapply(
  mae_bias_map[imaps], 
  function(x) atan_trans(x - mae_bias_map[[icombi_ref]])
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

# MAE gradient
pdf(paste0("figures/", var, "_mae_gradient_map_era5.pdf"), height = 8)
colors <- pal_whiteredblack
ncolors <- length(colors)
mae <- mae_gradient[imaps]
maintitle <- sprintf(
  "(%s) %s",
  letters[seq.int(length(imaps))],
  names(mae)
)
subtitle <- sprintf("MAE_g = %.2f ", mae)
lvalues <- lapply(mae_gradient_map[imaps], atan_trans)
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


pdf(paste0("figures/", var, "_mae_gradient_hist_era5.pdf"), height = 8)
xrange <- range(sapply(mae_gradient_map[imaps], range))
breaks <- seq(xrange[1], xrange[2], length.out = 41)
lhist <- lapply(mae_gradient_map[imaps], hist, breaks = breaks, plot = FALSE)
ymax <- max(sapply(lhist, function(h) max(h$count))) * 1.1
par(mfrow = c(3, 2))
for (i in seq_along(imaps)) {
    hist(
        mae_gradient_map[[imaps[i]]],
        breaks = breaks, ylim = c(0, ymax), col = "grey",
        xlab = "", main = maintitle[i], sub = subtitle[i]
    )
}
dev.off()


# MAE gradient Skill Score
pdf(paste0("figures/", var, "_mae_gradient_diff_map_era5.pdf"), height = 8)
colors <- rev(pal_violetwhitegreen)
ncolors <- length(colors)
lvalues <- lapply(
  mae_gradient_map[imaps], 
  function(x) atan_trans(x - mae_gradient_map[[icombi_ref]])
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
