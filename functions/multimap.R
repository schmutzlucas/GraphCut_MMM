  multimap <- function(
    lon, lat, lvalues,
    main, sub,
    zlim, colors,
    breaks, labels,
    outersplit = rbind(c(0, 1, 1/8, 1), c(0, 1, 0, 1/8)), 
    innersplit = c(3,2),
    mai = c(0.6732, 0.5412, 0.5412, 0.2772),
    legend_par = list(mar=c(8, 1, 1, 1), plt = c(0.1, 0.9, 0.9, 1))
  ){
    fields::set.panel() # reset plotting device
    iouter <- split.screen(outersplit)
    iinner <- split.screen(innersplit, screen = iouter[1])
    ncolors <- length(colors)
    for(icombi in seq_along(lvalues)){
      screen(iinner[icombi])
      par(mai = mai)
      print(sub[icombi])
      image(
        lon, lat, lvalues[[icombi]],
        xlab = '', ylab = '',
        main = main[icombi],
        zlim = zlim,
        col = colors,
        breaks = breaks
        # axis.args = list(at=breaks, labels = sprintf("%.2f", atan_inv(breaks, zmax)), cex.axis = 0.5)
      )
      mtext(sub[icombi], side = 1, line = 2)
      maps::map("world2", add=T)
    }
    screen(iouter[2])
    par(legend_par)
    at <- seq.int(ncolors) - 0.5 
    breaks_legend = seq.int(0, ncolors)
    if(any(is.na(as.numeric(labels)))){
      ticks <- at
    } else {
      ticks <- breaks_legend
    }
    image(
      x = at,
      y = 1, 
      z = matrix(at, ncol = 1),
      col = colors,
      breaks = breaks_legend,
      xaxt = "n", yaxt = "n", xlab = "", ylab = ""
    )
    axis(1, at = ticks, labels = labels, las = 3)
    close.screen(all=TRUE)
  }
