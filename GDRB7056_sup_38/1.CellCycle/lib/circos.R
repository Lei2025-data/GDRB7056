circ.plot <- function(data.use, color.use = NULL, scale.text.size = 1, show.axes = FALSE){
  g.name <- unlist(lapply(data.use[,1:2], as.character))
  max.name <- g.name[which.max(nchar(g.name))]
  mar <- strwidth(s = max.name, cex = scale.text.size, units = "inches") * 1.2 / (7/2)
  
  library(circlize)
  pdf() ## ensure recordPlot only capture current plot
  dev.control(displaylist="enable")
  
  if ( show.axes) {
    mar <- mar + 0.1 * 1.2 / 3.5
    circos.par("canvas.xlim" = c(-1 - mar, 1 + mar),
               "canvas.ylim" = c(-1 - mar, 1 + mar),
               points.overflow.warning = FALSE)
    
    chordDiagram(data.use, annotationTrack = "grid",
                 grid.col = color.use, order = names(color.use),
                 link.sort = TRUE, link.decreasing = TRUE)
    circos.track(track.index = 1, bg.border = NA,
                 panel.fun = function(x, y) {
                   circos.axis(labels.cex = scale.text.size * 0.3)
                   
                   xlim = get.cell.meta.data("xlim")
                   ylim = get.cell.meta.data("ylim")
                   xplot = get.cell.meta.data("xplot")
                   by = ifelse(abs(xplot[2] - xplot[1]) > 30, 0.2, ifelse(abs(xplot[2] - xplot[1]) > 10, 0.5, 1))
                   
                   circos.lines(xlim, c(2, 2), lty = 3) # dotted line
                   for(p in seq(by, 1, by = by)) {
                     circos.text(x = p * (xlim[2] - xlim[1]) + xlim[1], y = 2 + 0.1,
                                 labels = paste0(p * 100, "%"), cex = scale.text.size * 0.3,
                                 adj = c(0.5, 0), niceFacing = TRUE)
                   }
                   
                   circos.text(x = CELL_META$xcenter,
                               y = 2.5, 
                               labels = CELL_META$sector.index,
                               facing = "clockwise",
                               cex = scale.text.size,
                               adj = c(0, 0.5),
                               niceFacing = TRUE)
                 }
    )	  
  } else {
    circos.par("canvas.xlim" = c(-1 - mar, 1 + mar),
               "canvas.ylim" = c(-1 - mar, 1 + mar),
               points.overflow.warning = FALSE)
    
    chordDiagram(data.use, annotationTrack = "grid",
                 grid.col = color.use,
                 order = names(color.use),
                 link.sort = TRUE, link.decreasing = TRUE
    )
    circos.track(track.index = 1, bg.border = NA,
                 panel.fun = function(x, y) {
                   circos.text(x = CELL_META$xcenter,
                               y = CELL_META$ylim[2] + 0.2,
                               labels = CELL_META$sector.index,
                               facing = "clockwise",
                               cex = scale.text.size,
                               adj = c(0, 0.5),
                               niceFacing = TRUE)
                 }
    )
  }
  circos.clear()
  
  p <- recordPlot()
  dev.off()
  return(p)
}

.plotCircos <- function(data, color.1 = NULL, color.2 = NULL, surfix.1 = NULL, surfix.2 = NULL, sep = ".", ...) {
  dt <- reshape2::melt(table(data[,c(1,2)]))
  dt[[1]] <- factor(dt[[1]], levels = if(is.factor(data[[1]])) levels(data[[1]]) else unique(data[[1]]))
  dt[[2]] <- factor(dt[[2]], levels = if(is.factor(data[[2]])) levels(data[[2]]) else unique(data[[2]]))

  if ( ! is.null(surfix.1) ) {
    levels(dt[[1]]) <- paste(surfix.1, levels(dt[[1]]), sep = sep)
  }
  if ( ! is.null(surfix.2) ) {
    levels(dt[[2]]) <- paste(surfix.2, levels(dt[[2]]), sep = sep)
  }
  
  if ( is.null(color.1) ) color.1 <- scales::hue_pal()(length(unique(dt[[1]])))
  if ( is.null(color.2) ) color.2 <- scales::hue_pal()(length(unique(dt[[2]])))
  if ( is.null(names(color.1)) ) names(color.1) <- levels(dt[[1]])
  if ( is.null(names(color.2)) ) names(color.2) <- levels(dt[[2]])
  color.rainbow <- c(color.1, color.2)
  
  p <- circ.plot(dt, color.use = color.rainbow, ...)
  return(p)
}

