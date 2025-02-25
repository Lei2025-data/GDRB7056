#library(ggplot2)
#library(grid)

findW <- function(grob, plot.width = 0, unit = "in", is.legend.overlap = TRUE, is.each = FALSE) {
  if ( ! grid::is.grob(grob) ) {
    grob <- ggplot2::ggplotGrob(grob)
  }
  plot.width <- unit(plot.width, unit)
  
  plot.index <- if ( "unit.list" %in% class(grob$widths) )
                  which(unlist(sapply(grob$widths, attr, "unit")) == "null")
                else if ( "unit_v2" %in% class(grob$widths) )  ## since R 4.0.0
                  which(grid::unitType(grob$widths) == "null")
                else
                  which(attr(grob$widths, "unit") == "null")
  plot.range <- range(plot.index)
  w.plot  <- grid::convertUnit(sum(grob$widths[plot.range[1]:plot.range[2]]), unit, valueOnly = T)
  w.plot  <- w.plot + grid::convertUnit(plot.width, unit, valueOnly = T) * length(plot.index)
  w.other <- grid::convertUnit(sum(grob$widths[-plot.range[1]:-plot.range[2]]), unit, valueOnly = T)
  
  w.lg <- NULL
  if ( is.legend.overlap ) {
    lg.index <- grep("guide-box", sapply(grob$grobs, function(x) x$name))
    if (length(lg.index)){
      lg <- grob$grobs[[lg.index]]
      w.lg <- grid::convertUnit(sum(lg$widths), unit, valueOnly = T)
    }else{
      w.lg <- 0
    }
    #   w.plot <- max(w.plot, w.lg)
  }
  
  if ( is.each ){
    return(c(w.plot, w.other, w.lg))
  }
  
  w.plot <- max(w.plot, w.lg)
  w <- w.other + w.plot
  return(w)
}

findH <- function(grob, plot.height = 0, unit = "in", is.legend.overlap = TRUE, is.each = FALSE) {
  if ( ! grid::is.grob(grob) ) {
    grob <- ggplot2::ggplotGrob(grob)
  }
  plot.height <- unit(plot.height, unit)
  
  plot.index <- if ( "unit.list" %in% class(grob$heights) )
                  which(unlist(sapply(grob$heights, attr, "unit")) == "null")
                else if ( "unit_v2" %in% class(grob$heights) )  ## since R 4.0.0
                  which(grid::unitType(grob$heights) == "null")
                else
                  which(attr(grob$heights, "unit") == "null")
  plot.range <- range(plot.index)
  h.plot <- grid::convertUnit(sum(grob$heights[plot.range[1]:plot.range[2]]), unit, valueOnly = T)
  h.plot <- h.plot + grid::convertUnit(plot.height, unit, valueOnly = T) * length(plot.index)
  h.other <- grid::convertUnit(sum(grob$heights[-plot.range[1]:-plot.range[2]]), unit, valueOnly = T)
  
  h.lg <- NULL
  if ( is.legend.overlap ) {
    lg.index <- grep("guide-box", sapply(grob$grobs, function(x) x$name))
    if ( length(lg.index) ){
      lg <- grob$grobs[[lg.index]]
      h.lg <- grid::convertUnit(sum(lg$heights), unit, valueOnly = T)
    }else{
      h.lg <- 0
    }
#    h.plot <- max(h.plot, h.lg)
  }
  
  if ( is.each ){
    return(c(h.plot, h.other, h.lg))
  }
  
  h.plot <- max(h.plot, h.lg)
  h <- h.other + h.plot
  return(h)
}


findWH <- function(grob, plot.size = 0,
                   plot.width = plot.size, plot.height = plot.size,
                   unit = "in", is.legend.overlap = TRUE, is.square = TRUE) {
  ws <- findW(grob, plot.width = plot.width, unit = unit, is.legend.overlap = is.legend.overlap, is.each = TRUE)
  hs <- findH(grob, plot.height = plot.height, unit = unit, is.legend.overlap = is.legend.overlap, is.each = TRUE)
  if ( is.legend.overlap ){
    if ( length(grob$theme$legend.position) == 2){
    }else if (sum(grob$theme$legend.position %in% c("top", "bottom"))){
      if ( ! is.null(ws[3]) || ! is.na(ws[3]) ) {
        if ( is.square ) hs[1] <- max(hs[1], ws[3])
        ws[1] <- max(ws[1], ws[3])
      }
    }else{
      if ( ! is.null(hs[3]) || ! is.na(hs[3]) ) {
        if ( is.square ) ws[1] <- max(ws[1], hs[3])
        hs[1] <- max(hs[1], hs[3])
      }
    }
  }
  w <- ws[1] + ws[2]
  h <- hs[1] + hs[2]
  return(c(w,h))
}

  

SavePlot <- function(plot=last_plot(), filename = NULL, plot.size = 4,
                     plot.width = plot.size, plot.height = plot.size, units = "in",
                     is.legend.overlap = TRUE, is.square = TRUE, ...) {
  wh <- findWH(plot, plot.width = plot.width, plot.height = plot.height, unit = units,
               is.legend.overlap = is.legend.overlap, is.square = is.square)
  w <- grid::convertUnit(grid::unit(wh[1], units), unitTo = "in", valueOnly = T)
  h <- grid::convertUnit(grid::unit(wh[2], units), unitTo = "in", valueOnly = T)
  ggsave(plot = plot, filename = filename, width = w, height = h, units = "in", limitsize = FALSE, ...)
}


