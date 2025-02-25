
bin_stat <- function(data, bin = 256){
		data <- as.numeric(data)
		b <- floor((data - min(data) ) / diff(range(data)) * bin)
		b[b == bin] <- bin - 1
		names(b) <- data
		return(b)
}

core_otsu <- function(dt){
  weight.bg <- sapply(seq(0, max(dt)), function(i) sum(dt <= i) / length(dt))
  weight.fg <- 1 - weight.bg
  
  mean.bg <- sapply(seq(0, max(dt)), function(i) mean(dt[dt <= i]))
  mean.fg <- sapply(seq(0, max(dt)), function(i) mean(dt[dt > i]))
  
  var.between <-  weight.bg * weight.fg * ( mean.bg - mean.fg ) ^ 2
  
  thres <- which.max(var.between)
  return(thres)
}

Otsu.bak <- function(data, digits = 0, is.return.range = TRUE){
  dt0 <- table(round(data, digits = digits))
  dt <- dt0
  dt[dt > 100] <- 100
  
  thres <- core_otsu(dt)
  
  if ( is.return.range )
    return(range(as.numeric(names(dt0)[dt0 >= thres])))
  else
    return(thres)
}

Otsu <- function(data, bin = 256, is.return.range = TRUE){
		dt0 <- bin_stat(data, bin = bin)
		dt  <- table(dt0)
		dt[dt > 100] <- 100

		thres <- core_otsu(dt)

		if ( is.return.range )
				return(range(as.numeric(names(dt0)[dt0 %in% names(dt)[dt >= thres]])))
		else
				return(as.numeric(names(dt)[dt == thres - 1]))
}

#Otsu <- Otsu3
#formals(Otsu)$bin <- 100


autothres <- function(data, name = c("default", "nCount_RNA", "nFeature_RNA", "percent.mito"), bin = 100, digits = 2){
		thres <- Otsu(data = data, bin = bin)
		name <- match.arg(name)
		if ( name == "nCount_RNA" ) {
				up   <- signif(thres[2], digits)
				down <- -Inf
		} else if ( name == "nFeature_RNA" ) {
				up   <- signif(thres[2], digits)
				down <- signif(thres[1], digits)
				if ( up > 200 ) {
						down <- max(down, 200)
				}
		} else if ( name == "percent.mito" ) {
				if ( all(data <= 1) )
						up <- if(thres[2] <= 0.1) 0.1 else if (thres[2] <= 0.25) 0.25 else signif(thres[2], digits)
				else
						up <- if(thres[2] <= 10) 10 else if (thres[2] <= 25) 25 else signif(thres[2], digits)
						down <- -Inf
		}else{
				up   <- signif(thres[2], digits)
				down <- signif(thres[1], digits)
		}
		return(c(down, up))
}


