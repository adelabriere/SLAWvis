#
# path <- "U:/users/Sammy/Alexis/Slaw processed/mzML_unfiltered/10_1.mzML"
# peaktable <- "U:/users/Sammy/Alexis/Slaw processed/CENTWAVE/peaktables/10_1.csv"
#
#
# path2 <- "U:/users/Sammy/Alexis/Slaw processed/mzML_unfiltered/16_1.mzML"
# peaktable2 <- "U:/users/Sammy/Alexis/Slaw processed/CENTWAVE/peaktables/16_1.csv"

SLAWsample <- function(path,peaktable){
  #Load an xraw in memory
  nsd <- xcmsRaw(path)
  peaks <- fread(peaktable)
  structure(list(peaks=peaks,raw=nsd),class="SLAWsample")
}

#Function to retrieve a chromatogram
#' @export
schromatogram.SLAWsample <- function(x,mzlim=NULL,rtlim=NULL){
  if(is.null(mzlim)){
    mzlim <- c(0,10000)
  }
  reic <- rawEIC(x$raw,mzrange=mzlim,rtrange=rtlim)
  reic[[1]] <- xraw@scantime[reic$scan]
  return(reic)
}

#Retrieve multiple chromatograms from a single samples
#' @export
schromatograms.SLAWsample <- function(x,mzlims,rtlims=NULL){
  if(is.null(rtlims)){
    vals <- apply(mzlims,1,function(x,xraw){
      rawEIC(xraw,mzrange=xl)
    },xraw=x$raw)
  }else{
    lims <- cbind(mzlims,rtlims)
    vals <- apply(lims,1,function(x,xraw){
      mzl <- x[c(1,2)]
      rtl <- x[c(3,4)]*60
      rawEIC(xraw,mzrange=mzl,rtrange=rtl)
    },xraw=x$raw)
  }
  vals <- lapply(vals,function(x,vtime){
    x1 <- vtime[x[[1]]]
    x2 <- x[[2]]
    list(x1,x2)
  },vtime=x$raw@scantime)
  vals
}


#We return binary vector with the peak detected
#' @export
peaks.SLAWsample <- function(x,mzlim=NULL,rtlim=NULL){
  if(is.null(rtlim)){
    rtlim <- c(-1,100000)
  }
  if(is.null(mzlim)){
    mzlim <- c(-1,100000)
  }
  p <- .extract_peaks(x$peaks,mzlim,rtlim)
  return(p)
}

#We return binary vector with the peak detected
#' @export
mpeaks.SLAWsample <- function(x,mzlims,rtlims=NULL){
  if(is.null(rtlims)){
    rtlims <- cbind(rep(-1,nrow(mzlims)),rep(100000,nrow(mzlims)))
  }
  return(apply(cbind(mzlims,rtlims),1,function(x,peaks){
    mzl <- x[c(1,2)]
    rtl <- x[c(3,4)]
    return(.extract_peaks(peaks,mzl,rtl))
  },peaks=x$peaks))
}

.extract_peaks <- function(peaks,mzlim,rtlim){
  sel_peaks <- which(peaks$rt_min<rtlim[2]&peaks$rt_max>rtlim[1]&
                       peaks$mz_min<mzlim[2]&peaks$mz_max>mzlim[1])
  peaks <- peaks[sel_peaks,]
  return(peaks)
}
#
#
# ss <- SLAWsample(path,peaktable)
# ss2 <-  SLAWsample(path2,peaktable2)
#
# mzlim <- c(250,260)
# rtlim <- c(1,60)
#
#
# schromatogram(ss,mzlim=mzzlim)



