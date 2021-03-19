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
  def_val <- list(seq(rtlims[1],rtlims[2],length=10),rep(0,10))
  if(is.null(rtlims)){
    vals <- apply(mzlims,1,function(x,xraw){
      rawEIC(xraw,mzrange=xl)
    },xraw=x$raw)
  }else{
    lims <- cbind(mzlims,rtlims)
    vals <- apply(lims,1,function(x,xraw){
      mzl <- x[c(1,2)]
      rtl <- x[c(3,4)]*60
      vals <- tryCatch(rawEIC(xraw,mzrange=mzl,rtrange=rtl),error=function(e){
        return(NA)
      })
      vals
    },xraw=x$raw)
  }
  vals <- lapply(vals,function(x,vtime){
    if(is.na(x)) return(NA)
    x1 <- vtime[x[[1]]]
    x2 <- x[[2]]
    list(x1,x2)
  },vtime=x$raw@scantime)

  missing <- which(sapply(vals,function(x){
    if(length(x)!=1) return(FALSE)
    return(TRUE)
  }))

  for(mm in missing){
    vals[[mm]] <- def_val
  }
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

#' @export
plot_map.SLAWsample <- function(x,mzlim,rtlim=NULL,mz_margin=0.1,rt_margin=1/60){
  p <- peaks(x,mzlim=mzlim,rtlim=rtlim)
  rmat <- rawMat(x$raw,mzlim,rtlim*60)
  rmat[,"time"] <- rmat[,"time"]/60
  p <- as.data.frame(p[,c("rt_min","rt_max","mz_min","mz_max")])
  p[["rt_min"]] <- p[["rt_min"]]-rt_margin
  p[["rt_max"]] <- p[["rt_max"]]+rt_margin
  p[["mz_max"]] <- p[["mz_max"]]+mz_margin
  p[["mz_min"]] <- p[["mz_min"]]-mz_margin
  p[["color"]] <- scales::hue_pal()(nrow(p))
  ggp <- ggplot(data=as.data.frame(rmat),aes(x=time,y=mz,color=intensity))+geom_point()+viridis::scale_colour_viridis(begin=1,end=0.1,trans="log",option="magma",name="Intensity")+
    geom_rect(data=p,mapping = aes(xmin=rt_min,xmax=rt_max,ymin=mz_min,ymax=mz_max,fill=color),alpha=0.2,inherit.aes = FALSE)+
    xlab("TIme(min)")+ylab("m/z")+scale_fill_discrete(guide = FALSE)
  # ggplot()+geom_rect(data=p,mapping = aes(xmin=rt_min,xmax=rt_max,ymin=mz_min,ymax=mz_max),inherit.aes = FALSE)
  plot(ggp)
  invisible(ggp)
}



#
# ss <- SLAWsample(path,peaktable)
# ss2 <-  SLAWsample(path2,peaktable2)
#
# mzlim <- c(250,260)
# rtlim <- c(1,60)
#
#
# schromatogram(ss,mzlim=mzzlim)



