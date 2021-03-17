select_samples <- function(x,...) {
  UseMethod("select_samples",x)
}

schromatogram <- function(x,...){
  UseMethod("schromatogram",x)
}

schromatograms <- function(x,...){
  UseMethod("schromatograms",x)
}

peaks <- function(x,...){
  UseMethod("peaks",x)
}

mpeaks <- function(x,...){
  UseMethod("mpeaks",x)
}

get_eics <- function(x,...){
  UseMethod("get_eics",x)
}


plot_peaks <- function(x,...){
  UseMethod("plot_peaks",x)
}

plot_features <- function(x,...){
  UseMethod("plot_features",x)
}


#' SLAWexplorer TO Visualise SLAW Results
#'
#' @param path The path of the SLAW output directory
#' @param raw_path The path of the rawe files for EICs visualisation
#'
#' @return A SLAWexplorer object
#' @export
#'
#' @examples
#' print("Examples to be put here")
SLAWexplorer <- function(path,raw_path){
  #We check if if is a valid SLAW directory
  if(!dir.exists(path)){
    stop("'path' be a valid directory.")
  }

  # This is the path fothe database
  path_db <- file.path(path,"processing_db.sqlite")
  query <- "SELECT processing.id,path,output_ms,output_ms2 FROM samples INNER JOIN processing on samples.id=processing.sample WHERE level='MS1' AND output_ms!='NOT PROCESSED' AND valid=1"
  query_dm <- "SELECT annotated_peaktable_full,fused_msms FROM peakpicking"
  # Extracting samples and tables informations
  dbi <- dbConnect(RSQLite:::SQLite(),path_db)
  samples <- dbGetQuery(dbi,query)
  res_dm <- dbGetQuery(dbi,query_dm)
  path_peaktable <- str_replace(res_dm[1],"/output",path)
  path_mgf <- str_replace(res_dm[2],"/output",path)

  #We replace the path by the necessary path
  samples$path <- str_replace(samples$path,"/input",raw_path)
  samples$output_ms <- str_replace(samples$output_ms,"/output",path)
  samples$output_ms2 <- str_replace(samples$output_ms2,"/output",path)
  dbDisconnect(dbi)

  selected_samples <- numeric()

  #CONVENTION spectrum is NULL is empty
  spectra <- NULL
  if(file.exists(path_mgf)) spectra <- Spectra(path_mgf, source = MsBackendMgf())
  #THe MGF and the datamatrix are always read
  dm <- fread(path_peaktable)
  raws <- list()
  structure(list(infos=samples,selected=selected_samples,datamatrix=dm,ms2=spectra,samples=raws),class="SLAWexplorer")
}



rload_file <- function(slexp,idx){
  path <- slexp$infos$path[idx]
  peaktable <- slexp$infos$output_ms[idx]
  if(!path %in% slexp$samples){
    return(SLAWsample(path,peaktable))
  }else{
    return(slexp[[path]])
  }
}

#' Select samples from a SLAWexplorer acquisition
#'
#' @param slexp A SLAWexplorer object
#' @param samples A set of samples id as integer (Their position in the columns)
#'
#' @return
#' @export
#'
#' @examples
#' print("Examples to be put here")
select_samples.SLAWexplorer<- function(slexp,samples){
  #We check if the ides makes sens
  if(!is.numeric(samples)&any(samples>nrow(slexp$samples))){
    stop("Invalid 'samples'.")
  }
  slexp$select <- samples
  paths <- slexp$infos$path[samples]
  slexp$samples <- lapply(samples,FUN = rload_file,slexp=slexp)
  slexp
}


schromatogram.SLAWexplorer <- function(slexp,mzlim,rtlim=NULL){
  lapply(slexp$samples,schromatogram,mzlim=mzlim,rtlim=rtlim)
}

schromatograms.SLAWexplorer <- function(slexp,mzlims,rtlims=NULL){
  lapply(slexp$samples,schromatograms,mzlims=mzlims,rtlims=rtlims)
}

peaks.SLAWexplorer <- function(slexp,mzlim,rtlim){
  lapply(slexp$samples,peaks,mzlim=mzlim,rtlim=rtlim)
}


mpeaks.SLAWexplorer <- function(slexp,mzlims,rtlims){
  lapply(slexp$samples,mpeaks,mzlim=mzlims,rtlim=rtlims)
}


ordered_or_equal <- function(x,decreasing=FALSE){
  if(decreasing){
    return(all(diff(x)<=0))
  }else{
    return(all(diff(x)>=0))
  }
}

#' Title
#'
#' @param slexp A SLAWexplorer object
#' @param mzlim An mz limit
#' @param rtlim A rtlimit
#'
#' @return None
#' @export
#'
#' @examples
#' print("Examples to be put here")
plot_peaks.SLAWexplorer <- function(slexp,mzlim,rtlim,colors,...){
  #extracting informat♦ion
  chroms <- schromatogram(slexp,mzlim,rtlim)
  speaks <- peaks(slexp,mzlim,rtlim)
  if(missing(colors)){
    colors <- rainbow(length(chroms))
  }
  .plot_peaks(chroms,speaks,colors,...)
}

.get_peaks_limits <- function(slexp,fids,mz_margin=0.01,rt_margin=0.1){
  subfeatures <-  slexp$datamatrix[fids]
  mzmin <- subfeatures[["min_mz"]]-mz_margin
  mzmax <- subfeatures[["max_mz"]]+mz_margin
  rtmin <- subfeatures[["min_rt"]]-rt_margin
  rtmax <- subfeatures[["max_rt"]]+rt_margin
  #Filtering negative values
  rtmin[rtmin<0] <- 0
  mzmin[mzmin<0] <- 0

  return(list(cbind(mzmin,mzmax),cbind(rtmin,rtmax)))
}


layoutMatrix <- function(n) {
  if (n == 1) {
    return(matrix(c(1)))
  }
  if (n == 2) {
    return(matrix(c(1, 2), nrow = (2)))
  }
  if (n == 3) {
    return(matrix(c(1, 2, 3), nrow = (3)))
  }
  if (n == 4) {
    return(matrix(c(1, 2, 3, 4), nrow = (2), byrow = TRUE))
  }
  if (n == 5) {
    return(matrix(c(1, 2, 3, 4, 5, 6), nrow = (2), byrow = TRUE))
  }
  if (n == 6) {
    return(matrix(c(1, 2, 3, 4, 5, 6), nrow = (2), byrow = TRUE))
  }
  if (n == 7) {
    return(matrix(
      c(1, 2, 3, 4, 5, 6, 7, 8, 9),
      nrow = (3),
      byrow = TRUE
    ))
  }
  if (n == 8) {
    return(matrix(
      c(1, 2, 3, 4, 5, 6, 7, 8),
      nrow = (2),
      byrow = TRUE
    ))
  }
  if (n == 9 | n > 12) {
    return(matrix(
      c(1, 2, 3, 4, 5, 6, 7, 8, 9),
      nrow = (3),
      byrow = TRUE
    ))
  }
  if (n == 10) {
    return(matrix(
      c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12),
      nrow = (4),
      byrow = TRUE
    ))
  }
  if (n == 11) {
    return(matrix(
      c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12),
      nrow = (4),
      byrow = TRUE
    ))
  }
  if (n == 12) {
    return(matrix(
      c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12),
      nrow = (4),
      byrow = TRUE
    ))
  }
}



#' Plot feature across the selected files
#'
#' @param slexp
#' @param features
#'
#' @return
#' @export
#'
#' @examples
#' print("Examples to be put here.")
plot_features.SLAWexplorer <- function(slexp,fids,mz_margin=0.005,rt_margin=0.05,layout_matrix=NULL,...){
  if(is.null(layout_matrix)){
    layout_matrix <- layoutMatrix(length(fids))
  }
  layout(layout_matrix)
  peaks_limit <- .get_peaks_limits(slexp,fids,mz_margin=mz_margin,rt_margin=rt_margin)
  chroms <- schromatograms(slexp,mzlims=peaks_limit[[1]],rtlims=peaks_limit[[2]],...)
  peaks <- mpeaks(slexp,mzlims=peaks_limit[[1]],rtlims=peaks_limit[[2]],...)
  for(fid in seq_along(fids)){
    sel_chroms <- lapply(chroms,"[[",i=fid)
    sel_peaks <- lapply(peaks,"[[",i=fid)
    .plot_peaks(sel_chroms,sel_peaks,...)
  }
}

get_eics.SLAWexplorer <- function(slexp,fids,mz_margin=0.005,rt_margin=0.05,...){
  peaks_limit <- .get_peaks_limits(slexp,fids,mz_margin=mz_margin,rt_margin=rt_margin)
  chroms <- schromatograms(slexp,mzlims=peaks_limit[[1]],rtlims=peaks_limit[[2]],...)
  peaks <- mpeaks(slexp,mzlims=peaks_limit[[1]],rtlims=peaks_limit[[2]],...)
  for(sid in seq_along(chroms)){
    temp <- mapply(chroms[[sid]],peaks[[sid]],FUN=.filter_chromatogram,SIMPLIFY = FALSE)
    chroms$in_peak <- temp
  }
  return(chroms)
}

.filter_chromatogram <- function(chrom,peaks){
  peaks <- peaks[order(peaks[,"rt_min"],decreasing = FALSE)]
  # We try to plot the peaks
  rtbins <- as.numeric(t(peaks[,c("rt_min","rt_max")]))
  if(!ordered_or_equal(rtbins)) return(rep(FALSE,length(chrom[[1]])))
  vb <- .bincode(chrom[[1]]/60,rtbins)
  return(!is.na(vb))
}


.plot_peaks <- function(chromatograms,peaks,colors,title="EICs",...){
  if(missing(colors)){
    colors <- rainbow(length(chromatograms))
  }

  #min and max
  tmin <- min(sapply(chromatograms,function(x){x[[1]][1]}))
  tmax <- max(sapply(chromatograms,function(x){x[[1]][length(x[[1]])]}))
  imax <- max(sapply(chromatograms,function(x){max(x[[2]])}))

  #Plotting the peaks
  plot(NULL,xlim=c(tmin,tmax),ylim=c(0,imax),xlab="Time",ylab="Count/Intensity",main=title,...)
  for(idx in seq_along(chromatograms)){
    lines(chromatograms[[idx]][[1]],chromatograms[[idx]][[2]],type="l",lwd=1,lty=1)
    #We bin the peak
    cpeaks <- peaks[[idx]]
    if(nrow(cpeaks)==0) next

    in_peak <- .filter_chromatogram(chromatograms[[idx]],cpeaks)
    cpeaks <- cpeaks[order(cpeaks[,"rt_min"],decreasing = FALSE)]
    # We try to plot the peaks
    rtbins <- as.numeric(t(cpeaks[,c("rt_min","rt_max")]))
    if(!ordered_or_equal(rtbins)) next
    vb <- .bincode(chromatograms[[idx]][[1]]/60,rtbins)
    pos_peaks <- which(!is.na(vb))
    if(length(pos_peaks)==0) next
    lines(chromatograms[[idx]][[1]][pos_peaks],chromatograms[[idx]][[2]][pos_peaks],
          type="l",lwd=1,lty=1,col=colors[idx])
  }
}


# Test on Sammy
#
# path <- "U:/users/Sammy/Alexis/Slaw processed"
# raw_path <- "U:/users/Sammy/Alexis/Slaw processed/mzML_unfiltered"
#
#
# sexplo <- SLAWexplorer(path,raw_path)
#
# sexplo <- select_samples(sexplo,c(1,2,5,3,8,10))
#
#
# mzlims <- matrix(c(200,210,400,410),nrow = 2)
#
#
#
# mzlims
# vbv <- schromatograms(sexplo,mzlims)
# length(vbv[[1]])
# length(sexplo$samples)
#
# plot_peaks(sexplot,mzlim,rtlim=c())
#
# peaks(sexplo,c(400,450),NULL)
#
#
# plot_features(sexplo,sample.int(700,size=200))
