#' Fits a gamma mixture model to mIF imaging data with multiple slides at one time
#'
#' Takes a cfGMM model and forces the right-most component
#' probabilities to be monotonic increasing with respect
#' to the input values.
#'
#' @param expressionMarkers An ncell by nmarker matrix of expression values. Assumed to be coarsely normalized and transformed.
#' @param slide Single vector of length n. Identifies slide number for each cell.
#' @param subBatch If there are multiple subBatch on a slide, subBatch can be used to return probability estimates independently for each region.
#' @param boundaryMarkers A nmarker list of 4x4 matrices giving the boundaries for the modes of the unexpressed and expressed cell distributions.
#' @param qboundaryMarkers A nmarker list of 4x4 matrices giving the qauntile boundaries for the modes for the unexpressed and expressed cell distributions.
#' @param n.cores Number of cores. Default to one less than the number of cores available on the machine.
#' @param ... Arguments passed to cfGMM function
#' @importFrom cfGMM cfGMM
#' @importFrom stats quantile
#' @importFrom reactable reactable
#' @importFrom parallel mclapply detectCores
#' @export
#' @details Fits cfGMM models to each marker channel in a matrix of marker channels for multiple slides
groupGluster <- function(expressionMarkers, slide, boundaryMarkers=NULL, qboundaryMarkers=NULL, subBatch=NULL, n.cores=NULL, ...){
  if(is.null(n.cores)){n.cores = detectCores()-1}
  if(any(is.na(slide))){expressionMarkers = expressionMarkers[!is.na(slide)]; slide = slide[!is.na(slide)]}
  expressionMarkers = split(expressionMarkers, slide)
  constrCfGMMbunch = mclapply(expressionMarkers, function(x) gluster(x, boundaryMarkers, qboundaryMarkers, subBatch),   mc.cores = n.cores)
  class(constrCfGMMbunch) = c('groupGluster', class(constrCfGMMbunch))
  return(constrCfGMMbunch)
}

#' Fits a gamma mixture model to mIF imaging data with multiple slides at one time
#'
#' Takes a cfGMM model and forces the right-most component
#' probabilities to be monotonic increasing with respect
#' to the input values.
#'
#' @param x A groupGluster object
#' @param marker Select which marker to plot. Can be a character or integer.
#' @param component Integer specifying which component to plot, 1 is unexpressed nonzero cells, 2 is expressed cells.
#' @param diagnostic logical indicating whether to create the diagnostic plot. Default value is TRUE.
#' @param histogram logical indicating whether to create the slide histograms.
#' @param subBatch If there are multiple subBatch on a slide, subBatch can be used to return probability estimates independently for each region.
#' @param title Title for the plot. Default is the marker name.
#' @param color color for points.
#' @param p a ggplot2 object to add to.
#' @param print logical whether to display the plot. Default value TRUE.
#' @param ... Arguments passed to XX
#' @importFrom ggplot2 aes ggplot ggtitle stat_function geom_vline unit annotation_custom geom_point xlab ylab
#' @importFrom hrbrthemes theme_ipsum
#' @importFrom stats na.omit
#' @export
#' @details Various diagnostic and QC plots for groupGluster fits.
plot.groupGluster <- function(x, marker=1, component=2, diagnostic=TRUE, histogram=FALSE, title=NULL, color='grey', p=NULL, print=TRUE, ...){
  if(is.numeric(marker)){ marker=colnames(x[[1]][["expressionX"]])[marker] }
  if(is.null(title)) title=marker
  if(diagnostic){
    plot.df = as.data.frame(do.call(rbind, lapply(x, function(y){
      do.call(rbind, lapply(y$params[[marker]], function(z) z[,paste0('comp',component)]) )
    } ) ))
    plot.df$mode <- (plot.df$alpha-1)*plot.df$beta

    p1 <- ggplot(plot.df, aes(x=mode,y=lambda)) +
      geom_point(alpha=0.4, color=color, size=3) +
      xlab("Mode of Expressed Component")+ylab("Lambda of Expressed")+
      theme_ipsum(axis_title_size = 12) + ggtitle(title[1])
    if(print){
      print(p1)
    } else {
      return(p1)
    }
  }
  if(histogram)
    invisible(capture.output(sapply(names(x), function(y) plot(x[[y]], marker = markerind, title = x ) )))
}
