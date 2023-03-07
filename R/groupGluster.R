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
  cells = split(expressionMarkers, slide)
  constrCfGMMbunch = mclapply(expressionMarkers, function(x) gluster(x[,nzNormedMarkers], boundaryMarkers, qboundaryMarkers, subBatch),   mc.cores = n.cores)
  trash = sapply(names(constrCfGMMbunch), function(x) plot.gluster(constrCfGMMbunch[[x]], marker = 1, boundary = quantile(constrGMMfit$expressionX[,1], probs=quantileBoundaries[[1]][2,1]), title = x ) )
  reactable(sapply(constrCfGMMbunch, function(x) sapply(x$fit, function(y){ y$convergence})))
  class(constrCfGMMbunch) = c('groupGluster', class(constrCfGMMbunch))
  return(constrCfGMMbunch)
}
