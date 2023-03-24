#' Fits a gamma mixture model to aid clustering of mIF imaging data
#'
#' Takes a cfGMM model and forces the right-most component
#' probabilities to be monotonic increasing with respect
#' to the input values.
#'
#' @param gluster.fit gluster or groupGluster object.
#' @param out.all Whether to print out the full convergence result for all markers and slides, or just summary of non-convergence.
#' @param return.result Whether to return the convergence result into an object.
#' @importFrom reactable reactable
#' @export
#' @details convergence check to gluster or batch gluster fits
convCheck <- function(gluster.fit, out.all=FALSE, return.result=FALSE){
  class.fit <- class(gluster.fit)[[1]]
  if(class.fit=="gluster"){
    conv.vec <- sapply(gluster.fit$fit, function (x) x$convergence)
    if(out.all){
      print(conv.vec)
    } else {
      if(!all(conv.vec)){
        print(conv.vec[which(!conv.vec)])
      } else {print("All converged")}
    }
  } else if (class.fit=="groupGluster"){
    conv.tab <- sapply(gluster.fit, function(x) sapply(x$fit, function(y){ y$convergence}))
    if(out.all){
      reactable(conv.tab)
    } else{
      if(!all(conv.tab)){
        uncov.ind <- which(!conv.tab)
        uncov.names <- expand.grid(dimnames(conv.tab))[uncov.ind,]
        colnames(uncov.names) <- c("marker", "slide")
        print(uncov.names, row.names = FALSE)
      } else {print("All converged.")}
    }
  } else {print("Wrong input type :o")}
  if(return.result){return(conv.tab)}
}

#' Diagnostic plot
#'
#' Takes a cfGMM model and forces the right-most component
#' probabilities to be monotonic increasing with respect
#' to the input values.
#'
#' @param gluster.fit gluster or groupGluster object.
#' @param out.all Whether to print out the full convergence result for all markers and slides, or just summary of non-convergence.
#' @param return.result Whether to return the convergence result into an object.
#' @importFrom reactable reactable
#' @export
#' @details convergence check to gluster or batch gluster fits
convCheck <- function(gluster.fit, out.all=FALSE, return.result=FALSE){
  class.fit <- class(gluster.fit)[[1]]
  if(class.fit=="gluster"){
    conv.vec <- sapply(gluster.fit$fit, function (x) x$convergence)
    if(out.all){
      print(conv.vec)
    } else {
      if(!all(conv.vec)){
        print(conv.vec[which(!conv.vec)])
      } else {print("All converged")}
    }
  } else if (class.fit=="groupGluster"){
    conv.tab <- sapply(gluster.fit, function(x) sapply(x$fit, function(y){ y$convergence}))
    if(out.all){
      reactable(conv.tab)
    } else{
      if(!all(conv.tab)){
        uncov.ind <- which(!conv.tab)
        uncov.names <- expand.grid(dimnames(conv.tab))[uncov.ind,]
        colnames(uncov.names) <- c("marker", "slide")
        print(uncov.names, row.names = FALSE)
      } else {print("All converged.")}
    }
  } else {print("Wrong input type :o")}
  if(return.result){return(conv.tab)}
}
