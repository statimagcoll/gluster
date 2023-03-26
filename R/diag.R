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

#' Diagnostic plot - contrast
#'
#' @param fit1 A groupGluster object.
#' @param fit2 Another groupGluster object, possibly that is a refit of fit1 with updated constraints.
#' @param marker Marker to be plotted, defaults to the first one.
#' @param subBatch subBatch to be plotted, defaults to the first one.
#' @param title Title of plot.
#' @param add.table logical indicating whether to add table of fit2 parameter values.
#' @param fit.names Names of fit objects. Default to c("fit1", "fit2")
#' @param ... Arguments passed to plot.gluster
#' @importFrom ggplot2 aes ggplot ggtitle geom_histogram after_stat stat_function geom_vline unit annotation_custom geom_segment scale_color_manual
#' @importFrom hrbrthemes theme_ipsum
#' @importFrom gridExtra tableGrob
#' @import patchwork
#' @details Plot a histogram of a gluster model. Plots one model.
#' Takes a gluster object and plots the histogram with the fitted model and parameter values.
diag_contrast <- function(fit1, fit2, marker=1, component=2, title="Contrastive Diagnostic Plot", fit.names=c("fit1", "fit2"),boundaries=NULL,...){
  plot.df.fit1 = as.data.frame(do.call(rbind, lapply(fit1, function(y){
    do.call(rbind, lapply(y$params[[marker]], function(z) z[,paste0('comp',component)]))})))
  plot.df.fit1$mode <- (plot.df.fit1$alpha-1)*plot.df.fit1$beta
  plot.df.fit1$fit <- fit.names[1]
  plot.df.fit2 = as.data.frame(do.call(rbind, lapply(fit2, function(y){
    do.call(rbind, lapply(y$params[[marker]], function(z) z[,paste0('comp',component)]))})))
  plot.df.fit2$mode <- (plot.df.fit2$alpha-1)*plot.df.fit2$beta
  plot.df.fit2$fit <- fit.names[2]
  plot.df <- rbind(plot.df.fit1, plot.df.fit2)
  connect.df <- data.frame(x1=plot.df.fit1[,"mode"],
                           y1=plot.df.fit1[,"lambda"],
                           x2=plot.df.fit2[,"mode"],
                           y2=plot.df.fit2[,"lambda"])

  p <-  ggplot()+
    geom_point(data=plot.df, aes(x=mode, y=lambda,color=fit),alpha=0.3, size=2.5)+
    theme_ipsum()+ scale_color_manual(values=c("red","forestgreen","burlywood4"))+
    geom_segment(data=connect.df, aes(x = x1, y = y1, xend = x2, yend = y2), linetype="dotted")
  if(!is.null(boundaries)){p <- geom_vline(xintercept = boundaries, color="red", linetype="dashed")}
  if(!is.null(title)){p <- p + ggtitle(title)}
  return(p)
}
