get.mode <- function(par.mat){
  modes <- c(NA, (par.mat[2,2]-1)*par.mat[3,2], (par.mat[2,3]-1)*par.mat[3,3])
  par.mat <- rbind(par.mat, modes)
  return(par.mat)
}


#' Fits a gamma mixture model to aid clustering of mIF imaging data
#'
#'
#' @param x A gluster object.
#' @param marker Marker to be plotted, defaults to the first one.
#' @param subBatch subBatch to be plotted, defaults to the first one.
#' @param zero_include Logical whether the histogram should include zeros. Default is FALSE.
#' @param breaks number of breaks for histogram
#' @param boundary Vertical line drawn on plot to indicate boundary constraint.
#' @param hist Logical indicating whether to plot histogram and base layer of plot.
#' @param dens Logical indicating whether to draw parametric density over histogram.
#' @param tabl Logical indicating whether to add table of parameter values.
#' @param col Color of density.
#' @param p ggplot object to add layers too. If hist=FALSE this must be specified.
#' @param main Title of plot.
#' @param ... Arguments passed to cfGMM function
#' @importFrom stats dgamma
#' @importFrom ggplot2 aes ggplot ggtitle geom_histogram after_stat stat_function geom_vline unit annotation_custom
#' @importFrom rlang UQ
#' @importFrom hrbrthemes theme_ipsum
#' @importFrom gridExtra tableGrob
#' @importFrom stats na.omit density
#' @import patchwork
#' @export
#' @details Plot a histogram of a gluster model. Plots one model.
#' Takes a gluster object and plots the histogram with the fitted model and parameter values.
plot.gluster <- function(x, marker=1, subBatch=1, zero_include=FALSE, breaks=40, main="histogram of x", boundary=NULL, hist=TRUE, dens=TRUE, tabl=FALSE, col='forestgreen', p=NULL, ...){
  if(is.numeric(marker)){ marker=colnames(x[["expressionX"]])[marker] }
  p.range <- c(0, max(x[["expressionX"]][[marker]], na.rm=TRUE))
  xs <- seq(0, p.range[2], length.out=100)
  pars <- x$params[[marker]][[subBatch]]

  fun1 <- function(xs){pars[1,2] * dgamma(xs, shape=pars[2,2],scale=pars[3,2])}
  fun2 <- function(xs){pars[1,3]  * dgamma(xs, shape=pars[2,3],scale=pars[3,3])}
  plot.x <- na.omit(x[["expressionX"]][marker])

  if(hist){
    p <- ggplot2::ggplot(plot.x, ggplot2::aes(UQ(as.name(marker)) ) ) +
      ggplot2::ylab("")+hrbrthemes::theme_ipsum()+ggplot2::ggtitle(main)+
      ggplot2::geom_histogram(ggplot2::aes(y=after_stat(density)),bins = breaks, alpha=0.2)
  }
  if(dens){
     p<- p +
      ggplot2::stat_function(fun = fun1, n = 101, alpha=0.3, color=col, xlim = p.range) +
      ggplot2::stat_function(fun = fun2, n = 101, alpha=0.3, color=col, xlim=p.range)
  }
  if(!is.null(boundary)){
    p <- p + ggplot2::geom_vline(xintercept=boundary, linetype="dashed", color = col)
  }
  if(tabl){
    p2 <- gridExtra::tableGrob(round(as.data.frame(pars),3))
    # Set widths/heights to 'fill whatever space I have'
    p2$widths <- unit(rep(1, ncol(p2)), "null")
    p2$heights <- unit(rep(1, nrow(p2)), "null")
    # Format table as plot
    p3 <- ggplot() +
      annotation_custom(p2)
    # Patchwork magic
    p <- p + p3 + plot_layout(ncol = 2)
  }
  print(p)
}


#' Fits a gamma mixture model to aid clustering of mIF imaging data
#'
#'
#' @param fit1 A gluster object.
#' @param fit2 Another gluster object, possibly that is a refit of fit1 with updated constraints.
#' @param marker Marker to be plotted, defaults to the first one.
#' @param subBatch subBatch to be plotted, defaults to the first one.
#' @param title Title of plot.
#' @param add.table logical indicating whether to add table of fit2 parameter values.
#' @param ... Arguments passed to plot.gluster
#' @importFrom ggplot2 aes ggplot ggtitle geom_histogram after_stat stat_function geom_vline unit annotation_custom
#' @importFrom hrbrthemes theme_ipsum
#' @importFrom gridExtra tableGrob
#' @import patchwork
#' @export
#' @details Plot a histogram of a gluster model. Plots one model.
#' Takes a gluster object and plots the histogram with the fitted model and parameter values.
hist_sr_constrast <- function(fit1, fit2, marker=1, subBatch=1, title=NULL, add.table=FALSE, ...){
  p <- plot.gluster(x=fit1, marker=marker, subBatch=subBatch, main=title, ...)
  p <- plot.gluster(x=fit2, marker=marker, subBatch=subBatch, main=title, hist=FALSE, dens=TRUE, tabl=FALSE, p=p, col='red', ...)

  pars = fit2$params[[marker]][[subBatch]]
  # p <- p + ggplot2::stat_function(fun = fun1, n = 101, alpha=0.3, color="red", xlim=p.range) +
  #   ggplot2::stat_function(fun = fun2, n = 101, alpha=0.3, color="red", xlim=p.range) +
  #   ggplot2::geom_vline(xintercept=afterthresh, linetype="dashed", color = "red")
  if(add.table){
    pars <- rbind(pars, mode=c(NA, (pars[2,2]-1)*pars[3,2], (pars[2,3]-1)*pars[3,3]))
    p <- p + annotation_custom(tableGrob(pars), xmin=35, xmax=50, ymin=-2.5, ymax=-1)
  }
  print(p)
}

