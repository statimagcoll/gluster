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
#' @param print logical whether to print the plot.
#' @param title Title of plot.
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
plot.gluster <- function(x, marker=1, subBatch=1, zero_include=FALSE, breaks=40, title="histogram of x", boundary=NULL, hist=TRUE, dens=TRUE, tabl=FALSE, col='forestgreen', p=NULL, print=TRUE, ...){
  if(!all(is.na(x[["expressionX"]][[marker]]))){
    if(is.numeric(marker)){ marker=colnames(x[["expressionX"]])[marker] }
    if(is.numeric(subBatch) & !is.null(names(x$params[[marker]]))[1]){ subBatch = names(x$params[[marker]])[subBatch] }
    expr = x[["expressionX"]][marker]
    pars <- x$params[[marker]][[subBatch]]
    if(!zero_include){
      x[["subBatch"]] = x[["subBatch"]][expr[,marker]>0]
      expr = expr[expr>0, , drop=FALSE]
      pars[1,2:3] = pars[1, 2:3]/(1-pars[1,1])
    }
    p.range <- c(0, max(expr, na.rm=TRUE))
    xs <- seq(0, p.range[2], length.out=100)

    fun1 <- function(xs){pars[1,2] * dgamma(xs, shape=pars[2,2],scale=pars[3,2])}
    fun2 <- function(xs){pars[1,3]  * dgamma(xs, shape=pars[2,3],scale=pars[3,3])}
    plot.x <- na.omit(expr[ which(as.character(x[["subBatch"]])==as.character(subBatch)), ,drop=FALSE ])

    if(hist){
      p <- ggplot2::ggplot(plot.x, ggplot2::aes(UQ(as.name(marker)) ) ) +
        ggplot2::ylab("")+hrbrthemes::theme_ipsum()+ggplot2::ggtitle(title)+
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
      partab = rbind(x$params[[marker]][[subBatch]], NA)
      rownames(partab)[4] = 'mode'
      partab['mode', c(2,3)] =  (partab['alpha', c(2,3)]-1) * partab['beta', c(2,3)]
      p2 <- gridExtra::tableGrob(round(as.data.frame(partab),3))
      # Set widths/heights to 'fill whatever space I have'
      p2$widths <- unit(rep(1, ncol(p2)), "null")
      p2$heights <- unit(rep(1, nrow(p2)), "null")
      # Format table as plot
      p3 <- ggplot() +
        annotation_custom(p2)
      # Patchwork magic
      p <- p + p3 + plot_layout(ncol = 2)
    }
    if(print) print(p) else p
  }
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
  p <- plot.gluster(x=fit1, marker=marker, subBatch=subBatch, title=title, print=FALSE, ...)
  p <- plot.gluster(x=fit2, marker=marker, subBatch=subBatch, title=title, hist=FALSE, dens=TRUE, tabl=FALSE, p=p, col='red', print=FALSE, ...)

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



#' Fits a gamma mixture model to aid clustering of mIF imaging data
#'
#'
#' @param ... Matrices to compare to standard labels.
#' @param standard A matrix indicating gold/silver standard positive cells.
#' @param batch Split kappa computation by this variable.
#' @importFrom stats dgamma
#' @importFrom ggplot2 aes ggplot ggtitle geom_histogram after_stat stat_function geom_vline unit annotation_custom
#' @importFrom rlang UQ
#' @importFrom hrbrthemes theme_ipsum
#' @importFrom gridExtra tableGrob
#' @importFrom stats na.omit density
#' @export
#' @details Plot a histogram of a gluster model. Plots one model.
#' Takes a gluster object and plots the histogram with the fitted model and parameter values.
kappaGroupGluster = function(..., standard, batch){
  standard = as.data.frame(standard)
  # cohen's kappa
  x = lapply(list(...), as.data.frame)
  # check all dimensions match
  markers = colnames(standard)

  # Compute kappa for each method, slide, and marker
  cohen.kappa.df = lapply(x, function(xmat, standard, batch){
    result2 = as.data.frame(t(mapply(function(sxmat, sstandard){
      result = mapply(function(xv, yv){
        po=table(xv, yv)
        po = po/sum(po)
        pe=sum(rowSums(po) * colSums(po))
        po=sum(po[c(1,4)])
        (po-pe)/(1-pe)  }, xv=sxmat, yv=sstandard)
      #result$method = if(!is.null(rownames(result))) rownames(result) else
      return(result)
    }, sxmat=split(xmat, batch), sstandard=split(standard, batch) ) ))
    result2$batch = rownames(result2)
    result2
  }, standard=standard, batch=batch)
  if(is.null(names(cohen.kappa.df))) names(cohen.kappa.df)=paste0('method', 1:length(cohen.kappa.df))
  cohen.kappa.df = lapply(names(cohen.kappa.df), function(xname){res=cohen.kappa.df[[xname]]; names(res) = c(markers, 'batch'); res$method=xname; res})
  cohen.kappa.df = do.call(rbind,  cohen.kappa.df)

  # reshape to long format?


  # data_breaks <- data.frame(start = c(0, 0.1, 0.21, 0.41, 0.61, 0.81),  # Create data with breaks
  #                           end = c(0.1, 0.21, 0.41, 0.61, 0.81, 1),
  #                           colors = factor(c("no agreement","slight agreement", "fair agreement","moderate agreement",
  #                                             "substantial agreement","near-perfect agreement")))
  # library(metR)
  # library(magrittr)
  # ggplot() + theme_ipsum(plot_title_size = 10,base_size = 8)  + coord_flip()+
  #   geom_boxplot(data=cohen.kappa.df, aes(x=marker, y=cohen.kappa, fill=method), width=0.7, alpha=0.3)+ggtitle("Cohen's Kappa compared\nto Silver Standard")+xlab("Marker")
  # geom_rect(data = data_breaks,
  #           aes(ymin = start,ymax = end,xmin = 0,xmax = 4,fill = colors), alpha = 0.5) +scale_fill_grey()
  return(cohen.kappa.df)
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
#' @param diag Logical indicating whether to plot interactive diagnostic plot.
#' @param dens Logical indicating whether to draw parametric density over histogram.
#' @param tabl Logical indicating whether to add table of parameter values.
#' @param col Color of density.
#' @param p ggplot object to add layers too. If hist=FALSE this must be specified.
#' @param print logical whether to print the plot.
#' @param title Title of plot.
#' @param ... Arguments passed to cfGMM function
#' @importFrom stats dgamma
#' @importFrom ggplot2 aes ggplot ggtitle geom_histogram after_stat stat_function geom_vline unit annotation_custom
#' @importFrom rlang UQ
#' @importFrom hrbrthemes theme_ipsum
#' @importFrom gridExtra tableGrob
#' @importFrom stats na.omit density
#' @importFrom plotly plot_ly layout
#' @import patchwork
#' @export
#' @details Plot a histogram of a gluster model. Plots one model.
#' Takes a gluster object and plots the histogram with the fitted model and parameter values.
plot.groupGluster <- function(x, marker=1, slide=1, subBatch=1, zero_include=FALSE, breaks=40, title="histogram of x", boundary=NULL, hist=TRUE, dens=TRUE, tabl=FALSE, diag=TRUE, col='forestgreen', p=NULL, print=TRUE, ...){
  x[[slide]] <- x.slide
  plot.gluster(x.slide, slide=slide, subBatch=subBatch, zero_include=FALSE, breaks=40, title=title, boundary=boundary, hist=hist, dens=dens, tabl=tabl, col=col, p=p, print=print)
  if(diag){
    pars <- lapply(x, function(x){x[["params"]]})
    pars <- lapply(pars, function(x){x[[marker]]})
    pars <- do.call(c, pars)
    pars <- lapply(pars, get.mode)
    slide.names <- names(pars)

    plot.df <- data.frame(lapply(pars, function(par.mat){par.mat[c(1,4),3]}))
    plot.df1 <- data.frame(t(plot.df))
    plot.df1$slide_id <- slide.names
    fig1 <- plot_ly( plot.df1, x = ~modes, y = ~lambda,
                     marker = list(size = 10,color = 'rgba(255, 182, 193, .9)',
                                   line = list(color = 'rgba(152, 0, 0, .8)',width = 2)),
                     text = ~slide_id, type = "scatter", mode = 'markers')
    fig1 <- fig1 %>% layout(title = marker,
                            yaxis = list(zeroline = FALSE),
                            xaxis = list(zeroline = FALSE))
    return(fig1)
    }
}

#' Get mode from parameters
#'
#'
#' @param par.mat parameter matrix from gluster or groupGluster object.
#' @details Get mode from a parameter matrix
get.mode <- function(par.mat){
  modes <- c(NA, (par.mat[2,2]-1)*par.mat[3,2], (par.mat[2,3]-1)*par.mat[3,3])
  par.mat <- rbind(par.mat, modes)
  return(par.mat)
}
