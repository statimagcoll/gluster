#' Histogram of data with density curve
#'
#' This plot funtion can plot the histogram of the data and with the addition of density
#' curve (optional) overlay on top the histogram.
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
#' @return This function gives a plot of histogram and optional fitted curve from a gluster object.
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


#' Histogram contrast (inner function, not shown)
#'
#'
#' @param fit1 A gluster object.
#' @param fit2 Another gluster object, possibly that is a refit of fit1 with updated constraints.
#' @param marker Marker to be plotted, defaults to the first one.
#' @param subBatch subBatch to be plotted, defaults to the first one.
#' @param title Title of plot.
#' @param ... Arguments passed to plot.gluster
#' @importFrom ggplot2 aes ggplot ggtitle geom_histogram after_stat stat_function geom_vline unit annotation_custom
#' @importFrom hrbrthemes theme_ipsum
#' @importFrom gridExtra tableGrob
#' @import patchwork
#' @details Plot a histogram of a gluster model. Plots one model.
#' Takes a gluster object and plots the histogram with the fitted model and parameter values.
hist_sr_constrast_inner <- function(fit1, fit2, marker=1, subBatch=1, title=NULL, ...){
  p <- plot(x=fit1, marker=marker, subBatch=subBatch, title=title, print=FALSE, ...)
  p <- plot(x=fit2, marker=marker, subBatch=subBatch, title=title, hist=FALSE, dens=TRUE, tabl=FALSE, p=p, col='red', print=FALSE, ...)

  pars = fit2$params[[marker]][[subBatch]]
  # p <- p + ggplot2::stat_function(fun = fun1, n = 101, alpha=0.3, color="red", xlim=p.range) +
  #   ggplot2::stat_function(fun = fun2, n = 101, alpha=0.3, color="red", xlim=p.range) +
  #   ggplot2::geom_vline(xintercept=afterthresh, linetype="dashed", color = "red")
  # add.table argument currently on hold
  # if(add.table){
  #   pars <- rbind(pars, mode=c(NA, (pars[2,2]-1)*pars[3,2], (pars[2,3]-1)*pars[3,3]))
  #   p <- p + annotation_custom(tableGrob(pars), xmin=35, xmax=50, ymin=-2.5, ymax=-1)
  # }
  print(p)
}

#' Contrastive density curve on histogram
#'
#' This function can compare two fitted model of a marker on a slide, from gluster/groupGluster object.
#' @param fit1 A gluster or groupGluster object.
#' @param fit2 Another gluster or groupGluster object, possibly that is a refit of fit1 with updated constraints. Has to be of the same class as fit1.
#' @param marker Marker to be plotted, defaults to the first one.
#' @param subBatch subBatch to be plotted, defaults to the first one.
#' @param slide Slide to be plotted, defaults to the first one.
#' @param title Title of plot.
#' @param ... Arguments passed to plot.gluster
#' @importFrom ggplot2 aes ggplot ggtitle geom_histogram after_stat stat_function geom_vline unit annotation_custom
#' @importFrom hrbrthemes theme_ipsum
#' @importFrom gridExtra tableGrob
#' @import patchwork
#' @export
#' @details Plot histogram constrast of two gluster or groupGluster model. Plots one model.
#' Takes a gluster object and plots the histogram with the fitted model and parameter values.
hist_sr_constrast <- function(fit1, fit2, marker=1, slide=1, subBatch=1, title=NULL, ...){
  if(any(class(fit1)!=class(fit2))){
    print("Two objects are not of same class! Please change input.")
    } else if (class(fit1)[1]=="gluster"){
      marker.name <- names(fit1[[1]])[marker]
      hist_sr_constrast_inner(fit1, fit2, marker=marker, subBatch=subBatch, title=title, ...)
    } else if (class(fit1)[1]=="groupGluster"){
      hist_sr_constrast_inner(fit1[[slide]],fit2[[slide]], marker = marker, subBatch=subBatch, title = title )
    } else {print("Invalid input class! Please change input.")}
}


#' Compare binary classification result using Cohen's Kappa
#'
#' This function can compare the agreement of classification result against a gold/silver standard.
#' The result can be returned as in either a data frame of Cohen's Kappa or a side-by-side box plot.
#' @param ... Matrices to compare to standard labels. Rows are cells, columns are markers.
#' @param standard A matrix indicating gold/silver standard positive cells.
#' @param batch Split kappa computation by this variable. Similar to slide_id.
#' @importFrom stats dgamma
#' @importFrom ggplot2 aes ggplot ggtitle geom_histogram after_stat stat_function geom_vline unit annotation_custom geom_boxplot coord_flip
#' @importFrom rlang UQ
#' @importFrom hrbrthemes theme_ipsum
#' @importFrom gridExtra tableGrob
#' @importFrom stats na.omit density
#' @importFrom reshape2 melt
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
  method.names <- names(list(...))
  if(is.null(method.names)) method.names=paste0('method', 1:length(cohen.kappa.df))
  names(cohen.kappa.df)=method.names
  cohen.kappa.df = lapply(names(cohen.kappa.df), function(xname){res=cohen.kappa.df[[xname]]; names(res) = c(markers, 'batch'); res$method=xname; res})
  cohen.kappa.df = do.call(rbind,  cohen.kappa.df)
  colnames(cohen.kappa.df)[1:length(markers)] <- markers
  cohen.kappa.plot <- melt(cohen.kappa.df, id.vars = c("batch", "method"), variable.name = "marker",
                           value.name = "Cohen.Kappa")


  p <- ggplot() + theme_ipsum(plot_title_size = 10,base_size = 8)  + coord_flip()+
     geom_boxplot(data=cohen.kappa.plot, aes(x=marker, y=Cohen.Kappa, fill=method), width=0.7, alpha=0.3)+ggtitle("Cohen's Kappa compared\nto Silver Standard")+xlab("Marker")
  print(p)
  return(cohen.kappa.df)
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
