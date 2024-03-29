% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/groupGluster.R
\name{plotGroupGluster}
\alias{plotGroupGluster}
\title{Inner: Plot function of fitted model in groupGluster object}
\usage{
plotGroupGluster(
  x,
  marker = 1,
  slide = 1,
  component = 2,
  diagnostic = TRUE,
  interactive = FALSE,
  histogram = FALSE,
  title = NULL,
  boundary = NULL,
  color = "grey",
  print = TRUE,
  tabl = FALSE,
  ...
)
}
\arguments{
\item{x}{A groupGluster object}

\item{marker}{Select which markers to plot. Can be a vector of character or integer. If not specified, will print all.}

\item{slide}{Select which slides to plot. Can be a vector of character or integer. If not specified, will print all.}

\item{component}{Integer specifying which component to plot, 1 is unexpressed nonzero cells, 2 is expressed cells.}

\item{diagnostic}{logical indicating whether to create the diagnostic plot. Default value is TRUE.}

\item{interactive}{logical indicating whether diagnostic plot should be interactive.}

\item{histogram}{logical indicating whether to create the slide histograms.}

\item{title}{Title for the plot. Default is the marker name.}

\item{boundary}{Boundary (vertial dashed line) to be plotted on the histogram.}

\item{color}{color for points.}

\item{print}{logical whether to display the plot. Default value TRUE.}

\item{...}{Arguments passed to XX}
}
\description{
Takes a slide and marker input, and returns the corresponding histogram of expression value,
along with fitted density curve.
}
\details{
Various diagnostic and QC plots for groupGluster fits.
}
