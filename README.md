Guide for using Gluster in semi-automated cell gating
================
Ruby Xiong and Simon Vandekar
2023-03-27

# gluster

Gamma Mixture Models for Multiplexed Immunofluorescence Imaging Data

Install Package

``` r
# if( !require(devtools) ) install.packages("devtools")
# devtools::install_github( "statimagcoll/gluster" )
library(gluster)
```

## Loading mIF data from `VectraPolarisData`

Install instructions can be found at
<https://github.com/julia-wrobel/VectraPolarisData>. After installation,
run the following:

``` r
library(VectraPolarisData)
spe_lung <- HumanLungCancerV3()
# get the cell level imaging data with spatial coordinates
cell = cbind(as.data.frame(colData(spe_lung)), as.data.frame(spatialCoords(spe_lung)) )
```

## Prepare the data for fitting gamma mixture model

Before fitting the model, the expression levels are normalized to reduce
batch effect and variation among slides. Quantile boundaries and fixed
value boundaries are picked based on prior knowledge.

``` r
allMarkers = grep('^entire_cell_.*total$', names(cell), value=TRUE)
# normalize the data 
nzNormedMarkers = paste0('nzNorm_',allMarkers)
cell = do.call(rbind, lapply(split(cell, cell$slide_id), function(df){df[,nzNormedMarkers] = log10(1+sweep(df[,allMarkers], 2, colMeans(replace(df[,allMarkers], df[,allMarkers]==0, NA), na.rm = TRUE ), FUN = '/' ) ); df }) )
# split the data frame
cells = split(cell, cell$slide_id)
# set up boundary for all markers
quantileBoundaries = list(matrix(c(0,.8, .9, 1), byrow = TRUE, nrow=2),
                          matrix(c(0,.7, .8, 1), byrow = TRUE, nrow=2),
                          matrix(c(0,.6, .8, 1), byrow = TRUE, nrow=2),
                          matrix(c(0,.6, .8, 1), byrow = TRUE, nrow=2),
                          matrix(c(0,.7, .8, 1), byrow = TRUE, nrow=2),
                          NULL,
                          NULL,
                          NULL)
boundaries = list(matrix(c(0,.3, .5, 1), byrow = TRUE, nrow=2),
                  NULL,
                  NULL,
                  NULL,
                  matrix(c(0,.3, .45, Inf), byrow = TRUE, nrow=2),
                  matrix(c(0,.3, .5, Inf), byrow = TRUE, nrow=2),
                  matrix(c(0,.4, .45, Inf), byrow = TRUE, nrow=2),
                  NULL)
names(quantileBoundaries) = names(boundaries) = nzNormedMarkers
```

## Running gluster on a single slide

Here the fitting is run both with and without constraints. Fitting run
without constraints are susceptible to statistically optimal results
(local maximum) that are not reasonable biologically.

``` r
unconstrGMMfit = gluster(cells[[1]][,nzNormedMarkers])
constrGMMfit = gluster(cells[[1]][,nzNormedMarkers], boundaryMarkers = boundaries,
                 qboundaryMarkers=quantileBoundaries)
```

### Quality checking model fit

``` r
# check convergence
convCheck(constrGMMfit)
# visualize constrained fits with quantile drawn on histogram
plots_hist <- sapply(nzNormedMarkers, function(x) plot(constrGMMfit, marker = x, boundary = quantile(constrGMMfit$expressionX[,x], probs=quantileBoundaries[[x]][2,1]), title=x))
```

![](README_files/figure-gfm/visualization-1.png)<!-- -->![](README_files/figure-gfm/visualization-2.png)<!-- -->![](README_files/figure-gfm/visualization-3.png)<!-- -->![](README_files/figure-gfm/visualization-4.png)<!-- -->![](README_files/figure-gfm/visualization-5.png)<!-- -->![](README_files/figure-gfm/visualization-6.png)<!-- -->![](README_files/figure-gfm/visualization-7.png)<!-- -->![](README_files/figure-gfm/visualization-8.png)<!-- -->

``` r
# comparing constrained and unconstrained fits
contrast_hist = sapply(nzNormedMarkers, function(x) hist_sr_constrast(constrGMMfit, unconstrGMMfit, marker=x, title=x) )
```

![](README_files/figure-gfm/visualization-9.png)<!-- -->![](README_files/figure-gfm/visualization-10.png)<!-- -->![](README_files/figure-gfm/visualization-11.png)<!-- -->![](README_files/figure-gfm/visualization-12.png)<!-- -->![](README_files/figure-gfm/visualization-13.png)<!-- -->![](README_files/figure-gfm/visualization-14.png)<!-- -->![](README_files/figure-gfm/visualization-15.png)<!-- -->![](README_files/figure-gfm/visualization-16.png)<!-- -->

## Running with sub-batches

Sub-batches can be used if there are multiple regions per slide to
assess the fit within each region. The model is fit at the slide level
and then model fit is assessed by batch.

``` r
# fitting on the same slide and getting model fit results by batch
batchConstrGMMfit = gluster(cells[[1]][,nzNormedMarkers], boundaryMarkers = boundaries,
                 qboundaryMarkers=quantileBoundaries, subBatch = cells[[1]]$sample_id)
# visualize by subBatch for the first markers
trash = sapply(unique(cells[[1]]$sample_id), function(x) plot(batchConstrGMMfit, subBatch = x, title=x) )
```

![](README_files/figure-gfm/subBatch-1.png)<!-- -->![](README_files/figure-gfm/subBatch-2.png)<!-- -->![](README_files/figure-gfm/subBatch-3.png)<!-- -->![](README_files/figure-gfm/subBatch-4.png)<!-- -->![](README_files/figure-gfm/subBatch-5.png)<!-- -->

## Running gluster on multiple slides

We can use mclapply from the parallel package to run gluster on multiple
slides. The example line below runs on a subset of 3 slides from the
`cell` data frame.

Alternatively, use groupGluster:

``` r
# subset of cell, first 3 slides
slide3 <- names(cells)[1:3]
cells3 <- cell[cell$slide_id %in% slide3,]
# fit cfGMM with groupGluster
constrCfGMMbunch <- groupGluster(cells3[,nzNormedMarkers], slide = cells3$slide_id,
                                 boundaryMarkers = boundaries,qboundaryMarkers=quantileBoundaries, n.cores = 5)
```

### Quality check on model fitting

Check for convergence:

``` r
convCheck(constrCfGMMbunch)
```

    ## [1] "All converged."

Before plotting and looking over the (# slides $\times$ \# markers)
histograms, it might be easier to identify the problematic ones through
diagnostic plot first, then look at all histograms:

``` r
for(i in nzNormedMarkers){
  plot(constrCfGMMbunch, diagnostic=TRUE,  interactive=FALSE, histogram=FALSE, marker=i, title=i)
}
```

![](README_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-4-2.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-4-3.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-4-4.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-4-5.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-4-6.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-4-7.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-4-8.png)<!-- -->

``` r
slides.fit <- names(constrCfGMMbunch)
for(i in slides.fit){
  cat('\n#### ', i, " \n")
  for(j in nzNormedMarkers){
    plot(constrCfGMMbunch, marker=j, slide=i, title=paste0(i,"\n",j), diagnostic=FALSE, histogram=TRUE, 
         boundary = quantile(constrCfGMMbunch[[i]]$expressionX[,j], probs=quantileBoundaries[[j]][2,1]))
  }
  cat('\n\n')
}
```

#### \#01 0-889-121

![](README_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-5-2.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-5-3.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-5-4.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-5-5.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-5-6.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-5-7.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-5-8.png)<!-- -->

#### \#02 1-037-393

![](README_files/figure-gfm/unnamed-chunk-5-9.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-5-10.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-5-11.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-5-12.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-5-13.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-5-14.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-5-15.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-5-16.png)<!-- -->

#### \#03 2-080-378

![](README_files/figure-gfm/unnamed-chunk-5-17.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-5-18.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-5-19.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-5-20.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-5-21.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-5-22.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-5-23.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-5-24.png)<!-- -->

### Quality check on model fitting: compare

Again, we can compare it fitting with no boundaries.

``` r
# same as above , but fitting without boundaries
CfGMMbunch = groupGluster(cells3[,nzNormedMarkers], slide = cells3$slide_id,n.cores = 5)
convCheck(CfGMMbunch)
```

\[1\] “All converged.”

``` r
for(i in slides.fit){
  cat('\n#### ', i, " \n")
  for(j in nzNormedMarkers){
    hist_sr_constrast(constrCfGMMbunch, CfGMMbunch, marker=j, slide=i, title=paste0(i,"\n",j), diagnostic=FALSE, histogram=TRUE,boundary = quantile(constrGMMbunch[[i]]$expressionX[,j], probs=quantileBoundaries[[j]][2,1]))
  }
  cat('\n\n')
}
```

#### \#01 0-889-121

![](README_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-6-2.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-6-3.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-6-4.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-6-5.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-6-6.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-6-7.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-6-8.png)<!-- -->

#### \#02 1-037-393

![](README_files/figure-gfm/unnamed-chunk-6-9.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-6-10.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-6-11.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-6-12.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-6-13.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-6-14.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-6-15.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-6-16.png)<!-- -->

#### \#03 2-080-378

![](README_files/figure-gfm/unnamed-chunk-6-17.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-6-18.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-6-19.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-6-20.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-6-21.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-6-22.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-6-23.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-6-24.png)<!-- -->

### Constrastive diagnostic plot

We can overlay the diagnostic plot to further check what are the
difference between the two fitted models.

``` r
for(i in nzNormedMarkers){
  print(diag_contrast(constrCfGMMbunch, CfGMMbunch, marker=i, title=i, fit.names = c("Constrained", "Un-constrained")))
}
```

![](README_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-7-2.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-7-3.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-7-4.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-7-5.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-7-6.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-7-7.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-7-8.png)<!-- -->

### Cohen’s Kappa

To compare agreement between two binomial variables, we can use cohen’s
kappa. Here is the usage of the function with simulated data:

``` r
test1 <- matrix(sample(c(0,1),400, replace = TRUE), nrow=200)
test2 <- matrix(sample(c(0,1),400, replace = TRUE), nrow=200)
standard.mat <- matrix(sample(c(0,1),400, replace = TRUE), nrow=200)
colnames(test1) <- colnames(test2) <- colnames(standard.mat) <- c("a", "b")
kappaGroupGluster(test1, test2, standard=standard.mat, batch=rep(1:5, each=10))
```

![](README_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

    ##              a           b batch  method
    ## 1  -0.25000000 -0.19047619     1 method1
    ## 2   0.04040404  0.26108374     2 method1
    ## 3   0.01234568 -0.18518519     3 method1
    ## 4  -0.10000000  0.15422886     4 method1
    ## 5   0.05472637 -0.15662651     5 method1
    ## 11 -0.10000000  0.18987342     1 method2
    ## 21 -0.06060606  0.10000000     2 method2
    ## 31  0.13253012 -0.25000000     3 method2
    ## 41 -0.26903553  0.19799499     4 method2
    ## 51 -0.22448980  0.07317073     5 method2
