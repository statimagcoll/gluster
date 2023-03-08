# gluster
Gamma Mixture Models for Multiplexed Immunofluorescence Imaging Data

Install Package

```
if( !require(devtools) ) install.packages("devtools")
devtools::install_github( "statimagcoll/gluster" )
```

## Loading mIF data from `VectraPolarisData`

Install instructions can be found at https://github.com/julia-wrobel/VectraPolarisData

```{r, message=FALSE}
library(VectraPolarisData)
spe_lung <- HumanLungCancerV3()
# get the cell level imaging data with spatial coordinates
cell = cbind(as.data.frame(colData(spe_lung)), as.data.frame(spatialCoords(spe_lung)) )
# get the clinical data for each sample
#cd = as.data.frame(metadata(spe_lung)$clinical_data)
```

## Running gluster on a single slide

```{r}
library(gluster)
#devtools::install_github('statimagcoll/gluster')
allMarkers = grep('^entire_cell_.*total$', names(cell), value=TRUE)
# normalize the data 
nzNormedMarkers = paste0('nzNorm_',allMarkers)
cell = do.call(rbind, lapply(split(cell, cell$slide_id), function(df){df[,nzNormedMarkers] = log10(1+sweep(df[,allMarkers], 2, colMeans(replace(df[,allMarkers], df[,allMarkers]==0, NA), na.rm = TRUE ), FUN = '/' ) ); df }) )
# split the data frame
cells = split(cell, cell$slide_id)
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
unconstrGMMfit = gluster(cells[[1]][,nzNormedMarkers])
constrGMMfit = gluster(cells[[1]][,nzNormedMarkers], boundaryMarkers = boundaries,
                 qboundaryMarkers=quantileBoundaries)
```

### Quality checking model fit

```{r, visualization}
# check convergence
sapply(constrGMMfit$fit, function (x) x$convergence)
# visualize constrained fits with quantile drawn on histogram
trash = sapply(nzNormedMarkers, function(x) plot.gluster(constrGMMfit, marker = x, boundary = quantile(constrGMMfit$expressionX[,x], probs=quantileBoundaries[[x]][2,1], title=x) ) )
# comparing constrained and unconstrained fits
trash = sapply(nzNormedMarkers, function(x) hist_sr_constrast(constrGMMfit, unconstrGMMfit, marker=x, title=x) )
```

## Running with sub-batches

Sub-batches can be used if there are multiple regions per slide to assess the fit within each region. The model is fit at the slide level and then model fit is assessed by batch.

```{r subBatch}
# fitting on the same slide and getting model fit results by batch
batchConstrGMMfit = gluster(cells[[1]][,nzNormedMarkers], boundaryMarkers = boundaries,
                 qboundaryMarkers=quantileBoundaries, subBatch = cells[[1]]$sample_id)
# visualize by subBatch for the first markers
trash = sapply(unique(cells[[1]]$sample_id), function(x) plot.gluster(batchConstrGMMfit, subBatch = x, title=x) )
```

## Running gluster on multiple slides

We can use mclapply from the parallel package to run gluster on multiple slides.

```{r, batchFit}
library(parallel)
library(reactable)
constrCfGMMbunch = mclapply(cells[1:5], function(x) gluster(x[,nzNormedMarkers], boundaryMarkers = boundaries, qboundaryMarkers=quantileBoundaries),  mc.cores = 5)
trash = sapply(names(constrCfGMMbunch), function(x) plot(constrCfGMMbunch[[x]], marker = 1, boundary = quantile(constrGMMfit$expressionX[,1], probs=quantileBoundaries[[1]][2,1]), title = x ) )
reactable(sapply(constrCfGMMbunch, function(x) sapply(x$fit, function(y){ y$convergence})))
```
