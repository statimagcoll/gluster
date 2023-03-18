library(VectraPolarisData)
spe_lung <- HumanLungCancerV3()

# get the cell level imaging data with spatial coordinates
cell = cbind(as.data.frame(colData(spe_lung)), as.data.frame(spatialCoords(spe_lung)) )

#devtools::install_github('statimagcoll/gluster')
allMarkers = grep('^entire_cell_.*total$', names(cell), value=TRUE)
nzNormedMarkers = paste0('nzNorm_',allMarkers)
cell = do.call(rbind, lapply(split(cell, cell$slide_id), function(df){df[,nzNormedMarkers] = log10(1+sweep(df[,allMarkers], 2, colMeans(df[,allMarkers], na.rm = TRUE ), FUN = '/' ) ); df }) )

boundry <- matrix(c(0,.3, .5, 1), byrow = TRUE, nrow=2)

fitx <- cell$nzNorm_entire_cell_cd19_opal_650_total[cell$slide_id=="#01 0-889-121"]
fit.cfgmm <- cfGMM::cfGMM(fitx, k=2, constraint = boundry)
fit.cfgmm$gamma.pars
