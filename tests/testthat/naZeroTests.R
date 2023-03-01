devtools::load_all()
x = rgamma(1000, shape=1)
x[1:10 ] = NA
x[ 30:40] = 0
subClust = rep(c(1,2,3), c(400, 100, 500))
test = glusterX(x, constraints=matrix(c(-Inf, 0, 0.5, Inf), byrow = TRUE, nrow=2), subBatch = subClust)
#test = glusterX(x, constraints=matrix(c(-Inf, 0, 0.5, Inf), byrow = TRUE, nrow=2))
#test = glusterX(x, constraints=matrix(c(-Inf, 0, 0.5, Inf), byrow = TRUE, nrow=2), subBatch = subClust, max.iter=5)
#test = glusterX(x, constraints=matrix(c(-Inf, 0, 0.5, Inf), byrow = TRUE, nrow=2), max.iter=5)

# test of when there are too many zeros
xzeros = x; xzeros[1:900] = 0
test = glusterX(xzeros, constraints=matrix(c(-Inf, 0, 0.5, Inf), byrow = TRUE, nrow=2), subBatch = subClust)

# test of when there are all NAs
test = glusterX(rep(NA, 1000), constraints=matrix(c(-Inf, 0, 0.5, Inf), byrow = TRUE, nrow=2), subBatch = subClust)


# test of when there are fewer than 200 values
test = glusterX(x[1:199], constraints=matrix(c(-Inf, 0, 0.5, Inf), byrow = TRUE, nrow=2), subBatch = subClust[1:199])

y = rgamma(1000, shape=1)
expMat = cbind(x, y)
test = gluster(expressionMarkers = expMat, boundaryMarkers = rep(list(matrix(c(-Inf, 0, 0.5, Inf), byrow = TRUE, nrow=2)), 2), subBatch = subClust )
test = gluster(expressionMarkers = expMat,
               qboundaryMarkers = rep(list(matrix(c(0, 0.25, 0.5, 1), byrow = TRUE, nrow=2)), 2) )
test = gluster(expressionMarkers = expMat,
               boundaryMarkers = rep(list(matrix(c(-Inf, 0, 0.5, Inf), byrow = TRUE, nrow=2)), 2),
               qboundaryMarkers = rep(list(matrix(c(0, 0.25, 0.5, 1), byrow = TRUE, nrow=2)), 2) )

plot.gluster(test)
#debug(hist_sr_constrast)#; debug(plot.gluster)
hist_sr_constrast(test, test, add.table=TRUE)
