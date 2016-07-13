belong.to.disk  <- function(x, center, halfDiameter.power2) t(dif <- x - center) %*% dif  < halfDiameter.power2

expand.to.fullGrid  <- function(grid.size=256, z, missing.sites) {
	## it should be checked that length(z)+ (nb of true booleans) = n1*n1
  if(length(grid.size)==1) grid.size <- rep(grid.size,2)
	zz <- rep(0,prod(grid.size))
	zz[!missing.sites]  <- z
	matrix(zz, nrow=grid.size[1],ncol=grid.size[2])
}
