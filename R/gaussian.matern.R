gaussian.matern <- function(range=1/10,smoothness=.5,grid.size=c(256,256),dim=2,factor=1) {
	obj <- new.env()
	obj$dim <- dim
	if(length(grid.size)==1) grid.size <- rep(grid.size,obj$dim)
	obj$grid.size <- grid.size
	obj$n1 <- grid.size[1]
	obj$factor <- factor
	obj$range <- range
	obj$smoothness <- smoothness
	class(obj) <- "gaussian.matern"
	obj
}

simulate.gaussian.matern <- function(obj) {
  library(fields)
	#Simulate a Gaussian random field with an exponential covariance function, 
	#range parameter =1./10 and the domain is [0, factor]X [0, factor] evaluating the 
	#field at a sizeGridsimul x sizeGridsimul grid. 
	#one simulates 1 realization using fctions of fields -package
	sizeGridsimul<-obj$factor*obj$n1
	grid <- list( x= seq( 0, obj$factor, length=sizeGridsimul), y= seq(0, obj$factor,length=sizeGridsimul)) 
	#XXXXXX attention en 2016: il faut utiliser theta au lieu de range et matern ds le nom de la fct
	mic<-matern.image.cov( setup=TRUE, grid=grid, theta=obj$range,smoothness=obj$smoothness)
	obj$look <- sim.rf(mic)
}

########### plot of this realization
plot.gaussian.matern <- function(obj) {
	x1 <- x2 <- seq(1, obj$n1, 1)
	trueZ <-obj$look[1:obj$n1,1:obj$n1]
	image(x1,x2,trueZ)
}
