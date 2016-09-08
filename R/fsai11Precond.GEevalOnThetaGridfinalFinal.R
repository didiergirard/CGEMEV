

fsai11Precond.GEevalOnThetaGrid <- function(z,candidateThetas1DGrid,nu=0.5,grid.domain,tolPGC=1e-03) {
  n1 <- grid.domain$n1
  sparseG <- grid.domain$sparseG
  transOfSparseG <- t(sparseG)

  GEvaluesOnthegrid<- matrix(NA,length(candidateThetas1DGrid))
  #trace.term<- matrix(NA,length(candidateThetas1DGrid))
  niter <- matrix(NA,length(candidateThetas1DGrid),2)
  if(grid.domain$non.missing.number != length(z)) stop("Not proper length of z!")
  n <- length(z)
  DeltaTy <-   sparseG  %*%  z
  ##########################
  # the matrice-vector product required by conjugate.gradient:
  bEV<-sum(z**2)/n
  viaFFTwithMissings.prod.DeltaTCorrelDelta.Timesx <- function(x) {
    x <- transOfSparseG   %*%  x
    xprovFull<- expand.to.fullGrid(n1,  x, grid.domain$missing.sites)
    result<- c( stationary.image.cov(Y= xprovFull,cov.obj=cov.obj) )
    result<-  sparseG %*%  result[! grid.domain$missing.sites]
    result
  }
  ##########################
  coefProvFory<- z
  # loop over ranges  :
  for(k in 1:length(candidateThetas1DGrid)) {
     tethaCand <- candidateThetas1DGrid[k]
     #XXXXXX attention en 2016: il faut utiliser theta au lieu de range et matern
     #XXXXXX ds le nom de la fct
     cov.obj<- matern.image.cov( setup=TRUE,
              grid=list(x=1:n1,y=1:n1),
              theta= (1/tethaCand)*(n1-1),smoothness=nu)
     #
     out2<- conjugate.gradient(DeltaTy,  multAx=
              viaFFTwithMissings.prod.DeltaTCorrelDelta.Timesx,
              start= coefProvFory,
              tol= tolPGC,    kmax=200,verbose=FALSE)
     niter[k,1]<-out2$conv$niter
     coefProvFory<- out2$x  #save the sol for "start" for the next theta
     #
     GEvaluesOnthegrid[k] <- sum(DeltaTy *  coefProvFory) /n
  }
  list(values = GEvaluesOnthegrid,
              niterForY=niter[,1]
              #, niterForW=niter[,2],trace.term=trace.term #noisy case
              )
}
########### end of fsai11Precond.GEevalOnThetaGrid ##################
#####################################################################


fsai11Precond.GEevalOnThetaGridNEW <- function(z,candidateThetas1DGrid,nu=0.5,grid.domain,tolPGC=1e-03) {
  n1 <- grid.domain$n1
  sparseG <- grid.domain$sparseG
  transOfSparseG <- t(sparseG)
  #
  sparseG2 <- grid.domain$sparseGusingRange10smaller
  transOfSparseG2 <- t(sparseG2)
  #
  sparseG3 <- grid.domain$sparseGusingRange100smaller
  transOfSparseG3 <- t(sparseG3)
  #
  GEvaluesOnthegrid<- matrix(NA,length(candidateThetas1DGrid))
  #trace.term<- matrix(NA,length(candidateThetas1DGrid))
  niter <- matrix(NA,length(candidateThetas1DGrid),2)
  if(grid.domain$non.missing.number != length(z)) stop("Not proper length of z!")
  n <- length(z)
  #DeltaTy <-   sparseG  %*%  z
  ##########################
  # the matrice-vector product required by conjugate.gradient:
  bEV<-sum(z**2)/n
  ##########################
	# the (3 possible) matrice-vector product(s) required by conjugate.gradient:
		viaFFTwithMissings.prod.firstChoiceForG.DeltaTCorrelDelta.Timesx <- function(x) {
		    x <- transOfSparseG   %*%  x
		    xprovFull<- expand.to.fullGrid(n1,  x, grid.domain$missing.sites)
		    result<- c( stationary.image.cov(Y= xprovFull,cov.obj=cov.obj) )
		    result<-  sparseG %*%  result[! grid.domain$missing.sites]
		    result
	    }
		viaFFTwithMissings.prod.secondChoiceForG.DeltaTCorrelDelta.Timesx <- function(x) {
		    x <- transOfSparseG2   %*%  x
		    xprovFull<- expand.to.fullGrid(n1,  x, grid.domain$missing.sites)
		    result<- c( stationary.image.cov(Y= xprovFull,cov.obj=cov.obj) )
		    result<-  sparseG2 %*%  result[! grid.domain$missing.sites]
		    result
	    }
		viaFFTwithMissings.prod.thirdChoiceForG.DeltaTCorrelDelta.Timesx <- function(x) {
		    x <- transOfSparseG3   %*%  x
		    xprovFull<- expand.to.fullGrid(n1,  x, grid.domain$missing.sites)
		    result<- c( stationary.image.cov(Y= xprovFull,cov.obj=cov.obj) )
		    result<-  sparseG3 %*%  result[! grid.domain$missing.sites]
		    result
	    }
  ##########################
  coefProvFory<- z
  # loop over ranges  :
  for(k in 1:length(candidateThetas1DGrid)) {
     tethaCand <- candidateThetas1DGrid[k]
     candidateTheta <- tethaCand
     #XXXXXX attention en 2016: il faut utiliser theta au lieu de range et matern
     #XXXXXX ds le nom de la fct
     cov.obj<- matern.image.cov( setup=TRUE,
              grid=list(x=1:n1,y=1:n1),
              theta= (1/tethaCand)*(n1-1),smoothness=nu)
     #
     if (log(candidateTheta,10)< 0.5) {
		    DeltaTy<-  sparseG  %*%  z
		    out2<- conjugate.gradient(DeltaTy,
	            viaFFTwithMissings.prod.firstChoiceForG.DeltaTCorrelDelta.Timesx,  start= coefProvFory,
	            tol= tolPGC,    kmax=200,verbose=FALSE)
	  } 
	  else {
			    	if(log(candidateTheta,10)< 1.5) {
				    	DeltaTy<-  sparseG2 %*%  z
				    	out2<- conjugate.gradient(DeltaTy,
				            viaFFTwithMissings.prod.secondChoiceForG.DeltaTCorrelDelta.Timesx,  start= coefProvFory,
				            tol= tolPGC,    kmax=200,verbose=FALSE)
				   } 
			      else {
				    	DeltaTy<-  sparseG3 %*%  z
				    	out2<- conjugate.gradient(DeltaTy,
				            viaFFTwithMissings.prod.thirdChoiceForG.DeltaTCorrelDelta.Timesx,  start= coefProvFory,
				            tol= tolPGC,    kmax=200,verbose=FALSE)
	
			    	}
	    }
     niter[k,1]<-out2$conv$niter
     coefProvFory<- out2$x  #save the sol for "start" for the next theta
     #
     GEvaluesOnthegrid[k] <- sum(DeltaTy *  coefProvFory) /n
  }
  list(values = GEvaluesOnthegrid,
              niterForY=niter[,1]
              #, niterForW=niter[,2],trace.term=trace.term #noisy case
              )
}
########### end of fsai11Precond.GEevalOnThetaGrid ##################
#####################################################################
