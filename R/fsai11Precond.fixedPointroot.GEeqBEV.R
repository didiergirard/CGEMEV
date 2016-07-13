#################################################################
## definition of fsai11Precond.fixedPointroot.GEeqBEV ###############
#################################################################
fsai11Precond.fixedPointroot.GEeqBEV <- function(z,
#w,
nu=0.5,
grid.domain,tolPGC=1e-03,
candidateTethas.LBound,tol)
{
  n1 <- grid.domain$n1
  sparseG <- grid.domain$sparseG
  transOfSparseG <- t(sparseG)

  if(grid.domain$non.missing.number != length(z)) stop("Not proper length of z!")
  n <- length(z)
  DeltaTy <-   sparseG  %*%  z
  bEV<-sum(z**2)/n	
  #
#############################################
	fsai11Precond.GEeval <- function(n1,z,
	startForz, 
	candidateTheta
	#,nu=0.5, listOf.belongTo.OneOfTheDisks , sparseG,transOfSparseG
	,tolPGC= tolPGC)
	{
		##########################
		# the matrice-vector product required by conjugate.gradient:
		viaFFTwithMissings.prod.DeltaTCorrelDelta.Timesx <- function(x) {
	    x <- transOfSparseG   %*%  x
	    xprovFull<- expand.to.fullGrid(n1,  x, grid.domain$missing.sites)
	    result<- c( stationary.image.cov(Y= xprovFull,cov.obj=cov.obj) )
	    result<-  sparseG %*%  result[! grid.domain$missing.sites]
	    result
  }
	coefProvFory<-startForz
	cov.obj<- matern.image.cov( setup=TRUE, 
            grid=list(x=1:n1,y=1:n1), 
            theta= (1/candidateTheta)*(n1-1),smoothness=nu)
    #
    out2<- conjugate.gradient(DeltaTy, multAx=
            viaFFTwithMissings.prod.DeltaTCorrelDelta.Timesx,  start= coefProvFory,
            tol= tolPGC,    kmax=200,verbose=FALSE)   
    coefProvFory <- out2$x
    # coefProv<- solve( Correl.mat,y)
    GEvalue <- sum(DeltaTy *  coefProvFory) /n 
    list(value = GEvalue,  niterForY=out2$conv$niter,  
                 coefForY=coefProvFory)
	}
	## end of CG.eval ####################
#

	
	
	
	
}
