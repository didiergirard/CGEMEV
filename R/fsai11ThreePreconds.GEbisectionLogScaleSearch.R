#################################################################
## definition of fsai11ThreePrecond.GEbisectionLogScaleSearchNEW ###############
#################################################################
fsai11ThreePreconds.GEbisectionLogScaleSearchNEW <- function(z,
#w,
nu=0.5,
grid.domain,tolPGC=1e-03,
candidateThetas.LBound,candidateThetas.UBound,tolBisec=1e-04)
{
  n1 <- grid.domain$n1
  sparseG <- grid.domain$sparseG
  transOfSparseG <- t(sparseG)
  #
  sparseG2 <- grid.domain$sparseGusingRange10smaller
  transOfSparseG2 <- t(sparseG2)
  #
  sparseG3 <- grid.domain$sparseGusingRange100smaller
  transOfSparseG3 <- t(sparseG3)

	if(grid.domain$non.missing.number != length(z)) stop("Not proper length of z!")
	n=length(z)
	#DeltaTy <-   sparseG  %*%  z
	
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

	coefProvFory<-startForz
	#coefProvFory<-0 * startForz
	#
	cov.obj<- matern.image.cov( setup=TRUE,
            grid=list(x=1:n1,y=1:n1),
            theta= (1/candidateTheta)*(n1-1),smoothness=nu)
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
       
    coefProvFory <- out2$x
    # coefProv<- solve( Correl.mat,y)
    GEvalue <- sum(DeltaTy *  coefProvFory) /n
    list(value = GEvalue,  niterForY=out2$conv$niter,
                 coefForY=coefProvFory)
	}
	## end of CG.eval ####################
#
niterBisectionMax <- 22
niterCGiterationsHistory<-matrix(0., niterBisectionMax, 2)
candidateThetasHistory<-rep(NA, niterBisectionMax)
#
x1 <- candidateThetas.LBound
x2 <- candidateThetas.UBound
startForz  <- 0 *z
out1 <- fsai11Precond.GEeval(n1,z, startForz, x1,tolPGC=tolPGC)
out2 <- fsai11Precond.GEeval(n1,z, startForz, x2,tolPGC=tolPGC)
startForz  <- out2$coefForY
#
f1 <-    bEV - out1$value
f2 <-    bEV - out2$value
niterCGiterationsHistory[1,1] <- out1$niterForY
niterCGiterationsHistory[2,1] <- out2$niterForY
candidateThetasHistory[1] <- x1
candidateThetasHistory[2] <- x2
#listSuccessiveValues<-c(x1,f1,x2,f2)
if (f1 > f2)
        stop(" f1 must be < f2 ")
        #
        for (k in 1:niterBisectionMax) {
        xm <- sqrt(x1 * x2)     ################ instead of the mean
        outm <- fsai11Precond.GEeval(n1,z,  startForz, xm,tolPGC=tolPGC)
        startForz  <- outm$coefForY
        #
        niterCGiterationsHistory[k+2,1] <- outm$niterForY
				candidateThetasHistory[k+2] <- xm
        #
        fm <-    bEV - outm$value
        if (fm < 0) {
            x1 <- xm
            f1 <- fm
        }
        else {
            x2 <- xm
            f2 <- fm
        }
        if (abs(1- x2/x1) < tolBisec) {
            break
        }
       #listSuccessiveValues<-append(listSuccessiveValues,xm,fm)
    }
    xm <- sqrt(x1 * x2)         ################ instead of the mean
#
list(root=xm, candidateThetasHistory = candidateThetasHistory,
	niterCGiterationsHistory = c(niterCGiterationsHistory[,1]))
}
## end of fsai11Precond.GEbisectionLogScaleSearch ####################
######################################################################
