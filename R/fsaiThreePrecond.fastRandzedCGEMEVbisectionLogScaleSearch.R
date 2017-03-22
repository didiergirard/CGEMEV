#################################################################
## definition of fsaiThreePrecond.fastRandzedCGEMEVbisectionLogScaleSearch ###############
#################################################################
fsaiThreePrecond.fastRandzedCGEMEVbisectionLogScaleSearch <- function(sigma2Noise,y,
w,
nu=0.5,
grid.domain,tolPGC=1e-03,
candidateThetas.LBound,candidateThetas.UBound,tolBisec=1e-04)
{ yStandardized  <- y/sqrt(sigma2Noise)
  n1 <- grid.domain$n1
  sparseG <- grid.domain$sparseG
  transOfSparseG <- t(sparseG)
  #
  sparseG2 <- grid.domain$sparseGusingRange10smaller
  transOfSparseG2 <- t(sparseG2)
  #
  sparseG3 <- grid.domain$sparseGusingRange100smaller
  transOfSparseG3 <- t(sparseG3)

	if(grid.domain$non.missing.number != length(y)) stop("Not proper length of y!")
	n=length(y)
	#DeltaTy <-   sparseG  %*%  z
	
	bEV<-sum(yStandardized**2)/n -1.0
	#
	#############################################
	fsai11Precond.CGEMeval <- function(n1, yStandardized,
	startForz, w,  startForw,
	candidateTheta
	#,nu=0.5, listOf.belongTo.OneOfTheDisks , sparseG,transOfSparseG
	,tolPGC= tolPGC)
	{
	##########################   lambda is global
	# the (3 possible) matrice-vector product(s) required by conjugate.gradient:
		viaFFTwithMissings.prod.firstChoiceForG.DeltaTCorrelPlusLambdaIdentDelta.Timesx <- function(x) {
		    x <- transOfSparseG   %*%  x
		    xprovFull<- expand.to.fullGrid(n1,  x, grid.domain$missing.sites)
		    result<- c( stationary.image.cov(Y= xprovFull,cov.obj=cov.obj) )
		    result<-  sparseG %*%  (result[! grid.domain$missing.sites]+ lambda* x )
		    result
	    }
		viaFFTwithMissings.prod.secondChoiceForG.DeltaTCorrelPlusLambdaIdentDelta.Timesx <- function(x) {
		    x <- transOfSparseG2   %*%  x
		    xprovFull<- expand.to.fullGrid(n1,  x, grid.domain$missing.sites)
		    result<- c( stationary.image.cov(Y= xprovFull,cov.obj=cov.obj) )
		    result<-  sparseG2 %*%  (result[! grid.domain$missing.sites]+ lambda* x )
		    result
	    }
		viaFFTwithMissings.prod.thirdChoiceForG.DeltaTCorrelPlusLambdaIdentDelta.Timesx <- function(x) {
		    x <- transOfSparseG3   %*%  x
		    xprovFull<- expand.to.fullGrid(n1,  x, grid.domain$missing.sites)
		    result<- c( stationary.image.cov(Y= xprovFull,cov.obj=cov.obj) )
		    result<-  sparseG3 %*%  (result[! grid.domain$missing.sites]+ lambda* x )
		    result
	    }

	coefProvFory<-startForz
	#coefProvFory<-0 * startForz
	coefProvForw <- startForw 
	#
	cov.obj<- matern.image.cov( setup=TRUE,
            grid=list(x=1:n1,y=1:n1),
            theta= (1/candidateTheta)*(n1-1),smoothness=nu)
    #
    
    if (log(candidateTheta,10)< 0.5) {
		    DeltaTy<-  sparseG %*%  yStandardized
				    	lambda   <- 1/ bEV
				    	out2<- conjugate.gradient(DeltaTy,
				            viaFFTwithMissings.prod.firstChoiceForG.DeltaTCorrelPlusLambdaIdentDelta.Timesx,  start= coefProvFory,
				            tol= tolPGC,    kmax=200,verbose=FALSE)
				        coefProvFory <- out2$x
				        lambda   <- 0.
	                    coefProvPROVFory <- viaFFTwithMissings.prod.firstChoiceForG.DeltaTCorrelPlusLambdaIdentDelta.Timesx (out2$x)
	                    QuadraticTermInConditionalGEmean <- sum(coefProvPROVFory *  coefProvFory) /n 
	                    DeltaTw<-  sparseG %*%  w
				    	lambda   <- 1/ bEV
				    	out2W<- conjugate.gradient(DeltaTw,
				            viaFFTwithMissings.prod.firstChoiceForG.DeltaTCorrelPlusLambdaIdentDelta.Timesx,  start= coefProvForw,
				            tol= tolPGC,    kmax=200,verbose=FALSE)
				        coefProvForw <- out2W$x
				        lambda   <- 0
	                    traceTermInConditionalGEmean <- sigma2Noise * sum(DeltaTw *  coefProvForw) /  sum( w *  w) 
	  } 
	  else {
			    	if(log(candidateTheta,10)< 1.5) {
				    	DeltaTy<-  sparseG2 %*%  yStandardized
				    	lambda   <- 1/ bEV
				    	out2<- conjugate.gradient(DeltaTy,
				            viaFFTwithMissings.prod.secondChoiceForG.DeltaTCorrelPlusLambdaIdentDelta.Timesx,  start= coefProvFory,
				            tol= tolPGC,    kmax=200,verbose=FALSE)
				        coefProvFory <- out2$x
				        lambda   <- 0.
	                    coefProvPROVFory <- viaFFTwithMissings.prod.secondChoiceForG.DeltaTCorrelPlusLambdaIdentDelta.Timesx (out2$x)
	                    QuadraticTermInConditionalGEmean <- sum(coefProvPROVFory *  coefProvFory) /n 
	                    DeltaTw<-  sparseG2 %*%  w
				    	lambda   <- 1/ bEV
				    	out2W<- conjugate.gradient(DeltaTw,
				            viaFFTwithMissings.prod.secondChoiceForG.DeltaTCorrelPlusLambdaIdentDelta.Timesx,  start= coefProvForw,
				            tol= tolPGC,    kmax=200,verbose=FALSE)
				        coefProvForw <- out2W$x
				        lambda   <- 0
	                    traceTermInConditionalGEmean <- sigma2Noise * sum(DeltaTw *  coefProvForw) /  sum( w *  w) 
				   } 
			      else {
				    	DeltaTy<-  sparseG3 %*%  yStandardized
				    	lambda   <- 1/ bEV
				    	out2<- conjugate.gradient(DeltaTy,
				            viaFFTwithMissings.prod.thirdChoiceForG.DeltaTCorrelPlusLambdaIdentDelta.Timesx,  start= coefProvFory,
				            tol= tolPGC,    kmax=200,verbose=FALSE)
				        coefProvFory <- out2$x
				        lambda   <- 0.
	                    coefProvPROVFory <- viaFFTwithMissings.prod.thirdChoiceForG.DeltaTCorrelPlusLambdaIdentDelta.Timesx (out2$x)
	                    QuadraticTermInConditionalGEmean <- sum(coefProvPROVFory *  coefProvFory) /n 
	                    DeltaTw<-  sparseG3 %*%  w
				    	lambda   <- 1/ bEV
				    	out2W<- conjugate.gradient(DeltaTw,
				            viaFFTwithMissings.prod.thirdChoiceForG.DeltaTCorrelPlusLambdaIdentDelta.Timesx,  start= coefProvForw,
				            tol= tolPGC,    kmax=200,verbose=FALSE)
				        coefProvForw <- out2W$x
				        lambda   <- 0
	                    traceTermInConditionalGEmean <- sigma2Noise * sum(DeltaTw *  coefProvForw) /  sum( w *  w) 
			    	}
	    }
       
    # coefProvFory <- out2$x
    # # coefProv<- solve( Correl.mat,y)
    # GEvalue <- sum(DeltaTy *  coefProvFory) /n
    CGEMvalue <- QuadraticTermInConditionalGEmean + traceTermInConditionalGEmean
    list(value = CGEMvalue,  niterForY=out2$conv$niter,niterForW=out2W$conv$niter,
                 coefForY= coefProvFory,coefForW= coefProvForw)
	}
	## end of CG.eval ####################
#
niterBisectionMax <- 22
niterCGiterationsHistory<-matrix(0., niterBisectionMax, 2)
candidateThetasHistory<-rep(NA, niterBisectionMax)
#
x1 <- candidateThetas.LBound
x2 <- candidateThetas.UBound
startForz  <- 0 *yStandardized
startForw  <- 0 *yStandardized
out1 <- fsai11Precond.CGEMeval(n1,yStandardized, startForz, w,  startForw, x1,tolPGC=tolPGC)
out2 <- fsai11Precond.CGEMeval(n1,yStandardized, startForz, w,  startForw, x2,tolPGC=tolPGC)
startForz  <- out2$coefForY
startForw  <- out2$coefForW
#
f1 <-    bEV - out1$value
f2 <-    bEV - out2$value
niterCGiterationsHistory[1,1] <- out1$niterForY
niterCGiterationsHistory[2,1] <- out2$niterForY
niterCGiterationsHistory[1,2] <- out1$niterForW
niterCGiterationsHistory[2,2] <- out2$niterForW
#
candidateThetasHistory[1] <- x1
candidateThetasHistory[2] <- x2
#listSuccessiveValues<-c(x1,f1,x2,f2)
if (f1 > f2)
        stop(" f1 must be < f2 ")
        #
        for (k in 1:niterBisectionMax) {
        xm <- sqrt(x1 * x2)     ################ instead of the mean
        outm <- fsai11Precond.CGEMeval(n1,yStandardized,  startForz, w,  startForw, xm,tolPGC=tolPGC)
        startForz  <- outm$coefForY
        startForw  <- outm$coefForW
        #
        niterCGiterationsHistory[k+2,1] <- outm$niterForY
        niterCGiterationsHistory[k+2,2] <- outm$niterForW
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
	niterCGiterationsHistory = c(niterCGiterationsHistory[,1]),
	niterCGiterationsHistoryForw = c(niterCGiterationsHistory[,2]))
}
## end of fsai11Precond.GEbisectionLogScaleSearch ####################
######################################################################
