#################################################################
## definition of fsai11Precond.GEbisectionLogScaleSearch ###############
#################################################################
fsai11Precond.GEbisectionLogScaleSearch <- function(z,
#w,
nu=0.5,
grid.domain,tolPGC=1e-03,
candidateTethas.LBound,candidateTethas.UBound,tolBisec=1e-04)
{
	n1 <- grid.domain$n1
  sparseG <- grid.domain$sparseG
  transOfSparseG <- t(sparseG)
if(grid.domain$non.missing.number != length(z)) stop("Not proper length of z!")
n=length(z)
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
    out2<- conjugate.gradient(DeltaTy,
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
niterBisectionMax <- 22
niterCGiterationsHistory<-matrix(0., niterBisectionMax, 2)
candidateTethasHistory<-rep(NA, niterBisectionMax)
#
x1 <- candidateTethas.LBound
x2 <- candidateTethas.UBound
out1 <- fsai11Precond.GEeval(n1,z,z, x1,tolPGC=tolPGC)
out2 <- fsai11Precond.GEeval(n1,z,z, x2,tolPGC=tolPGC)
startForz  <- out2$coefForY
#
f1 <-    bEV - out1$value
f2 <-    bEV - out2$value
niterCGiterationsHistory[1,1] <- out1$niterForY
niterCGiterationsHistory[2,1] <- out2$niterForY
candidateTethasHistory[1] <- x1
candidateTethasHistory[2] <- x2
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
				candidateTethasHistory[k+2] <- xm
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
list(root=xm, niterCGiterationsHistory = niterCGiterationsHistory)
}
## end of fsai11Precond.GEbisectionLogScaleSearch ####################
##########################################################
