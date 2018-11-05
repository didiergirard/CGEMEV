#################################################################
## definition of sig2noiseByVariogExtrapolAtZero ###############
#################################################################
sig2noiseByVariogExtrapolAtZero <- function(y,
grid.domain)
{
  n1grid <- grid.domain$n1
	tmp <- array(y,c(n1grid,n1grid))
	imageOfyWithMissingData<-tmp
	imageOfyWithMissingData[grid.domain$missing.sites!=FALSE]<-0
  #
	outWithMissing<-structurogram.matrix(imageOfyWithMissingData)
	variog<-outWithMissing$vgram
	lag<-outWithMissing$d
	lag2<-lag^2
	quadratic.model <-lm(variog ~ lag + lag2)
	#
	sig2noiseVEZ <-quadratic.model$coefficients[1]
#
sig2noiseVEZ
}
## end of sig2noiseByVariogExtrapolAtZero ####################
##########################################################
