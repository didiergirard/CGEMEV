library(CGEMEV)

# md=missing.domains
ex1.md <- list(
	list(center=c(0.67841,0.67841),radius=0.17841),
	list(center=c(0.278412, 0.228412),radius=0.071365)
)

ex2.md  <- c(ex1.md,list(center=c(0.8, 0.25),radius=0.03568))

# gd=grid.domain
print(system.time(ex1.gd <- grid.domain(ex1.md,ngrid<-64)))
#
#
cat("gaussian matern creation\n")
gm <- gaussian.matern(grid.size=ngrid)
cat("-> done\n")


cat("simulation and plot\n")
set.seed(321))  # so that it is reproducible #
simulate(gm)
#plot(gm)


# only observed outside the disks:
z <- gm$look[1:gm$n1,1:gm$n1][!ex1.gd$missing.sites]
cat("number of observations\n",length(z),"\n")
plotUncompleteGrid.gaussian.matern(gm,ex1.gd)
#
sum(z**2)/ex1.gd$non.missing.number
cat("-> done\n")

candidateThetas.Grid <- 1/gm$range * 10**seq(-1.1,1.1,,15)

print("iciiiiii")
system.time(print(out <-
	 fsai11Precond.GEevalOnThetaGrid(z, candidateThetas.Grid, gm$smoothness, ex1.gd ,tolPGC=1e-03)))

plot( candidateThetas.Grid, out$values, type="l")

system.time(print(out <-
	 fsai11Precond.GEbisectionLogScaleSearch(z,
		  gm$smoothness, ex1.gd ,tolPGC=1e-03, 1, 100,  1e-04)))


#########################################################
			####### generation of 100 replicates            #########
			#######  and  CGEM-EV estimates for each one   ##########
			bHatEV<-matrix(NA,100)
			cHatGEEVLogScaleSearch<-matrix(NA,100)
			nCGiterationsMaxForYLogScaleSearch<-matrix(NA,100)
			#nCGiterationsMaxForWLogScaleSearch<-matrix(NA,500)
			ut<-system.time(
			for(indexReplcitate in 1:100){
			set.seed(indexReplcitate)  # so that it is reproducible #
			simulate(gm)
			z <- gm$look[1:gm$n1,1:gm$n1][!ex1.gd$missing.sites]
			bEV<-sum(z**2)/ex1.gd$non.missing.number # -1.
			#
			# w <- rnorm(n)
			# w<- c( w)
			# w <- c( w)*sqrt((n/sum( w *w)))
			# #
			out <- fsai11Precond.GEbisectionLogScaleSearch(z,
	 		  gm$smoothness, ex1.gd ,tolPGC=1e-03, 1, 100,  1e-04)
			cHatGEEVLogScaleSearch[indexReplcitate]<-
							bEV*(out$root)**(2*gm$smoothness)
			nCGiterationsMaxForYLogScaleSearch[indexReplcitate]<-
				 max(out$niterCGiterationsHistory[,1])
			# nCGiterationsMaxForWLogScaleSearch[indexReplcitate]<-
			# 	max(out$niterCGiterationsHistory[,2])
			#
			}
			)
			ut
			ctrue<-1.*(1/gm$range)**(2*gm$smoothness)
			summary(cHatGEEVLogScaleSearch[1:100]/ctrue)
			sd(cHatGEEVLogScaleSearch[1:100]/ctrue)
			nCGiterationsMaxForYLogScaleSearch
