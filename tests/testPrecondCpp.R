library(CGEMEV)

# md=missing.domains
ex1.md <- list(
	list(center=c(0.67841,0.67841),radius=0.17841),
	list(center=c(0.278412, 0.228412),radius=0.071365)
)

ex2.md  <- c(ex1.md,list(center=c(0.8, 0.25),radius=0.03568))

# gd=grid.domain
print(system.time(ex1.gd <- grid.domain(ex1.md,ngrid<-512)))
#
#
cat("gaussian matern creation\n")
gm <- gaussian.matern(grid.size=ngrid)
cat("-> done\n")
set.seed(321)  # so that it is reproducible #

cat("simulation and plot\n")
simulate(gm)
plot(gm)
cat("-> done\n")


# only observed outside the disks:
z <- gm$look[1:gm$n1,1:gm$n1][!ex1.gd$missing.sites]
length(z)
#
sum(z**2)/ex1.gd$non.missing.number


candidateThetas.Grid <- 1/gm$range * 10**seq(-1.1,1.1,,15)

print("iciiiiii")
system.time(print(out <-
	 fsai11Precond.GEevalOnThetaGrid(z, candidateThetas.Grid, gm$smoothness, ex1.gd ,tolPGC=1e-03)))

plot( candidateThetas.Grid, out@values, type="l")
