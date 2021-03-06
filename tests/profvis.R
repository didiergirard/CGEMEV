library(CGEMEV)
require(profvis)

# md=missing.domains
ex1.md <- list(
	list(center=c(0.67841,0.67841),radius=0.17841),
	list(center=c(0.278412, 0.228412),radius=0.071365)
)

ex2.md  <- c(ex1.md,list(center=c(0.8, 0.25),radius=0.03568))

# gd=grid.domain
print(profvis({
  ex1.gd <- grid.domain(ex1.md,ngrid <- 256)
  gm <- gaussian.matern(grid.size=ngrid)

  set.seed(321)  # so that it is reproducible #


  simulate(gm)
  z <- gm$look[1:gm$n1,1:gm$n1][!ex1.gd$missing.sites]
  length(z)
  #
  sum(z**2)/ex1.gd$non.missing.number

  candidateThetas.Grid <- 1/gm$range * 10**seq(-1.1,1.1,,15)

  out <- fsai11Precond.GEevalOnThetaGrid(z, candidateThetas.Grid, gm$smoothness, ex1.gd ,tolPGC=1e-03)

}))
