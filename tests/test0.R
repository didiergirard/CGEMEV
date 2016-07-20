library(CGEMEV)

# md=missing.domains
ex1.md <- list(
	list(center=c(0.67841,0.67841),radius=0.17841),
	list(center=c(0.278412, 0.228412),radius=0.071365)
)

ex2.md  <- c(ex1.md,list(center=c(0.8, 0.25),radius=0.03568))

# gd=grid.domain
cat("n=",ngrid<-1024,"\n")
print(system.time(ex1.gd <- grid.domain(ex1.md,ngrid)))
#tmp <- preconditioning2.grid.domain(ex1.gd)
