# sites
grid <- function(dim) {
  # TODO
}

# matrix (only called when size od$f sites is reasonable)
cor.matern.matrix <- function(sites,theta,b) { #notation ref à l'article de Didier

}

# model
matern.cov.field <- function(nu,theta,b=1) { # contructor
  obj <- list(nu=nu,theta=theta,b=b)
  class(obj) <- "matern.cov.field"
  obj
}

# simulation

simulate.matern.cov.field <- function(n,sites=grid(c(100,100)),missing.sites) {
## z=b^(.5)L^Tw (ou version sans calcul de R ni de L => circulante embedding but only for grid !!!!)
## return z
}


################# precomputations : ##########



#### INPUTS
## arrayOfIndicesOfObservationInVectorZ -> obj$z
## array.listOfPixels.in.OneOfTheDisks -> obj$missing.sites
## array.listOfDistantEnoughPixels.from.BoundaryAndTheDisks -> obj$distant.sites


cgemev <- function(z,missing.sites,tol.bissection=1e-4,tol.pcg=1e-4) {
## theta_r+1=theta_r (z^T R(theta_r)^(-1) z/ z^Tz)^(1/2*nu)
## preconditionning
## G^T R(theta) Gw=G^T z => x=Gw
## with FSAI: see algo
}
