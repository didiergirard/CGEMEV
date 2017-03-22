Untitled
================

Estimating a centered isotropic Matérn field from a (possibly incomplete and noisy) lattice observation
=======================================================================================================

The `CGEMEV` R package (dependent on the `field`, `spam` and `RCPP` packages) provides

-   a simple function for simulating realizations of a stationary isotropic Gaussian process when the correlation belongs to the common Matérn family with smoothness index *nu &gt; 0*
-   tools for estimating, when *nu* and the level *sigma<sub>NOISE</sub>* of the possible \`\`measurement errors'' are known, the correlation range also called ''decorrelation length'' (or its inverse that will be denoted *theta*) from one realization on a (possibly incomplete) lattice. The variance of the field is simply estimated by the empirical variance that will be denoted *b<sub>EV</sub>*.

These tools implement the estimation method proposed in Girard, D.A., 2016, Asymptotic near-efficiency of the Gibbs-energy and empirical-variance estimating functions for fitting Matérn models I: Densely sampled processes. Statistics and Probability Letters 110, 191-197, and in the preprint <https://hal.archives-ouvertes.fr/hal-00515832> .

In this first version of `CGEMEV`, the region of missing observations is defined as a union of disks.

Four main functions are used: `gaussian.matern()`, `simulate()`, `fsai11Precond.GEevalOnThetaGrid()` and `fsaiThreePrecond.fastRandzedCGEMEVbisectionLogScaleSearch()`, and a fifth function `grid.domain()` is required to precompute preconditioning sparse matrices.

These R-fonctions can be applied to a quite large grid even on a laptop (for example 512x512, provided the ''extension factor'' required for simulation, see below, is not too big). Indeed quite fast computation of the quadratic form which occurs in the estimating equation is possible by using a conjugate-gradient (CG) solver preconditioned by a classical factored sparse approximate inverse (FSAI) preconditioning, since the matrix-vector product, required in each CG iteration, can be obtained via FFT from the standard embedding of the correlation matrix in a circulant matrix.

Contents
--------

-   [Setting the probabilistic model](#Setting-the-probabilistic-model)
-   [Simulating one realization](#Simulating-one-realisation)
-   [Plotting (and saving) several realizations](#Plotting-(and-saving)-several-realizations)
-   [Setting the uncomplete lattice](#Setting-the-uncomplete-lattice)
-   [Plotting data](#Plotting-data)
-   [Computing (and plotting) the estimating function at log-equispaced ranges](#Computing-(and-plotting)-the-estimating-function-at-log-equispaced-ranges)
-   [Estimating theta and the micro-ergodic parameter](#Estimating-theta-and-the-micro-ergodic-parameter)

Setting the probabilistic model
-------------------------------

In this first version of the `CGEMEV` package, for simplicity, we restrict the spatial domain to be the unit square (0,1)X(0,1). For the example here, the simulations and the choice of observed sites are done using a 128x128 grid which partitions this domain. We consider the example `nu =1` (recall that is a widely used correlation following the seminal paper by Peter Whittle in Biometrika, 41(3–4), 1954 pp. 434–449). Let us choose the correlation range such that its inverse, denoted \(\theta\), satisfies \(\sqrt{2 nu} \theta^{-1} =0.2\). The resulting correlation function can then be considered as one with an ‘’effective range’‘ equal to 0.2 (the formulae for Matern correlations often use this quantity denoted \(\rho\), see <https://en.wikipedia.org/wiki/Mat%C3%A9rn_covariance_function> ).

``` r
library(CGEMEV)
n1grid <- 128
# gaussian matern creation
nu <-0.5
effectiveRange <- 0.2
gm <- gaussian.matern(grid.size=n1grid,smoothness=nu,
                      range=effectiveRange/sqrt(2*nu),factor=2)
```

NB: in the previous setting, `factor=2` specifies the required extension factor (assumed, for simplicity, to be an integer) of the observation domain. Indeed for this example the choice `factor=1` would entail (when calling `simulate(gm)`) the message “FFT of covariance has negative values” which means that generating a realization via the classical embedding method (which doubles each length of the considered rectangular domain) would not work.

Simulating one realization
--------------------------

This is simply:

``` r
set.seed(321)  # so that it is reproducible #
simulate(gm)
```

Plotting (and saving) 6 realizations *Z<sub>1</sub>*, ...,*Z<sub>6</sub>*
-------------------------------------------------------------------------

We can plot (and save), the previous realization and, for example, 5 further realizations:

``` r
fullLattice.sixZs<- array(NA,c(n1grid*n1grid,6))
set.panel(2,3)
```

    ## plot window will lay out plots in a 2 by 3 matrix

``` r
plot(gm)
fullLattice.sixZs[,1]<-gm$look[1:gm$n1,1:gm$n1]
ut <- system.time(
for (indexReplcitate in 2:6){
  set.seed(320+indexReplcitate)
  simulate(gm)
  plot(gm)
  fullLattice.sixZs[,indexReplcitate]<-gm$look[1:gm$n1,1:gm$n1]
})
```

![](TESTJUNK_files/figure-markdown_github/unnamed-chunk-3-1.png)

The following timing (in seconds) is for a MacMini (late2012) I7 2.3GHz :

``` r
ut   # for the simulation of 5 realizations :
```

    ##    user  system elapsed 
    ##   2.374   0.042   2.439

Setting the uncomplete lattice
------------------------------

Let us now define the regions (actually 5 disks) where the observations will be missing, and precompute the preconditioning matrix:

``` r
# md=missing.domains
listOfHalfDiameters <- c(1,0.4,0.2,0.1,0.05)*sqrt(.1/ pi)

listOfCenters <- t(matrix(c(c(0.5+sqrt(.1/ pi),0.5+sqrt(.1/ pi)),c(0.1+sqrt(.1/ pi),
0.05+sqrt(.1/ pi)),c(0.8,0.25),c(0.2,0.8),c(0.3,0.7)),2,length(listOfHalfDiameters)))

ex1.md <- list(
  list(center=listOfCenters[1,],radius= listOfHalfDiameters[1]),
list(center=listOfCenters[2,],radius= listOfHalfDiameters[2]),
list(center=listOfCenters[3,],radius= listOfHalfDiameters[3]),
list(center=listOfCenters[4,],radius= listOfHalfDiameters[4]),
list(center=listOfCenters[5,],radius= listOfHalfDiameters[5])
)
# gd=grid.domain
print(system.time(ex1WithN1eq128And5missingDisks.gd <- grid.domain(missing.domains=ex1.md,grid.size=n1grid,
            smoothness=nu)))
```

    ##    user  system elapsed 
    ##  10.136   0.894  10.668

Size of the data set:

``` r
(nObs <-dim(ex1WithN1eq128And5missingDisks.gd$sparseG)[1])
```

    ## [1] 14431

Average width of the sparse preconditioner:

``` r
length(ex1WithN1eq128And5missingDisks.gd$sparseG) / dim(ex1WithN1eq128And5missingDisks.gd$sparseG)[1]
```

    ## [1] 72.06645

Plotting observation sites
--------------------------

Each dataset will be made up of the values of one realization *Z* of the previous type at the following points:

``` r
xFull<- (1./n1grid)*matrix( c(rep(1: n1grid, n1grid),rep(1: n1grid,each= n1grid)), ncol=2)
x <- xFull[!ex1WithN1eq128And5missingDisks.gd$missing.sites,]
plot(x,asp=1, xlim=c(0,1), ylim=c(0,1), pch=".")
```

![](TESTJUNK_files/figure-markdown_github/unnamed-chunk-8-1.png)

Computing (and plotting) the estimating function at log-equispaced ranges
-------------------------------------------------------------------------

Choose a grid of candidates for the range parameter (more precisely, for the inverse-range parameter denoted `theta`) at which the estimating function is computed:

``` r
(thetaTrue  <- 1/gm$range)
```

    ## [1] 5

``` r
candidateThetas1DGrid <- thetaTrue * 10**seq(-0.8,0.8,,15)
```

Consider the first one of the above realizations, and the naive variance estimator *b<sub>EV</sub>*:

``` r
# only observed outside the disks:
#z <- gm$look[1:gm$n1,1:gm$n1][!ex1WithN1eq128And5missingDisks.gd$missing.sites]
indexReplcitate <- 1
z <- fullLattice.sixZs[,indexReplcitate][!ex1WithN1eq128And5missingDisks.gd$missing.sites]
#
(bEV  <- mean(z**2))
```

    ## [1] 1.558528

For this dataset `z`, let us give the whole output of the function `fsai11Precond.GEevalOnThetaGrid()` :

``` r
(out <- fsai11Precond.GEevalOnThetaGridNEW(z,candidateThetas1DGrid,nu=gm$smoothness,                          
grid.domain=ex1WithN1eq128And5missingDisks.gd,tolPGC=1e-04)
)
```

    ## $values
    ##            [,1]
    ##  [1,] 6.1280797
    ##  [2,] 4.7103954
    ##  [3,] 3.6207859
    ##  [4,] 2.7833609
    ##  [5,] 2.1397962
    ##  [6,] 1.6452701
    ##  [7,] 1.2653415
    ##  [8,] 0.9735654
    ##  [9,] 0.7496387
    ## [10,] 0.5780086
    ## [11,] 0.4468007
    ## [12,] 0.3470098
    ## [13,] 0.2719343
    ## [14,] 0.2167601
    ## [15,] 0.1783596
    ## 
    ## $niterForY
    ##  [1] 17 14 14 14 13 13 13 11  9  7  6  7  9 13 19

Let us repeat this computation for the five next data sets obtained from the above realizations, and plot the results:

``` r
set.panel(2,3)
```

    ## plot window will lay out plots in a 2 by 3 matrix

``` r
{plot(candidateThetas1DGrid, out$values, , type="l",
                 col=1,lty= "dashed",log="xy")
    abline(h= bEV, lty= "dotted" )
}
ut <- system.time(
for (indexReplcitate in 2:6){
  z <- fullLattice.sixZs[,indexReplcitate][!ex1WithN1eq128And5missingDisks.gd$missing.sites]
  bEV  <- mean(z**2)
  out <-     fsai11Precond.GEevalOnThetaGridNEW(z, candidateThetas1DGrid,
          gm$smoothness, ex1WithN1eq128And5missingDisks.gd ,tolPGC=1e-04)
#  
  plot(candidateThetas1DGrid, out$values, type="l",
                 col=1, lty= "dashed",log="xy")
    abline(h= bEV, lty= "dotted")
})
```

![](TESTJUNK_files/figure-markdown_github/unnamed-chunk-13-1.png)

Timing for a MacMini (late2012) I7 2.3GHz :

``` r
ut   # for computing the estimating equation for 5 realizations :
```

    ##    user  system elapsed 
    ##  22.427   5.935  28.454

Estimating theta and the micro-ergodic parameter from noisy observations with (unknown) SNR=1000
------------------------------------------------------------------------------------------------

For solving the CGEMEV estimating equation in theta, a specific bisection method has been elaborated, since it is useful to use "historical" results (more precisely the starting point of the PCG iterations at a given theta is chosen as the solution from the linear-solve associated with the previously considered theta); indeed a classic use of the `uniroot` R-fonction would be much slower. Let us analyse N=200 datasets *Y<sub>k</sub>*, *k=1,...,200* each one being *Y<sub>k</sub> = Z<sub>k</sub> +epsilon<sub>k</sub> *, where *epsilon<sub>k</sub>* is a field of white normal noise; they are simulated such that the true SNR (i.e. var(Z)/var(noise)) = 1000. Only *nu* and the variance of the noise ("measurement errors") are assumed known.

``` r
sigma2Process<- 1000
sigma2Noise <- 1.
bTrue <- sigma2Process/sigma2Noise
  
nbReplicates <- 200
bHatEV<-matrix(NA,nbReplicates)
cHatCGEMEV<-matrix(NA, nbReplicates)
thetaHatCGEMEV<-matrix(NA, nbReplicates)
nCGiterationsMaxForY<-matrix(NA, nbReplicates)
#
ut <- system.time(
for (indexReplcitate in 1: nbReplicates){
  set.seed(720+indexReplcitate)
  simulate(gm)
  z <- gm$look[1:gm$n1,1:gm$n1][!ex1WithN1eq128And5missingDisks.gd$missing.sites]
  y <-  sqrt(bTrue)* z +  c(rnorm(nObs))
#
  bEV  <- ( mean(y**2)-sigma2Noise ) / sigma2Noise
  if (bEV <0) {
            break}
  w <-  c(rnorm(nObs))
  out <-     fsaiThreePrecond.fastRandzedCGEMEVbisectionLogScaleSearch(
          sigma2Noise,y,
          w,
          gm$smoothness, ex1WithN1eq128And5missingDisks.gd ,tolPGC=1e-04,
        0.2,50, tolBis=1e-05)

#  
  thetaHatCGEMEV[indexReplcitate]<- (out$root)
  cHatCGEMEV[indexReplcitate]<-
                            bEV*(out$root)**(2*gm$smoothness)
  nCGiterationsMaxForY[indexReplcitate]<-
                 max(out$niterCGiterationsHistory)
})
ut
```

    ##     user   system  elapsed 
    ## 1580.074  452.702 2040.524

``` r
print(c("log10 of true range =",log(sqrt(2*nu)*gm$range,10)))
```

    ## [1] "log10 of true range =" "-0.698970004336019"

``` r
summary(log(sqrt(2*nu)/thetaHatCGEMEV,10))
```

    ##        V1         
    ##  Min.   :-1.0262  
    ##  1st Qu.:-0.8029  
    ##  Median :-0.7169  
    ##  Mean   :-0.7122  
    ##  3rd Qu.:-0.6298  
    ##  Max.   :-0.3396

``` r
cTrue <-  bTrue*(1/gm$range)**(2*gm$smoothness)
summary(cHatCGEMEV/cTrue)
```

    ##        V1        
    ##  Min.   :0.9687  
    ##  1st Qu.:0.9917  
    ##  Median :0.9994  
    ##  Mean   :1.0008  
    ##  3rd Qu.:1.0086  
    ##  Max.   :1.0421

``` r
sd(cHatCGEMEV/cTrue)
```

    ## [1] 0.01298273

This sd can be compared to the CR lower bound assuming theta is known and no noise in the data :

``` r
sqrt(2/length(z))
```

    ## [1] 0.01177245

Let us plot an estimate of the density of the CGEMEV estimates of the effective range, and of the density of the relative errors in the CGEMEV estimates of the microergodic-parameter cTrue:

``` r
set.panel(2,2)
```

    ## plot window will lay out plots in a 2 by 2 matrix

``` r
den <- density(log(sqrt(2*nu)/thetaHatCGEMEV,10))
den$x <- 10**(den$x)
plot(den, log="x",main="")
title("CGEMEVestimates of the range (=0.2)")
plot(density((cHatCGEMEV-cTrue)/cTrue),main="")
title(" N errors (cHatCGEMEV-cTrue)/cTrue")
```

![](TESTJUNK_files/figure-markdown_github/unnamed-chunk-18-1.png)
