

grid.domain  <- function(missing.domains,grid.size=256,oversampling.size=2.5,dim=2) {
  require(spam)
	obj <- new.env()
	obj$missing.domains <- missing.domains # disk only now! TODO: extend to any form of domain
	obj$n.domains <- length(missing.domains)

	obj$dim <- dim
	if(length(grid.size)==1) grid.size <- rep(grid.size,dim)
	obj$grid.size <- grid.size
  obj$sites.number <- prod(grid.size)

  ## Temporary!!
	obj$n1 <- obj$grid.size[1] #TODO: extend to other domain

  ## OLD: xFull => TODO: pourquoi "(1./(obj$n1-1)) *" ?
	obj$coords <- (1./(obj$n1-1)) * matrix( c(rep(1:obj$n1,obj$n1),rep(1:obj$n1,each=obj$n1)), ncol=2) # dim or 2

  obj$oversampling.size <- oversampling.size

  class(obj) <- "grid.domain"

  ## everything done here!!!
	update(obj)

	obj
}

update.grid.domain <-function(obj) {
	missing.sites.grid.domain(obj)
	distant.sites.grid.domain(obj)
	non.missing.coords.grid.domain(obj)
	non.missing.number.grid.domain(obj)
  expand.to.full.grid.domain(obj)
  preconditionning.grid.domain(obj)
}

missing.sites.grid.domain <- function(obj) {
	n1Xn1  <- obj$n1*obj$n1
	missing.sites <- rep(FALSE,  n1Xn1)
	for(indexPixel in 1:n1Xn1) for(id in 1:obj$n.domains) missing.sites[indexPixel] <- missing.sites[indexPixel] | any(belong.to.disk(obj$coords[indexPixel,] ,  obj$missing.domains[[id]]$center,obj$missing.domains[[id]]$radius^2 ))
	# OLD: listOfbelong.to.OneOfTheDisks. Rmk: matrix in R can be used as a vector since it is actually a vector
	obj$missing.sites <- matrix(missing.sites,nrow=obj$grid.size[1],ncol=obj$grid.size[2])
}

distant.sites.grid.domain <- function(obj) {
	n1Xn1  <- obj$n1*obj$n1


	halfWidthNeighborood  <-  obj$oversampling.size/obj$n1
	augmentedRadius  <-  sapply(obj$missing.domains,function(md) md$radius + halfWidthNeighborood)

  roundHalfWidth  <-  round(obj$oversampling.size)
	distant.sites <- rep(FALSE,  n1Xn1)
	for(i in (roundHalfWidth+1):(obj$n1-roundHalfWidth)) for(j in (roundHalfWidth+1):(obj$n1-roundHalfWidth))  {
		indexPixel <- j + (i-1)*obj$n1
    tmp <- TRUE
		for(id in 1:obj$n.domains) tmp <- tmp & !any(belong.to.disk(obj$coords[indexPixel,] ,  obj$missing.domains[[id]]$center,augmentedRadius[id]^2 ))
    distant.sites[indexPixel] <- tmp
	}
  
	# OLD: listOfDistantEnoughPixels.from.BoundaryAndTheDisks
	obj$distant.sites <- matrix(distant.sites,nrow=obj$grid.size[1],ncol=obj$grid.size[2])
}

# OLD: observationSitesWithMissingDisk
non.missing.coords.grid.domain <- function(obj) obj$non.missing.coords <- obj$coords[!obj$missing.sites,]

# OLD: n
non.missing.number.grid.domain <- function(obj) obj$non.missing.number <- sum(!obj$missing.sites) #if full grid then ==obj$n1*obj$n1


expand.to.full.grid.domain  <- function(obj) {
	z  <- rep(0,obj$sites.number)
	z[!obj$missing.sites]  <- 1:obj$non.missing.number
  # arrayOfIndicesOfObservationInVectorZ
  obj$z <- matrix(z,nrow=obj$grid.size[1],ncol=obj$grid.size[2])
}

## create the sparse matrix (spam class)
preconditionning.grid.domain <- function(obj,nu=.5,range=1,precond.bandwidth=2.5)  {
  library(pdist)
  library(fields)

  # function with no arguments since all variables inside inherited from the parent environment
  # only for code more readable
  first.preconditioning <- function()  {
    iz <- obj$z[irow, icol]
  	listPi <- as.vector(pdist(	obj$non.missing.coords[iz,], obj$non.missing.coords)@dist)<(precond.bandwidth/obj$grid.size[1])    #!!!!!!! magic !!!!
    if (iz < obj$non.missing.number)  listPi[(iz +1):obj$non.missing.number] <- FALSE
    j <-  which(listPi)
    #subR <-  Matern(rdist(j),nu=0.5,range=1)
    #
    length.listPi <- length(j)
    if(length.listPi==1) {
      i <- iz
    	gi <-  1
    } else {
      #subR <-  exp(- rdist(obj$non.missing.coords[j,]))
    	subR <-  Matern(rdist(obj$non.missing.coords[j,]), nu=nu, range=range)
      vectorEmi <- c(rep(0,length.listPi -1),1)
      gTildai <-  solve(subR, vectorEmi)
      gi <-   gTildai/sqrt(gTildai[length.listPi])
    }
    list(i=iz,j=j,g=gi,l=length.listPi,irow=irow,icol=icol)
  }

  ## init sparse matrix
  sparseG <-  spam(0, obj$non.missing.number, obj$non.missing.number)
  entriesRaw<-c()
  colindicesRaw<-c()
  rowpointersRaw <-c(1)

  first.precond <- NULL #found when first condition satisfied

  for(icol in 1:obj$grid.size[2]) for(irow in 1:obj$grid.size[1])  {
  	iz <- obj$z[irow, icol]

    if  (! obj$missing.sites[ irow, icol ]) {
      #	(! iz ==0)
      #(! listOfPixels.in.OneOfTheDisks[  ])
      if (obj$distant.sites[ irow, icol ]) {
        #print(c("row=",irow,"comumn=", icol))

        if(is.null(first.precond)) first.precond <- first.preconditioning()

        j  <- obj$z[first.precond$j +  rep(( irow  - first.precond$irow) + obj$n1*(icol - first.precond$icol),  first.precond$l)]
    	  length.listPi <- length(j)
    	  i  <- rep(iz, length.listPi)
    	  gi <- first.precond$g
      } else { 	#print(iz)
  				listPi <- as.vector(pdist(	obj$non.missing.coords[iz,],
  	#			obj$non.missing.coords[max(1, iz - 2*obj$n1 - 2):min(n, iz + 2*obj$n1 + 2)])@dist)<(2.5/obj$n1)
  				obj$non.missing.coords[max(1, iz - 2*obj$n1 - 2): iz,])@dist)<(precond.bandwidth/obj$n1)     #!!!!!!! magic !!!!
  				#
  			  #if (iz <n)  listPi[(iz +1):n] <- FALSE
  			  j <-  which(listPi)
  			  j <-  j  +  max(1, iz - 2*obj$n1 - 2) -1
  			  if   (  length(j)==1  ) {
            i <- iz
  				  gi <-  1
  				  length.listPi <- 1
          } else {  #subR <-  exp(- rdist(obj$non.missing.coords[j,]))
  				      subR <-  Matern(rdist(obj$non.missing.coords[j,]), nu=nu, range=range)
                length.listPi <- length(j)
                vectorEmi <- c(rep(0,length.listPi -1),1)
                gTildai <-  solve(subR, vectorEmi)
                gi <-   gTildai/sqrt(gTildai[length.listPi])
          }
  		}
      entriesRaw <-c(entriesRaw,gi)
      colindicesRaw <-c(colindicesRaw,j)
      rowpointersRaw <-c(rowpointersRaw, rowpointersRaw[length(rowpointersRaw)] +length.listPi)
    }
  }

  sparseG@entries <- entriesRaw
  sparseG@colindices <- as.integer(colindicesRaw)
  sparseG@rowpointers <- as.integer(rowpointersRaw)

  # Alternative:
  # If rowindicesRaw instead of rowpointersRaw
  # sparseG[cbind(as.integer(colindicesRaw), as.integer(rowindicesRaw))] <-  entriesRaw
  obj$sparseG <- sparseG
}
