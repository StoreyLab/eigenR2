eigenR2 <- function(dat, model, null.model=NULL, adjust=FALSE, eigen.sig=NULL, mod.fit.func=NULL){  

    if (!is.matrix(dat)) stop("Please input a data matrix!")
    
	nn <- ncol(dat)
	if (!is.null(mod.fit.func)){
		null.model <- NULL
	} else {
	  df <- ncol(model)
      df0 <- max(ncol(null.model), 1)
	}
	
	if (!is.null(null.model)) {
	    ## check if the model is nested in the null.model
	    idx.col <- NULL
		if (sum(null.model[, 1]==1) == nn) {
            sidx <- 2 
        } else {
            sidx <- 1
        }
		for (i in sidx:ncol(null.model)){
		    oo <- apply(model, 2, function(x) sum(x == null.model[, i]))
	        ooi <- which((oo[sidx:ncol(model)]) == nn)
            if (length(ooi) == 0) {
                stop("The model has to be nested in the null model")
            } else {
                idx.col <- c(idx.col, ooi)
            }
	    }
		rmodel <- matrix(model[, ((sidx: df)[-idx.col])], nrow= nn)
		H0 <- null.model %*% solve(t(null.model) %*% null.model) %*% t(null.model) 
	    rmodel <- rmodel - H0 %*% rmodel
		if (sidx == 2) rmodel <- model.matrix(~1+rmodel)
		model <- rmodel
	} else{
	    null.model <- matrix(rep(1, nn), ncol=1)
		H0 <- null.model %*% solve(t(null.model) %*% null.model) %*% t(null.model) 
	}
	dat <- dat - t(H0 %*% t(dat))


    svd.t <- fast.svd(dat, tol=0)
    eigenVec <- t(svd.t$v)  ## each column represents a right eigenvector
    N <- nrow(eigenVec)
    eigenV <- eigenVec[-N, ]
    ds <- svd.t$d[-N]
    weights <- ds^2 / sum(ds^2)
    ww <- c(weights, 0)

    if (!is.null(eigen.sig)) {
        oov <- sva.id(dat, mod=NULL, eigen.sig=eigen.sig)
		n.sv <- oov$n.sv
        if (n.sv>0) {
            eigenV <- matrix(eigenV[1:n.sv, ], nrow=n.sv)
            weights <- weights[1:n.sv]
        } else {
            return(list(eigenR2=0, weights=ww, eg.R2s=NULL, p.eg=oov$p.sv))
        }
    }


    if (!is.null(mod.fit.func)){
        R2s <- apply(eigenV, 1, mod.fit.func)
        v <- sum(weights*R2s)
    } else {  
        H <- model %*% solve(t(model) %*% model) %*% t(model) 
        res <- eigenV - t(H %*% t(eigenV))
        SSR <- apply(res, 1, function(x) sum((x-mean(x))^2))
        SST <- apply(eigenV, 1, function(x) sum((x-mean(x))^2))
        R2s <- 1- SSR/SST
        v <- sum(R2s*weights)
		if (adjust == TRUE) {
			v <- 1-(1-v)*(N-df0)/(N-df)        
		}
    }


    if (!is.null(eigen.sig)) {
        return(list(eigenR2=v, weights=ww, eg.R2s=R2s, p.eg=oov$p.sv))
    } else {
        return(list(eigenR2=v, weights=ww, eg.R2s=R2s, p.eg=NULL))
    }
}


## the function sva.id is a modified function from "sva" package by Leek and Storey.
sva.id <- function(dat, mod=NULL, B=20, eigen.sig=0.1, seed=NULL) {

    if(!is.null(seed)){set.seed(seed)}
    warn <- NULL
    n <- ncol(dat)
    m <- nrow(dat)
    if (is.null(mod)){
        res <- t(scale(t(dat),scale=F))
        ndf <- n-1
    } else {
        H <- mod %*% solve(t(mod) %*% mod) %*% t(mod) 
        res <- dat - t(H %*% t(dat))
        ndf <- n - ceiling(sum(diag(H)))
    }
    uu <- fast.svd(res,tol=0)
    dstat <- uu$d[1:ndf]^2/sum(uu$d[1:ndf]^2)
	dstat0 <- matrix(0,nrow=B,ncol=ndf)
	
	for(i in 1:B){
        res0 <- t(apply(res, 1, sample, replace=FALSE))
        if(is.null(mod)){
            res0 <- t(scale(t(res0),scale=F))
        } else {
            res0 <- res0 - t(H %*% t(res0))
        }
        uu0 <- fast.svd(res0, tol=0)
        dstat0[i,] <- uu0$d[1:ndf]^2/sum(uu0$d[1:ndf]^2)
	}
	psv <- rep(1,n)
	for(i in 1:ndf){
	    psv[i] <- mean(dstat0[,i] >= dstat[i])
	}
	for(i in 2:ndf){
	    psv[i] <- max(psv[(i-1)],psv[i]) 
	}
	
    nsv <- sum(psv <= eigen.sig)
    return(list(n.sv = nsv,p.sv=psv))
}



eigenR2plot <- function(eigenR2obj){
    N <- length(eigenR2obj$weights)
    nplot <- 1
    if (!is.null(eigenR2obj$p.eg)) nplot <- nplot+1
    if (!is.null(eigenR2obj$eg.R2s)) nplot <- nplot+1
    par(mfrow=c(round(nplot/2),2))
    plot(1:N, eigenR2obj$weights, xlab="Eigenvectors", ylab="Proportion of Variation Explained", main="Variation Dissection by Eigenvector")
    
    if (!is.null(eigenR2obj$eg.R2s)) {
        plot(1:N, 1:N, ylim=c(0, max(eigenR2obj$eg.R2s)), "n", xlab="Eigenvectors", ylab="R-square", main="R-square for Significant Eigenvectors")
        n.sv <- length(eigenR2obj$eg.R2s)
        points(1:n.sv, eigenR2obj$eg.R2s)
    }
    if (!is.null(eigenR2obj$p.eg)){
        plot(1:N, eigenR2obj$p.eg, xlab="Eigenvectors", ylab="p-values", main="P-value Plot for Eigenvectors")
    }
    par(mfrow=c(1,1))
}

plot.eigenR2 <- function(x,...){eigenR2plot(x,...)}

attr(plot.eigenR2, "source") <- NULL
attr(eigenR2, "source") <- NULL
attr(sva.id, "source") <- NULL


