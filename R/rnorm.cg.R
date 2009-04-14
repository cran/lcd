`rnorm.cg` <-
function(n, amat, Bstar)
{
    fcheck <- is.chaingraph(amat)
    if(!fcheck$result)
        stop("The input should be a chain graph!")
    p <- nrow(amat)
    corder <- rev(fcheck$vert.order)
    csize <- rev(fcheck$chain.size)
    diag(amat) <- 1
    if(!all(rownames(Bstar) == rownames(amat)[corder]) ||
       !all((t(Bstar) == 0) == (amat[corder, corder] == 0)))
        stop("The Bstar and amat are not compatible!")
    b <- cumsum(csize)
    W <- matrix(nrow = p, ncol = n)
    for(k in 1:length(csize)){
        idx <- (max(b[k-1],0)+1):b[k]
        if(csize[k] > 1)
            W[idx, ] <- t(mvrnorm(n=n, mu=rep(0,csize[k]), Sigma=Bstar[idx, idx]))
        if(csize[k] == 1)
            W[idx, ] <- rnorm(n=n, mean=0, sd=sqrt(Bstar[idx, idx]))
    }
    t(solve(Bstar, W))
}

