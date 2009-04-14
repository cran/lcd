`get.normal.dist` <- function(amat)
{
    fcheck <- is.chaingraph(amat)
    if(!fcheck$result)
        stop("The input should be a chain graph!")
    p <- nrow(amat)
    corder <- rev(fcheck$vert.order)
    csize <- rev(fcheck$chain.size)
    b <- cumsum(csize)
    if(length( rownames(amat) ) != p)
        rownames(amat) <- colnames(amat) <- as.character(1:nrow(amat))
    if(p>1)
        amat <- amat[corder,corder]
    Bstar <- matrix(nrow = p, ncol = p)
    rownames(Bstar) <- colnames(Bstar) <- rownames(amat)
    admissible <- FALSE
    niter <- 0
    ## repeat till we get the block diagonal of Bstar all pos-def!
    ## modified: @@ Tue, Apr 29, 2008, 17:17 @@
    while (!admissible) {
        admissible <- TRUE
        niter <- niter + 1
        for(k in 1:length(csize)){
            aa <- (p - max(b[k-1], 0))
            bb <- b[k] - max(b[k-1],0)
            randdist <- matrix(runif((aa*aa), 0.5/bb, 1.5/bb)*
                               sample(c(-1,1), aa*aa, replace = TRUE),
                               nrow = aa, ncol = aa)
            randdist[upper.tri(randdist)] <- t(randdist)[upper.tri(t(randdist))]
            diag(randdist) <- 1
            Bstar[(max(b[k-1],0)+1):p, (max(b[k-1],0)+1):p]  <-randdist
            if(k != length(csize))
                Bstar[(b[k]+1):p, (max(b[k-1],0)+1):b[k]] <- 0
        }
        diag(amat) <- rep(1, p)
        Bstar <- Bstar*t(amat)
        for (k in 1:length(csize)) {
            idx <- (max(b[k-1],0)+1):b[k]
            if (!all(Bstar[idx, idx] == t(Bstar[idx,idx])))
                warning("Block diagonal not sysmetric!")
            if (length(idx) > 1 && eigen(Bstar[idx,idx],
                                         only.values = TRUE)$values[length(idx)] < 0) {
                admissible <- FALSE
                break
            }
        }
    }
    cat("Number of iteration = ", niter, "\n")
    Bstar
}
