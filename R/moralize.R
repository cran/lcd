`moralize` <- function(amat)
{
    p <- nrow(amat)
    vset <- rownames(amat)
    pos <- which(amat == 1 & t(amat) == 0)
    s <- d <- c()
    n.arrow <- length(pos)
    if (n.arrow == 0) 
        return(amat)
    for (i in 1:n.arrow) {
        d <- c(d, ceiling(pos[i]/p))
        s <- c(s, pos[i] - (ceiling(pos[i]/p) - 1) * p)
    }
    if (n.arrow == 1) {
        amat[d[1], s[1]] <- amat[s[1], d[1]] <- 1
        return(amat)
    }
    a <- outer(rep(1, length(s)),1:length(s))
    b <- t(a)
    a <- as.numeric(a[upper.tri(a)])
    b <- as.numeric(b[upper.tri(b)])
    pa <- cbind(a,b)
    pair <- cbind(s[pa[,1]], d[pa[,1]], s[pa[,2]], d[pa[,2]])
    mat <- cbind(pair,
                 matrix(amat[pair[,1],], nrow = nrow(pair)),
                 matrix(amat[pair[,3],], nrow = nrow(pair)))
    rownames(mat) <- NULL
    l <- (ncol(mat)-4)/2
    cidx <- ncol(mat)*(0:(nrow(mat)-1))
    idx1 <- 4 + mat[,2] + cidx
    idx2 <- 4 + l + mat[,4] + cidx
    mat <- t(mat)
    mat[idx1] <- mat[idx2] <- 0
    check1 <- mat[1,] != mat[2,] & mat[3,] != mat[4,]
    idx3 <- 4 + mat[4,] + cidx
    idx4 <- 4 + l + mat[2,] + cidx
    check2 <- !(mat[idx3] | mat[idx4])
    idx5 <- mat[1,] + p*(mat[3,]-1)
    idx6 <- mat[3,] + p*(mat[1,]-1)
    check3 <- !(amat[idx5] | amat[idx6])
    .is.complex <- function(x){
        f <- !(x[5:(4+l)] | x[(5+l):(4+2*l)])
        .is.linked(amat[f,f], vset[x[2]], vset[x[4]])
    }
    status <- apply(mat, 2, .is.complex) & check1 & check2 & check3
##     idx7 <- mat[1,status] + p*(mat[3,status]-1)
##     idx8 <- mat[3,status] + p*(mat[1,status]-1)
##     amat[idx7] <- 1
##     amat[idx8] <- 1
    idx <- pair[status, c(1,3)]
    if (any(nrow(idx)))
        apply(idx, 1, function(x) amat[x[1],x[2]] <- amat[x[2],x[1]] <- 1)
    else
        amat[idx[1], idx[2]] <- amat[idx[2], idx[1]] <- 1
    skeleton(amat)
}
