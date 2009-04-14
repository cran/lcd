.get.localug.m.ic <- function(freq.tb, p.value)
{
    p <- length(freq.tb@levels)
    vset <- colnames(freq.tb@table)
    vset <- vset[-length(vset)]
    sep.pairs <- c()
    if (p > 1)
        for(i in 1:(p-1))
            for(j in (i+1):p){
                cand <- vset[-c(i,j)]
                res <- .get.sep.m(freq.tb, p.value, vset[i], vset[j], cand)
                if(res$seped)
                    sep.pairs <- append(sep.pairs, res$sep)
            }
    sep.pairs
}


.get.localug.m.pc <- function(freq.tb, p.value)
{
    p <- length(freq.tb@levels)
    vset <- colnames(freq.tb@table)
    vset <- vset[-length(vset)]
    sep.pairs <- c()
        G <- matrix(rep(TRUE, p*p), p, p)
    diag(G) <- FALSE
    seq_p <- 1:p
    done <- FALSE
    ord <- 0
    n.edgetests <- numeric(1)
    while (!done && any(G)) {
        n.edgetests[ord + 1] <- 0
        done <- TRUE
        ind <- which(G, arr.ind = TRUE)
        ind <- ind[order(ind[ ,1]), ]
        remainingEdgeTests <- nrow(ind)
        for(i in 1:remainingEdgeTests) {
            x <- ind[i,1]
            y <- ind[i,2]
            if (G[y, x]) {
                nbrsBool <- G[, x]
                nbrsBool[y] <- FALSE
                nbrs <- seq_p[nbrsBool]
                length_nbrs <- length(nbrs)
                if (length_nbrs >= ord) {
                    if (length_nbrs > ord)
                        done <- FALSE
                    S <- seq(length = ord)
                    repeat {
                        n.edgetests[ord + 1] <- n.edgetests[ord + 1] + 1
                        p.val <- multinom.ci.test(freq.tb, vset[x], vset[y],
                                                  vset[nbrs[S]])$p.value
                        if (p.val > p.value) {
                            G[x, y] <- G[y, x] <- FALSE
                            pair <- new("sep.pair", u = vset[x],
                                        v = vset[y], s = vset[nbrs[S]])
                            sep.pairs <- append(sep.pairs, pair)
                            break
                        }
                        else {
                            nextSet <- .getNextSet(length_nbrs, ord, S)
                            if (nextSet$wasLast)
                                break
                            S <- nextSet$nextSet
                        }
                    }
                }
            }
        }
        ord <- ord + 1
    }
    sep.pairs
}



.get.sep.m <- function(freq.tb, p.value, u, v, cand)
{
    seped <- FALSE
    p <- length(cand)
    if (p < 8) {
        mat <- as.matrix(t(sapply(1:(2^p),
                                  function(x) .binary(x-1, max(1,p))$dicotomy)))
        if(nrow(mat) < ncol(mat)) mat <- t(mat)
        p.val <- apply(mat, 1, function(x) multinom.ci.test(freq.tb, u, v, cand[x])$p.value)
        p.val.max <- max(p.val)
        idx <- which(p.val == p.val.max)
        sep <- new("sep.pair",
                   u=u,
                   v=v,
                   s=character(0))
        if(p.val.max >= p.value){
            sep@s <- cand[mat[idx[1],]]
            seped <- TRUE
        }
        return(list(seped=seped, sep=sep))
    } else {
        .get.sep.m.stepwise(freq.tb, p.value, u, v, cand)
    }
}


.get.sep.m.stepwise <- function(freq.tb, p.value, u, v, cand){
     .get.sep.m.step(freq.tb, p.value, u, v, cand, c())
 }

.get.sep.m.step <- function(freq.tb, p.value, u, v, current, rest)
{
    modified <- FALSE
    sep <- new("sep.pair",
               u = u,
               v = v,
               s = character(0))
    pp <- multinom.ci.test(freq.tb, u, v, current)$p.value
    if(pp >= p.value) {
        sep@s = current
        return(list(seped = TRUE, sep = sep))
    }
    if(length(current) == 0)
        return(list(seped = FALSE, sep = sep))
    n.curr <- length(current)
    n.rest <- length(rest)
    mat <- matrix(rep(current, n.curr), n.curr, byrow = TRUE)
    diag(mat) <- NA
    pval <- apply(mat, 1, function(x) multinom.ci.test(freq.tb, u, v, x[!is.na(x)])$p.value)
    todel <- which(pval == max(pval))
    if(length(todel) > 1)
        todel <- todel[sample(length(todel),1)]
    if(pval[todel] >= pp){
        todel <- current[todel]
        current <- setdiff(current, todel)
        rest <- union(rest, todel)
        pp <- max(pval)
        modified <- TRUE
    }
    if(pp >= p.value){
        sep@s = current
        return(list(seped = TRUE, sep = sep))
    }
    if(length(rest) == 0)
        return(list(seped = FALSE, sep = sep))
    if(modified){
        n.curr <- length(current)
        n.rest <- length(rest)
    }
    mat <- matrix(NA, n.rest, n.curr + 1)
    if(n.curr != 0)
        mat[,1:n.curr] <- matrix(rep(current, n.rest), n.rest, byrow = TRUE)
    mat[,n.curr + 1] <- rest
    pval <- apply(mat, 1, function(x) multinom.ci.test(freq.tb, u, v, x)$p.value)
    toadd <- which(pval == max(pval))
    if(length(toadd) > 1)
        toadd <- toadd[sample(length(toadd),1)]
    if(pval[toadd] > pp){
        toadd <- rest[toadd]
        current <- union(current, toadd)
        rest <- setdiff(rest, toadd)
        modified <- TRUE
    }
    if(modified)
        return(.get.sep.m.step(freq.tb, p.value, u, v, current, rest))
    else 
        return(list(seped = FALSE, sep = sep))
}


`learn.skeleton.multinom` <- function(tree, freq.tb, p.value, drop = TRUE)
{
    validObject(tree)
    validObject(freq.tb)
    local.ug <- c()
    vset <- colnames(freq.tb@table)
    vset <- vset[-length(vset)]
    n.clique <- length(tree@cliques)
    for(i in 1:n.clique){
        idx <- tree@cliques[[i]]$vset
#        if (length(idx) >= 10)
            new.ug <- .get.localug.m.pc(compress.freq.tb(freq.tb, idx), p.value)
#        else
#            new.ug <- .get.localug.m.ic(compress.freq.tb(freq.tb, idx), p.value)
        local.ug <- append(local.ug, new.ug)
    }
    p <- length(vset)
    amat <- matrix(0, p, p)
    rownames(amat) <- colnames(amat) <- vset
    n.clique <- length(tree@cliques)
    for(i in 1:n.clique){
        idx <- tree@cliques[[i]]$vset
        amat[idx, idx] <- 1
    }
    diag(amat) <- 0
    sep.pairs <- c()
    n.loc.sep <- length(local.ug)
    if(n.loc.sep>0)
        for(i in 1:n.loc.sep){
            u <- local.ug[[i]]@u
            v <- local.ug[[i]]@v
            if(amat[u,v] == 1){
                amat[u,v] <- amat[v,u] <- 0
                sep.pairs <- append(sep.pairs, local.ug[i])
            }
        }


    ## code partially adapted from the "pcAlgo" function from the
    ## "pcalg" package
    
    if (drop) {
        ind <- .get.exed.cand1(tree, amat)
#        if (any(ind) && nrow(ind) > 50) {
        if (any(ind)) {
            ind <- ind[order(ind[,1]),]
            ord <- 0
            seq_p <- 1:p
            done <- FALSE
            remainingEdgeTests <- nrow(ind)
            while (!done && any(as.logical(amat))) {
                done <- TRUE
                for (i in 1:remainingEdgeTests) {
                    x <- ind[i, 1]
                    y <- ind[i, 2]
                    if (amat[y, x]) {
                        nbrsBool <- amat[, x] == 1
                        nbrsBool[y] <- FALSE
                        nbrs <- seq_p[nbrsBool]
                        length_nbrs <- length(nbrs)
                        if (length_nbrs >= ord) {
                            if (length_nbrs > ord)
                                done <- FALSE
                            S <- seq(length = ord)
                            repeat {
                                p.val <- multinom.ci.test(freq.tb, vset[x], vset[y],
                                                      vset[nbrs[S]])$p.value
                                if (p.val > p.value) {
                                    amat[x, y] <- amat[y, x] <- 0
                                    pair <- new("sep.pair", u = vset[x],
                                                v = vset[y], s = vset[nbrs[S]])
                                    sep.pairs <- append(sep.pairs, pair)
                                    break
                                }
                                else {
                                    nextSet <- .getNextSet(length_nbrs, ord, S)
                                    if (nextSet$wasLast)
                                        break
                                    S <- nextSet$nextSet
                                }
                            }
                        }
                    }
                }
                ord <- ord + 1
            }
        }
##         } else {
##             if (any(ind)) {
##                 for(i in 1:nrow(ind)){
##                     pair <- new("sep.pair", u = vset[ind[i,1]],
##                                 v = vset[ind[i,2]], s = character(0))
##                     cand <- setdiff(vset[amat[pair@u,]==1], pair@v)
##                     idx <- c(pair@u, pair@v, cand)
##                     res <- .get.sep.m(compress.freq.tb(freq.tb, idx), p.value, pair@u, pair@v, cand)
##                     if(res$seped){
##                         amat[pair@u, pair@v] <- amat[pair@v, pair@u] <- 0
##                         sep.pairs <- append(sep.pairs, res$sep)
##                     } 
##                 }
##             }
##         }
    }
    return(list(amat=amat, sep.pairs=sep.pairs))    
}

