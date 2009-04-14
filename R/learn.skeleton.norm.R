`learn.skeleton.norm` <- function(tree, cov, n, p.value, drop = TRUE)
{
    validObject(tree)
    local.ug <- c()
    vset <- rownames(cov)
    n.clique <- length(tree@cliques)
    for(i in 1:n.clique){
        idx <- tree@cliques[[i]]$vset
#        if (length(idx) >= 10)
            new.ug <- .get.localug.pc(cov[idx, idx], n, p.value)
#        else
#            new.ug <- .get.localug.ic(cov[idx, idx], n, p.value)
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
    
    ## the following code is partially adapted from the "pcAlgo" function
    ## from "pcalg" package in R
    
    if (drop) {
        ind <- .get.exed.cand1(tree, amat)
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
                                p.val <- norm.ci.test(cov, n, vset[x], vset[y],
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
##                     res <- .get.sep(cov[idx, idx], n, p.value, pair@u, pair@v, cand)
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

