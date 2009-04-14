 
## `learn.complex.multinom` <- function(skel, tree, freq.tb, p.value)
## {
##     validObject(tree)
##     validObject(freq.tb)
##     wmat <- skel$amat
##     vset <- rownames(wmat)
##     p <- nrow(wmat)
##     sep.pairs <- skel$sep.pairs
##     cliques <- tree@cliques
##     n.clique <- length(tree@cliques)
##     for(i in 1:p)
##         for(j in 1:p){
##             if(wmat[i,j] != 0){
##                 u <- vset[i]
##                 w <- vset[j]
##                 nb.u <- vset[skel$amat[u,]==1]
##                 v.cand <- c()
##                 for(k in 1:n.clique){
##                     if(is.element(u, cliques[[k]]$vset))
##                         v.cand <- append(v.cand, cliques[[k]]$vset)
##                 }
##                 v.cand <- unique(v.cand)
##                 v.cand <- setdiff(v.cand, c(nb.u, u))
##                 n.cand <- length(v.cand)
##                 if(n.cand > 0)
##                     for(l in 1:n.cand){
##                         pair <- new("sep.pair", u=u, v=v.cand[l])
##                         a <- lapply(sep.pairs, function(x) all.equal(x, pair))
##                         sep <- sep.pairs[[which(a==TRUE)[1]]]@s
##                         sep <- unique(append(sep, w))
##                         idx <- c(u, v.cand[l], sep)
##                         res <- multinom.ci.test(compress.freq.tb(freq.tb, idx), u, v.cand[l], sep)
##                         ## cat(idx, val, '\n')
##                         if(res$p.value < p.value && (-1-res$deviance) < wmat[w,u]){
##                             wmat[w,u] <- -1 - res$deviance
##                         }
##                     }
##             }
##         }
##     idx <- which(wmat - t(wmat) < 0)
##     wmat[idx] <- 0
##     wmat <- 0 + (wmat != 0)
##     pattern(wmat)
## }


## function for complex recovery in chain graphs
## last modified: @@ Fri, Apr 25, 2008, 10:53 @@


`learn.complex.multinom` <- function(skel, freq.tb, p.value)
{
    validObject(freq.tb)
    wmat <- skel$amat
    vset <- rownames(wmat)
    sep.pairs <- skel$sep.pairs
    n.sep <- length(sep.pairs)
    if (n.sep == 0) return(wmat)
    for (i in 1:n.sep) {
        pair <- sep.pairs[[i]]
        for (turn in 1:2) {
            u <- if(turn == 1) pair@u else pair@v
            v <- if(turn == 1) pair@v else pair@u
            sep <- pair@s
            nb.u <- vset[skel$amat[u,] == 1]
            nb.u.size <- length(nb.u)
            if (nb.u.size > 0) {
                for (j in 1:nb.u.size) {
                    w <- nb.u[j]
                    newsep <- unique(append(sep, w))
                    idx <- c(u, v, newsep)
                    res <- multinom.ci.test(freq.tb,
                                            u, v, newsep)
                    if (res$p.value < p.value &&
                        (-1 - res$deviance) < wmat[w,u]) {
                        wmat[w, u] <-  -1 - res$deviance
                    }
                }
            }
        }
    }
    idx <- which(wmat - t(wmat) < 0)
    wmat[idx] <- 0
    wmat <- 0 + (wmat != 0)
    pattern(wmat)
}

