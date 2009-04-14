## function for complex recovery in chain graphs
## Last modified: @@ Fri, Apr 25, 2008, 10:18 @@

`learn.complex.norm` <- function(skel, cov, n, p.value)
{
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
                    res <- norm.ci.test(cov[idx, idx], n, u, v, newsep)
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
