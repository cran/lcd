is.separated <- function(u, v, sep, amat)
{
    vset <- rownames(amat)
    curr <- c(u, v, sep)
    s <- length(curr)
    q <- curr
    if(s > 0) {
        for(i in 1:s)
            q <- union(q,
                       vset[sapply(vset, function(x)
                                   .is.linked(amat, curr[i], x))])
    }
    curr <- union(curr, q)
    while(TRUE){
        cc <- pa(curr, amat)
        cc <- setdiff(cc, curr)
        scc <- length(cc)
        if (scc > 0) {
            cccliq <- cc
            for(j in 1:scc)
                cccliq <- union(cccliq,
                                vset[sapply(vset, function(x)
                                            .is.linked(amat, cc[j], x))])
            curr <- union(curr, cccliq)
        } else {
            break
        }
    }
    anamat <- amat[curr, curr]
    anamat <- moralize(anamat)
    idx <- match(sep, rownames(anamat))
    !.is.linked(anamat[-idx, -idx], u, v)
}


