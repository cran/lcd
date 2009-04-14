`learn.v` <-
    function(skel, tree)
{
    validObject(tree)
    wmat <- skel$amat
    vset <- rownames(wmat)
    sep.pairs <- skel$sep.pairs
    cliques <- tree@cliques
    n.clique <- length(tree@cliques)
    for(i in 1:n.clique){
        cvset <- cliques[[i]]$vset
        p <- length(cvset)
        if(p > 2){
            for(j in 1:(p-1))
                for(k in (j+1):p)
                    for(l in 1:p){
                        u <- cvset[j]
                        v <- cvset[k]
                        w <- cvset[l]
                        if(skel$amat[u,v] == 0 &&
                           wmat[u,w] == 1 && wmat[v,w] == 1 &&
                           wmat[w,u] + wmat[w,v] != 0){
                            pair <- new("sep.pair", u=u, v=v)
                            a <- lapply(sep.pairs, function(x) all.equal(x,pair))
                            sep <- sep.pairs[[which(a==TRUE)[1]]]@s
                            if(!is.element(w, sep))
                                wmat[w,u] <- wmat[w,v] <- 0
                        }
                    }
        }
    }
    wmat
}
