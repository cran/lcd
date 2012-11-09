`is.chaingraph` <-
function(amat)
{
    wmat <- matrix(as.integer((amat + t(amat)) > 1), nrow = nrow(amat))
    wg <- graph.adjacency(wmat, mode = "undirected")
    cc <- clusters(wg)
    neworder <- order(cc$membership)
    a <- matrix(0, nrow = length(cc$csize), ncol = length(cc$csize))
    b <- cumsum(cc$csize)
    wmat <- amat[neworder, neworder]
    for(i in 1: length(cc$csize)){
        for(j in 1: length(cc$csize)){
            if(j != i){
                a[i,j] <- as.integer(sum(wmat[(max(b[i-1],0)+1):b[i],
                                              (max(b[j-1],0)+1):b[j]]) > 0)
            }
        }
    }
    rownames(a) <- colnames(a) <- as.character(1:length(b))
    output <- isAcyclic(a)
    for(i in 1:length(b)){
        temp <- wmat[(max(b[i-1],0)+1):b[i], (max(b[i-1],0)+1):b[i]]
        if (!all(temp == t(temp))) {
            output <- FALSE
            break
        }
    }
    chainorder <- topOrder(a)
    vertorder<-c()
    chainsize<-c()
    if(output == TRUE){
        for(k in 1:length(b)){
            ## vertorder <- c(vertorder, which(cc$membership == chainorder[k]-1))
            vertorder <- c(vertorder, which(cc$membership == chainorder[k]))
            chainsize <- c(chainsize, cc$csize[chainorder[k]])
        }
    }
    return(list(result = output,
                vert.order = vertorder,
                chain.size = chainsize))
}

