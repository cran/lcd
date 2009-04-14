`draw` <-
    function(amat, plain = TRUE){
    e.total <- length(which(skeleton(amat) == 1))/2
    if(e.total == 0)
        return("The graph has no edge!")
    el <- matrix("", nrow = e.total, ncol = 2)
    vset <- sort(rownames(amat))
    amat <- amat[vset,vset]
    p <- nrow(amat)
    arrow <- c()
    k <- 1
    for(i in 1:(p-1))
        for(j in (i+1):p){
            if(amat[i,j] == 1){
                el[k,1] <- i-1
                el[k,2] <- j-1
                k <- k + 1
                if(amat[j,i] == 0){
                    arrow <- c(arrow, ">")
                } else {
                    arrow <- c(arrow, "-")
                }
            } else {
                if(amat[j,i] == 1){
                    el[k,1] <- j-1
                    el[k,2] <- i-1
                    k <- k + 1
                    arrow <- c(arrow, ">")
                }
            }
        }
    g <- graph.empty()
    g <- add.vertices(g, p)
    V(g)$name <- vset
    g <- add.edges(g, t(el))
    if(plain){
        plot.igraph(g, layout = layout.reingold.tilford,
                    edge.arrow.mode = arrow,
                    vertex.size = min(15, 500/p),
                    edge.arrow.size = min(1, 30/p))
    } else{
        tkplot(g, layout = layout.reingold.tilford,
               edge.arrow.mode = arrow,
               vertex.size = min(15, 500/p),
               edge.arrow.size = min(1, 30/p))
    }
    V(g)
}
