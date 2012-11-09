setClass("sep.tree",
         representation(tree.struct = "matrix",
                        cliques = "list",
                        separators = "list"))

.is.tree <- function(amat){
    if(nrow(amat) <= 1){
        return(TRUE)
     } else {
        deg <- apply(amat, 1, sum)
        a <- deg == 1
        if(!any(a)) return(FALSE)
        else{
            leaf <- which(a)[1]
            return(.is.tree(as.matrix(amat[-leaf,-leaf])))
        }
    }
}

.valid.tree <- function(object){
    n.clique <- length(object@cliques)
    n.sep <- length(object@separators)
    if(n.clique != n.sep + 1) return("Cliques and separators don't match!")
    if(n.clique > 0)
        for(i in 1:n.clique){
            if(!all(names(object@cliques[[i]]) == c("name", "vset")) ||
               length(names(object@cliques[[i]])) != 2)
                return(paste("The names in clique", i, "is invalid!"))
        }
    if(n.sep > 0)
        for(i in 1:n.sep){
            if(!all(names(object@separators[[i]]) == c("separator", "edge")) ||
               length(names(object@separators[[i]])) != 2)
                return(paste("The name in separator", i, "is invalid!"))
        }
    if(n.sep > 0)
        for(i in 1:n.sep){
            sep <- object@separators[[i]]
            if(length(sep$edge) != 2)
                return(paste("The edge attribute in separator", i, "is invalid!"))
            u <- sep$edge[1]
            v <- sep$edge[2]
            if(object@tree.struct[u,v] != 1)
                return(paste("Separator", i, "doesn't match the tree structure!"))
            sepp <- sep$separator
            btr <- object@tree.struct
            btr[u,v] <- btr[v,u] <- 0
            bt <- graph.adjacency(btr, "undirected")
            cc <- clusters(bt)
            l.vset <- r.vset <- c()
            l <- r <- c()
            for(j in 1:(n.sep+1)){
                ## if(cc$membership[j]==0){
              if(cc$membership[j]==1){
                    l.vset <- c(l.vset, object@cliques[[j]]$vset)
                    l <- c(l, j)
                } else {
                    r.vset <- c(r.vset, object@cliques[[j]]$vset)
                    r <- c(r, j)
                }
            }
            l.vset <- unique(l.vset)
            r.vset <- unique(r.vset)
            if(!setequal(intersect(l.vset, r.vset), sepp))
                return("Separator definition doesn't match that of separation tree!")
        }
    TRUE
}

setValidity("sep.tree", .valid.tree)
