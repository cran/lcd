`rmultinom.cg` <-
    function(n, amat, distn)
{
    vset <- rownames(amat)
    p <- length(vset)
    data <- matrix(NA, n, p)
    colnames(data) <- vset
    n.c <- length(distn)
    for(i in 1:n.c){
        temp <- distn[[i]]
        V.idx <- match(colnames(temp$Vtable), vset)
        if(ncol(temp$cond.dist) > 1){
            B.idx <- match(colnames(temp$Btable), vset)
            if (length(B.idx) > 1) {
              config <- apply(as.matrix(data[,B.idx]), 1, function(x)
                              which(apply(temp$Btable, 1, function(y) identical(x,y))))
            } else {
              config <- apply(as.matrix(data[,B.idx]), 1, function(x)
                              which(apply(temp$Btable, 1, function(y) y==x)))
            }
            row <- sapply(config, function(x)
                          sample(nrow(temp$Vtable), 1, prob = temp$cond.dist[,x]))
            data[,V.idx] <- temp$Vtable[row,]
        } else {
            data[,V.idx] <- temp$Vtable[sample(nrow(temp$Vtable), n, TRUE,
                                              prob = temp$cond.dist[,1]),]
        }
    }
    data
}

