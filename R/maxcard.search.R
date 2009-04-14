`maxcard.search` <-
function(amat)
{
    output <- TRUE
    k <- nrow(amat)
    L <- c()
    V <- 1:k
    C <- rep(0,k)
    card <- pi_record <- c()
    for(j in 1:k){
        U <- setdiff(V, L)
        pos <- which(C[U] == max(C[U]))
        if(length(pos) > 1){
            v <- U[sample(pos,1)]
        } else {
            v <- U[pos]
        }
        nb <- which(amat[v,]==1)
        pi_v <- intersect(nb, L)
        m <- length(pi_v)
        card <- c(card, m)
        pi_record <- c(pi_record, pi_v)
        ne <- sum(amat[pi_v, pi_v])/2
        if(ne < m*(m-1)/2){
            output <- FALSE
            break
        }
        nbU <- intersect(nb, U)
        C[nbU] <- C[nbU] + 1
        L <- c(L, v)
    }
    return(list(is.triangulated = output,
                perfect.numbering = L,
                card = card,
                pi.record = pi_record))
}
