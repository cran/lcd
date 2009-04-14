`tri.ug` <-
function(amat)
{
    k <- nrow(amat)
    unlab <- rep(TRUE, k)
    numorder <- c()
    counter <- k
    while(counter > 0){
        cc <- sample(which(unlab))
        minToFill <- k^2
        toLabel <- k+1
        Nb <- c()
        nNb <- 0
        for(j in 1:counter){
            cNb <- which(amat[cc[j],]==1 & unlab)
            cNb <- union(cNb, cc[j])
            ncNb <- length(cNb)
            toFill <- ncNb*(ncNb-1)/2 -
                sum(amat[cNb, cNb])/2
            if(toFill < minToFill){
                toLabel <- cc[j]
                Nb <- cNb
                minToFill <- toFill
                nNb <- ncNb
            }
        }
        if(nNb > 0 && minToFill > 0){
            newm <- matrix(1, nrow=nNb, ncol=nNb)
            diag(newm) <- 0
            amat[Nb, Nb] <- newm
        }
        numorder <- c(numorder, toLabel)
        unlab[toLabel] <- FALSE
        counter <- counter - 1
    }
    amat
}

