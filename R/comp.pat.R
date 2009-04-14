`comp.pat` <- function(truepat, pat)
{
    vset <- rownames(truepat)
    truearr <- which(truepat - t(truepat) == 1)
    pat <- pat[vset, vset]
    arr <- which(pat - t(pat) == 1)
    a.missing <- length(truearr)-length(which(match(truearr, arr)>0))
    a.extra <- length(arr)-length(which(match(arr, truearr)>0))
    
    ## computing structural Hamming distance

    idx <- which(skeleton(truepat) - skeleton(pat) != 0)
    pat[idx] <- skeleton(truepat)[idx]
    shd <- length(idx) / 2
    c <- abs(truepat-pat)
    shd <- shd + length(which(c+t(c) != 0)) / 2
    return(list(a.total = length(truearr),
                a.missing = a.missing,
                a.extra = a.extra,
                shd = shd))
}
