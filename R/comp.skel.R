`comp.skel` <- function(trueskel, skel)
{
    vset <- rownames(trueskel)
    a <- trueskel - skel[vset, vset]
    e.missing <- length(which(a == 1))/2
    e.extra <- length(which(a == -1))/2
    return(list(e.total = length(which(trueskel == 1))/2,
                e.missing = e.missing,
                e.extra = e.extra))
}
