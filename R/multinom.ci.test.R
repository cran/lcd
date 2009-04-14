 `multinom.ci.test` <- function(tb, u, v, cond = c())
{
    if(class(tb) != "freq.tb")
        stop("tb should be an object of class freq.tb!")
    if(!all(is.na(match(c(u,v), cond))))
        stop("The conditional set should not contain u or v!")
    if (nrow(tb@table) == 1) wtb <- tb
    else wtb <- compress.freq.tb(tb, c(u, v, cond))
    res <- .Call("cond_ind_test",
                 as.matrix(wtb@table),
                 as.vector(wtb@levels),
                 as.integer(0),
                 as.integer(1),
                 PACKAGE = "lcd")
    return(list(deviance = res[1], df = res[2], p.value = res[3]))
}
