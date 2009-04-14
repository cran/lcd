`norm.ci.test` <- function(cov, n, u, v, cond = c())
{
    vset <- rownames(cov)
    if(is.character(u)) u <- which(vset == u)
    if(is.character(v)) v <- which(vset == v)
    if(is.character(cond)) cond <- match(cond, vset)
    if(!is.numeric(u) || !is.numeric(v) || (!is.numeric(cond) && length(cond) > 0))
        stop("Invalid vertex!")
    res <- .Call("g_ci_test", cov,
                 as.integer(n),
                 as.integer(u), as.integer(v), as.integer(cond),
                 PACKAGE = "lcd")
    if (is.na(res[1])) {
        res[3] <- 0
        res[1] <- qchisq(1 - 1e-16, 1)
    }
    return(list(deviance = res[1], df = res[2], p.value = res[3]))
}
