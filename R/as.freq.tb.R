`as.freq.tb` <- function(mat)
{
    n <- nrow(mat)
    p <- ncol(mat)
    vnames <- colnames(mat)
    mat <- matrix(as.integer(mat), n, p)
    levels <- apply(mat, 2, function(x) length(unique(x)))
    tb <- new("freq.tb",
              table = cbind(mat, as.integer(rep(1, n))),
              levels = levels
              )
    colnames(tb@table) <- c(vnames, "Freq")
    compress.freq.tb(tb)
}
