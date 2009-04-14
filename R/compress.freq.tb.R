 
`compress.freq.tb` <- function(tb, subset = colnames(tb@table)[-ncol(tb@table)])
{
    vnames <- colnames(tb@table)
    last <- length(vnames)
    idx <- match(subset, vnames[-last])
    if(any(is.na(idx)))
        stop("Invalid subset!")
    newlev <- tb@levels[idx]
    idx <- c(idx, last)
    newtb <- .Call("compress_freq_table",
                   tb@table[,idx],
                   PACKAGE = "lcd")
    colnames(newtb) <- vnames[idx]
    new("freq.tb", table = newtb, levels = newlev)
}
