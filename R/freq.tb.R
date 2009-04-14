setClass("freq.tb",
         representation(table = "matrix",
                        levels = "integer"))

.valid.freq.tb <- function(object){
    n.var <- ncol(object@table) - 1
    if(n.var < 1)
        return("Table should contain at least one variable!")
    if(n.var != length(object@levels))
        return("Table and levels do not match!")
    if(any(object@table[,n.var + 1] < 0))
        return("The Frequencies should be non-negative integers!")
    TRUE
}

setValidity("freq.tb", .valid.freq.tb)
