
# some helper functions

# parse.graph
# ===========
# construct a formula for log-linear model from an undirected graph

.parse.graph <- function(amat)
{
    if(any(skeleton(amat) != amat))
        stop("You should input an undirected graph!")
    cliq <- .lcd.cliques
    cl <- cliq(amat)
    ele <- sapply(cl, function(x) paste(x, collapse=":"))
    as.formula(paste("Freq ~", paste(ele, collapse="+")))
}

# 
# fit: original fit
# amat: orginal graph structure
# freq.tb: frequency table as data
# edge: a number of where to add edge
# no intended to be called by user!
#

.update.fit <- function(fit, amat, edge, freq.tb)
{
    amat[upper.tri(amat)][edge] <- 1
    amat <- skeleton(amat)
    form <- .parse.graph(amat)
    newfit <- loglm(form, data = as.data.frame(freq.tb@table))[c(3,8)]
    res <- c(fit$deviance - newfit$deviance,
             fit$df - newfit$df)
    c(newfit$deviance, newfit$df, abs(res[1]), abs(res[2]),
      1 - pchisq(abs(res[1]), abs(res[2])))
}


#
# amat: original graph
# fit: original log-linear model
# freq.tb: frequency table as data
# p.value: thresholding p-value for picking edge
# 

.try.addedge <- function(amat, fit, freq.tb, p.value)
{
    change <- FALSE
    if(any(skeleton(amat) != amat))
        stop("You should input an undirected graph!")
    candlist <- which(!amat[upper.tri(amat)] == 1)
    if (length(candlist) > 0) {
        res <- t(sapply(candlist, function(x)
                        .update.fit(fit, amat, x, freq.tb)))
        p.val <- res[,5]
        p.val.min <- min(p.val)
        idx <- which(p.val == p.val.min)[1]
        if(p.val.min < p.value){
            change <- TRUE
            amat[upper.tri(amat)][candlist[idx]] <- 1
            amat <- skeleton(amat)
            fit <- list(df = res[idx,2], deviance = res[idx,1])
        }
    }
    return(list(amat = amat,
                change = change,
                fit = fit))
}


#
# if data is sparse, do forward selection
# else do simultaneous testing (or backward selection? may be slow...)
#

naive.getug.multinom <- function(freq.tb, p.value, method = "mkb")
{
    mm <- switch(method,
                 mkb = 1,
                 simple = 2,
                 fwd = 3,
                 0)
    if (mm == 0) stop("Invalid method!")
    vnames <- colnames(freq.tb@table)
    vnames <- vnames[-length(vnames)]
    p <- length(vnames)
    amat <- matrix(0, p, p)
    rownames(amat) <- colnames(amat) <- vnames
    if (mm == 1)
        return(.naive.getug.mkb(freq.tb, p.value))
    if (mm == 2) {
        for(i in 1:(p-1)) for(j in (i+1):p)
            amat[i,j] <- multinom.ci.test(freq.tb,
                                          vnames[i],
                                          vnames[j],
                                          vnames[-c(i,j)])$p.value < p.value + 0
        return(skeleton(amat))
    }
    if (mm == 3){
        initform <- .parse.graph(amat)
        fit <- loglm(initform, data = as.data.frame(freq.tb@table))[c(3,8)]
        while (TRUE) {
            res <- .try.addedge(amat, fit, freq.tb, p.value)
            if (res$change) {
                amat <- res$amat
                fit <- res$fit
            } else {
                break
            }
        }
    }
    amat
}


#
# alternative method: markov blanket selection
# ============================================
# a grow-shrink markov blanket selection procedure
# last modified: @@ Mon, Apr 28, 2008, 13:35 @@
#

.learn.mkvblkt <- function(freq.tb, var, curr = c(), p.value)
{
    forbid <- curr
    vnames <- colnames(freq.tb@table)
    vnames <- vnames[-length(vnames)]
    rest <- setdiff(vnames, c(curr,var))
    continue <- TRUE
    while (continue) {          # grow the Markov blanket first
        p.val <- sapply(rest, function(x)
                        multinom.ci.test(freq.tb, var, x, curr)$p.value)
        p.val.min <- min(p.val)
        idx <- which(p.val == p.val.min)[1]
        if(p.val.min < p.value) {
            curr <- c(curr, rest[idx])
            rest <- rest[-idx]
        } else {
            continue <- FALSE
        }
    }
    continue <- TRUE
    freq.tb <- compress.freq.tb(freq.tb, c(curr,var))
    delcand <- setdiff(curr, forbid)    # only those added later is allowed to be deleted
    if (length(delcand) == 0)
        continue <- FALSE
    while (continue) {       # shrink the Markov blanket
        p.val <- sapply(delcand, function(x)
                        multinom.ci.test(freq.tb, var, x, setdiff(curr,x))$p.value)
        ## this step could be speeded up significantly!!!
        p.val.max <- max(p.val)
        idx <- which(p.val == p.val.max)[1]
        if(p.val.max > p.value) {
            curr <- setdiff(curr, delcand[idx])
            delcand <- delcand[-idx]
        } else {
            continue <- FALSE
        }
    }
    curr
}

.naive.getug.mkb <- function(freq.tb, p.value)
{
    vnames <- colnames(freq.tb@table)
    vnames <- vnames[-length(vnames)]
    p <- length(vnames)
    amat <- matrix(0, p, p)
    vset <- vnames[sample(p)]
    rownames(amat) <- colnames(amat) <- vset
    for(i in 1:p) {
        curr <- vset[amat[vset[i],] == 1]
        res <- .learn.mkvblkt(freq.tb, vset[i], curr, p.value)
        amat[vset[i], res] <- 1
        amat[res, vset[i]] <- 1
    }
    amat[vnames, vnames]
}
