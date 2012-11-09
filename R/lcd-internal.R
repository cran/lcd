`.find.clique` <-
function(input)
{
    if(!input$is.triangulated)
        stop("Invalid Input!")
    n <- input$perfect.numbering
    card <- input$card
    rcd <- input$pi.record
    lad <- c()
    kep <- c()
    for(i in 1:length(n)){
        if (i == length(n) | (card[i]+1)>card[i+1]){
            lad <- c(lad, n[i])
            kep <- c(kep, i)
        }
    }
    clique <- vector("list", length(lad))
    for(j in 1:length(lad)){
        if(card[kep[j]] > 0){
            a <- sum(card[1:(kep[j]-1)])
            clique[[j]] <- c(lad[j], rcd[(a+1):(a+card[kep[j]])])
        } else {
            clique[[j]] <- c(lad[j])
        }
    }
    return(list(n.clique = length(lad),
                cliques = clique))
}

`.junc.tree` <-
function(c, vset)
{
    n <- c$n.clique
    clique <- vector("list", n)
    cnames <- paste("C", 1:n, sep="")
    for(i in 1:n) clique[[i]] <- list(name = cnames[i],
                                      vset = vset[c$cliques[[i]]])
    adjmat <- matrix(0, nrow=n, ncol=n)
    rownames(adjmat) <- colnames(adjmat) <- cnames
    sep <- vector("list", n-1)
    temp <- c()
    k <- 1
    if(n > 1){
        for(i in 2:n) {
            temp <- c(temp, clique[[i-1]]$vset)
            test <- intersect(clique[[i]]$vset, temp)
            for(j in 1:(i-1)){
                if(all(is.element(test, clique[[j]]$vset))){
                    adjmat[i,j] <- 1
                    sep[[k]] <- list(separator = intersect(clique[[j]]$vset,
                                     clique[[i]]$vset),
                                     edge = cnames[c(j,i)])
                    k <- k+1
                    break
                }
            }
        }
    }
    adjmat <- adjmat + t(adjmat)
    return(list(tree.struct = adjmat,
                cliques = clique,
                separators = sep))
#    return(list(Structure = adjmat,
#                Cliques = clique,
#                Separators = sep))
}


`.binary` <- function (x, dim)         # quoted from the "wle" package
{
    if (x == 0) {
        pos <- 1
    }
    else {
        pos <- floor(log(x, 2)) + 1
    }
    if (!missing(dim)) {
        if (pos <= dim) {
            pos <- dim
        }
        else {
            warning("the value of `dim` is too small")
        }
    }
    bin <- rep(0, pos)
    dicotomy <- rep(FALSE, pos)
    for (i in pos:1) {
        bin[i] <- floor(x/2^(i - 1))
        dicotomy[i] <- bin[i] == 1
        x <- x - ((2^(i - 1)) * bin[i])
    }
    return(list(binary = bin, dicotomy = dicotomy))
}

# need to be accelerated!!!

#`norm.ci.test` <- function(S, N, i, j, Q = c())
#{
#    q <- length(Q)
#    Mmar <- S[c(i, j, Q), c(i, j, Q)]
#    S11 <- Mmar[1, 1]
#    S12 <- Mmar[1, -1]
#    S21 <- Mmar[-1, 1]
#    S22 <- Mmar[-1, -1]
#    S22inv <- solve(S22)
#    betahat <- S12 %*% S22inv[, 1]
#    sigma <- sqrt((S11 - S12 %*% S22inv %*% S21) * (N - 1)/(N - 
#        q - 2))
#    se <- sigma * sqrt(S22inv[1, 1]/(N - 1))
#    t.value <- betahat/se
#    p.value <- 2 * (1 - pt(abs(t.value), N - q - 2))
#    return(list(t.value = t.value, p.value = p.value))
#}

#`norm.ci.test` <- function(cov, n, u, v, sep = c())
#{
#    vset <- rownames(cov)
#    if(is.character(u)) u <- which(vset == u)
#    if(is.character(v)) v <- which(vset == v)
#    if(is.character(sep)) sep <- match(sep, vset)
#    if(!is.numeric(u) || !is.numeric(v) || (!is.numeric(sep) && length(sep) > 0))
#        stop("Invalid vertex!")
#    res <- .Call("g_ci_test", cov,
#                 as.integer(n),
#                 as.integer(u), as.integer(v), as.integer(sep),
#                 PACKAGE = "lcd")
#    return(list(t.value = res[1], dof = res[2], p.value = res[3]))
#}


# Function for learning the undirected graph in each tree node
# Last modified: March-06-2008
# Just return the separation pairs seems to be a good strategy?
# At least easy for future implementation!

`.get.localug.ic` <- function(cov, n, p.value){
    p <- nrow(as.matrix(cov))
    vset <- rownames(cov)
    sep.pairs <- c()
    if(p > 1)
        for(i in 1:(p-1))
            for(j in (i+1):p){
                cand <- vset[-c(i,j)]
                res <- .get.sep(cov, n, p.value, vset[i], vset[j], cand)
                if(res$seped)
                    sep.pairs <- append(sep.pairs, res$sep)
            }
    sep.pairs
}

## adapt code from the "pcAlgo" from "pcalg" package in R

`.get.localug.pc` <- function(cov, n, p.value){
    p <- nrow(as.matrix(cov))
    vset <- rownames(cov)
    sep.pairs <- c()
    G <- matrix(rep(TRUE, p*p), p, p)
    diag(G) <- FALSE
    seq_p <- 1:p
    done <- FALSE
    ord <- 0
    while (!done && any(G)) {
        done <- TRUE
        ind <- which(G, arr.ind = TRUE)
        ind <- ind[order(ind[ ,1]), ]
        remainingEdgeTests <- nrow(ind)
        for(i in 1:remainingEdgeTests) {
            x <- ind[i,1]
            y <- ind[i,2]
            if (G[y, x]) {
                nbrsBool <- G[, x]
                nbrsBool[y] <- FALSE
                nbrs <- seq_p[nbrsBool]
                length_nbrs <- length(nbrs)
                if (length_nbrs >= ord) {
                    if (length_nbrs > ord)
                        done <- FALSE
                    S <- seq(length = ord)
                    repeat {
                        p.val <- norm.ci.test(cov, n, vset[x], vset[y],
                                              vset[nbrs[S]])$p.value
                        if (p.val > p.value) {
                            G[x, y] <- G[y, x] <- FALSE
                            pair <- new("sep.pair", u = vset[x],
                                        v = vset[y], s = vset[nbrs[S]])
                            sep.pairs <- append(sep.pairs, pair)
                            break
                        }
                        else {
                            nextSet <- .getNextSet(length_nbrs, ord, S)
                            if (nextSet$wasLast)
                                break
                            S <- nextSet$nextSet
                        }
                    }
                }
            }
        }
        ord <- ord + 1
    }
    sep.pairs
}


`.getNextSet` <- function (n, k, set) 
{
    chInd <- k - (zeros <- sum((seq(n - k + 1, n) - set) == 0))
    wasLast <- (chInd == 0)
    if (!wasLast) {
        set[chInd] <- set[chInd] + 1
        if (chInd < k) 
            set[(chInd + 1):k] <- seq(set[chInd] + 1, set[chInd] + 
                zeros)
    }
    list(nextSet = set, wasLast = wasLast)
}



#`.get.localug` <- function(cov, n, p.value){
#    p <- nrow(as.matrix(cov))
#    vset <- rownames(cov)
#    if(p>1){
#         a <- outer(rep(1,p),1:p)
#         b <- t(a)
#         a <- as.numeric(a[upper.tri(a)])
#         b <- as.numeric(b[upper.tri(b)])
#         pair <- cbind(b, a)
#         mat <- t(apply(pair, 1, function(x) vset[-x]))
#         mat <- cbind(vset[b], vset[a], mat)
#         res <- apply(mat, 1, function(x)
#                      .get.sep(cov, n, p.value, x[1], x[2], x[-c(1,2)]))
#    }
#    idx <- which(lapply(res, function(x) x$seped) == TRUE)
#    lapply(res[idx], function(x) x$sep)
#}

# Function for getting separator for a vertex pair

## `.get.sep` <- function(cov, n, p.value, u, v, cand)
## {
##     seped <- FALSE
##     p <- length(cand)
##     sep <- new("sep.pair",
##                u=u,
##                v=v,
##                s=character(0))
##     for(k in 1:(2^p)){
##         idx <- .binary(k-1, max(1, p))$dicotomy
##         if(any(idx)) tset <- cand[idx] else tset <- c()
##         if(norm.ci.test(cov, n, u, v, tset)$p.value >= p.value){
##             sep@s=cand[idx]
##             seped <- TRUE
##             break
##         }
##     }
##     return(list(seped=seped, sep=sep))
## }


.get.sep <- function(cov, n, p.value, u, v, cand)
{
    seped <- FALSE
    p <- length(cand)
    if (p < 8) {
        mat <- as.matrix(t(sapply(1:(2^p),
                                  function(x) .binary(x-1, max(1,p))$dicotomy)))
        if(nrow(mat) < ncol(mat)) mat <- t(mat)
        p.val <- apply(mat, 1, function(x) norm.ci.test(cov, n, u, v, cand[x])$p.value)
        p.val.max <- max(p.val)
        idx <- which(p.val == p.val.max)
        sep <- new("sep.pair",
                   u=u,
                   v=v,
                   s=character(0))
        if(p.val.max >= p.value){
            sep@s <- cand[mat[idx[1],]]
            seped <- TRUE
        }
        return(list(seped=seped, sep=sep))
    } else {
        .get.sep.stepwise(cov, n, p.value, u, v, cand)
    }
}


.get.sep.stepwise <- function(cov, n, p.value, u, v, cand){
     .get.sep.step(cov, n, p.value, u, v, cand, c())
 }

.get.sep.step <- function(cov, n, p.value, u, v, current, rest)
{
    modified <- FALSE
    sep <- new("sep.pair",
               u = u,
               v = v,
               s = character(0))
    pp <- norm.ci.test(cov, n, u, v, current)$p.value
    if(pp >= p.value) {
        sep@s = current
        return(list(seped = TRUE, sep = sep))
    }
    if(length(current) == 0) 
        return(list(seped = FALSE, sep = sep))
    n.curr <- length(current)
    n.rest <- length(rest)
    mat <- matrix(rep(current, n.curr), n.curr, byrow = TRUE)
    diag(mat) <- NA
    pval <- apply(mat, 1, function(x) norm.ci.test(cov, n, u, v, x[!is.na(x)])$p.value)
    todel <- which(pval == max(pval))
    if(length(todel) > 1)
        todel <- todel[sample(length(todel),1)]
    if(pval[todel] >= pp){
        todel <- current[todel]
        current <- setdiff(current, todel)
        rest <- union(rest, todel)
        pp <- max(pval)
        modified <- TRUE
    }
    if(pp >= p.value) {
        sep@s = current
        return(list(seped = TRUE, sep = sep))
    }
    if(length(rest) == 0)
        return(list(seped = FALSE, sep = sep))
    if(modified){
        n.curr <- length(current)
        n.rest <- length(rest)
    }
    mat <- matrix(NA, n.rest, n.curr + 1)
    if(n.curr != 0)
        mat[,1:n.curr] <- matrix(rep(current, n.rest), n.rest, byrow = TRUE)
    mat[,n.curr + 1] <- rest
    pval <- apply(mat, 1, function(x) norm.ci.test(cov, n, u, v, x)$p.value)
    toadd <- which(pval == max(pval))
    if(length(toadd) > 1)
        toadd <- toadd[sample(length(toadd),1)]
    if(pval[toadd] > pp){
        toadd <- rest[toadd]
        current <- union(current, toadd)
        rest <- setdiff(rest, toadd)
        modified <- TRUE
    }
    if(modified)
        return(.get.sep.step(cov, n, p.value, u, v, current, rest))
    else 
        return(list(seped = FALSE, sep = sep))
}



#`.get.sep` <- function(cov, n, p.value, u, v, cand)
#{
#    seped <- FALSE
#    p <- length(cand)
#    mat <- as.matrix(t(sapply(1:(2^p),
#                              function(x) .binary(x-1, max(1,p))$dicotomy)))
#    res <- apply(mat, 1, function(x) norm.ci.test(cov, n, u, v, cand[x])$p.value)
#    idx <- which(res >= p.value)
#    sep <- new("sep.pair",
#               u=u,
#               v=v,
#               s=character(0))
#    if(any(idx)){
#        sep@s <- cand[mat[idx[1],]]
#        seped <- TRUE
#    }
#    return(list(seped=seped, sep=sep))
#}

# Function for reading separation info from tree
# Last modified: March-06-2008

`.break.tree` <- function(tree)
{
    set.seed(as.numeric(Sys.time()))
    sep.pairs <- c()
    n.sep <- length(tree@separators)
    if(n.sep == 0)
        return(sep.pairs)
    b <- sample(n.sep, 1)
    e <- tree@separators[[b]]$edge
    sep <- tree@separators[[b]]$separator
    btr <- tree@tree.struct
    btr[e[1], e[2]] <- btr[e[2], e[1]] <- 0
    bt <- graph.adjacency(btr, "undirected")
    cc <- clusters(bt)
    l.vset <- r.vset <- c()
    l <- r <- c()
    for(i in 1:(n.sep+1)){
        if(cc$membership[i]==0){
            l.vset <- c(l.vset, tree@cliques[[i]]$vset)
            l <- c(l, i)
        } else {
            r.vset <- c(r.vset, tree@cliques[[i]]$vset)
            r <- c(r, i)
        }
    }
    l.vset <- unique(l.vset)
    r.vset <- unique(r.vset)
    for(i in 1:length(l.vset))
        for(j in 1:length(r.vset)){
            if(!is.element(l.vset[i], sep) &&
               !is.element(r.vset[j], sep)){
                sepp <- new("sep.pair",
                           u=l.vset[i],
                           v=r.vset[j],
                           s=sep)
                sep.pairs <- append(sep.pairs, sepp)
            }
        }
    l.cliques <- tree@cliques[l]
    r.cliques <- tree@cliques[r]
    l.struct <- tree@tree.struct[l,l]
    r.struct <- tree@tree.struct[r,r]
    l.separators <- vector("list", length(l)-1)
    r.separators <- vector("list", length(r)-1)
    ll <- rr <- 1
    for(i in 1:n.sep){
        if(i != b){
            if(is.element(tree@separators[[i]]$edge[1], rownames(l.struct))){
                l.separators[ll] <- tree@separators[i]
                ll <- ll+1
            } else {
                r.separators[rr] <- tree@separators[i]
                rr <- rr+1
            }
        }
    }
    l.tree <- list(tree.struct = l.struct,
                   cliques = l.cliques,
                   separators = l.separators)
    r.tree <- list(tree.struct = r.struct,
                   cliques = r.cliques,
                   separators = r.separators)
    sep.pairs <- append(sep.pairs, .break.tree(l.tree))
    sep.pairs <- append(sep.pairs, .break.tree(r.tree))
}


# Fucntion for getting the possible extra edges
# Last modified: March-06-2008

`.get.exed.cand` <- function(tree, amat){
    n.clique <- length(tree@cliques)
    cand.pairs <- c()
    if(n.clique == 1) return(cand.pairs)
    for(i in 1:(n.clique-1))
    {
        sepset <- tree@separators[[i]]$separator
        sep.amat <- amat[sepset, sepset]
        if(length(sepset)>=2){
            for(j in 2:length(sepset))
                for(k in 1:(j-1)){
                    if(sep.amat[j,k]!=0){
                        newpair <- new("sep.pair",
                                       u = sepset[j],
                                       v = sepset[k],
                                       s = character(0))
                        cand.pairs <- append(cand.pairs, newpair)
                    }
                }
        }
    }
    cand.pairs
}


`.get.exed.cand1` <- function(tree, amat){
    vset <- rownames(amat)
    n.clique <- length(tree@cliques)
    cand.pairs <- c()
    if(n.clique == 1) return(cand.pairs)
    for (i in 1:(n.clique - 1)) {
        sepset <- tree@separators[[i]]$separator
        idx <- match(sepset, vset)
        G <- amat[idx, idx]
        ind <- which(G == 1, arr.ind = TRUE)
        if (any(ind)) {
            ind[,1] <- idx[ind[,1]]
            ind[,2] <- idx[ind[,2]]
            cand.pairs <- rbind(cand.pairs, ind)
        }
    }
    unique(cand.pairs)
}


# Function for checking whether two vertices are (undirected) linked in a graph
# Last modified: March-09-2008

#`.is.linked` <- function(amat, u, v)
#{
#    amat <- as.matrix(amat)
#    if(nrow(amat) == 0) return(FALSE)
#    q <- c()
#    q <- c(q, u)
#    vset <- rownames(amat)
#    if(!all(is.element(c(u,v), vset))){
#        #warning("Invalid vertex!")
#        return(FALSE)
#    }
#    un <- rep(TRUE, length(vset))
#    while(length(q) > 0){
#        if(is.element(v, q))
#            return(TRUE)
#        toadd <- c()
#        for(i in 1:length(q)){
#            idx <- amat[q[i],] & amat[,q[i]] & un
#            toadd <- c(toadd, vset[idx])
#            un[idx] <- FALSE
#        }
#        toadd <- unique(toadd)
#        q <- toadd
#    }
#    return(FALSE)
#}

`.is.linked` <- function(amat, u, v){
    amat <- as.matrix(amat)
    if(nrow(amat) == 0) return(FALSE)
    vset <- rownames(amat)
    if(!all(is.element(c(u,v), vset))){
        return(FALSE)
    }
    i <- which(vset == u)
    j <- which(vset == v)
    as.logical(.Call("is_linked", amat,
                     as.integer(i), as.integer(j),
                     PACKAGE = "lcd"))
}


# helper function to get.multinom.dist

.make.table <- function(vset, n.state, subset, alpha, beta, freq = TRUE){
    idx <- match(subset, vset)
    phi <- matrix(0, prod(n.state[idx]), length(idx))
    colnames(phi) <- c(subset)
    for(k in 1:length(idx)){
        a <- n.state[idx[k]]
        b <- ifelse(k == length(idx), 1, prod(n.state[idx[(k+1):length(idx)]]))
        c <- ifelse(k == 1, 1, prod(n.state[idx[1:(k-1)]]))
        phi[,k] <- rep(rep(1:a, rep(b, a)), c)
    }
    if(freq){
        n <- nrow(phi)
        Freq <- rbeta(n, alpha, beta)
        ## Freq <- rbeta(n, 1/3, 1/3)
        ## key step in generating a possibly FAITHFUL distribution
        ## the current trick is to make each
        ## potential close to either 0 or 1
        phi <- cbind(phi, Freq)
    }
    phi
}


.lcd.cliques <- function(amat)
  {
    vertices <- rownames(amat)
    amat.igraph <- graph.adjacency(amat, mode = "max")
    cliq.list <- cliques(amat.igraph)
    lapply(cliq.list, function(x) vertices[x])
  }
