`naive.getug.norm` <-
function(data, p.value)
{
    conmat <- solve(cor(data))
    p <- ncol(data)
    n <- nrow(data)
    parcor <- cov2cor(conmat)
    t.stat <- parcor*sqrt((n-p-2)/(1-parcor^2))
    thres <- -qt(p.value/2, df=n-p-2)
    amat <- matrix(as.integer(abs(t.stat)>thres), nrow=p, ncol=p)
    diag(amat) <- rep(0,p)
    rownames(amat) <- colnames(amat) <- rownames(conmat)
    amat
}

