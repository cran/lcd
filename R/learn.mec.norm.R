learn.mec.norm <- function(tree, cov, n, p.value, method = "CG")
{
    if(!is.element(method, c("CG", "DAG", "UG"))){
        cat("Invalid method!\nMethod must be one of: \"CG\", \"DAG\" and \"UG\"!\n")
        return(invisible())
    }
    skel <- learn.skeleton.norm(tree, cov, n, p.value, method != "DAG")
    switch(method,
           CG = learn.complex.norm(skel, cov, n, p.value),
           DAG = learn.v(skel, tree),
           UG = skel$amat)
}
