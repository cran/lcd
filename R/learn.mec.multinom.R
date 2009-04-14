learn.mec.multinom <- function(tree, freq.tb, p.value, method = "CG")
{
    if(!is.element(method, c("CG", "DAG", "UG"))){
        cat("Invalid method!\nMethod must be one of: \"CG\", \"DAG\" and \"UG\"!\n")
        return(invisible())
    }
    skel <- learn.skeleton.multinom(tree, freq.tb, p.value, method != "DAG")
    switch(method,
           CG = learn.complex.multinom(skel, freq.tb, p.value),
           DAG = learn.v(skel, tree),
           UG = skel$amat)
}
 
