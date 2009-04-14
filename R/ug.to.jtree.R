`ug.to.jtree` <-
function(ugamat)
{
    vset <- rownames(ugamat)
    a <- maxcard.search(tri.ug(ugamat))
    tree <- .junc.tree(.find.clique(a), vset)
    new("sep.tree",
        tree.struct = tree$tree.struct,
        cliques = tree$cliques,
        separators = tree$separators)
}
