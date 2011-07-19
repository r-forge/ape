bionjs <- function(X,fs=15)
{   
   
    if (is.matrix(X)) X <- as.dist(X)
    X[is.na(X)]=-1
    X[X<0]=-1	
    X[is.nan(X)]=-1
    N <- attr(X, "Size")
    labels <- attr(X, "Labels")
    if (is.null(labels)) labels <- as.character(1:N)
    ans <- .C("bionjs", as.double(X), as.integer(N), integer(2*N - 3),
              integer(2*N - 3), double(2*N - 3),as.integer(fs),
              DUP = FALSE, NAOK = TRUE, PACKAGE = "ape")
    obj <- list(edge = cbind(ans[[3]], ans[[4]]), edge.length = ans[[5]],
                tip.label = labels, Nnode = N - 2L)
    class(obj) <- "phylo"
    reorder(obj)
}