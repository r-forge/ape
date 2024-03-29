mvrs <- function(X,V,fs=15)
{      
    if (is.matrix(X)) X <- as.dist(X)
    X
    if (is.matrix(V)) V <- as.dist(V)
    V
    
    X[is.na(X)]=-1
    X[X<0]=-1	
    X[is.nan(X)]=-1    
    
    N <- attr(X, "Size")
    labels <- attr(X, "Labels")
    if (is.null(labels)) labels <- as.character(1:N)
    ans <- .C("mvrs", as.double(X), as.double(V),as.integer(N), integer(2*N - 3),
              integer(2*N - 3), double(2*N - 3),as.integer(fs),
              DUP = FALSE, NAOK = TRUE, PACKAGE = "ape")
    obj <- list(edge = cbind(ans[[4]], ans[[5]]), edge.length = ans[[6]],
                tip.label = labels, Nnode = N - 2L)
    class(obj) <- "phylo"
    reorder(obj)
}
