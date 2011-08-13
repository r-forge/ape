mvr <- function(X,V)
{
    if (is.matrix(X)) X <- as.dist(X)
    X
    if (is.matrix(V)) V <- as.dist(V)
    V
    if (any(is.na(X)))
        stop("missing values are not allowed in the distance matrix")
    if (any(is.na(V)))
        stop("missing values are not allowed in the variance matrix")
    N <- attr(X, "Size")
    labels <- attr(X, "Labels")
    if (is.null(labels)) labels <- as.character(1:N)
    ans <- .C("mvr", as.double(X), as.double(V),as.integer(N), integer(2*N - 3),
              integer(2*N - 3), double(2*N - 3),
              DUP = FALSE, NAOK = TRUE, PACKAGE = "ape")
    obj <- list(edge = cbind(ans[[4]], ans[[5]]), edge.length = ans[[6]],
                tip.label = labels, Nnode = N - 2L)
    class(obj) <- "phylo"
    reorder(obj)
}
