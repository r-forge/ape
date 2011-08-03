additive <- function(X)
{
    if (is.matrix(X)) X <- as.dist(X)
    X[is.na(X)]=-1
    X[X<0]=-1	
    X[is.nan(X)]=-1
    N <- attr(X, "Size")
    labels <- attr(X, "Labels")
    if (is.null(labels)) labels <- as.character(1:N)
    m<-length(which(is.na(d)))
    ans<- .C("additive",as.double(X),as.integer(N),as.integer(m),double(N*N))
    ret<-matrix(ans[[4]],N,N)   
  ret 
}