treePop <- function(obj)
{
  mf=obj$x
  n=obj$n
  labels=obj$labels
  imf=as.integer(mf)
  mimf=matrix(imf,nrow(mf),ncol(mf))
  ans<- .C("treePop",mimf,as.integer(ncol(mf)),as.integer(n),integer(2*n-3),integer(2*n-3),double(2*n - 3)
           ,DUP = FALSE, NAOK = TRUE, PACKAGE = "ape")
  obj <- list(edge = cbind(ans[[3]], ans[[4]]), edge.length = ans[[5]],
                tip.label = labels, Nnode = n - 2L)
  class(obj) <- "phylo"
  reorder(obj)

}