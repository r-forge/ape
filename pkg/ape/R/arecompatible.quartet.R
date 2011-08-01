arecompatible.quartet <-function(x,...)
{
  n=x$n #number of taxa
  mf=x$mf #quartet label matrix
  lab=x$labels
  w=x$w
  ret=TRUE
  intq=matrix(0,nrow(mf),4)
  for(i in c(1:nrow(mf)))
   for(j in c(1:4))
    {
     intq[i,j]=which(lab==mf[i,j])
    }
  print(intq)
  print("c")
  ret=0
  ans<-.C("quartComp",as.integer(intq),as.double(w),as.integer(nrow(mf)),as.integer(n),as.integer(1))
  ret=ans[[5]]
  print("ret") 
  print(ret)
 ret 
}