iscompatible <-function(obj,...)
{
  m=obj$x
  n=ncol(m)
  ret=TRUE
  for(i in c(1:(n-1)))
   for(j in c((i+1):n))
    {
      if(!arecompatible(m[,i],m[,j],obj$n))
        {
          ret=FALSE
        }
    }
  ret
}