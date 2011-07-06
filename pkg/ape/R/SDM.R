SDM <- function(...)
{ 
  st<-list(...)#first half contains matrices, second half s_p
  k<-length(st)/2 
  labels <- attr(as.dist(st[[1]]),"Labels")
  for(i in c(2:k))
   { 
     labels <- union(labels,attr(as.dist(st[[i]]),"Labels"))    
   }
  sp<-mat.or.vec(1,k)
  for(i in c(k+1:k))
   { print(i)
     sp[i-k] <- st[[i]] 
   } 
  n <- length(labels)
  X <- mat.or.vec(n,n)
  V <- mat.or.vec(n,n)
  rownames(X) <- labels
  colnames(X) <- labels
  
  rownames(V) <- labels
  colnames(V) <- labels
  
  for(i in c(1:(n-1)))
   for(j in c(i+1:(n-i)))
    { wij=0
      for(p in c(1:k))
        {d=st[[p]]
          if(is.element(rownames(X)[i],rownames(d)) & is.element(colnames(X)[j],colnames(d)))
            {
              wij=wij+sp[p]
            } 
        } 
      sumij=0 
      vsumij=0
      for(p in c(1:k))
        {d=st[[p]]
          if(is.element(rownames(X)[i],rownames(d)) & is.element(colnames(X)[j],colnames(d)))
            {
             sumij=sumij+(sp[p]*d[rownames(d)==rownames(X)[i],colnames(d)==colnames(X)[j]])
             vsumij=vsumij+(sp[p]*d[rownames(d)==rownames(X)[i],colnames(d)==colnames(X)[j]]*d[rownames(d)==rownames(X)[i],colnames(d)==colnames(X)[j]])
            } 
        }
     X[j,i]=X[i,j]=(1/wij)*sumij 
     V[j,i]=V[i,j]=(1/(wij*wij))*vsumij
    }
  list(X,V)
}