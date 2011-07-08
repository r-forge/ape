SDM <- function(...)
{ 
  st<-list(...)#first half contains matrices, second half s_p
  k<-length(st)/2 
  labels <- attr(as.dist(st[[1]]),"Labels")
  tot=length(rownames(st[[1]]))
  for(i in c(2:k))
   { 
     labels <- union(labels,attr(as.dist(st[[i]]),"Labels"))    
     tot=tot+length(rownames(st[[i]]))
   }
  sp<-mat.or.vec(1,k)
  for(i in c(k+1:k))
   { 
     sp[i-k] <- st[[i]] 
   } 
  astart <- mat.or.vec(1,tot)#start of aip, astart[p] is start of aip  
  astart[1]=k
  for(i in c(2:k))
   {
     astart[i]=astart[i-1]+length(rownames(st[[i-1]]))
   }
  a <- mat.or.vec(1,k+tot+k+length(labels))#first k are alphas, subseqeunt ones aip
 			        #each matrix p starting at astart[p], next are 
			        #Lagrange multipliers
  n <- length(labels)
  X <- mat.or.vec(n,n)
  V <- mat.or.vec(n,n)
  w <- mat.or.vec(n,n) 
  rownames(X) <- labels
  colnames(X) <- labels
  
  rownames(V) <- labels
  colnames(V) <- labels

  rownames(w) <- labels
  colnames(w) <- labels
  
  for(i in c(1:(n-1)))
   for(j in c(i+1:(n-i)))
    { wij=0
      for(p in c(1:k))
        {d=st[[p]]
          if(is.element(rownames(X)[i],rownames(d)) & is.element(colnames(X)[j],colnames(d)))
            {
              w[i,j]=w[j,i]=(w[i,j]+sp[p])
            } 
        }  
    }
 Q=mat.or.vec(length(labels)+k+k+tot,length(labels)+k+k+tot) 
 for(p in c(1:k)) 
  { d_p=st[[p]]    
    for(l in c(1:k))#first compute coeficients of alphas
     { sum=0	 
       if(l==p)#calculate alpha_p
        { 
          for(i in c(1:n))
          {
           for(j in c(1:n))
            { #check if {i,j}\subset L_l
		  d=st[[l]]             		  	
              if(i!=j & is.element(labels[i],rownames(d)) & is.element(labels[j],colnames(d)))
               {
                sum=sum+((d[rownames(d)==labels[i],colnames(d)==labels[j]]*d[rownames(d)==labels[i],colnames(d)==labels[j]])-sp[l]*d[rownames(d)==labels[i],colnames(d)==labels[j]]/w[i,j])  
		   }
            }
          }
	   }else
          {
            for(i in c(1:n))
             {
         	  for(j in c(1:n))
               { #check if {i,j}\subset L_l
	   	     d=st[[l]]             		  	
                 if(i!=j & is.element(labels[i],rownames(d)) & is.element(labels[j],colnames(d)) & is.element(labels[i],rownames(d_p)) & is.element(labels[j],colnames(d_p)))
                  {
                   sum=sum-sp[l]*d[rownames(d)==labels[i],colnames(d)==labels[j]]/w[i,j] 
		      }
               }
             }
          }
      Q[p,l]=sum
     }
  } 
 print(Q)
 list(X,V)
}