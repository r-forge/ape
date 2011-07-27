arecompatible <-function(x,y,n)
{
  res=vector("raw")
  msk=as.raw(2^(8-(n%%8))-1)
  msk=!msk
    
  ret=FALSE
  nE=0
  
  res=x&y   
  bt=res[ceiling(n/8)]
  bt=bt&msk
  res[ceiling(n/8)]=bt
  if(all(res==as.raw(0)))
    {
      nE=nE+1
    }
  
  res=x&(!y)
  bt=res[ceiling(n/8)]
  bt=bt&msk
  res[ceiling(n/8)]=bt
  if(all(res==as.raw(0)))
    {
      nE=nE+1
    }

  res=(!x)&y
  bt=res[ceiling(n/8)]
  bt=bt&msk
  res[ceiling(n/8)]=bt
  if(all(res==as.raw(0)))
    {
      nE=nE+1
    }
 
  res=(!x)&(!y)
  bt=res[ceiling(n/8)]
  bt=bt&msk
  res[ceiling(n/8)]=bt
  if(all(res==as.raw(0)))
    {
      nE=nE+1
    }
  
  if(nE==1)ret=TRUE
  
  ret

}