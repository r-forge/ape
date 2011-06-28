#include <R.h>


void triangMtdS(double* d, int* np, int* ed1,int* ed2, double* edLen)
 {
    int n=*np;
    int k=0;
    int i=0;
    int j=0;
    int ij=-1;
    int m[n+1];
    int c[n+1];
    int s[n+1];
    int o[n+1];
    for(i=1;i<=n;i++)
    {m[i]=0;
     c[i]=i;
     s[i]=0;
     for(j=1;j<=n;j++)
      {
        if(i==j){m[i]++;continue;}
        if(d[give_index(i,j,n)]==-1)continue;
        m[i]++;
      }
    }

   for(i=1;i<n;i++)
    for(j=i+1;j<=n;j++)
      {
        if(m[i]<m[j])
         {
            int aux=m[i];
            m[i]=m[j];
            m[j]=aux;
            aux=c[i];
            c[i]=c[j];
            c[j]=aux;
         }
      }
   ij=1;
   while(m[ij]==n){ij++;}
   for(i=ij;i<n;i++)
     {if(s[c[i]]==1)continue;
      for(j=i+1;j<=n;j++)
      {
       if(s[c[j]]==1)continue;      
       if(d[give_index(c[i],c[j],n)]==-1)
         {s[c[j]]=1;
         }
      }
     }
   
   for(i=1;i<=n;i++)
     {
       if(s[i]==0)
         {
          k++;
          o[k]=i;
         }
     }
   double* sub_d=(double*)R_alloc(k*(k - 1)/2, sizeof(double));
   ij=0;
   for(i=1;i<n;i++)
   {if(s[i]==1)continue;
    for(j=i+1;j<=n;j++)
     {
      if(s[j]==1)continue;
      sub_d[ij++]=d[give_index(i,j,n)];
     }
   }
   
   triangMtd(sub_d,&k,ed1,ed2,edLen);
   for(i=0;i<2*k-3;i++)
     {
       if(ed1[i]>k)
        {
          ed1[i]+=(n-k);
        }
       if(ed2[i]>k)
        {
          ed2[i]+=(n-k);
        }
     }
   for(i=0;i<2*k-3;i++)
     {
      if(ed2[i]<=k)
       {
          ed2[i]=o[ed2[i]];
       }
     }
   for(i=1;i<=n;i++)
    {
      if(s[i]==0)continue;//take only leaves not in Y
      m[i]=0;
      for(j=1;j<=n;j++)
       {
         if(s[j]==1)continue;//take only leaves already in Y
         if(d[give_index(i,j,n)]==-1)continue;//igonore if distance unknown
         m[i]++;         
       }
    }
 while(k<n)
 {
   //now find element in X\Y such that it has most known distances to already
   //built tree until all leaves are added or we can not find a place to attach
   //the new element
   //s[i]=1 => i not added to tree
    int max=-1e50;
    int maxPos=-1;
    for(i=1;i<n;i++)
    {
     if(s[i]==0)continue;
     if(m[i]>max)
      {
        max=m[i];
        maxPos=i;
      }
    }
   s[maxPos]=0;//mark maxPos as added
   //calculate new m values for leaves not added, i.e we just increment any 
   //already present value by 1 if we know the distance between i and maxPos
   for(i=1;i<=n;i++)
    {
      if(s[i]==0)continue;
      if(d[give_index(i,maxPos,n)]==-1)continue;
      m[i]++;
    }

   //find path to attach maxPos to, grow tree
   //TO DO 
   k++;
 }
 }
