#define DINDEX(i, j) n*(i - 1) - i*(i - 1)/2 + j - i - 1

void ultrametric(double *dd, int* np,int* mp,double *ret)//d received as dist object, -1 for missing entries
{
    int n=*np;
    int m=*mp;
    int i=0,j=0;
    double max=dd[0];
    double d[n][n];
    for(i=1;i<n;i++)
    {d[i-1][i-1]=0;
     for(j=i+1;j<=n;j++)
      {
         d[i-1][j-1]=d[j-1][i-1]=dd[give_index(i,j,n)];
         if(dd[give_index(i,j,n)]>max)
          {
            max=dd[give_index(i,j,n)];
          }
      }
    }
    d[n-1][n-1]=0;
   
  int entrCh=0;
   do{
    entrCh=0;
    for(i=0;i<n-1;i++)
     for(j=i+1;j<n;j++)
      {  
         if(d[i][j]!=-1)continue;
         double minimax=max;
         int k=0;
         int sw=0;
         for(k=0;k<n;k++)
          {  
             if(d[i][k]==-1 || d[j][k]==-1)continue;
             sw=1;
             double mx=d[i][k]>d[j][k]?d[i][k]:d[j][k];
             if(mx<minimax){minimax=mx;}
          }
        if(sw==1)
          {
            d[i][j]=d[j][i]=minimax;
            m--;
            entrCh=1;
          }
      }
   }while(entrCh==1);
  int ij=0;
  for(i=0;i<n;i++)
   for(j=0;j<n;j++)
    {
       ret[ij++]=d[i][j];
    }
}
