#include <math.h>
#include <R.h>

short* setdiff(short* x,short *y,int nrow)//x-y
  { int i=0;
    short* ret=(short*)R_alloc(nrow , sizeof(short));
    for(i=0;i<nrow;i++)
      {
        short tmp=(~y[i]);
        
        Rprintf("x[%i]=%i\n",i,x[i]);
        Rprintf("y[%i]=%i\n",i,y[i]);
        Rprintf("tmp=%i\n",tmp);

        ret[i]=(x[i]&tmp);
      }
    return ret;
  }

void treePop(int* splits, int* ncolp,int* np, int* ed1, int* ed2, double* edLen)
  {
    int n=*np;
    int ncol=*ncolp;
    int nrow=ceil(n/(double)8);
    short split[nrow][ncol];
    int i=0,j=0;

    Rprintf("n=%i nrow=%i ncol=%i\n",n,nrow,ncol);
    Rprintf("got\n");
    for(i=0;i<ncol;i++)
    {
     for(j=0;j<nrow;j++)
       {
          Rprintf("%i ",splits[i*nrow+j]);
       }
     Rprintf("\n");
    }

    for(i=0;i<ncol;i++)
      for(j=0;j<nrow;j++)
       { 
         split[j][i]=(short)splits[i*nrow+j];
       }

    Rprintf("short-ed\n");
    for(i=0;i<nrow;i++)
    {
      for(j=0;j<ncol;j++)
       {
         Rprintf("%i ",split[i][j]);
       }
      Rprintf("\n");
    }

   short vlabs[2*n-1][nrow];
   for(i=0;i<nrow-1;i++)
    {
       vlabs[1][i]=255;
    }
   vlabs[1][nrow-1]=~(pow(2,8-(n%8))-1);
  }
