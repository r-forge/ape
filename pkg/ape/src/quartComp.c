#include <math.h>
#include <R.h>
#include <stdint.h>

int isElement(int a,int* q)
{ int j;
  for(j=0;j<4;j++)
   {
     if(q[j]==a)
      {
       return 1;  
      }
   }
 return 0;
}

void quartComp(int* quarts, double* w,int*mp, int*np,int* ret)
{ *ret=1;
  int i=0,j=0, m=*mp, n=*np,a,b,c,d,e;
  int qt[m][4];
  for(i=0;i<4;i++)
  {
   for(j=0;j<m;j++)
    {
      qt[j][i]=quarts[i*m+j];
    }  
  }

  for(i=0;i<m;i++)
  {
   for(j=0;j<4;j++)
    {
      Rprintf("%i ",qt[i][j]);
    }
   Rprintf("\n");
  }

 for(a=1;a<=n;a++)
  {
   for(b=1;b<=n;b++)
    {if(a!=b)
      for(c=1;c<=n;c++)
       {
        if(a!=c && b!=c)
         for(d=1;d<=n;d++)
	    {
           if(a!=d && b!=d && c!=d)
            {		
              int n0=0;
              for(i=0;i<m;i++)
  		   { int q[4];
                     int ii=0;
                     for(ii=0;ii<4;ii++)
                      {
                        q[ii]=qt[i][ii];
                      }
 		     if(isElement(a,q) && isElement(b,q) && isElement(c,q) && isElement(d,q))
                       {
			  if(w[i]==0)
  		            {
                              n0=n0+1;
                            }
                       }
                  }
              if(n0<2){Rprintf("first condition\n");*ret=0;return;}
             for(e=1;e<=n;e++)
              {
                if(a!=e && b!=e && c!=e && d!=e)
                 {
                   int p1=-1; //ab|cd pos
                   int p2=-1; //bc|de pos
                   int p3=-1; //ab|de pos
                   int p4=-1; //ab|ce pos
                   int p5=-1; //ae|cd pos
                   int p6=-1; //be|cd pos
 			 for(i=0;i<m;i++)
 			  {
		      if(((qt[i][0]==a && qt[i][1]==b)|| (qt[i][1]==a && qt[i][0]==b)) && ((qt[i][2]==d && qt[i][3]==c)||(qt[i][3]==d && qt[i][2]==c)))
			      {
				  p1=i;
			      }
                      if(((qt[i][0]==c && qt[i][1]==b)||(qt[i][1]==c && qt[i][0]==b))&& ((qt[i][2]==d && qt[i][3]==e)||(qt[i][3]==d && qt[i][2]==e)))
			      {
				  p2=i;
			      }
                      if(((qt[i][0]==a && qt[i][1]==b)||(qt[i][1]==a && qt[i][0]==b))&& ((qt[i][2]==d && qt[i][3]==e)||(qt[i][3]==d && qt[i][2]==e)))
			      {
				  p3=i;
		              }
 	             if(((qt[i][0]==a && qt[i][1]==b)||(qt[i][1]==a && qt[i][0]==b))&& ((qt[i][2]==c && qt[i][3]==e)||(qt[i][3]==c && qt[i][2]==e)))
			      {
				  p4=i;
			      }
  	             if(((qt[i][0]==a && qt[i][1]==e)||(qt[i][1]==a && qt[i][0]==e))&& ((qt[i][2]==d && qt[i][3]==c)||(qt[i][3]==d && qt[i][2]==c)))
			      {
				  p5=i;
			      }
                     if(((qt[i][0]==b && qt[i][1]==e)||(qt[i][1]==b && qt[i][0]==e))&& ((qt[i][2]==d && qt[i][3]==c)||(qt[i][3]==d && qt[i][2]==c)))
			      {
				  p6=i;
			      }

                    }
                  if(p1!=-1 && p2!=-1 && p3!=-1 && w[p1]>0 && w[p2]>0)
                    {
                      if(w[p3]!=(w[p1]+w[p2]))
                        { Rprintf("second condition for taxa a=%i b=%i c=%i d=%i e=%i p1=%i p2=%i p3=%i\n",a,b,c,d,e,p1,p2,p3);
			  *ret=0;
                          return;
                        }
                    }
                  if(p1!=-1 && p4!=-1 && p5!=-1 && w[p1]>w[p4] && w[p4]>0)
                    {
                      if(w[p5]!=(w[p1]-w[p4]))
                        {  Rprintf("third condition for taxa a=%i b=%i c=%i d=%i e=%i\n",a,b,c,d,e);
                           *ret=0;
                           return;
                        }
                    }
                   if(p1!=-1 && p3!=-1 && p4!=-1 && p5!=-1 && p6!=-1 && w[p1]>0)
                    {
                     if(((w[p4]<0 || w[p3]<0) && (w[p5]<0 || w[p6]<0)))
                       {Rprintf("fourth condition for taxa a=%i b=%i c=%i d=%i e=%i p4=%i p3=%i p1=%i p5=%i p6=%i\n",a,b,c,d,e,p4,p3,p1,p5,p6);
                        Rprintf("w[p4]=%f w[p3]=%f w[p5]=%f w[p6]=%f\n",w[p4],w[p3],w[p5],w[p6]);
                        *ret=0;
                        return;
                       }
                    }
                 }
                
              }
 		}
          }
       }
    }
  }
}
