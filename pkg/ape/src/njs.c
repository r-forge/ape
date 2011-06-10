/* nj.c       2010-04-21 */

/* Copyright 2006-2010 Emmanuel Paradis

/* This file is part of the R-package `ape'. */
/* See the file ../COPYING for licensing issues. */

#include <R.h>

#define DINDEX(i, j) n*(i - 1) - i*(i - 1)/2 + j - i - 1
/* works if i < j strictly, and i = 1, ...;
   see give_index() below */

void njs(double *D, int *N, int *edge1, int *edge2, double *edge_length)
{       //assume missing values are denoted by -1
	double *R, Sdist, Ndist, *new_dist, A, B, smallest_S, x, y;
	int n, i, j, k, ij, smallest, OTU1, OTU2, cur_nod, o_l, *otu_label;

        int *s;//s contains |Sxy|, which is all we need for agglomeration
        double *newR;
        int *newS;

	R = &Sdist;
	new_dist = &Ndist;
	otu_label = &o_l;
        n = *N;
	cur_nod = 2*n - 2;

	R = (double*)R_alloc(n*(n - 1)/2, sizeof(double));
        newR = (double*)R_alloc(n*(n - 1)/2, sizeof(double));
	new_dist = (double*)R_alloc(n*(n - 1)/2, sizeof(double));
	otu_label = (int*)R_alloc(n + 1, sizeof(int));
        s = (int*)R_alloc(n*(n - 1)/2, sizeof(int));
        newS = (int*)R_alloc(n*(n - 1)/2, sizeof(int));

	for (i = 1; i <= n; i++) otu_label[i] = i; /* otu_label[0] is not used */

	k = 0;

        //compute Sxy and Rxy
        for(i=0;i<n*(n-1)/2;i++)
          {newR[i]=0;
           newS[i]=0;
           s[i]=0;
           R[i]=0;
          }

        for(i=1;i<n;i++)
         for(j=i+1;j<=n;j++)
         {//algorithm assumes i,j /in Sij, so skip pair if it is not known
          if(D[give_index(i,j,n)]==-1)
            {              
              continue;
            }
          for(k=1;k<=n;k++)
           {//ij is the pair for which we compute
            //skip k if we do not know the distances between it and i AND j
            
              if(D[give_index(i,k,n)]==-1 || D[give_index(j,k,n)]==-1)continue;             
              s[give_index(i,j,n)]++;
              if(i!=k)
              R[give_index(i,j,n)]+=D[give_index(i,k,n)];
              if(j!=k)
              R[give_index(i,j,n)]+=D[give_index(j,k,n)];
           }
         }



        k=0;
	while (n > 3) {

		

		ij = 0;
		smallest_S = -1e50;
                for(i=0;i<n*(n-1)/2;i++)
                  {newR[i]=0;
                   newS[i]=0;
                  }
		B = n - 2;
		for (i = 1; i < n; i++) {
			for (j = i + 1; j <= n; j++) {
                            
                                A=((R[give_index(i,j,n)]*100)/(s[give_index(i,j,n)]-2))-D[give_index(i,j,n)]*100;
                                //Rprintf("%f/%i - %f",R[give_index(i,j,n)],s[give_index(i,j,n)]-2,D[give_index(i,j,n)]);
                                //Rprintf("Q[%i,%i]=%f\n",i,j,A);
				if (A > smallest_S) {
					OTU1 = i;
					OTU2 = j;
					smallest_S = A;
					smallest = ij;
				}
				ij++;
			}
		}

                //update Rxy and Sxy
               /* Rprintf("agglomerating %i and %i, Q=%f \n",OTU1,OTU2,smallest_S);
                
                for(i=1;i<n;i++)
                  {
                    for(j=i+1;j<=n;j++)
                      {
                        Rprintf("R[%i,%i]=%f ",i,j,R[give_index(i,j,n)]);
                      }
                    Rprintf("\n");
                  }

                for(i=1;i<n;i++)
                  {
                    for(j=i+1;j<=n;j++)
                      {
                        Rprintf("s[%i,%i]=%i ",i,j,s[give_index(i,j,n)]);
                      }
                    Rprintf("\n");
                  }*/

                for(i=1;i<n;i++)
                {if(i==OTU1 || i==OTU2)continue;
                 for(j=i+1;j<=n;j++)
                  {if(j==OTU1 || j==OTU2)continue;
                     if(D[give_index(i,OTU1,n)]!=-1 && D[give_index(j,OTU1,n)]!=-1)
                      {//OTU1 was considered for Rij, so now subtract
                       R[give_index(i,j,n)]-=(D[give_index(i,OTU1,n)]+D[give_index(j,OTU1,n)]);
                       s[give_index(i,j,n)]--;
                      }
                     if(D[give_index(i,OTU2,n)]!=-1 && D[give_index(j,OTU2,n)]!=-1)
                      {//OTU2 was considered for Rij, so now subtract
                       R[give_index(i,j,n)]-=(D[give_index(i,OTU2,n)]+D[give_index(j,OTU2,n)]);
                       s[give_index(i,j,n)]--;
                      }
                  }
                }
		edge2[k] = otu_label[OTU1];
		edge2[k + 1] = otu_label[OTU2];
		edge1[k] = edge1[k + 1] = cur_nod;

		/* get the distances between all OTUs but the 2 selected ones
		   and the latter:
		   a) get the sum for both
		   b) compute the distances for the new OTU */

                double sum=0;

                for(i=1;i<=n;i++)
                 {if(i==OTU1 || i==OTU2)continue;
                 if(D[give_index(OTU1,i,n)]==-1 || D[give_index(OTU2,i,n)]==-1)continue;
                 sum+=(D[give_index(OTU1,i,n)]-D[give_index(OTU2,i,n)]);
                 }
                sum*=(1.0/(2*(s[give_index(OTU1,OTU2,n)]-2)));
                double dxy=D[give_index(OTU1,OTU2,n)]/2;

                //Rprintf("R[%i,%i]:%f \n",OTU1,OTU2,sum);
		edge_length[k] = dxy+sum;//OTU1
               // Rprintf("l1:%f \n",edge_length[k]);
		edge_length[k + 1] = dxy-sum;//OTU2
               // Rprintf("l2:%f \n",edge_length[k+1]);

		A = D[give_index(OTU1,OTU2,n)];
		ij = 0;
		for (i = 1; i <= n; i++) {
			if (i == OTU1 || i == OTU2) continue;
                        if(D[give_index(OTU1,i,n)]!=-1 && D[give_index(OTU2,i,n)]!=-1)
                         {
                            new_dist[ij]=0.5*(D[give_index(OTU1,i,n)]-edge_length[k]+D[give_index(OTU2,i,n)]-edge_length[k+1]);
                         }else{
                         if(D[give_index(OTU1,i,n)]!=-1)
                                {
                                 new_dist[ij]=D[give_index(OTU1,i,n)]-edge_length[k];
                                }else{
                                      if(D[give_index(OTU2,i,n)]!=-1)
                                        {
                                            new_dist[ij]=D[give_index(OTU2,i,n)]-edge_length[k+1];
                                        }else{new_dist[ij]=-1;}
                                     }
                              }

			ij++;
		}

                for (i = 1; i < n; i++) {
			if (i == OTU1 || i == OTU2) continue;
			for (j = i + 1; j <= n; j++) {
				if (j == OTU1 || j == OTU2) continue;
				new_dist[ij] = D[DINDEX(i, j)];
				ij++;
			}
		}

                /*for(i=1;i<n-1;i++)
                {
                  for(j=i+1;j<=n-1;j++)
                   {Rprintf("%f ",new_dist[give_index(i,j,n-1)]);
                   }
                  Rprintf("\n");
                }*/
                //compute Rui
                ij=0;
                for(i=2;i<n;i++)
                  {
                   ij++;
                   if(new_dist[give_index(i,1,n-1)]==-1)continue;
                   
                   for(j=1;j<n;j++)
                     {
                       if(new_dist[give_index(i,j,n-1)]!=-1 && new_dist[give_index(1,j,n-1)]!=-1)
                        {                        
                          newS[give_index(1,i,n-1)]++;
                          if(i!=j)
                          newR[give_index(1,i,n-1)]+=new_dist[give_index(i,j,n-1)];
                          if(1!=j)
                          newR[give_index(1,i,n-1)]+=new_dist[give_index(1,j,n-1)];
                        }
                     }
                  }
                //fill in the rest of R and S
                for(i=1;i<n;i++)
                {if(i==OTU1 || i==OTU2)continue;
                 for(j=i+1;j<=n;j++)
                  {if(j==OTU1 || j==OTU2)continue;                                    
                   newR[ij]=R[give_index(i,j,n)];
                   newS[ij]=s[give_index(i,j,n)];
                   ij++;
                  }
                }
                //update newR and newS with the new taxa

                for(i=2;i<n-1;i++)
                {if(new_dist[give_index(1,i,n-1)]==-1)continue;
                 for(j=i+1;j<=n-1;j++)
                  {if(new_dist[give_index(1,j,n-1)]==-1)continue;
                   newR[give_index(i,j,n-1)]+=(new_dist[give_index(1,i,n-1)]+new_dist[give_index(1,j,n-1)]);
                   newS[give_index(i,j,n-1)]++;
                  }
                }
		/* compute the branch lengths */
		 
                
                
		/* update before the next loop
		   (we are sure that OTU1 < OTU2) */
		if (OTU1 != 1)
			for (i = OTU1; i > 1; i--)
				otu_label[i] = otu_label[i - 1];
		if (OTU2 != n)
			for (i = OTU2; i < n; i++)
				otu_label[i] = otu_label[i + 1];
		otu_label[1] = cur_nod;

		

		n--;
		for (i = 0; i < n*(n - 1)/2; i++)
                  {
                    D[i] = new_dist[i];
                    R[i] = newR[i];
                    s[i] = newS[i];
                  }
		cur_nod--;
		k = k + 2;
	}
        int dK=0;//number of known distances in final distance matrix
        int iUK=-1;//index of unkown distance, if we have one missing distance
        int iK=-1;//index of only known distance, only needed if dK==1
        for (i = 0; i < 3; i++) {
		edge1[*N*2 - 4 - i] = cur_nod;
		edge2[*N*2 - 4 - i] = otu_label[i + 1];
                if(D[i]!=-1){dK++;iK=i;}else{iUK=i;}
	}        
        if(dK==2)
         {//if two distances are known: assume our leaves are x,y,z, d(x,z) unknown
          //and edge weights of three edges are a,b,c, then any b,c>0 that
          //satisfy c-b=d(y,z)-d(x,y) a+c=d(y,z) are good edge weights, but for
          //simplicity we assume a=c if d(yz)<d(xy) a=b otherwise, and after some
          //algebra we get that we can set the missing distance equal to the
          //maximum of the already present distances
            int max=-1e50;
          for(i=0;i<3;i++)
            {if(i==iUK)continue;
             if(D[i]>max)max=D[i];
            }
          D[iUK]=max;          
         }
        if(dK==1)
         {//through similar motivation as above, if we have just one known distance
          //we set the other two distances equal to it
          for(i=0;i<3;i++)
            {if(i==iK)continue;
             D[i]=D[iK];
            }
         }
        if(dK==0)
         {//no distances are known, we just set them to 1
          for(i=0;i<3;i++)
           {D[i]=1;
           }
         }
        edge_length[*N*2 - 4] = (D[0] + D[1] - D[2])/2;
	edge_length[*N*2 - 5] = (D[0] + D[2] - D[1])/2;
	edge_length[*N*2 - 6] = (D[2] + D[1] - D[0])/2;
}

