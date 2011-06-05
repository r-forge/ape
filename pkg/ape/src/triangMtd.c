#include<R.h>

/*
 * leafs labelled 1 to n. root labelled n+1. other nodes labelled n+1 to m
 */

int pred(int k, int* ed1, int* ed2, int numEdges)
//find the predecesor of vertex k
 {  int i=0;
    for(i=0;i<=numEdges;i++)
       {
        if(ed2[i]==k){return ed1[i];}
       }
   return -1;
 }
int* getPathBetween(int x,int y, int n, int* ed1, int* ed2, int numEdges)
//get the path between vertices x and y in an array ord.
//ord[i]=j means that we go between i and j on the path between x and y
 {  int i=0;
    int k=x;
    int ch[2*n-2];//ch[i]==1 implies {k,pred(k)} on path between x and y
    for(i=0;i<2*n-2;i++)
      {ch[i]=0;
      }

    while(k!=n)
      {
        ch[k]=1;
        k=pred(k,ed1,ed2,numEdges);
      }
    k=y;
    while(k!=n)
      {
        ch[k]++;
        k=pred(k,ed1,ed2,numEdges);
      }

    int co=0;

        int *ord=(int*)malloc(sizeof(int)*(2*n-2));
        //starting from x, fill ord


        int p=x;

        while(ch[p]==1)
          {
            int aux=p;
            p=pred(p,ed1,ed2,numEdges);
            ord[aux]=p;
          }
        p=y;
        while(ch[p]==1)
          {
            int aux=p;
            p=pred(p,ed1,ed2,numEdges);
            ord[p]=aux;//other way
          }

      return ord;
 }
int getLength(int x, int y, int* ed1, int* ed2, int numEdges,int* edLen)
//get length of edge {x,y}, -1 if edge does not exist
 {  int i=0;
    for(i=0;i<=numEdges;i++)
      {
        if((ed1[i]==x && ed2[i]==y) || (ed1[i]==y && ed2[i]==x)){return edLen[i];}
      }
    return -1;
 }
void triangMtd(double* dd, int* np, int* ed1,int* ed2, double* edLen)
 {
    int n=*np;   
    int i=0;
    int j=0;  
    int ij=-1;

    double d[n*n];
    for(i=0;i<n;i++)
    {
     for(j=0;j<i;j++)
       {ij++;
         d[i*n+j]=dd[ij];
	 Rprintf("%f ",dd[ij]);
       }
     Rprintf("\n"); 	
    }

    for(i=0;i<n;i++)
    {d[i*n+i]=0;
     for(j=i;j<n;j++)
       {d[i*n+j]=d[j*n+i];
       }
    }
    for(i=0;i<n;i++)
    {
     for(j=0;j<n;j++)
       {Rprintf("%f ",d[i*n+j]);
       }
     Rprintf("\n");
    }

    double minW=0;
    int x=-1,y=-1,numEdges=0;
    int st[n];//array holding at position i the starting point of attachment path of leaf i
    int ed[n];//array holding at position i the ending point of attachment path of leaf i
    double l[n];//array holding at position i the distance of leaf i from the partially constructed tree
    int w[n];//presence array:w[i]=1 if leaf i is added to the tree, 0 otherwise
    int wSize=0;//number of leaves added so far

    for(i=0;i<n;i++){w[i]=0;}
    //find the minimum distance in our matrix
    
    int sw=0;
    //choose two leaves x,y such that d[x][y] is minimal
    for(i=0;i<n;i++)
    {
     for(j=0;j<n;j++)
       {if(d[i*n+j]<=0)continue;
        minW=d[i*n+j];
        x=i;
        y=j;
        sw=1;
        break;
       }
     if(sw)break;
    }
    for(i=0;i<n;i++)
     for(j=0;j<n;j++)
       { if(i==j)continue;
         if(d[i*n+j]<minW)
           {
             minW=d[i*n+j];
             x=i;
             y=j;
           }
       }

    w[x]=1; w[y]=1;//mark x and y as added

    //calculate the distance between leaves not in w and leaves in w
    for(i=0;i<n;i++)
      {
        if(w[i]){continue;}
        st[i]=x;ed[i]=y;
        l[i]=0.5*(d[i*n+x]+d[i*n+y]-d[x*n+y]);//distance of leaf i from the initial tree
      }


    wSize+=3;//since we construct a star tree on three leaves
    int nv=n;//since first n numbers are reserved for leaves
    int minDist=1000;
    int x3=0;
    //search for x3 that is closest to the tree
    for(i=0;i<n;i++)
          {if(w[i])continue;//skip if leaf already added
            if(l[i]<minDist)
              {
                minDist=l[i];
                x3=i;
              }
          }
    //construct initial star tree on three leaves:x,y,x3. nv is the interior
    //vertex and also the root of the tree
    ed1[numEdges]=nv;
    ed2[numEdges]=x;
    edLen[numEdges]=0.5*(d[x*n+x3]+d[x*n+y]-d[x3*n+y]);
    numEdges++;
    ed1[numEdges]=nv;
    ed2[numEdges]=y;
    edLen[numEdges]=0.5*(d[x3*n+y]+d[x*n+y]-d[x*n+x3]);
    numEdges++;
    ed1[numEdges]=nv;
    ed2[numEdges]=x3;
    edLen[numEdges]=0.5*(d[x3*n+x]+d[x3*n+y]-d[x*n+y]);
    w[x3]=1;

    //calculate distance of leaves not yet added to the star tree
    int s;
    for(s=0;s<n;s++)
         {if(w[s])continue;
           for(i=0;i<n;i++)
            {
               if(!w[i])continue;
               double newL=0.5*(d[i*n+s]+d[x3*n+s]-d[i*n+x3]);
               if(newL<l[s])
                  {
                   l[s]=newL;
                   st[s]=i;
                   ed[s]=x3;
                  }
            }
         }

    //partial tree construction begins here

    while(wSize<n)
      { int minDist=1000;
        int z=0;
        int i=0;
        //search for leaf z which is closest to partial tree
        for(i=0;i<n;i++)
          {if(w[i])continue;//skip if leaf already added
            if(l[i]<minDist)
              {
                minDist=l[i];
                z=i;
              }
          }
        //x and y are the leaves on the path between which z is attached
        x=st[z];y=ed[z];


        int* ord=getPathBetween(x,y,n,ed1,ed2,numEdges);


        //order the path from x to y, in an array ord, ord[i]=j means vertex i comes before j
        //first count number of edges on path (i.e count all i s.t ch[i]==1)


        //look for the edge on the path x to y to subdivide
        int p=x;
        double sum=0;
        double prevSum=0;
        int aux=0;
        int subdiv=-1;//index of edge to subdivide
        double lx=0.5*(d[x*n+y]+d[x*n+z]-d[z*n+y]);//distance of attachment point from x
        //printf("Adding leaf %i, between %i and %i\n",z,x,y);
        //printf("%i situated at a distance of %d from tree",z,minDist);
        int sw=0;
        while(p!=y && sum<lx)
          { aux=p;
          //printf("%i\n",p);
            p=ord[p];
            prevSum=sum;
            for(i=0;i<=numEdges;i++)
              {
                if((ed1[i]==aux && ed2[i]==p)||(ed2[i]==aux && ed1[i]==p))
                  {
                    if(ed2[i]==aux && ed1[i]==p){sw=1;}
                    subdiv=i;
                    sum+=edLen[i];
                  }
              }

          }


        nv++;
        //subdivide subdiv with a node labelled nv
        //length calculation
        //multifurcating vertices
        if((sw==1 && sum-lx == 0) || (sw==0 && lx-prevSum == 0))
        {
            numEdges++;
            ed1[numEdges]=ed1[subdiv];
            ed2[numEdges]=z;
            edLen[numEdges]=minDist;
        }else{
        int edd=ed2[subdiv];
        ed2[subdiv]=nv;
        edLen[subdiv]= (sw==1)?(lx-prevSum):(sum-lx);//check which 'half' of the
                                                     //path the leaf belongs to
                                                     //and updates accordingly
        numEdges++;
        ed1[numEdges]=nv;
        ed2[numEdges]=edd;
        edLen[numEdges]= (sw==1)?(sum-lx):(lx-prevSum);
        numEdges++;
        edLen[numEdges]=minDist;
        ed1[numEdges]=nv;
        ed2[numEdges]=z;
        }



        wSize++;
        w[z]=1;



        /*update distance matrix, only needed for incomplete distances
        int ii;
      for(ii=0;ii<n;ii++)
      {if(!w[ii])continue;
       printf("path between %i and %i\n",ii,z);
        int* ord=getPathBetween(ii,z,n,ed1,ed2,numEdges);

        p=ii;
        int newDist=0;
        printf("path");
        for(i=0;i<2*n-2;i++)
          {printf("ord[%i]=%i ",i,ord[i]);
          }
       printf("\n");
       printf("innn\n");
        while(p!=z)
          { //printf("p=%i\n",p);
            int aux=p;
            p=ord[p];
            newDist+=getLength(ii,z,ed1,ed2,numEdges,edLen);
          }

        printf("outt\n");

        d[ii][z]=d[z][ii]=newDist;
       }*/
        //update l[s] for all s not yet added
        int s;
        for(s=0;s<n;s++)
         {if(w[s])continue;
           for(i=0;i<n;i++)
            {
               if(!w[i])continue;
               int newL=0.5*(d[i*n+s]+d[z*n+s]-d[i*n+z]);//one of leaves is z, since
                                                      //all pairs not cotaining z
                                                      //will remain unchanged
               if(newL<l[s])
                  {
                   l[s]=newL;
                   st[s]=i;
                   ed[s]=z;
                  }
            }
         }
        free(ord);
        /*printf("new tree\n");
        for(i=0;i<2*n-3;i++)
        {
        printf("%i->%i of length:%d \n",ed1[i],ed2[i],edLen[i]);
        }
        printf("end new tree\n");*/
      }
   for(i=0;i<=numEdges;i++){ed1[i]++;ed2[i]++;}
 }