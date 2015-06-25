#include <mex.h>

/* MEX function [s,coef]=mhg([MAX,K],alpha,p,q,x,y)
 computes the truncated hypergeometric function pFq ^alpha(p;q;x;y)
 of one or two matrix arguments.
 
 The sum is only over those partitions kappa such that 
 a) |kappa|<=MAX 
 b) kappa has no more than n=length(x) parts
 c) kappa_1<=K (only if K is specified, if K is omitted, this restriction 
    is not applied)
 
 p and q are arrays, so mhg(30,9,[3 4],[5 6 7],[0.5 0.6],[0.8,0.9]) is 
 2F3^9([3 4],[5 6 7];[0.5 0.6], [0.8 0.9]) summed over all kappa with 
 |kappa|<=30

 K and y may be omitted. If y is omitted, the hypergeometric function of one 
 matrix argument is computed.

 Copyright, Plamen Koev, Massachusetts Institute of Technology
 Written: May 2004. 
 
 Updated:
 
 December 2005: introduced coef -- array from 0 to MAX, such that 
                coef[i]=sum over all partitions of size i. 
                Then c=sum_{i=0}^MAX coef[i].
 
 May 2007:      The Jacks are now Schur normalized leading to 
                   easier updates, i.e., Sx_kappa=Jx_kappa/H^*_kappa
                   where H^_kappa is the product of upper hook lengths
                Introduced Raymond Kan's faster Q_kappa and 
                   beta_{lambda mu} updates; 
                Removed all recursions; 
                Introduced a faster update for partitions with 
                   exactly n parts from Kaneko (1993) eq. (5.3);

 June 2007 ideas to consider: 
        1. malloc not more than choose(n+K,n) space if K is specified
        2. can we get J^(1/alpha) from J^(alpha), Stanley 89, Cor. 3.5<
*/


#include <stdio.h>
#include <stdlib.h>
#include <memory.h>

void mexFunction(int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[] )
{  int    h,sl,nmu,i,j,k,*f;
   double *x, *Sx, *xn, *prodx;
   double *y, *Sy, *yn, *prody;
   double c, *s, alpha, *p, *q, *z, *coef, *kt, *mt, *blm, zn, dn,
          t, q1, q2, cc;

   int lg, n, np, nq, MAX, *D, heap, *l, w, *mu, *d, *g, *ww, *lmd,
        slm, nhstrip, gz;
   
  if ((nrhs<5)||(nrhs>6))
    mexErrMsgTxt("Must have five or six input arguments.");
  if (nlhs>2)
    mexErrMsgTxt("Too many output arguments.");
  
  
  if (!mxIsDouble(prhs[0]) || !mxIsDouble(prhs[1]) || !mxIsDouble(prhs[2]) || 
      !mxIsDouble(prhs[3]) || !mxIsDouble(prhs[4]))
    mexErrMsgTxt("All inputs must be double arrays.");

  if (nrhs==6) if (!mxIsDouble(prhs[5])) 
    mexErrMsgTxt("All inputs must be double arrays.");

  if (mxGetM(prhs[0])!=1 || mxGetN(prhs[0])<1 || mxGetN(prhs[0])>2)
    mexErrMsgTxt("First input must be a scalar = max size of partition.");
  if (mxGetM(prhs[1])!=1 || mxGetN(prhs[1])!=1)
    mexErrMsgTxt("Second input must be a scalar = alpha.");
  if (mxGetM(prhs[2])>1 && mxGetN(prhs[2])>1)
    mexErrMsgTxt("Third input must be a vector = (a_1,...,a_p).");
  if (mxGetM(prhs[3])>1 && mxGetN(prhs[3])>1)
    mexErrMsgTxt("Forth input must be a vector = (b_1,...,b_q).");
  if (mxGetM(prhs[4])!=1 && mxGetN(prhs[4])!=1)
    mexErrMsgTxt("Fifth input must be a vector = x.");

  if (nrhs==6)
  if (mxGetM(prhs[5])!=1 && mxGetN(prhs[5])!=1)
    mexErrMsgTxt("Sixth input must be a vector = y.");

  plhs[0]=mxDuplicateArray(prhs[1]);   /* output is scalar */
  
  s=mxGetPr(plhs[0]); 
  
  x=mxGetPr(prhs[4]);
  n=mxGetNumberOfElements(prhs[4]);  /* size of x */
  
  if (nrhs==6) y=mxGetPr(prhs[5]);
  else y=NULL;
  p=mxGetPr(prhs[2]);
  np=mxGetNumberOfElements(prhs[2]);  /* size of p */
  q=mxGetPr(prhs[3]);
  nq=mxGetNumberOfElements(prhs[3]);  /* size of q */


  MAX=(int) *mxGetPr(prhs[0]);
  alpha=*mxGetPr(prhs[1]);
  
  if (nlhs==2) {
     plhs[1]=mxCreateDoubleMatrix(1,MAX+1,mxREAL); 
     coef=mxGetPr(plhs[1]); 
  }
  else coef=(double*)mxMalloc(sizeof(double)*(MAX+1)); 
  
  for (i=1;i<=MAX;i++) coef[i]=0;  /* set to zero, these are the coefficients of the polynomial */ 
  coef[0]=1; /* free term equals one */
  
  w = 0; /* index of the zero partition, currently l*/

  /* figure out the number of partitions |kappa|<= MAX with at most n parts */

  f=(int*) mxCalloc (MAX+2,sizeof(int));
  for (i=1;i<=MAX+1;i++) f[i]=i;
	
  for (i=2;i<n;i++) for (j=i+1;j<=MAX+1;j++) f[j]+=f[j-i];

  w=f[MAX+1];
  mxFree(f);
  
  D     = (int*) mxCalloc(w+1,sizeof(int));
  Sx    = (double*) mxCalloc(n*(w+1),sizeof(double));
  xn    = (double*) mxMalloc(sizeof(double)*(n+1)*(MAX+2));
  prodx = (double*) mxMalloc(sizeof(double)*(n+1));
  prodx[1]=x[0];
  for (i=2;i<=n;i++) prodx[i]=prodx[i-1]*x[i-1];
  for (i=1; i<=n; i++) {
    Sx[n+i-1]=1;
    xn[(MAX+2)*i+1]=1;
    for (j=2;j<=MAX+1;j++) xn[(MAX+2)*i+j]=xn[(MAX+2)*i+j-1]*x[i-1];
  }

  if (y!=NULL) {
     Sy    = (double*) mxCalloc(n*(w+1),sizeof(double));
     yn    = (double*) mxMalloc(sizeof(double)*(n+1)*(MAX+2));
     prody = (double*) mxMalloc(sizeof(double)*(n+1));
     prody[1]=y[0];
     for (i=2;i<=n;i++) prody[i]=prody[i-1]*y[i-1];
     
     for (i=1; i<=n; i++) {
        Sy[n+i-1]=1;
        yn[(MAX+2)*i+1]=1;
        for (j=2;j<=MAX+1;j++) yn[(MAX+2)*i+j]=yn[(MAX+2)*i+j-1]*y[i-1];
     }
  }
  
  l     = (int*)mxCalloc(n+1,sizeof(int));
  
  
  if (mxGetNumberOfElements(prhs[0])==2) l[0]=(int) *(mxGetPr(prhs[0])+1);
     /* l[0] is the second element of MAX, when MAX has 2 elements */
  else l[0]=MAX; 
  /* this is what limits l[1] by the second element of MAX if needed and 
     allows for the check l[i]<l[i-1] to be OK even for i=1 */
  
  z     = (double*)mxMalloc((n+1) * sizeof(double));
  for (i=1;i<=n;i++) z[i]=1;
  
  mu    = (int*)mxCalloc(n+1,sizeof(int));
  kt    = (double*)mxCalloc(n+1,sizeof(double));
  for (i=1;i<=n;i++) kt[i]=-i;
  
  ww    = (int*)mxMalloc((n+1) * sizeof(int));
  for (i=1;i<=n;i++) ww[i]=1;
  
  d     = (int*)mxCalloc(n,sizeof(int));
  g     = (int*)mxCalloc(n+1,sizeof(int));
  mt    = (double*)mxCalloc(n+1,sizeof(double));
  blm   = (double*)mxCalloc(n+1,sizeof(double));
  lmd   = (int*) mxCalloc(n+1,sizeof(int));
  
  heap  = MAX+2;
  cc=1;
  h=1;
  sl=1;  /* sl= sum(l) */
  
  while (h>0) {
      if ((l[h]<l[h-1]) && (MAX>=sl) && (z[h]!=0)) {
          
          l[h]++;
          
          if ((l[h]==1) && (h>1) && (h<n)) {
              D[ww[h]]=heap;
              ww[h]=heap;
              k=MAX-sl+l[h];
              if (k>l[h-1]) k=l[h-1];
              heap+=k;
          }
          else ww[h]++;
          w=ww[h];
              
          /* Update Q */
          c=(1-h)/alpha+l[h]-1;
          zn=alpha;
          dn=kt[h]+h+1;
          for (j=0;j<np;j++)  zn*=p[j]+c;
          for (j=0;j<nq;j++)  dn*=q[j]+c;
          if (y!=NULL) {
              zn*=alpha*l[h];
              dn*=n+alpha*c;
              for (j=1;j<h;j++) {
                  t=kt[j]-kt[h];
                  zn*=t;
                  dn*=t-1;
              }
              zn/=dn;
              dn=1; /* trying to prevent overflow */
          }
          kt[h]+=alpha;
          for (j=1;j<h;j++) {
              t=kt[j]-kt[h];
              zn*=t;
              dn*=t+1;
          }
          z[h]*=zn/dn;
          
          /* Working hard only when l has less than n parts */
          
          
          if (h<n) {
              t=h+1-alpha; cc=1; for (j=1; j<=h;j++) cc*=(t+kt[j])/(h+kt[j]); 
              
              /* computing the index of l-ones(1,h) */
              nmu=l[1]; k=2; while ((k<=h)&&(l[k]>1)) nmu=D[nmu]+l[k++]-2;
                
              Sx[w*n+h-1]=cc*prodx[h]*Sx[nmu*n+h-1];

              if (y!=NULL) Sy[w*n+h-1]=cc*prody[h]*Sy[nmu*n+h-1];
              cc=1; /* this way we can just update from 1 in the h=n case*/

              d[h-1]--; /* technically this has to execute only when h>1 
                           but is OK if it is always executed; d[0] will 
                           end up being -MAX at the end of the code */

              d[h]=l[h];  /* for (k=1;k<h;k++) d[k]=l[k]-l[k+1]; 
                             this happens automatically now via updates */
              
              lg=0; for (k=1;k<=h;k++) if (d[k]>0) {lg++; g[lg]=k;}
              slm=1; /* this is sum(l-mu) */
              nhstrip=1; for (k=1;k<=lg;k++) nhstrip*=d[g[k]]+1; nhstrip--;

              memcpy(&mu[1],&l[1],sizeof(int)*h);
              memcpy(&mt[1],&kt[1],sizeof(double)*h);
              for (k=1;k<=lg;k++) { blm[k]=1; lmd[k]=l[g[k]]-d[g[k]]; }
              
              for (i=1;i<=nhstrip;i++) {
                  j=lg;
                  gz=g[lg];
                  while (mu[gz]==lmd[j]) {
                      mu[gz]=l[gz];
                      mt[gz]=kt[gz];
                      slm-=d[gz];
                      j--;
                      gz=g[j];
                  }
                  t=kt[gz]-mt[gz];
                  
                  zn=1+t;
                  dn=t+alpha;
                  for (k=1; k<gz; k++) {
                      q1=mt[k]-mt[gz];
                      q2=kt[k]-mt[gz];
                      zn*=(alpha-1+q1)*(1+q2);
                      dn*=q1*(alpha+q2); 
                  }
                  blm[j]*=zn/dn;
                                    
                  mu[gz]--;
                  mt[gz]-=alpha;
                  slm++;

                  for (k=j+1;k<=lg;k++) blm[k]=blm[j];    
                  
                  /* next, find the index of mu */
                  nmu=mu[1]+1; for (k=2;k<=h-(mu[h]==0);k++) nmu=D[nmu]+mu[k]-1;
                  
                  for (k=h+1; k<=n;k++) 
                      Sx[w*n+k-1]+=blm[j]*Sx[nmu*n+k-2]*xn[k*(MAX+2)+slm];
                  
                  if (y!=NULL) for (k=h+1; k<=n;k++)
                      Sy[w*n+k-1]+=blm[j]*Sy[nmu*n+k-2]*yn[k*(MAX+2)+slm];
             }
             
             for (k=h; k<n; k++) Sx[w*n+k]+=Sx[w*n+k-1];
             if (y!=NULL) {
                  for (k=h; k<n; k++) Sy[w*n+k]+=Sy[w*n+k-1];
                  coef[sl]+=z[h]*Sx[w*n+n-1]*Sy[w*n+n-1];
             }
             else coef[sl]+=z[h]*Sx[w*n+n-1];

          } /* of "if h<n" */
          else {
              /* computing the index of the partition l-l[n]*ones(1,n) */
              nmu=l[1]-l[n]+1;
              k=2; while ((k<n)&&(l[k]>l[n])) nmu=D[nmu]+l[k++]-1-l[n];
              /* cc is 1 if l[n]==1, (guaranteed by the h<n case); 
                 we then update from the previous */ 
            
              if (y!=NULL) {
                  t=(1/alpha+l[n]-1)/l[n];
                  for (k=1;k<n;k++) t*=(1+kt[k]-kt[n])/(alpha+kt[k]-kt[n]);
                  cc*=t*t*prodx[n]*prody[n];
                  coef[sl]+=z[n]*cc*Sx[nmu*n+n-1]*Sy[nmu*n+n-1];
              }
              else {
                  cc*=(1/alpha+l[n]-1)*prodx[n]/l[n];
                  for (k=1;k<n;k++) cc*=(1+kt[k]-kt[n])/(alpha+kt[k]-kt[n]);
                  coef[sl]+=z[n]*cc*Sx[nmu*n+n-1];
              }
          }
          if (h<n) {
              z[h+1]=z[h];
              h++;
              ww[h]=w;
          }
          sl++;
      }
      else { 
          sl-=l[h];
          l[h]=0;
          kt[h]=-h;
          h--;
      }
  } /* of while h>0 */
  
  *s=0; for (i=0;i<MAX+1;i++) (*s)+=coef[i];  

  mxFree(lmd);
  mxFree(blm);
  mxFree(mt);
  mxFree(g);
  mxFree(d);
  mxFree(ww);
  mxFree(kt);
  mxFree(mu);
  mxFree(z);
  mxFree(l);
  if (y!=NULL) { mxFree(prody); mxFree(yn); mxFree(Sy); }
  mxFree(prodx);
  mxFree(xn);
  mxFree(Sx);
  mxFree(D);
  if (nlhs!=2) mxFree(coef); 
}

