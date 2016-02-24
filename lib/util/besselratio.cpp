#include <mex.h> 
#include <matrix.h> 
#include <math.h>
   
// Number of iterations
#define N 10

/*
 *  Computes the Ratio of Bessel-I Functions
 *  The Ratio I_{v+1}(x)/I_v(x) is computed, where
 *  I_v(x) is the Bessel-I function of order v evaluated at x.
 *  The algorithm is based on Amos(1974)
 */

void mexFunction ( 
                   int nlhs,       mxArray *plhs[], 
                   int nrhs, const mxArray *prhs[] 
                   ) 
{ 
	double *xp, *vp, *result;
    double x;
    int v;
    int ord, i, k;
    double r[N+1];
   
	/* Check for proper number of arguments */ 
	if (nrhs != 2 ) { 
		mexErrMsgTxt ("requires 2 input arguments."); 
	} else if (nlhs > 1) { 
		mexErrMsgTxt ("returns just one value."); 
	} 
	
	/* Get arguments from MATLAB */ 
	vp = mxGetPr(prhs[0]); 
	if(mxGetM(prhs[0])!= 1 || mxGetN(prhs[0]) !=1){
        mexErrMsgTxt ("v must be 1x1."); 
    }
    xp = mxGetPr(prhs[1]); 
	if(mxGetM(prhs[1])!= 1 || mxGetN(prhs[1]) !=1){
        mexErrMsgTxt ("x must be 1x1."); 
    }
   
    plhs[0] = mxCreateDoubleMatrix (1, 1, mxREAL);
	result = mxGetPr(plhs[0]);  
    *result = 0;
    
    v=(int)(*vp);
    x=*xp;
    
    if (x==0){ //special case, ratio is zero
        return;
    }
    
    //Ensure faster convergence
    if (v<10){
        ord=10;
    }
    else{
        ord=v;
    }

    // Initial calulation of r^0_a(x) for a=ord, ..., ord + N
    // See Amos(1974), eq. (20a).
    //r=zeros(1,N+1);
    for(i=0; i<=N; i++){
        r[i]=x/(ord + i + 0.5 + sqrt((ord + i +1.5)*(ord + i +1.5) + x*x));
	}

    // Calculation of r^N_v(x)
    // See Amos(1974), eq. (20b)
    for(i=1; i<=N; i++){
        for(k=0; k<=(N-i); k++){
            r[k]=x/(ord + k + 1 + sqrt((ord + k +1)*(ord + k +1) + x*x*r[k+1]/r[k]));
        }
    }

    *result = r[0];

    // Perform order reduction, see Amos(1974) eq. (2).
    for(i=ord; i>v; i--){
        *result = 1/(2*i/x + *result);
    }

}