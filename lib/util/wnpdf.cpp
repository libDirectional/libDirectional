#include <mex.h> 
#include <matrix.h> 

#define _USE_MATH_DEFINES 
#include <math.h>

/* Probability density function of a wrapped normal distribution.
 *
 * Parameters:
 *  xa (1 x n row vector)
 *      where to evaluate the pdf
 *  mu (scalar)
 *      location parameter in [0,2*pi)
 *  sigma (scalar)
 *      uncertainty parameter >0
 * Returns:
 *  result (1 x n  row vector)
 *      value of the pdf at each point in xa
 */

void mexFunction ( 
                   int nlhs,       mxArray *plhs[], 
                   int nrhs, const mxArray *prhs[] 
                   ) 
{ 
	double *xa, *x, *result, *mu, *sigma;
    size_t cols; 
   
	/* Check for proper number of arguments */ 
	if (nrhs != 3 ) { 
		mexErrMsgTxt ("requires 3 input arguments (xa, mu, sigma)."); 
	} else if (nlhs > 1) { 
		mexErrMsgTxt ("returns just one value."); 
	} 
	
	/* Get arguments from MATLAB */ 
	xa = mxGetPr(prhs[0]); 
	if(mxGetM(prhs[0])!= 1){
        mexErrMsgTxt ("xa must have one row."); 
    }
    cols = mxGetN(prhs[0]);
       
    mu = mxGetPr(prhs[1]); 
	if(mxGetM(prhs[1])!= 1 || mxGetN(prhs[1]) !=1){
        mexErrMsgTxt ("mu must be 1x1."); 
    }
    
    sigma = mxGetPr(prhs[2]); 
	if(mxGetM(prhs[2])!= 1 || mxGetN(prhs[2]) !=1){
        mexErrMsgTxt ("sigma must be 1x1."); 
    }
    
    if(*sigma<=0){
        mexErrMsgTxt ("sigma must be >0."); 
    }
        
       
    plhs[0] = mxCreateDoubleMatrix (1, (int)cols, mxREAL);
	result = mxGetPr(plhs[0]);  
    
    // check if sigma is large and return uniform distribution in this case
    if (*sigma > 10){
        for(int i=0;i<(int)cols;i++){
            result[i] = 1.0/(2*M_PI);
        }
        return;
    }
    
    x = (double*)malloc(cols*sizeof(double));
    
    const int maxIterations = 1000; 
    
    const double tmp = -1.0/2/(*sigma)/(*sigma);
    const double nc = 1/sqrt(2*M_PI)/(*sigma);
        
    #pragma omp parallel for
    for(int i=0;i<(int)cols;i++){ //iterate over columns of xa
        x[i] = fmod(xa[i], 2*M_PI);
        if(x[i]<0){
            x[i] += 2*M_PI;
        }
        x[i] -=  *mu;
        
        double oldResult;
        
        result[i]=exp(x[i]*x[i]*tmp);
        for(int k=1;k<=maxIterations;k++){
            double xp = x[i] + 2*M_PI*k;
            double xm = x[i] - 2*M_PI*k;
            double tp = xp*xp*tmp;
            double tm = xm*xm*tmp;
            oldResult = result[i];
            result[i] += exp(tp) + exp(tm);
            if(result[i]==oldResult){ //break if result is accurate up to double precision
                break;
            }
        }
        //normalize
        result[i] *= nc;
    }
    
    free(x);

}