#include <mex.h> 
#include <matrix.h> 

#define _USE_MATH_DEFINES 
#include <math.h>

/* Probability density function of a bivariate wrapped normal distribution.
 *
 * Parameters:
 *  xa (2 x m matrix)
 *      where to evaluate the pdf
 *  mu (2 x 1 column vector)
 *      location parameter in [0,2*pi)
 *  sigma (2 x 2 matrix)
 *      uncertainty parameter, symmetric pos. def. matrix
 *  n (scalar integer >= 0)
 *      determines number of summands to consider (per dimension, -n to n)
 * Returns:
 *  result (1 x m  row vector)
 *      value of the pdf at each column of xa
 */

void mexFunction ( 
                   int nlhs,       mxArray *plhs[], 
                   int nrhs, const mxArray *prhs[] 
                   ) 
{ 
	double *xa, *result, *mu, *C;
    int cols;
    int n;
    double Cinv[4];
   
	/* Check for proper number of arguments */ 
	if (nrhs != 4 ) { 
		mexErrMsgTxt ("requires 4 input arguments (xa, n, mu, C)."); 
	} else if (nlhs > 1) { 
		mexErrMsgTxt ("returns just one value."); 
	} 
	
	/* Get arguments from MATLAB */ 
	xa = mxGetPr(prhs[0]); 
	if(mxGetM(prhs[0])!= 2){
        mexErrMsgTxt ("v must have two rows."); 
    }
    cols = (int)mxGetN(prhs[0]);
    
    n = int(*mxGetPr(prhs[1])); 
	if(mxGetM(prhs[1])!= 1 || mxGetN(prhs[1]) !=1){
        mexErrMsgTxt ("n must be 1x1."); 
    }
    
    mu = mxGetPr(prhs[2]); 
	if(mxGetM(prhs[2])!= 2 || mxGetN(prhs[2]) !=1){
        mexErrMsgTxt ("mu must be 2x1."); 
    }
    
    C = mxGetPr(prhs[3]); 
	if(mxGetM(prhs[3])!= 2 || mxGetN(prhs[3]) !=2){
        mexErrMsgTxt ("C must be 2x2."); 
    }
    
    const double detC = C[3]*C[0]-C[1]*C[2];
    const double sqrtDetC = sqrt(detC);
    if(detC==0){
        mexErrMsgTxt ("C must be invertible");
    }
    
    Cinv[0] = C[3]/detC;
    Cinv[1] = -C[1]/detC;
    Cinv[2] = -C[2]/detC;
    Cinv[3] = C[0]/detC;
   
    plhs[0] = mxCreateDoubleMatrix (1, mxGetN(prhs[0]), mxREAL);
	result = mxGetPr(plhs[0]);  
    
    #pragma omp parallel for
    for(int i=0;i<cols;i++){ //iterate over columns of xa
        result[i]=0;
        for(int j=-n;j<=n;j++){
            for(int k=-n;k<=n;k++){
                double x = xa[0+2*i] + 2*M_PI*j - mu [0];
                double y = xa[1+2*i] + 2*M_PI*k - mu [1];
                double t = -0.5*(x*x*Cinv[0] + 2*x*y*Cinv[1] + y*y*Cinv[3]);
                result[i] += exp(t);
            }
        }
        //normalize
        result[i] /= (2*M_PI*sqrtDetC);
    }

}