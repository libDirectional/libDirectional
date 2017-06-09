#include "mex.h"
#include "hypergeometricRatioInverseLookUp.h"
#define EPSILON 1e-15
#define MAXKAPPA 1000.0

double hypergeometricRatioInverse(double lambda, unsigned int dim)
{
    unsigned int dimLambdaOffset;
    unsigned int dimCoeffOffset;
    unsigned int spline;
    double kappa = -1;
    if ((dim < INVUPSILON__maxDimension+1) && (dim>1)) {
        dimLambdaOffset=(dim-2)*INVUPSILON__splines;
        dimCoeffOffset=(dim-2)*INVUPSILON__splines*4;
        //mexPrintf("___ %u\n", dimLambdaOffset+INVUPSILON__splines-1);
        //mexPrintf("___ %f\n", INVUPSILON__slicesLambda[dimLambdaOffset+INVUPSILON__splines-1]);
        if (lambda<1.0/dim) {
            return 0.0;
        } else if ((1.0-lambda)<EPSILON) {
            return MAXKAPPA;
        } else if (lambda>INVUPSILON__slicesLambda[dimLambdaOffset+INVUPSILON__splines-1]) {
            return (dim-1.0)/((1.0-lambda));
        } else {    
            for(spline=0; spline<INVUPSILON__splines; spline++) {
                if(lambda<INVUPSILON__slicesLambda[dimLambdaOffset+spline+1]) {
                    kappa= ((INVUPSILON__slicesCoeff[dimCoeffOffset+spline*4]*lambda+
                            INVUPSILON__slicesCoeff[dimCoeffOffset+spline*4+1])*lambda+
                            INVUPSILON__slicesCoeff[dimCoeffOffset+spline*4+2])*lambda+
                            INVUPSILON__slicesCoeff[dimCoeffOffset+spline*4+3];
                    break;
                }
            }       
            if(kappa<0.0) {
                mexErrMsgTxt("kappa<0! Something is wrong with your InvUpsilon lookup table.\n");
                kappa=0.0;
            }     
            if(kappa>MAXKAPPA) {
                mexErrMsgTxt("kappa>0! Something is wrong with your InvUpsilon lookup table.\n");
                kappa=MAXKAPPA;
            }     
            return kappa;
        }
    } else {
        mexErrMsgTxt("Invalid number of dimensions.\n");
        return -1.0;
    }            
}

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
    //Define 
    const double *Lambda;    
    double MaxKappa;    
    size_t NumberOfElements;
    unsigned int nb;
    unsigned int Dimensions;
    double *Kappa;
    
    //Check
    if(nlhs!=1) {
        mexErrMsgTxt("Wrong number of output arguments");
    }
    if (nrhs!=3) {
        mexErrMsgTxt("Wrong number of input arguments.");
    }    
    if(mxIsComplex(prhs[0])) {
            mexErrMsgTxt("X should not be a complex.");
    }                       
    
    //Get System Param
    NumberOfElements=mxGetNumberOfElements(prhs[0]);    
    MaxKappa=mxGetScalar(prhs[2]);
             
    //Allocate        
    plhs[0]=mxCreateDoubleMatrix(NumberOfElements, 1, mxREAL);
    plhs[0]=mxCreateNumericArray(mxGetNumberOfDimensions(prhs[0]), mxGetDimensions(prhs[0]), mxDOUBLE_CLASS,  mxREAL);
     
    //GetPtr 
    Lambda=mxGetPr(prhs[0]);
    Dimensions=(int) mxGetScalar(prhs[1]);
    Kappa=mxGetPr(plhs[0]);     
    
    //mexPrintf("___ %u\n", NumberOfElements);
    for(nb=0; nb<NumberOfElements; nb++) {        
        Kappa[nb]=hypergeometricRatioInverse(Lambda[nb], Dimensions);
        if (Kappa[nb]>MaxKappa) {
           Kappa[nb]=MaxKappa;
        }
    }
}
