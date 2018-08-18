#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <mex.h> 
#include <matrix.h> 
#include <iostream>
#include "../external/Eigen/Dense"

#define _USE_MATH_DEFINES
#include <math.h> 

#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

#include "binghamNormalizationConstant.h"


using namespace Eigen;

/*
 * Igor Gilitschenski, Gerhard Kurz, Simon J. Julier, Uwe D. Hanebeck,
 * Efficient Bingham Filtering based on Saddlepoint Approximations
 * Proceedings of the 2014 IEEE International Conference on Multisensor Fusion 
 * and Information Integration (MFI 2014), Beijing, China, September 2014.
*/


int binghamMLE(int dim, double *in, double *res, double *initval);

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    double *om; // Diagonal entries of Covariance
    double *res; // Result of computation.
    size_t m,n;
    
    /* Check argument count */
    if (nrhs != 1 ) { 
		mexErrMsgTxt("requires 1 input argument."); 
	} else if (nlhs > 1) { 
		mexErrMsgTxt("returns just one value."); 
	}
    
    /* Get arguments from MATLAB */
	om = mxGetPr(prhs[0]); // Vector containing om.
    m = mxGetM(prhs[0]); //rows of la
    n = mxGetN(prhs[0]); //columns of la
    
	/* Check argument dimensionality */
    if(n != 1 || m==0) {
        mexErrMsgTxt("column vector expected.");
    }

    plhs[0] = mxCreateDoubleMatrix(m,1,mxREAL);
    res = mxGetPr(plhs[0]);
    
    binghamMLE((int)m, om, res, 0);
}



/**
 * binghamMLE
 * Implements a maximum likelihood estimator of the bingham distribution
 * based on Newton's Method. We try to find z_i such that
 *      N(Z)^-1 * (d N(Z)/d res[i]) = in[i] ,
 * where Z=diag(res[0], ..., res[dim-1]). This implementation assumes a 
 * correct length / dimensionality of parameters and does not perform any 
 * checks.
 *
 * Parameters:
 * dim      - dimension of considered distribution.
 * in       - values to be fitted (basically eigenvalues of covariance)
 * res      - memory for storing the result.
 * initval  - initial value for optimization.
 */
int binghamMLE(int dim, double *in, double *res, double *initval) {
    
    int i,j;
    
    double *ncApprox = new double[3]; // SP approx of Norm Const.
    double *derivApprox = new double [3*dim]; // SP approx. of derrivatives
    
    // Norm of objective function (shall be minimized)
    double oldNorm;
    double normObjFun=0.0;
    
    // Temporary variables for derivative computation.
    double *tmpEig=new double[dim+2]; 
    double *tmpNc=new double[3]; 
    double *tmpDeriv=new double[3*(dim+2)]; 
    
    VectorXd ncDeriv(dim); // Partial derivatives
    MatrixXd ncDeriv2(dim, dim); // Second order partial derivatives
    VectorXd objFun(dim); // Vector holding value of function which is searched for roots
    MatrixXd objFunJacobian(dim,dim);
    
    // Make handling of predefined arrays as Eigen object possible.
    Map<VectorXd,0,InnerStride<3> > eTmpDeriv(tmpDeriv+2,dim+2); 
    Map<VectorXd> eInput(in,dim);
    Map<VectorXd> eRes(res, dim);
    
    // Copy initialization values.
    if(initval == 0) {
        memset(res, 0, dim*sizeof(double));
    
        // Heuristic choice of starting values for the optimization.
        for(i=0; i<dim; i++)
            res[i] = 1/(2*in[i]);
        
    } else {
        memcpy(res, initval, dim*sizeof(double));
        eRes = -eRes;
    }
    
    // Main optimization loop.
    for(i=0; i<=1000; i++)
    {
        // Set maximum entry to be zero.
        eRes.array()=eRes.array()+(-eRes.minCoeff());
        
        //std::cout << eRes.array()+(-eRes.minCoeff()) << "\n";
        
        // Compute the normalization constant.
        binghamNormalizationConstant(dim, res, ncApprox, 0);

        // Compute its gradient and Hessian
        for(j=0; j<dim; j++) {
            memcpy(tmpEig, res, j*sizeof(double));
            tmpEig[j]=res[j];
            tmpEig[j+1]=res[j];
            tmpEig[j+2]=res[j];
            memcpy(tmpEig+j+3, res+j+1, (dim-j-1)*sizeof(double*) );


            binghamNormalizationConstant(dim+2, tmpEig, tmpNc, tmpDeriv);

            // Gradient of normalization constant.
            ncDeriv(j) = -tmpNc[2]/(2*M_PI);

            // Hessian of normalization constant.
            ncDeriv2.row(j).head(j) = -eTmpDeriv.head(j)/(2*M_PI);
            ncDeriv2(j,j) = - 3*tmpDeriv[3*j+2]/(2*M_PI);
            ncDeriv2.row(j).tail(dim-j-1) = -eTmpDeriv.tail(dim-j-1)/(2*M_PI);
        }

        // Compute objective function and its jacobian.
        objFun = ncDeriv / ncApprox[2] - eInput;
        objFunJacobian = 
                (ncDeriv2*ncApprox[2] - ncDeriv*ncDeriv.transpose())
                / (ncApprox[2]*ncApprox[2]);


        oldNorm = normObjFun;
        normObjFun = objFun.norm();

        // Gaus-Newton Step
        eRes = eRes - 
                (objFunJacobian.transpose() * objFunJacobian).inverse() 
                * objFunJacobian.transpose() * objFun;
        
        if(fabs(oldNorm-normObjFun)<=1e-10) {
            std::cout << "Stopped after "<< i+1 << " iterations.\n";
            break;
        }
    }

    
    eRes.array()=-(eRes.array()+(-eRes.minCoeff()));
    
    delete[] ncApprox;
    delete[] derivApprox;
    delete[] tmpEig;
    delete[] tmpNc;
    delete[] tmpDeriv;
    
    return 0;
}