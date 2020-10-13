//#define EIGEN_NO_DEBUG
//#define NDEBUG 
//#define EIGEN_UNROLLING_LIMIT 1000

#define _USE_MATH_DEFINES

#include <mex.h> 

#include "Mex/Matrix.h"

#include <cmath>
#include <iostream>

#include "../external/fmath.hpp"

using namespace Mex;

/*
 * Efficient implementation of the probability density function of a multivariate Gaussian.
 *
 * Call with mvnpdf(x, mu, C).
 * Parameters:
 *  x (n x d matrix)
 *      locations where the pdf should be evaluated (n points with d dimensions)
 *  mu (1 x d vector) 
 *      mean vector
 *  C (d x d matrix)
 *      covariance matrix, must be symmetric positive definite
 *
 * This function can serve as a drop-in replacement for MATLAB's mvnpdf, if the parameters are given as above. 
 * It does not support all allowed parameter combinations found in MATLAB's mvnpdf. 
 * For efficiency reasons, we do not verify that C is really s.p.d. and the results will be incorrect in case it is not.
 *
 * Compilation:
 * on Windows: 
 * mex('OPTIMFLAGS=/openmp /Ox /arch:AVX','mvnpdffast.cpp')
 * (remove /arch:AVX if your CPU does not support AVX)
 *
 * on Linux: 
 * mex('CXXFLAGS=\$CXXFLAGS -std=c++0x -fopenmp -O3 -march=native','LDFLAGS=\$LDFLAGS -fopenmp','mvnpdffast.cpp')
 */


template<int N> void calculateGaussPdf(Eigen::Matrix<double,Eigen::Dynamic,N> X, Eigen::Matrix<double,1,N> mu, Eigen::Matrix<double,N,N> C, double* result){
    Eigen::Matrix<double,N,N> L = C.llt().matrixL().transpose(); // cholesky decomposition
    Eigen::Matrix<double,N,N> Linv = L.inverse();
    
    double det = L.diagonal().prod(); //determinant of L is equal to square rooot of determinant of C
	double lognormconst = -log(2 * M_PI)*X.cols()/2 - log(fabs(det));

    Eigen::Matrix<double,1,N> x = mu;
    Eigen::Matrix<double,1,N> tmp = x;
    for (int i=0; i<X.rows(); i++){
        x.noalias() = X.row(i) - mu;
        tmp.noalias() = x*Linv;
        double exponent = -0.5 * (tmp.cwiseProduct(tmp)).sum();
        result[i] = lognormconst+exponent;
    }
    
    /*
	Eigen::Matrix<double,Eigen::Dynamic,N> X0 = (X.rowwise() - mu)*Linv;
	Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,1> > resultMap(result, X.rows());
	resultMap = (X0.rowwise().squaredNorm()).array() * (-0.5) + lognormconst;
    */
    
    fmath::expd_v(result, X.rows());
}

void mexFunction (int numOutputs, mxArray *outputs[],
                 int numInputs, const mxArray *inputs[])
{   
    if(numInputs == 3 && numOutputs == 1){
        //Eigen::initParallel();

        ConstMatrix<double, Eigen::Dynamic, Eigen::Dynamic> X(inputs[0]);

        // Get size of data. 
        ConstMatrix<double, Eigen::Dynamic, Eigen::Dynamic>::Index samples = X.rows(); 
        ConstMatrix<double, Eigen::Dynamic, Eigen::Dynamic>::Index dimension = X.cols();

        if (dimension<1){
            mexErrMsgTxt("stats:mvnpdf:TooFewDimensions");
        }    
   
        ConstMatrix<double, 1, Eigen::Dynamic> mu(inputs[1]);
        ConstMatrix<double, Eigen::Dynamic, Eigen::Dynamic> C(inputs[2]);
        
        if (mu.cols()!=dimension || C.rows()!=dimension || C.cols()!=dimension ){
            mexErrMsgTxt("stats:mvnpdf:invalidparameters");
        }

        outputs[0] = mxCreateDoubleMatrix (samples, 1, mxREAL);
        double* result = mxGetPr(outputs[0]);  
        
        double lognormconst;
        switch(dimension){
            case 1:
                {
                    lognormconst = log(1.0/sqrt(2 * M_PI*C(0,0)));
                    double factor =  -0.5 /C(0,0);
                    //#pragma omp parallel for
                    for (int i=0; i<samples; i++){
                        double x = X(i,0) - mu(0,0);
                        double exponent =x*x*factor;
                        result[i] = lognormconst + exponent;
                    }
                    fmath::expd_v(result, X.rows());
                }
                break;
            case 2:
                calculateGaussPdf<2>(X, mu, C, result);
                break;
            case 3:
                calculateGaussPdf<3>(X, mu, C, result);
                break;
            case 4:
                calculateGaussPdf<4>(X, mu, C, result);
                break;
                /*
            case 5:
                calculateGaussPdf<5>(X, mu, C, result);
                break;                
            case 6:
                calculateGaussPdf<6>(X, mu, C, result);
                break;                
            case 7:
                calculateGaussPdf<7>(X, mu, C, result);
                break;       */         
            default:
                //calculateGaussPdf<Eigen::Dynamic>(X, mu, C, result);
                //break;
                
                Eigen::MatrixXd L = C.llt().matrixL().transpose(); // cholesky decomposition
                Eigen::MatrixXd Linv = L.inverse();

                double det = L.diagonal().prod(); //determinant of L is equal to square rooot of determinant of C
                lognormconst = -log(2 * M_PI)*dimension/2 - log(fabs(det));
                
                
                Eigen::MatrixXd X0 = (X.rowwise() - mu)*Linv;
                Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,1> > resultMap(result, samples);
                resultMap = (X0.rowwise().squaredNorm()).array() * (-0.5) + lognormconst;

                /*
                Eigen::Matrix<double,1,Eigen::Dynamic> x = mu;
                Eigen::Matrix<double,1,Eigen::Dynamic> tmp = x;

                for (int i=0; i< samples; i++){
                    x = X.row(i) - mu;
                    tmp.noalias() = x*Linv;
                    double exponent = -0.5 * (tmp.cwiseProduct(tmp)).sum();
                    result[i] = lognormconst+exponent;
                }*/
                
                fmath::expd_v(result, X.rows());
                break;
        }
    }
    else{
        mexErrMsgTxt("stats:mvnpdf:invalidparameters");
    }
}
