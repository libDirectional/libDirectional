
/* 
 * This is an adaptation of Glover's method for calculation of the Bingham 
 * normalization constant to a MATLAB mex file. The same method was previously
 * used in libbingham to create a precomputed lookup table. For some reason, it
 * does not work if Z contains zero entries except the last entry.
 * 
 * See 
 * https://code.google.com/p/bingham/
 * and
 * 
 * Glover, J. & Kaelbling, L. P. 
 * Tracking 3-D Rotations with the Quaternion Bingham Filter 
 * MIT, 2013
 * http://dspace.mit.edu/handle/1721.1/78248
 *
 *
 * Parameters:
 *      Z (d x 1 vector)
 *          concentration matrix with ascending entries (for d=2,3,4)
 * Returns:
 *      normalization constant of Bingham distribution
 */


#include <mex.h> 
#include <matrix.h> 

#define _USE_MATH_DEFINES
#include <math.h> 

#include "Gamma.h"

#define EPSILON 1e-8
#define ITERATION_MULT 10
#define MIN_ITERATIONS 10

#define MAXFACT 10000
#define MAX(x,y) ((x) > (y) ? (x) : (y))

// computes the log factorial of x
double lfact(int x)
{
	static double logf[MAXFACT];
	static int first = 1;
	int i;

	if (first) {
		first = 0;
		logf[0] = 0;
		for (i = 1; i < MAXFACT; i++)
			logf[i] = log((double)i) + logf[i-1];
	}

	return logf[x];
}

// computes the surface area of a unit sphere with dimension d
double surface_area_sphere(int d)
{
	switch(d) {
	case 0:
		return 2;
	case 1:
		return 2*M_PI;
	case 2:
		return 4*M_PI;
	case 3:
		return 2*M_PI*M_PI;
	}

	return (2*M_PI/((double)d-1))*surface_area_sphere(d-2);
}

/*
* Computes the hypergeometric function 1F1(a;b;z) in canonical form (z > 0)
*/
static double compute_1F1_1d_canon(double a, double b, double z, int iter)
{
	//printf("compute_1F1_1d_canon(%f, %f, %f, %d)\n", a, b, z, iter);

	int i;
	double g, F = 0.0, logz = log(z);

	for (i = 0; i < iter; i++) {
		g = LogGamma(i+a) - LogGamma(i+b) + i*logz - lfact(i);
		if (i > z && exp(g) < EPSILON * F)  // exp(g) < 1e-8 * F
			break;
		F += exp(g);
	}

	return 2*sqrt(M_PI)*F;
}

/*
* Computes the hypergeometric function 1F1(a;b;z) with a = 1/2, b = (dim+1)/2
*/
static double compute_1F1_1d(int dim, double z, int iter)
{
	if (fabs(z) < EPSILON)  // uniform
		return surface_area_sphere(dim);

	if (z < 0)
		return exp(z)*compute_1F1_1d(dim, -z, iter);

	return compute_1F1_1d_canon(.5, .5*(dim+1), z, iter);
}

/*
* Computes the hypergeometric function 1F1(a;b;z1,z2) in canonical form (z1 > z2 > 0)
*/
static double compute_1F1_2d_canon(double a, double b, double z1, double z2, int iter)
{
	int i, j;
	double g, F = 0.0, logz1 = log(z1), logz2 = log(z2);

	for (i = 0; i < iter; i++) {
		for (j = 0; j < iter; j++) {
			g = LogGamma(i+a) + LogGamma(j+a) - LogGamma(i+j+b) + i*logz1 + j*logz2 - lfact(i) - lfact(j);
			if ((i > z1 || j > z2) && exp(g) < EPSILON * F)  // exp(g) < 2e-9
				break;
			F += exp(g);
		}
	}

	return 2*sqrt(M_PI)*F;
}

/*
* Computes the hypergeometric function 1F1(a;b;z1,z2) with a = 1/2, b = (dim+1)/2
*/
static double compute_1F1_2d(int dim, double z1, double z2, int iter)
{
	if (fabs(z1) < EPSILON)  // uniform
		return surface_area_sphere(dim);

	if (fabs(z2) < EPSILON)  // z2 = 0
		return sqrt(M_PI)*compute_1F1_1d(dim, z1, iter);

	if (z1 < 0)
		return exp(z1)*compute_1F1_2d(dim, -z1, z2-z1, iter);

	return compute_1F1_2d_canon(.5, .5*(dim+1), z1, z2, iter);
}


/*
* Computes the hypergeometric function 1F1(a;b;z1,z2,z3) in canonical form (z1 > z2 > z3 > 0)
*/
static double compute_1F1_3d_canon(double a, double b, double z1, double z2, double z3, int iter)
{
	int i, j, k;
	double g, F = 0.0, logz1 = log(z1), logz2 = log(z2), logz3 = log(z3);

	for (i = 0; i < iter; i++) {
		for (j = 0; j < iter; j++) {
			for (k = 0; k < iter; k++) {
				g = LogGamma(i+a) + LogGamma(j+a) + LogGamma(k+a) - LogGamma(i+j+k+b) + i*logz1 + j*logz2 + k*logz3 - lfact(i) - lfact(j) - lfact(k);
				if ((i > z1 || j > z2 || k > z3) && exp(g) < EPSILON * F)  // exp(g) < 2e-9
					break;
				F += exp(g);
			}
		}
	}

	return 2*sqrt(M_PI)*F;
}

/*
* Computes the hypergeometric function 1F1(a;b;z1,z2,z3) with a = 1/2, b = (dim+1)/2
*/
static double compute_1F1_3d(int dim, double z1, double z2, double z3, int iter)
{
	if (fabs(z1) < EPSILON)  // uniform
		return surface_area_sphere(dim);

	if (fabs(z3) < EPSILON)  // z3 = 0
		return sqrt(M_PI)*compute_1F1_2d(dim, z1, z2, iter);

	if (z1 < 0)
		return exp(z1)*compute_1F1_3d(dim, -z1, z3-z1, z2-z1, iter);

	return compute_1F1_3d_canon(.5, .5*(dim+1), z1, z2, z3, iter);
}

inline double bingham_F_1d(double z)
{
	int iter = MAX((int)fabs(z)*ITERATION_MULT, MIN_ITERATIONS);
	return compute_1F1_1d(1, z, iter);
}

inline double bingham_F_2d(double z1, double z2)
{
	int iter = MAX((int)MAX(fabs(z1), fabs(z2))*ITERATION_MULT, MIN_ITERATIONS);
	return compute_1F1_2d(2, z1, z2, iter);
}

inline double bingham_F_3d(double z1, double z2, double z3)
{
	int iter = MAX((int)MAX(MAX(fabs(z1), fabs(z2)), fabs(z3))*ITERATION_MULT, MIN_ITERATIONS);
	return compute_1F1_3d(3, z1, z2, z3, iter);
}


void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double *la; // Eigenvalues of the matrix.
	size_t m,n;
	double *result;
	int i;


	/* Check argument count */
	if (nrhs != 1 ) { 
		mexErrMsgTxt("Requires 1 input argument."); 
	} else if (nlhs > 1) { 
		mexErrMsgTxt("Returns just one value."); 
	} 

	/* Get arguments from MATLAB */ 
	la = mxGetPr(prhs[0]); // vector containing la.
	m = mxGetM(prhs[0]); //rows of la
	n = mxGetN(prhs[0]); //columns of la

	/* Check argument dimensionality */
	if(n != 1 || m<2) {
		mexErrMsgTxt("Column vector expected.");
	}      

	for(i=1; i<(int)m; i++) {
		if(la[i]<la[i-1]) mexErrMsgTxt("Entries. must be given in increasing order");
	}
	if (la[m-1] != 0){
		mexErrMsgTxt("Last entry must be zero.");
	}

	plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
	result = mxGetPr(plhs[0]);

	switch(m){
		case 2:
			result[0] = bingham_F_1d(la[0]); 
			break;
		case 3:
			result[0] = bingham_F_2d(la[0],la[1]); 
			break;
		case 4:
			result[0] = bingham_F_3d(la[0],la[1],la[2]); 
			break;
		default:
			mexErrMsgTxt("Only 2D, 3D, and 4D case supported.");        
	}
}
