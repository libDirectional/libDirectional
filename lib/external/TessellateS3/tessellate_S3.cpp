#include <stdio.h>
#include "hypersphere.h"
#include "util.h"
#include "mex.h"

/*
int main(int argc, char *argv[])
{
  if (argc < 2) {
    printf("usage: %s <num points>\n", argv[0]);
    return 1;
  }

  hypersphere_tessellation_t *T = tessellate_S3(atoi(argv[1]));

  int i,j;

  printf("%d %d\n", T->n, T->d + 1);
  for (i = 0; i < T->n; i++) {
    for (j = 0; j < T->d; j++)
      printf("%f ", T->centroids[i][j]);
    printf("%f\n", T->volumes[i]);
  }

  return 0;
}*/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, 
  const mxArray *prhs[]){
    if(nrhs != 1){
        mexErrMsgIdAndTxt("tessellate_S3:nrhs",
                      "One input required.");
    }
    if( !mxIsDouble(prhs[0]) || 
        mxIsComplex(prhs[0]) ||
        mxGetNumberOfElements(prhs[0]) != 1 ) {
        mexErrMsgIdAndTxt("tessellate_S3:notScalar",
                      "Input must be a scalar.");
    }
    if(nlhs != 1) {
        mexErrMsgIdAndTxt("tessellate_S3:nlhs",
                          "One output required.");
    }
    hypersphere_tessellation_t *T = tessellate_S3((int)(*mxGetPr(prhs[0])));
    
    plhs[0] = mxCreateDoubleMatrix(T->d,T->n,mxREAL);
    double* p = mxGetPr(plhs[0]);
    for (int i = 0; i < T->n; i++) {
        for (int j = 0; j < T->d; j++)
            p[i*T->d + j] = T->centroids[i][j];
    }
}
