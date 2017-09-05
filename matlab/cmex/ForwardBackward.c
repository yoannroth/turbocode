#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#include "mex.h"
#include <stdio.h>

void mexFunction(int nlhs, mxArray *plhs[], 
             int nrhs, const mxArray *prhs[]) 
{
/* variable declarations here */
	double  *B,                         /* input matrix*/
            *alpha, *beta;              /* output*/
    int *nextState, *prevState;         /* input matrix*/
	mwSize N, nStates, i, j, mlmap, infty;
    double alpha0, alpha1, beta0, beta1, norm;

/* Check number of inputs and outputs: */
	if (nrhs != 6) {
       mexErrMsgIdAndTxt("ForwardBackward:BadNInput",
                         "6 inputs required.");
    }
    if(nlhs != 2) {
        mexErrMsgIdAndTxt("ForwardBackward:BadNOutput",
                          "2 outputs required.");
    }
    
 if( !mxIsInt32(prhs[1]) ) {
    mexErrMsgTxt("Input argument is not int32 class.");
}

/* Read inputs: */
    B           = mxGetPr(prhs[0]);
    nextState   = (int*)mxGetData(prhs[1]);
    prevState	= (int*)mxGetData(prhs[2]);
    nStates     = (mwSize) mxGetScalar(prhs[3]);
    N           = (mwSize) mxGetScalar(prhs[4]);
    mlmap       = (mwSize) mxGetScalar(prhs[5]);
    
/* Create output: */
    plhs[0]     = mxCreateDoubleMatrix(nStates, N+1, mxREAL);
    alpha       = mxGetPr(plhs[0]);
    plhs[1]     = mxCreateDoubleMatrix(nStates, N+1, mxREAL);
    beta        = mxGetPr(plhs[1]);
    
/* Computations: */
    infty = (1<<10);
    
    /* MAP */
    if(!mlmap){
        /* Forward */
        /* Init */
        alpha[0] = 1;
        for (j = 1; j < nStates; j++) {
            alpha[j] = 0;
        }
        /* Compute probabilities */
        for (i = 1; i < N+1; i++) {     /* For each step*/ 
            norm = 0;
            for (j=0; j < nStates; j++) {     /* For each state*/
                alpha[j + nStates*i] = alpha[prevState[j] + nStates*(i-1)] * B[prevState[j] + (2*nStates)*(i-1)] 
                                  + alpha[prevState[j + nStates] + nStates*(i-1)] * B[prevState[j + nStates] + nStates + (2*nStates)*(i-1)];
                norm += alpha[j + nStates*i];   
            }
            /* normalisation */
            for (j=0; j < nStates; j++) {
                alpha[j + nStates*i] = alpha[j + nStates*i]/norm;
            }
        }
        
        /* Backward */
        /* Init */
        beta[nStates*N] = 1;
        for (j = 1; j < nStates; j++) {
            beta[nStates*N + j] = 0;
        }
        
        /* Compute probabilities (log) */
        for (i = N-1; i > -1; i--) {     /* For each step*/ 
            norm = 0;
            for (j=0; j < nStates; j++) {     /* For each state*/
                beta[j + nStates*i] = beta[nextState[j] + nStates*(i+1)] * B[j + (2*nStates)*(i)]
                              + beta[nextState[j + nStates] + nStates*(i+1)] * B[j + nStates + (2*nStates)*(i)];
                                
                norm += beta[j + nStates*i]; 
            }
            
            /* normalisation */
            for (j=0; j < nStates; j++) {
                beta[j + nStates*i] = beta[j + nStates*i]/norm;
            }
            
        }
        
    } else{
        /* max log-MAP */
        
        /* Forward */
        /* Init */
        alpha[0] = 0;
        for (j = 1; j < nStates; j++) {
            alpha[j] = -infty;
        }
        /* Compute probabilities */
        for (i = 1; i < N+1; i++) {     /* For each step*/ 
            norm = -infty;
            for (j=0; j < nStates; j++) {     /* For each state*/
                alpha[j + nStates*i] = MAX(alpha[prevState[j] + nStates*(i-1)] + B[prevState[j] + (2*nStates)*(i-1)], 
                                  alpha[prevState[j + nStates] + nStates*(i-1)] + B[prevState[j + nStates] + nStates + (2*nStates)*(i-1)]);
                if(alpha[j + nStates*i] > norm){
                    norm = alpha[j + nStates*i];      /* max search */
                } 
            }
            /* normalisation */
            for (j=0; j < nStates; j++) {
                alpha[j + nStates*i] = alpha[j + nStates*i] - norm;
            }
        }
        
        /* Backward */
        /* Init */
        beta[nStates*N] = 0;
        for (j = 1; j < nStates; j++) {
            beta[nStates*N + j] = -infty;
        }
        
        /* Compute probabilities (log) */
        for (i = N-1; i > -1; i--) {     /* For each step*/ 
            norm = -infty;
            for (j=0; j < nStates; j++) {     /* For each state*/
                beta[j + nStates*i] = MAX(beta[nextState[j] + nStates*(i+1)] + B[j + (2*nStates)*(i)],
                               beta[nextState[j + nStates] + nStates*(i+1)] + B[j + nStates + (2*nStates)*(i)]);
               
                if(beta[j + nStates*i] > norm){
                    norm = beta[j + nStates*i];      /* max search */
                } 
            }
            
            /* normalisation */
            for (j=0; j < nStates; j++) {
                beta[j + nStates*i] = beta[j + nStates*i] - norm;
            }
            
        }
        
    }
    
    return;
}

/* Debug :
  
 printf("b[%d,%d] = b[%d,%d]*gamma[%d,%d] + b[%d,%d]*gamma[%d,%d] = %f\n",
                     j, i, nextState[j], (i+1), j, (i), nextState[j+nStates], (i+1), j+nStates, (i), beta[j + nStates*i]);
 printf("%f = %f*%f + %f*%f\n",
                  beta[j + nStates*i], beta[nextState[j] + nStates*(i+1)], B[j + (2*nStates)*(i)], beta[nextState[j + nStates] + nStates*(i+1)], B[j + nStates + (2*nStates)*(i)]);
                */
