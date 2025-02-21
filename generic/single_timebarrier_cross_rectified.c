#include "mex.h"
/* single_timebarrier_cross_rectified(ev,Bound,lower_bound)
 Inputs: a matrix [Time x Trials], time-varying bound [Time], and the
 *minimum value that the cumsum can take.
 * Computes the cumsum sum along 
 * the rows of the matrix. When the barrier is reached, the rest of the column
 * is set to zero. Only a single positive bound is permitted. The second output is the
 * sample (linear interpolation) at which the bound was reached; it's set to Nan when the
 * barrier was not reached.
 
 */
void
        mexFunction(int nlhs, mxArray* plhs[], int nrhs,
        const mxArray* prhs[]) {
    mwIndex i, j, k;
    mwSize n;
    double *vri, *vro, *out,*bound,*reflective_bound_input;
    double element,Sprev,Scurr;
    int ini, cont, bReached;
    double reflective_bound;
    
    const mwSize *dim_array;
    mwSize mrows, ncols;
    
    if (nrhs < 3 || ! mxIsNumeric(prhs[0]))
        mexErrMsgTxt("wrong inputs");
    
    n = mxGetNumberOfElements(prhs[0]);
    
    /* Get the number of dimensions in the input argument. */
    dim_array = mxGetDimensions(prhs[0]);
    
    plhs[0] = (mxArray *) mxCreateNumericArray
            (mxGetNumberOfDimensions(prhs[0]),
            dim_array, mxGetClassID(prhs[0]),
            mxIsComplex(prhs[0]));
    
    
    vri = mxGetPr(prhs[0]);
    bound = mxGetPr(prhs[1]);
    
    reflective_bound_input = mxGetPr(prhs[2]);
    reflective_bound = reflective_bound_input[0];
    
    
    vro = mxGetPr(plhs[0]);
    
    mrows = mxGetM(prhs[0]);
    ncols = mxGetN(prhs[0]);
    
    plhs[1] = mxCreateDoubleMatrix(ncols, 1, mxREAL);
    out = mxGetPr(plhs[1]); 
    
    cont = 0;
    for (j = 0; j < ncols; j++) {/*trials*/
        bReached = 0;
        out[j] = mxGetNaN();/*set out of trial to nan*/
        for (i = 0; i < mrows; i++) {/*time*/
            if (bReached == 0){
                element = vri [mrows*j + i];
                if (mxIsNaN(element))
                    element = 0;
                if (i==0) 
                    vro [cont] = element;
                else
                    vro[cont] = vro[cont-1] + element;
                    
                if (vro[cont]<reflective_bound)
                    vro[cont] = reflective_bound;
                
                if (vro[cont]>=bound[i]){
                    bReached = 1;
                    /* if bounf reached, interpolate decision-time */
                    Sprev = vro[cont] - element;/*val before last sample*/
                    Scurr = vro[cont];
                    
                    /*interpolate linearly:*/
                    if (i>0)
                        out[j] = i + (bound[i-1]-Sprev)/(bound[i-1]-Sprev+Scurr-bound[i]);
                    else /*crossed on first sample*/
                        out[j] = i + (bound[i]-Sprev)/(bound[i]-Sprev+Scurr-bound[i]);
                }
            }
            else {
                vro[cont] = mxGetNaN();
            }
            cont++;
        }
        
    }
}

