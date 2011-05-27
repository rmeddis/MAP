/* intpow - Nick Clark - 3/5/2011
 *
 * Function that raises a double precision vector to an integer power 
 * vector. The following call to this mex function... 
 *
 * >>      z = intpow(x,y)
 * 
 * ... does the equivalent of the following in Matlab...
 *
 * >>      z = x .^ ceil(y)
 * 
 * NOTE: Under most circumstances, this function is slower than MATLAB's 
 * built in ^ operator, but for the (generally) small integer powers used 
 * in MAP we observe massive performance boosts using this C function.
 */

/*************************************************************************/
/* Header(s)                                                             */
/*************************************************************************/
#include "mex.h"

/*************************************************************************/
/* Input vars                                                            */
/*************************************************************************/
#define IN_x 	prhs[0]
#define IN_y 	prhs[1]

/*************************************************************************/
/* Output vars                                                           */
/*************************************************************************/
#define OUT_z 	plhs[0]

/*************************************************************************/
/* Gateway function and error checking                                   */
/*************************************************************************/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* variable declarations */
    int numElements, xM, xN, yM, yN; 
    int nn, kk; 
    double *x, *y, *z;
     
    /*  check the number of input and output parameters  */  
    if(nrhs!=2)
        mexErrMsgTxt("intpow : Two input args expected");
    if(nlhs > 1)
        mexErrMsgTxt("intpow : Too many outputs");
    
    /* check that x and y have equal size */ 
    x = mxGetPr(IN_x);
    xM = mxGetM(IN_x);
    xN = mxGetN(IN_x);
    
    y = mxGetPr(IN_y);
    yM = mxGetM(IN_y);
    yN = mxGetN(IN_y);
    
    if  (xM != yM || xN != yN)
        mexErrMsgTxt("intpow : x and y must have equal size");
    
    /* find upper loop boundary */
    numElements = xM * xN; 
    
    /*  allocate memory and pointer for the output array */ 
    OUT_z = mxCreateDoubleMatrix(xM,xN,mxREAL);
    z = mxGetPr(OUT_z);
    
    /*  do stuff */ 
    for( nn = 0;  nn<numElements; ++nn )
    {
        z[nn] = 1.0;
        for (kk = 0; kk<y[nn]; ++kk)
        {
            z[nn] = z[nn]*x[nn];
        }
    }            
}
