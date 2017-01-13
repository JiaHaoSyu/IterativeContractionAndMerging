#include "mex.h" /* Always include this */
#include <math.h>
#include <stdlib.h>

#ifndef max
#define max(a,b) (((a) > (b)) ? (a) : (b))
#endif

#ifndef min
#define min(a,b) (((a) < (b)) ? (a) : (b))
#endif


using namespace std;

/* nlhs: number of output variables 
   plhs : Array of mxArray pointers to the output variables
   nrhs : number of input variables
   prhs : Array of mxArray pointers to the input variables
*/

// [Output] = SparseMatrix(RowInd,ColInd,Value,len,tempHeight,tempWidth);

void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{
    int tempHeight,tempWidth,len;
	double *Output;
	double *RowInd,*ColInd,*Value;

    len =  mxGetScalar(prhs[3]);
	tempHeight = mxGetScalar(prhs[4]);
	tempWidth = mxGetScalar(prhs[5]);
			
	// output array
	plhs[0] = mxCreateDoubleMatrix(tempHeight,tempWidth,mxREAL);
	Output = mxGetPr(plhs[0]);
    
    for(int i=0;i<tempHeight;i++)
       for(int j=0;j<tempWidth;j++)
            Output[j*tempHeight+ i] = 0; 
	// input array
	RowInd = mxGetPr(prhs[0]);
	ColInd = mxGetPr(prhs[1]);
    Value = mxGetPr(prhs[2]);
    int tempRow,tempCol;
    double tempVal;
    
    for(int i=0;i<len;i++){
            tempRow = RowInd[i]-1;
            tempCol = ColInd[i]-1;
            tempVal = Value[i];
            Output[tempCol*tempHeight + tempRow] = Output[tempCol*tempHeight + tempRow] + tempVal;
    }
  
  return;  
  
}

