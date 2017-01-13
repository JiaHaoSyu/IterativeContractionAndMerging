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

// [Output] = SparseVector(tempInd,MappingV,tempHeight,tempWidth);

void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{
    int tempHeight,tempWidth;
	double *Output;
	double *tempInd,*MappingV;
    if (nlhs>1){
		mexErrMsgTxt("one return values");}
    if (nrhs!=4){
		mexErrMsgTxt("Require Four inputs");}

	tempHeight = mxGetScalar(prhs[2]);
	tempWidth = mxGetScalar(prhs[3]);
			
	// output array
	plhs[0] = mxCreateDoubleMatrix(tempHeight,tempWidth,mxREAL);
	
	Output = mxGetPr(plhs[0]);
	// input array
	tempInd = mxGetPr(prhs[0]);
	MappingV = mxGetPr(prhs[1]);
    int tempI,tempC;
    for(int i=0;i<tempWidth;i++){
        for(int j=0;j<tempHeight;j++){
            tempI = tempInd[i*tempHeight + j]-1;
            tempC = MappingV[tempI];
            Output[i*tempHeight + j] = tempC;
        }
    }
  
  return;  
  
}

