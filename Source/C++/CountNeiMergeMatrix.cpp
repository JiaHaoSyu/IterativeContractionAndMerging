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

// [CountNew] = CountNeiMergeMatrix(CountOld,LabelTC,Height,Width,row,col);

void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{
    int Height,Width,row,col;
	double *CountNew;
	double *CountOld,*LabelTC;
    if (nlhs>1){
		mexErrMsgTxt("one return values");}
    if (nrhs!=6){
		mexErrMsgTxt("Require six inputs");}

	Height = mxGetScalar(prhs[2]);
	Width = mxGetScalar(prhs[3]);
	row = mxGetScalar(prhs[4]);    
	col = mxGetScalar(prhs[5]);    	
	// output array
	plhs[0] = mxCreateDoubleMatrix(row,col,mxREAL);
	
	CountNew = mxGetPr(plhs[0]);
	// input array
	CountOld = mxGetPr(prhs[0]);
	LabelTC = mxGetPr(prhs[1]);
	   
	int indi,indj;
    for(int j=0;j<Width;j++){
//         indj = LabelTC[j]-1;
		for(int i=0;i<Height;i++){ 
			indi = LabelTC[i]-1;
            CountNew[j*(row)+indi] = CountNew[j*(row)+indi] + CountOld[j*Height+i];
		}
	}
  
  return;  
  
}

