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

// [CountNew] = CountNeiMerge(CountOld,LabelTC,Height,Width,LabelNum);

void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{
    int Height,Width,LabelNum;
	double *CountNew;
	double *CountOld,*LabelTC;
    if (nlhs>1){
		mexErrMsgTxt("one return values");}
    if (nrhs!=5){
		mexErrMsgTxt("Require five inputs");}

	Height = mxGetScalar(prhs[2]);
	Width = mxGetScalar(prhs[3]);
	LabelNum = mxGetScalar(prhs[4]);    
		
	// output array
	plhs[0] = mxCreateDoubleMatrix(LabelNum,LabelNum,mxREAL);
	
	CountNew = mxGetPr(plhs[0]);
	// input array
	CountOld = mxGetPr(prhs[0]);
	LabelTC = mxGetPr(prhs[1]);
	   
	int indi,indj;
	for(int j=0;j<Width;j++){
        indj = LabelTC[j]-1;
		for(int i=0;i<Height;i++){ 
			indi = LabelTC[i]-1;
            if(indj==indi){
                CountNew[indj*(LabelNum)+indi] = 0;
            }
            else{
                CountNew[indj*(LabelNum)+indi] = CountNew[indj*(LabelNum)+indi] + CountOld[j*Height+i];
            }
		}
	}
  
  return;  
  
}

