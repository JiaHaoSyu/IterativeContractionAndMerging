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

// [CountNew] = AccumulateRangeofTexture(CountOld,Label,R,Height,Width,MaxLabel,Dimension);

void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{
    int R,Height,Width,MaxLabel,Dimension;
	double *CountNew;
	double *CountOld,*Label;
    if (nlhs>1){
		mexErrMsgTxt("one return values");}
    if (nrhs!=7){
		mexErrMsgTxt("Require seven inputs");}

	Height = mxGetScalar(prhs[3]);
	Width = mxGetScalar(prhs[4]);
    R = mxGetScalar(prhs[2]);
	MaxLabel = mxGetScalar(prhs[5]);
    Dimension = mxGetScalar(prhs[6]);
	
	// output array
	plhs[0] = mxCreateDoubleMatrix(MaxLabel,Dimension,mxREAL);
    CountNew = mxGetPr(plhs[0]);
    
    for(int j=0;j<Dimension;j++)
        for(int i=0;i<MaxLabel;i++) 
            CountNew[j*(MaxLabel)+i] = 0;
    
	// input array
	CountOld = mxGetPr(prhs[0]);
    int Leftrows = mxGetM(prhs[0]);
    int Leftcols = mxGetN(prhs[0]);
    int maxL = Leftrows*Leftcols;
    
	Label = mxGetPr(prhs[1]);
	   
	int indj,indi;
    for(int i=0;i<maxL;i++){
//     for(int i=0;i<Height*Width*R*R;i++){
        indi = Label[i]-1;
        indj = CountOld[i]-1;
        CountNew[indj*(MaxLabel)+indi] = CountNew[indj*(MaxLabel)+indi] + 1;
    }
    	  
  return;  
  
}

