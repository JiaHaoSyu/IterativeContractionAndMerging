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

// [CountNew] = MeanToLabel(Label,MeanInput,Height,Width,MaxLabel);

void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{
    int Height,Width,MaxLabel;
	double *CountNew;
	double  *Label,*MeanInput;
    if (nlhs>1){
		mexErrMsgTxt("one return values");}
    if (nrhs!=5){
		mexErrMsgTxt("Require five inputs");}

	Height = mxGetScalar(prhs[2]);
	Width = mxGetScalar(prhs[3]);
   	MaxLabel = mxGetScalar(prhs[4]);
   	//LabelNum = mxGetScalar(prhs[4]);    
	//mexPrintf("Height = %d\n",Height);
	//mexPrintf("Width = %d\n",Width);
		
	// output array
	plhs[0] = mxCreateDoubleMatrix(Height,Width,mxREAL);
	CountNew = mxGetPr(plhs[0]);
    
    for(int i=0;i<Height;i++)
        for(int j=0;j<Width;j++)
            CountNew[j*Height+i] = 0;
    
	// input array
	Label = mxGetPr(prhs[0]);
	MeanInput = mxGetPr(prhs[1]);
	   
	int indj,indi;
    for(int i=0;i<Height*Width;i++){
        indi = Label[i]-1;
        CountNew[i] = MeanInput[indi];
    }
    	  
  return;  
  
}

