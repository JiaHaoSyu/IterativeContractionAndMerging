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

// [CountNew] = LabelSum(Label,Input,Height,Width,MaxLabel);

void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{
    int Height,Width,MaxLabel;
	double *CountNew;
	double  *Label,*Input;

	Height = mxGetScalar(prhs[2]);
	Width = mxGetScalar(prhs[3]);
   	MaxLabel = mxGetScalar(prhs[4]);
		
	// output array
	plhs[0] = mxCreateDoubleMatrix(MaxLabel,1,mxREAL);
	
	CountNew = mxGetPr(plhs[0]);
    for(int j=0;j<MaxLabel;j++)
         CountNew[j] = 0;
    
	// input array
	Label = mxGetPr(prhs[0]);
	Input = mxGetPr(prhs[1]);
	   
	int indj,indi;
    for(int i=0;i<Height*Width;i++){
        indi = Label[i]-1;
        CountNew[indi] = CountNew[indi] + Input[i];
    }
    	  
  return;  
  
}

