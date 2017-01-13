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

// [IMatrix] = IndexToMatrix(row_indsN,col_indsN,Dist1,length,MaxL);

void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{
    int MaxL,length;
	double *IMatrix;
	double *row_indsN,*col_indsN,*Dist1;

	length = mxGetScalar(prhs[3]);    
    MaxL = mxGetScalar(prhs[4]);    
// output array
	plhs[0] = mxCreateDoubleMatrix(MaxL,MaxL,mxREAL);
	         
	IMatrix = mxGetPr(plhs[0]);
	
     for(int j=0;j<MaxL;j++){
		for(int i=0;i<MaxL;i++){ 
			IMatrix[i+MaxL*j] = 0;
		}
	}
    // input array
	row_indsN = mxGetPr(prhs[0]);
	col_indsN = mxGetPr(prhs[1]);
	Dist1 = mxGetPr(prhs[2]);
    int indi,indj;
	for (int i=0;i<length;i++)
    {
        indi = (int)row_indsN[i]-1;
        indj = (int)col_indsN[i]-1;
        IMatrix[indi+MaxL*indj] = IMatrix[indi+MaxL*indj] + Dist1[i];
    }
   
  return;  
  
}

