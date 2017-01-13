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

// [row_inds,col_inds] = GetNeighborInformation(ComputeTimes,Index,TotalLength,I_Height,I_Width,Para.SlidingWindows);

void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{
    int TotalLength,Height,Width,WinSize;
	double *row_inds,*col_inds,*label_inds;
	double *ComputeTimes,*Index;
   
    TotalLength = mxGetScalar(prhs[2]);
	Height = mxGetScalar(prhs[3]);
	Width = mxGetScalar(prhs[4]);
	WinSize = mxGetScalar(prhs[5]);
    
	int maxX,maxY,minX,minY;
	// output array
	plhs[0] = mxCreateDoubleMatrix(1,TotalLength,mxREAL);
	plhs[1] = mxCreateDoubleMatrix(1,TotalLength,mxREAL);
	row_inds = mxGetPr(plhs[0]);
    col_inds = mxGetPr(plhs[1]);
	
	// input array
	ComputeTimes = mxGetPr(prhs[0]);
	Index = mxGetPr(prhs[1]);
	   
	int len = 0;
	int count = 0;
	int label = 1;
	for(int j=0;j<Width;j++){ 
		for(int i=0;i<Height;i++){ 
		    minX = max(i-WinSize,0);
			maxX = min(i+WinSize,Height-1);
			minY = max(j-WinSize,0);
			maxY = min(j+WinSize,Width-1);
			for(int k=0;k<ComputeTimes[count];k++)
			row_inds[len+k] = Index[count];
			int counxy = 0;
			for(int v=minY;v<=maxY;v++){
				for(int u=minX;u<=maxX;u++){
					col_inds[len+counxy] = v*Height+u+1;
					counxy++;
				}
			}

			len = len + ComputeTimes[count];
			count++;
		}
	}
  
  return;  
  
}

