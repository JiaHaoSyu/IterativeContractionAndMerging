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

// [row_inds col_inds Dist BoundaryCount] = GetBorderDist(lab,Label,Mask,TotalLength,I_Height,I_Width,Windows,Length);

void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{
    int TotalLength,Height,Width,WinSize,Length;
	double *row_inds,*col_inds,*Dist,*MaxMatrix,*MinMatrix,*BoundaryCount;
	double *lab,*Label,*Mask,*ComputeTimes;
    
	//ComputeTimes = mxGetN(prhs[0]);
    TotalLength = mxGetScalar(prhs[3]);
	Height = mxGetScalar(prhs[4]);
	Width = mxGetScalar(prhs[5]);
	WinSize = mxGetScalar(prhs[6]);
    Length = mxGetScalar(prhs[7]);
	//mexPrintf("TotalLength = %d\n",TotalLength);
	//mexPrintf("Height = %d\n",Height);
	//mexPrintf("Width = %d\n",Width);
	//mexPrintf("WinSize = %d\n",WinSize);
		
	int maxX,maxY,minX,minY;
	// output array
	plhs[0] = mxCreateDoubleMatrix(1,TotalLength,mxREAL);
	plhs[1] = mxCreateDoubleMatrix(1,TotalLength,mxREAL);
	plhs[2] = mxCreateDoubleMatrix(1,TotalLength,mxREAL);
	plhs[3] = mxCreateDoubleMatrix(Length,Length,mxREAL);

    row_inds = mxGetPr(plhs[0]);
    col_inds = mxGetPr(plhs[1]);
	Dist = mxGetPr(plhs[2]);
// 	MaxMatrix = mxGetPr(plhs[3]);
// 	MinMatrix = mxGetPr(plhs[4]);
	BoundaryCount = mxGetPr(plhs[3]);
	// input array
	lab = mxGetPr(prhs[0]);
	Label = mxGetPr(prhs[1]);
	Mask = mxGetPr(prhs[2]);
	//ComputeTimes = mxGetPr(prhs[3]);   
	int Frame = Width*Height;
	int len = 0;
	int count = 0;
	int IndexT,IndexMid;
	int MaxMinIndex;

	for(int j=0;j<Width;j++){ //360
		for(int i=0;i<Height;i++){ //240
			if(Mask[count]==1){
		    minX = max(i-WinSize,0);
			maxX = min(i+WinSize,Height-1);
			minY = max(j-WinSize,0);
			maxY = min(j+WinSize,Width-1);
			
			IndexMid = j*Height+i;

			int counxy = 0;
			for(int v=minY;v<=maxY;v++){
				for(int u=minX;u<=maxX;u++){
					IndexT = v*Height+u;
					//IndexT = v*Height+u+1;
					row_inds[len+counxy] = Label[IndexMid];
					col_inds[len+counxy] = Label[IndexT];
                   
                    Dist[len+counxy] = sqrt((lab[IndexT]-lab[IndexMid])*(lab[IndexT]-lab[IndexMid])+(lab[IndexT+Frame]-lab[IndexMid+Frame])*(lab[IndexT+Frame]-lab[IndexMid+Frame])+(lab[IndexT+2*Frame]-lab[IndexMid+2*Frame])*(lab[IndexT+2*Frame]-lab[IndexMid+2*Frame]));

					MaxMinIndex =  (Label[IndexMid]-1)+Length*(Label[IndexT]-1);
                    
                    BoundaryCount[MaxMinIndex] = BoundaryCount[MaxMinIndex] + 1;
                    
                    counxy++;
				}
			}
			len = len + counxy;
			}
			count++;
		}
	}
  
  return;  
  
}

