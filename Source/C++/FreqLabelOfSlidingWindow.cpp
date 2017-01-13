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

// [Labelfreq] = FreqLabelOfSlidingWindow(Label,Windows,Height,Width);

void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{
    int Windows,Height,Width;
	double *Labelfreq;
	double *Label,*Index;
    if (nlhs>1){
		mexErrMsgTxt("one return values");}
    if (nrhs!=4){
		mexErrMsgTxt("Require four inputs");}

    Windows = mxGetScalar(prhs[1]);
	Height = mxGetScalar(prhs[2]);
	Width = mxGetScalar(prhs[3]);
    	
    int R = 2*Windows+1;
	// output array
	plhs[0] = mxCreateDoubleMatrix(Height*Width,R*R,mxREAL);
//     plhs[1] = mxCreateDoubleMatrix(Height*Width,R*R,mxREAL);
//     plhs[2] = mxCreateDoubleMatrix(Height*Width,R*R,mxREAL);
	
	Labelfreq = mxGetPr(plhs[0]);
//     Test2 = mxGetPr(plhs[1]);
//     Test3 = mxGetPr(plhs[2]);
    
    // input array
	Label = mxGetPr(prhs[0]);
// 	Mask2 = mxGetPr(prhs[1]);
//     Mask3 = mxGetPr(prhs[2]);
//     Label = mxGetPr(prhs[1]);
//     IndexT = mxGetPr(prhs[4]);
	
    int IndexMid,count,IndexNei;

    for(int j=0+Windows;j<Width-Windows;j++){
        	for(int i=0+Windows;i<Height-Windows;i++){
            IndexMid = j*(Height)+i;
            count = 0;
            for(int u=j-Windows;u<=j+Windows;u++){
                 for(int v=i-Windows;v<=i+Windows;v++){
                    IndexNei = u*(Height)+v;
                    Labelfreq[IndexMid+count*Height*Width] = Label[IndexNei];
//                     Test2[IndexMid+count*Height*Width] = Mask2[IndexNei];
//                     Test3[IndexMid+count*Height*Width] = Mask3[IndexNei];
                    count++;
                }
            }
               
        }
    }
   
  return;  
  
}

