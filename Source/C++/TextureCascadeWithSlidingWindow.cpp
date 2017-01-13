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

// [Test1T Test2T Test3T] = TextureCascadeWithSlidingWindow(Mask1,Mask2,Mask3,Index,IndexT,Windows,Height,Width);

void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{
    int Windows,Height,Width;
	double *Test1,*Test2,*Test3;
	double *Mask1,*Mask2,*Mask3,*Index,*IndexT;
    if (nlhs>3){
		mexErrMsgTxt("one return values");}
    if (nrhs!=8){
		mexErrMsgTxt("Require seven inputs");}

    Windows = mxGetScalar(prhs[5]);
	Height = mxGetScalar(prhs[6]);
	Width = mxGetScalar(prhs[7]);
    	
    int R = 2*Windows+1;
	// output array
	plhs[0] = mxCreateDoubleMatrix(Height*Width,R*R,mxREAL);
    plhs[1] = mxCreateDoubleMatrix(Height*Width,R*R,mxREAL);
    plhs[2] = mxCreateDoubleMatrix(Height*Width,R*R,mxREAL);
	
	Test1 = mxGetPr(plhs[0]);
    Test2 = mxGetPr(plhs[1]);
    Test3 = mxGetPr(plhs[2]);
    
    // input array
	Mask1 = mxGetPr(prhs[0]);
	Mask2 = mxGetPr(prhs[1]);
    Mask3 = mxGetPr(prhs[2]);
    Index = mxGetPr(prhs[3]);
    IndexT = mxGetPr(prhs[4]);
	
    int IndexMid,count,IndexNei;

    for(int j=0+Windows;j<Width-Windows;j++){
        	for(int i=0+Windows;i<Height-Windows;i++){
            IndexMid = Index[j*(Height)+i]-1;
            count = 0;
            for(int u=j-Windows;u<=j+Windows;u++){
                 for(int v=i-Windows;v<=i+Windows;v++){
                    IndexNei = Index[u*(Height)+v]-1;
                    Test1[IndexMid+count*Height*Width] = Mask1[IndexNei];
                    Test2[IndexMid+count*Height*Width] = Mask2[IndexNei];
                    Test3[IndexMid+count*Height*Width] = Mask3[IndexNei];
                    count++;
                }
            }
               
        }
    }
    int IndexTempMid;	
     for(int j=0;j<Width;j++){
            for(int i=0;i<Height;i++){
            IndexMid = Index[j*(Height)+i]-1;
            if((Index[IndexMid]-IndexT[IndexMid])!=0){
               count = 0;
               IndexTempMid = IndexT[IndexMid]-1;
               for(int k=0;k<R*R;k++){
                   Test1[IndexMid+count*Height*Width] = Test1[IndexTempMid+count*Height*Width];
                   Test2[IndexMid+count*Height*Width] = Test2[IndexTempMid+count*Height*Width]; 
                   Test3[IndexMid+count*Height*Width] = Test3[IndexTempMid+count*Height*Width]; 
                   count++;
               }
            }
        }
    }
    
  return;  
  
}

