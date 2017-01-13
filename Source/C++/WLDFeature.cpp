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

// [excitation orientation] = WLDFeature(Gray,f00,R,Height,Width);

void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{
    int R,Height,Width;
	double *excitation,*orientation;
	double *Gray,*f00;
    if (nlhs>2){
		mexErrMsgTxt("two return values");}
    if (nrhs!=5){
		mexErrMsgTxt("Require five inputs");}

	Height = mxGetScalar(prhs[3]);
	Width = mxGetScalar(prhs[4]);
	R = mxGetScalar(prhs[2]);    
		
    
    double BELTA=5; // to avoid that center pixture is equal to zero
    double ALPHA=3; // like a lens to magnify or shrink the difference between neighbors
    double EPSILON=0.0000001;
    double PI=3.141592653589;
    int bsizey=2*R+1;
    int bsizex=2*R+1;
	int dx = Width-bsizex+1;
	int dy = Height-bsizey+1;
	// output array
	plhs[0] = mxCreateDoubleMatrix(Height,Width,mxREAL);
    plhs[1] = mxCreateDoubleMatrix(Height,Width,mxREAL);
//     plhs[0] = mxCreateDoubleMatrix(dy+R,dx+R,mxREAL);
	//plhs[1] = mxCreateDoubleMatrix(dx+1,dy+1,mxREAL);
	excitation = mxGetPr(plhs[0]);
    orientation = mxGetPr(plhs[1]);
	// input array
	Gray = mxGetPr(prhs[0]);
	f00 = mxGetPr(prhs[1]);
    int maxX,maxY,minX,minY;
    int IndexT,Index;
    double v00,v01,v10,v11;
    double N1,N3,N5,N7;
    int countConv;
	for(int j=0+R;j<Width-R;j++){
        for(int i=0+R;i<Height-R;i++){ 
            minX = max(i-R,0);
			maxX = min(i+R,Height-R);
			minY = max(j-R,0);
			maxY = min(j+R,Width-R);
            Index = j*Height+i;
//          step 1 compute differential excitationt
            v01 = Gray[Index] + BELTA;
            countConv = 0;
            v00 = 0;
            for(int v=minY;v<=maxY;v++){
				for(int u=minX;u<=maxX;u++){
                    IndexT = v*Height+u;
                    v00 = v00 + Gray[IndexT]*f00[countConv];
                    countConv = countConv + 1;
                }
            }
            Index = j*Height+i;
            if(v01!=0){
                excitation[Index] = atan(ALPHA*(v00/v01));
            }
            else{
                excitation[Index] = 0.1;
            }
            IndexT = (j-R)*Height+i;
            N7 = Gray[IndexT];
            IndexT = (j+R)*Height+i;
            N3 = Gray[IndexT];
            IndexT = j*Height+(i-R);
            N1 = Gray[IndexT];
            IndexT = j*Height+(i+R);
            N5 = Gray[IndexT];
            // step 2 compute gradient orientation
            if(abs(N7-N3)< EPSILON){
                orientation[Index] = 0;
            }
            else{
                v10=N5-N1;v11=N7-N3;
                orientation[Index] = atan(v10/v11);
                orientation[Index] = orientation[Index]*180/PI;
                if(v11 > EPSILON && v10 >  EPSILON){
                    orientation[Index] = orientation[Index] + 0;
                } 
                else if(v11 < -EPSILON && v10 >  EPSILON){
                    orientation[Index] = orientation[Index] + 180;
                }
                else if(v11 < -EPSILON && v10 < -EPSILON) {
                    orientation[Index] = orientation[Index] + 180;
                }
                else if(v11 >  EPSILON && v10 < -EPSILON){
                    orientation[Index] = orientation[Index] + 360;
                }
            }
            
		}
	}
  
  return;  
  
}

