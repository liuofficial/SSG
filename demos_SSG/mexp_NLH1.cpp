/* *************************************************************************
 * mexp_NLH1.cpp
 * ||u-im0|| + mu||HW(u)||
 * *************************************************************************/
#include "mex.h"
#include <SDKDDKVer.h>
#include <stdio.h>
#include <tchar.h>
#include <ppl.h>

#include <stdlib.h>
#include <math.h>
#include <time.h>
using namespace Concurrency;

#define X(ix,iy) (ix)*iNy + (iy) // location of image
//#define X2(ix,iy) (ix)*im + (iy) // location of patch
//#define X3(ix,iy) (ix)*iw + (iy) // location of search window
#define X4(ix,iy,i) (i)*iNyNx + (iy)*iNx + (ix)
#define X5(iX,i) (i)*iNyNx + (iX)
#define Xb(ix,iy) (iy)*iNx + (ix)

//#define XX(ip,i) (i)*nPiece + (ip) // location of 3D image

#define ABS(x) ( x >= 0.0 ? x : -x )
#define SQR(x) (x)*(x)

typedef float DAT;

float SQRT(float number) 
{
    long i;
    float x, y;
    const float f = 1.5F;
    
    x = number * 0.5F;
    y  = number;
    i  = * ( long * ) &y;
    i  = 0x5f3759df - ( i >> 1 );
    y  = * ( float * ) &i;
    y  = y * ( f - ( x * y * y ) );
    y  = y * ( f - ( x * y * y ) );
    return number * y;
}

void alg
(int ip, int iNy, int iNx, DAT *pfuNew, DAT *pfu, 
int inner, int iNbNeigh, DAT mu, DAT *pfIm0, int *piY, 
DAT *pfW, int nPiece, int iNyNx)
{
	int iy, ix, iIterU, i, ixY, iyY, iX, iXb;
	// init u d b
	for (iy=0; iy<iNy; iy++)
	{
		for (ix=0; ix<iNx; ix++)
		{
			//int XX_ip_X_ixiy = XX(ip,X(ix,iy));
			pfuNew[X(ix,iy)] = pfu[X(ix,iy)];
		}
	}
    
	// u
	/* 0.5*||u-u0|| + 0.5*mu*||HW(u)|| */
	for (iIterU=0; iIterU<inner; iIterU++)
		for (iy=0; iy<iNy; iy++)
			for(ix=0; ix<iNx; ix++)
	{          
		iX = X(ix,iy); iXb = Xb(ix,iy);
		DAT fSum1 = 0.0; DAT fSum2 = 0.0;
        for (i=0; i<iNbNeigh; i++)
        {
            ixY = (int) piY[X5(iXb,1+2*i)];
            iyY = (int) piY[X5(iXb,1+2*i+1)];
            DAT fw = pfW[X5(iXb,2*i)];
            fSum1 += fw;
            //fSum2 += fw * pfuNew[XX(ip,X(ixY,iyY))];
			fSum2 += fw * pfuNew[X(ixY,iyY)];
        }
        DAT fDen = 1.0 + mu*fSum1; // modified
		//int XX_ip_iX = XX(ip,iX);
        pfuNew[iX] = (mu*fSum2 + pfIm0[iX]) / fDen;
    }
}

////////////////////////////////////////////////////////////////////////////////////
////////////////////// main function //////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	const mxArray *mx;
	int rhs = 0;
	mx = prhs[rhs++]; 
	//int nPiece = (int) mxGetM(mx); // band of image
	int nPiece = (int) mxGetN(mx); // band of image
	DAT *pfu = (DAT*)mxGetData(mx); // image
	mx = prhs[rhs++];
	DAT *pfIm0 = (DAT*)mxGetData(mx); // image
	mx = prhs[rhs++];
	DAT *pfW = (DAT*)mxGetData(mx); // weight
	mx = prhs[rhs++];
	int *piY = (int*)mxGetData(mx); // location of weights
	mx = prhs[rhs++];
	DAT *pfVecGeneralParameters = (DAT*)mxGetData(mx); // parameters
	mx = prhs[rhs++];
	int nPal = (int)mxGetScalar(mx);
	
	int r = 0;
	int iNy = (int) pfVecGeneralParameters[r++]; // rows
    int iNx = (int) pfVecGeneralParameters[r++]; // cols
    //int im = (int) pfVecGeneralParameters[2]; // patch
    //int iw = (int) pfVecGeneralParameters[3]; // search window
    int iNbNeigh = (int) pfVecGeneralParameters[r++]; // neighboors
    DAT mu = pfVecGeneralParameters[r++];
	int inner = (int) pfVecGeneralParameters[r++];
	
	int iNyNx = iNy * iNx;
	if ((int)mxGetM(prhs[0]) != iNyNx) mexErrMsgTxt("Wrong input parameters.");
	
	/* return value */	
	int iNdim = 2;
	mwSize  iDim[3];
    //iDim[0] = nPiece;
    //iDim[1] = iNyNx;   
	iDim[0] = iNyNx;
    iDim[1] = nPiece;	
    plhs[0] = mxCreateNumericArray(iNdim,(const mwSize*)iDim,mxSINGLE_CLASS,mxREAL);
	mx = plhs[0];
    DAT *pfuNew = (DAT*)mxGetData(mx);
	
	if (nPal == 1)
	{
		structured_task_group tasks;
		tasks.run_and_wait([&]
		{
			parallel_for(int(0), int(nPiece), [&](int n)
			{
				alg(n, iNy, iNx, pfuNew+n*iNyNx, pfu+n*iNyNx, inner, iNbNeigh, mu, pfIm0+n*iNyNx, piY, pfW, nPiece, iNyNx);
			});
		});
	}
	else
	{
		for (int ip=0; ip<nPiece; ip++)
		{
			alg(ip, iNy, iNx, pfuNew+ip*iNyNx, pfu+ip*iNyNx, inner, iNbNeigh, mu, pfIm0+ip*iNyNx, piY, pfW, nPiece, iNyNx);
		}
	}
	
}