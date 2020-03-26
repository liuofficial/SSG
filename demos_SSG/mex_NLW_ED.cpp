/* *************************************************************************
 * mex_nl_weight.cpp
 * *************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <mex.h>
#include <math.h>

#define X(ix,iy) (ix)*iNy + (iy) // location of image
#define X2(ix,iy) (ix)*im + (iy) // location of patch
#define X3(ix,iy) (ix)*iw + (iy) // location of search window
#define X4(ix,iy,i) (i)*iNyNx + (iy)*iNx + (ix)

#define XX(is,i) (i)*iBand + (is) // location of 3D image

#define ABS(x) ( x >= 0.0 ? x : -x )
#define SQR(x) (x)*(x)

static union
{
  double d;
  struct 
  {
    int j,i;
  } n;
} _eco;
#define EXP_A (1048576/0.69314718055994530942)
#define EXP_C 60801
#define DEXP(y) (_eco.n.i = EXP_A*(y) + (1072693248 - EXP_C), _eco.d)
#define EXP(y) (float)DEXP(y)

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

void v4NearNeig
(
int    ix,
int    iy,
int    ix3,
int    iy3,
int    iCptX,
float  *pfW,
int    *piY,
int    iNy,
int    iNx,
int    iNyNx,
float  *pfW2,
int    im,
int    iw,
int    iw2,
int    iNbNeigh
)
{ 
    pfW[X4(ix,iy,2*iCptX)] = pfW2[X3(ix3+iw2,iy3+iw2)];
    piY[X4(ix,iy,1+2*iCptX)] = ix+ix3;
    piY[X4(ix,iy,1+2*iCptX+1)] = iy+iy3;
    pfW2[X3(ix3+iw2,iy3+iw2)] = 0.0;
}

void vUpdate
(
int    ix,
int    iy,
int    ix3,
int    iy3,
int    iCptX,
float  *pfW,
int    *piY,
int    iNy,
int    iNx,
int    iNyNx,
float  *pfW2,
int    im,
int    iw,
int    iw2,
int    iNbNeigh
)
{
    int   ixY, iyY, ixY2, iyY2, iCptXc, iCptY, iCptYc, iTest, iTest2, i2;
        
    ixY = ix+ix3; iyY = iy+iy3;
    iTest = 0;
    for (i2=0; i2<iCptX; i2++)
    {
        ixY2 = piY[X4(ix,iy,1+2*i2)];
        iyY2 = piY[X4(ix,iy,1+2*i2+1)];
        if ( ixY==ixY2 && iyY==iyY2 ) iTest = 1;
    }
    if ( iTest==0 )
    {
        iCptXc = iCptX;
        if ( iCptX<iNbNeigh )
        {
            iCptX++; piY[X4(ix,iy,0)] = iCptX;
            piY[X4(ix,iy,1+2*iCptXc)] = ixY;
            piY[X4(ix,iy,1+2*iCptXc+1)] = iyY;
            pfW[X4(ix,iy,2*iCptXc)] = pfW2[X3(ix3+iw2,iy3+iw2)];
            iCptY = piY[X4(ixY,iyY,0)];
            if ( iCptY<iNbNeigh )
            {
                iTest2 = 0;
                for (i2=0; i2<iCptY; i2++)
                {
                    ixY2 = piY[X4(ixY,iyY,1+2*i2)];
                    iyY2 = piY[X4(ixY,iyY,1+2*i2+1)];
                    if ( ix==ixY2 && iy==iyY2 ) iTest2 = 1;
                }
                if ( iTest2==0 )
                {
                    iCptYc = iCptY; iCptY++; piY[X4(ixY,iyY,0)] = iCptY;
                    piY[X4(ixY,iyY,1+2*iCptYc)] = ix;
                    piY[X4(ixY,iyY,1+2*iCptYc+1)] = iy;
                    pfW[X4(ixY,iyY,2*iCptYc)] = pfW2[X3(ix3+iw2,iy3+iw2)];
                }
            }
        }
    }
    pfW2[X3(ix3+iw2,iy3+iw2)] = 0.0;    
}

float ifun(int ix, int iy, int ix3, int iy3, int iBand, int im2, int im, int iNy,
float *pfpx, float *pfIm0, float *pfG, float fh, int iym, int iyp, int ixm, int ixp, float *pfbuf, int iw, int iw2)
{
	int iYx = ix+ix3;
	int iYy = iy+iy3;
	float fSum = 0.0;
	for (int iy4=-iym; iy4<=iyp; iy4++)
	{
		for (int ix4=-ixm; ix4<=ixp; ix4++)
		{
			float fTemp = 0.0;
			int id = (X(ix + ix4, iy + iy4)) * iw * iw + ((ix3 + iw2) * iw + iy3 + iw2);
			if (pfbuf[id] > 0) 
			{
				fSum += pfbuf[id] * pfG[X2(ix4+im2,iy4+im2)];
				continue;
			}	
			int jd = (X(iYx + ix4, iYy + iy4)) * iw * iw + ((-ix3 + iw2) * iw - iy3 + iw2);
			if (pfbuf[jd] > 0) 
			{
				fSum += pfbuf[jd] * pfG[X2(ix4+im2,iy4+im2)];
				continue;
			}	
			for (int is=0; is<iBand; is++)
			{
				float fx = pfpx[XX(is, X2(ix4+im2,iy4+im2))]; float fy = pfIm0[XX(is, X(iYx+ix4,iYy+iy4))];
				fTemp += SQR(fx - fy);
			}
			float fval = fTemp;
			pfbuf[id] = fval; pfbuf[jd] = fval;
			fSum += fval * pfG[X2(ix4+im2,iy4+im2)];
		}
	}
	//float fexp = EXP(-fSum/fh);
	float fexp = exp(-fSum/fh);
	if (fexp < 1e-9) fexp = 1e-9;
	return fexp;
}

/* *************************************************************************
 *                          main function
 * *************************************************************************/
 void mexFunction(int iNbOut, mxArray *pmxOut[], int iNbIn, const mxArray *pmxIn[])
 {
	/* iNbOut: number of outputs
		pmxOut: array of pointers to output arguments */
    
	/* iNbIn: number of inputs
		pmxIn: array of pointers to input arguments */
	
	float   *pfIm0, *pfG, *pfVecGeneralParameters, *pfW;
    float   *pfpx, *pfW2, f, fSum, fh, fMax;
    float   fTemp, fNorm1, fNorm2, fx, fy;
    int     iNy, iNx, iNdim;
    int     iNbNeigh, iN3, im, im2, iw, iw2, ic1;
    int     iy, ix, iym, ixm, iyp, ixp, iy2, ix2, iy3, ix3, iy4, ix4, ixMax, iyMax, i;
    int     iB, ib, iys, iye, ixs, ixe;
    int     iCptX, ixY, iyY, iTest, i2, iYx, iYy;
    int     iNbNearNeigh, iNbBestNeigh, ixYX, iyYX;
    int     *piY, iNyNx, iYX, iTemp;
	int     is, iBand;
	mwSize  iDim[3];
	float   fexp;
	float   *pfbuf;
	
	pfIm0 = (float*)mxGetData(pmxIn[0]); // image
	pfG = (float*)mxGetData(pmxIn[1]); // gauss
	pfVecGeneralParameters = (float*)mxGetData(pmxIn[2]); // parameters
	
	iBand = (int) mxGetM(pmxIn[0]); // band of image
	
    iNy = (int) pfVecGeneralParameters[0]; // row
    iNx = (int) pfVecGeneralParameters[1]; // col
    im = (int) pfVecGeneralParameters[2]; // patch
    iw = (int) pfVecGeneralParameters[3]; // search window
    iNbNearNeigh = (int) pfVecGeneralParameters[4]; // nearest neighbors
    iNbBestNeigh = (int) pfVecGeneralParameters[5]; // best neighbors
    fh = pfVecGeneralParameters[6]; // variation
	
	iNbNeigh = iNbNearNeigh + iNbBestNeigh; // all neighbors
	
	// weights
	iN3 = iNbNeigh * 2;
    iNdim = 3;
    iDim[0] = iNx;
    iDim[1] = iNy;
    iDim[2] = iN3;
    pmxOut[0] = mxCreateNumericArray(iNdim,(const mwSize*)iDim,mxSINGLE_CLASS,mxREAL);
    pfW = (float*)mxGetData(pmxOut[0]);
    
    // location
    iN3 = 1 + iNbNeigh * 2;
    iNdim = 3;
    iDim[0] = iNx;
    iDim[1] = iNy;
    iDim[2] = iN3;
    pmxOut[1] = mxCreateNumericArray(iNdim,(const mwSize*)iDim,mxINT32_CLASS,mxREAL);
    piY = (int *)mxGetData(pmxOut[1]);
	
	im2 = (im-1)/2;
    iw2 = (iw-1)/2;
    ic1 = im2+iw2;
	
	pfpx = (float *) calloc( (unsigned)(im*im*iBand), sizeof(float) ); // buf of patch
    if (!pfpx)
        mexPrintf("Memory allocation failure\n");
    
    pfW2 = (float *) calloc( (unsigned)(iw*iw), sizeof(float) ); // W^2
    if (!pfW2)
        mexPrintf("Memory allocation failure\n");
	
	iNyNx = iNy * iNx; // size of image
	
	// weights buf
	pfbuf = (float *) calloc( (unsigned)(iw*iw*iNyNx), sizeof(float) );
	for (i=0; i<iw*iw*iNyNx; i++) pfbuf[i] = -1.0;
	
	// init(0)
	for (iy=0; iy<iNy; iy++)
	{
		for (ix=0; ix<iNx; ix++)
		{
			for (i=0; i<iNbNeigh*2; i++) pfW[X4(ix,iy,i)] = 0.0;
			for (i=1; i< 1+iNbNeigh*2; i++) piY[X4(ix,iy,i)] = 0;
			piY[X4(ix,iy,0)] = iNbNearNeigh;
		}
	}
	
	
	// CENTRAL 
    for (iy=ic1; iy<iNy-ic1; iy++)
	{
		for (ix=ic1; ix<iNx-ic1; ix++)
		{
			// save central patches
			for (iy2=-im2; iy2<=im2; iy2++)
			{
				for (ix2=-im2; ix2<= im2; ix2++)
				{
					for (is=0; is<iBand; is++)
						pfpx[XX(is, X2(ix2+im2,iy2+im2))] = pfIm0[XX(is, X(ix+ix2,iy+iy2))];
				}
			}   
			
			// traverse search window
			for (iy3=-iw2; iy3<=iw2; iy3++)
			{
				for (ix3=-iw2; ix3<=iw2; ix3++)
				{
					fexp = ifun(ix, iy, ix3, iy3, iBand, im2, im, iNy, pfpx, pfIm0, pfG, fh, im2, im2, im2, im2, 
					pfbuf, iw, iw2);
					pfW2[X3(ix3+iw2,iy3+iw2)] = fexp;
				}
			}
			
			// init
			pfW2[X3(iw2,iw2)] = 0.0;
			// 4 nearest neighbors
			if (4 == iNbNearNeigh)
			{
				ix3=1; iy3=0; iCptX=0; v4NearNeig(ix,iy,ix3,iy3,iCptX,pfW,piY,iNy,iNx,iNyNx,pfW2,im,iw,iw2,iNbNeigh);
				ix3=-1; iy3=0; iCptX=1; v4NearNeig(ix,iy,ix3,iy3,iCptX,pfW,piY,iNy,iNx,iNyNx,pfW2,im,iw,iw2,iNbNeigh);
				ix3=0; iy3=1; iCptX=2; v4NearNeig(ix,iy,ix3,iy3,iCptX,pfW,piY,iNy,iNx,iNyNx,pfW2,im,iw,iw2,iNbNeigh);
				ix3=0; iy3=-1; iCptX=3; v4NearNeig(ix,iy,ix3,iy3,iCptX,pfW,piY,iNy,iNx,iNyNx,pfW2,im,iw,iw2,iNbNeigh);
			}
			else if (8 == iNbNearNeigh)
			{
				ix3=1; iy3=0; iCptX=0; v4NearNeig(ix,iy,ix3,iy3,iCptX,pfW,piY,iNy,iNx,iNyNx,pfW2,im,iw,iw2,iNbNeigh);
				ix3=-1; iy3=0; iCptX=1; v4NearNeig(ix,iy,ix3,iy3,iCptX,pfW,piY,iNy,iNx,iNyNx,pfW2,im,iw,iw2,iNbNeigh);
				ix3=0; iy3=1; iCptX=2; v4NearNeig(ix,iy,ix3,iy3,iCptX,pfW,piY,iNy,iNx,iNyNx,pfW2,im,iw,iw2,iNbNeigh);
				ix3=0; iy3=-1; iCptX=3; v4NearNeig(ix,iy,ix3,iy3,iCptX,pfW,piY,iNy,iNx,iNyNx,pfW2,im,iw,iw2,iNbNeigh);
				ix3=1; iy3=1; iCptX=4; v4NearNeig(ix,iy,ix3,iy3,iCptX,pfW,piY,iNy,iNx,iNyNx,pfW2,im,iw,iw2,iNbNeigh);
				ix3=1; iy3=-1; iCptX=5; v4NearNeig(ix,iy,ix3,iy3,iCptX,pfW,piY,iNy,iNx,iNyNx,pfW2,im,iw,iw2,iNbNeigh);
				ix3=-1; iy3=1; iCptX=6; v4NearNeig(ix,iy,ix3,iy3,iCptX,pfW,piY,iNy,iNx,iNyNx,pfW2,im,iw,iw2,iNbNeigh);
				ix3=-1; iy3=-1; iCptX=7; v4NearNeig(ix,iy,ix3,iy3,iCptX,pfW,piY,iNy,iNx,iNyNx,pfW2,im,iw,iw2,iNbNeigh);
			}
			else if (iNbNearNeigh != 0)
				mexPrintf("Number of nearest neighbors must be 0 or 4 or 8!!!\n");
			// END 4 nearest neighbors
			// find best patches
			iCptX = piY[X4(ix,iy,0)];
			while (iCptX < iNbNeigh)
			{
				fMax = 0.0;
				for (iy3=-iw2; iy3<=iw2; iy3++)
				{
					for (ix3=-iw2; ix3<=iw2; ix3++)
					{
						f = pfW2[X3(ix3+iw2,iy3+iw2)];
						if (f > fMax)
						{
							fMax = f;
							ixMax = ix3; iyMax = iy3;
						}
					}
				}       
				if (fMax > 0.001)
				{
					vUpdate(ix,iy,ixMax,iyMax,iCptX,pfW,piY,iNy,iNx,iNyNx,pfW2,im,iw,iw2,iNbNeigh);
					iCptX = piY[X4(ix,iy,0)];
				}
				else
					iCptX = iNbNeigh+1;
			}
			// END find best patches
		} // end of inner loop of central
	}
	// END OF CENTRAL
	
	// BORDERS
    for (ib=0; ib<4; ib++)
	{
		// 4 kinds
		if (ib==0) { iys=0; iye=ic1; ixs=0; ixe=iNx; }
        if (ib==1) { iys=iNy-ic1; iye=iNy; ixs=0; ixe=iNx; }
        if (ib==2) { iys=ic1; iye=iNy-ic1; ixs=0; ixe=ic1; }
        if (ib==3) { iys=ic1; iye=iNy-ic1; ixs=iNx-ic1; ixe=iNx; }
		// traverse
		for (iy=iys; iy<iye; iy++)
		{
            for (ix=ixs; ix<ixe; ix++)
			{
				if (iy-im2>=0) iym=im2; else iym=iy;
				if (iy+im2<=iNy-1) iyp=im2; else iyp=iNy-1-iy;
				if (ix-im2>=0) ixm=im2; else ixm=ix;
				if (ix+im2<=iNx-1) ixp=im2; else ixp=iNx-1-ix;
				// save central patches
				for (iy2=-iym; iy2<=iyp; iy2++)
				{
					for (ix2=-ixm; ix2<=ixp; ix2++)
					{
						for (is=0; is<iBand; is++)
							pfpx[XX(is, X2(ix2+im2,iy2+im2))] = pfIm0[XX(is, X(ix+ix2,iy+iy2))];
					}
				}
				// traverse search window
				for (iy3=-iw2; iy3<=iw2; iy3++)
				{
					for(ix3=-iw2; ix3<=iw2; ix3++)
					{
						iB = 0;
						if (iy+iy3-iym<0) iB=1;
						if (iy+iy3+iyp>iNy-1) iB=1;
						if (ix+ix3-ixm<0) iB=1;
						if (ix+ix3+ixp>iNx-1) iB=1;
						if (iB==0)
						{
							fexp = ifun(ix, iy, ix3, iy3, iBand, im2, im, iNy, pfpx, pfIm0, pfG, fh, iym, iyp, ixm, ixp, 
							pfbuf, iw, iw2);
							pfW2[X3(ix3+iw2,iy3+iw2)] = fexp;
						}
						else
							pfW2[X3(ix3+iw2,iy3+iw2)] = 0.0;
					}
				}
				// init
				pfW2[X3(iw2,iw2)] = 0.0;
				// 4 nearest neighbors
				if (iNbNearNeigh == 4)
				{
					ix3=1; iy3=0;
					if ( ix+ix3>=0 && ix+ix3<=iNx-1 && iy+iy3>=0 && iy+iy3<=iNy-1 )
					{ iCptX=0; v4NearNeig(ix,iy,ix3,iy3,iCptX,pfW,piY,iNy,iNx,iNyNx,pfW2,im,iw,iw2,iNbNeigh); }                
					ix3=-1; iy3=0; 
					if ( ix+ix3>=0 && ix+ix3<=iNx-1 && iy+iy3>=0 && iy+iy3<=iNy-1 )
					{ iCptX=1; v4NearNeig(ix,iy,ix3,iy3,iCptX,pfW,piY,iNy,iNx,iNyNx,pfW2,im,iw,iw2,iNbNeigh); }                
					ix3=0; iy3=1; 
					if ( ix+ix3>=0 && ix+ix3<=iNx-1 && iy+iy3>=0 && iy+iy3<=iNy-1 )
					{ iCptX=2; v4NearNeig(ix,iy,ix3,iy3,iCptX,pfW,piY,iNy,iNx,iNyNx,pfW2,im,iw,iw2,iNbNeigh); }                
					ix3=0; iy3=-1; 
					if ( ix+ix3>=0 && ix+ix3<=iNx-1 && iy+iy3>=0 && iy+iy3<=iNy-1 )
					{ iCptX=3; v4NearNeig(ix,iy,ix3,iy3,iCptX,pfW,piY,iNy,iNx,iNyNx,pfW2,im,iw,iw2,iNbNeigh); }
				}
				else if (iNbNearNeigh == 8)
				{
					ix3=1; iy3=0;
					if ( ix+ix3>=0 && ix+ix3<=iNx-1 && iy+iy3>=0 && iy+iy3<=iNy-1 )
					{ iCptX=0; v4NearNeig(ix,iy,ix3,iy3,iCptX,pfW,piY,iNy,iNx,iNyNx,pfW2,im,iw,iw2,iNbNeigh); }                
					ix3=-1; iy3=0; 
					if ( ix+ix3>=0 && ix+ix3<=iNx-1 && iy+iy3>=0 && iy+iy3<=iNy-1 )
					{ iCptX=1; v4NearNeig(ix,iy,ix3,iy3,iCptX,pfW,piY,iNy,iNx,iNyNx,pfW2,im,iw,iw2,iNbNeigh); }                
					ix3=0; iy3=1; 
					if ( ix+ix3>=0 && ix+ix3<=iNx-1 && iy+iy3>=0 && iy+iy3<=iNy-1 )
					{ iCptX=2; v4NearNeig(ix,iy,ix3,iy3,iCptX,pfW,piY,iNy,iNx,iNyNx,pfW2,im,iw,iw2,iNbNeigh); }                
					ix3=0; iy3=-1; 
					if ( ix+ix3>=0 && ix+ix3<=iNx-1 && iy+iy3>=0 && iy+iy3<=iNy-1 )
					{ iCptX=3; v4NearNeig(ix,iy,ix3,iy3,iCptX,pfW,piY,iNy,iNx,iNyNx,pfW2,im,iw,iw2,iNbNeigh); }
					ix3=1; iy3=1;
					if ( ix+ix3>=0 && ix+ix3<=iNx-1 && iy+iy3>=0 && iy+iy3<=iNy-1 )
					{ iCptX=0; v4NearNeig(ix,iy,ix3,iy3,iCptX,pfW,piY,iNy,iNx,iNyNx,pfW2,im,iw,iw2,iNbNeigh); }                
					ix3=1; iy3=-1; 
					if ( ix+ix3>=0 && ix+ix3<=iNx-1 && iy+iy3>=0 && iy+iy3<=iNy-1 )
					{ iCptX=1; v4NearNeig(ix,iy,ix3,iy3,iCptX,pfW,piY,iNy,iNx,iNyNx,pfW2,im,iw,iw2,iNbNeigh); }                
					ix3=-1; iy3=1; 
					if ( ix+ix3>=0 && ix+ix3<=iNx-1 && iy+iy3>=0 && iy+iy3<=iNy-1 )
					{ iCptX=2; v4NearNeig(ix,iy,ix3,iy3,iCptX,pfW,piY,iNy,iNx,iNyNx,pfW2,im,iw,iw2,iNbNeigh); }                
					ix3=-1; iy3=-1; 
					if ( ix+ix3>=0 && ix+ix3<=iNx-1 && iy+iy3>=0 && iy+iy3<=iNy-1 )
					{ iCptX=3; v4NearNeig(ix,iy,ix3,iy3,iCptX,pfW,piY,iNy,iNx,iNyNx,pfW2,im,iw,iw2,iNbNeigh); }
				}
				else if ( iNbNearNeigh!=0 )
					mexPrintf("Number of nearest neighbors must be 0 or 4 or 8!!!\n");
				// END 4 nearest neighbors
				// find best patches
				iCptX = piY[X4(ix,iy,0)];
				while (iCptX < iNbNeigh)
				{
					fMax = 0.0;
					for (iy3=-iw2; iy3<=iw2; iy3++)
						for(ix3=-iw2; ix3<=iw2; ix3++)
					{
						f = pfW2[X3(ix3+iw2,iy3+iw2)];
						if (f>fMax)
						{
							fMax = f;
							ixMax = ix3; iyMax = iy3;
						}
					}               
					if (fMax>0.001)
					{
						vUpdate(ix,iy,ixMax,iyMax,iCptX,pfW,piY,iNy,iNx,iNyNx,pfW2,im,iw,iw2,iNbNeigh);
						iCptX = piY[X4(ix,iy,0)];
					}
					else
						iCptX = iNbNeigh+1;
				}
				// END find best patches
			}
		}
	}
	// END OF BORDERS

	// create symmetric W (NOT PERFECT!!!!)
    for (iy=0; iy<iNy; iy++)
        for(ix=0; ix<iNx; ix++)
    {       
        for (i=0; i<iNbNeigh; i++) // see x's neighbors
        {
            ixY = piY[X4(ix,iy,1+2*i)]; 
            iyY = piY[X4(ix,iy,1+2*i+1)]; // neighbor y of x
            // see y's neighbors to check if x is already a neighbor of y. if it is the case, do nothing.
            iTest = 0;
            for (i2=0; i2<iNbNeigh; i2++)
            {
                ixYX = piY[X4(ixY,iyY,1+2*i2)];
                iyYX = piY[X4(ixY,iyY,1+2*i2+1)];
                if ( ix==ixYX && iy==iyYX ) {iTest=1; iYX=i2; i2=iNbNeigh;} // iTest=1 then one neighbor of y is x.
            }
            if (iTest==0) // then remove y in x's neighbors
            {
                piY[X4(ix,iy,1+2*i)] = 0;
                piY[X4(ix,iy,1+2*i+1)] = 0;
                pfW[X4(ix,iy,2*i)] = 0.0;
                pfW[X4(ix,iy,2*i+1)] = 0.0;
            }
            if (iTest==1) // then put y in x's neighbors at the same location of x in y's neighbors
            {
                if (i!=iYX)
                {
                    iTemp=piY[X4(ixY,iyY,1+2*i)]; piY[X4(ixY,iyY,1+2*i)]=piY[X4(ixY,iyY,1+2*iYX)]; piY[X4(ixY,iyY,1+2*iYX)]=iTemp;
                    iTemp=piY[X4(ixY,iyY,1+2*i+1)]; piY[X4(ixY,iyY,1+2*i+1)]=piY[X4(ixY,iyY,1+2*iYX+1)]; piY[X4(ixY,iyY,1+2*iYX+1)]=iTemp;
                    fTemp=pfW[X4(ixY,iyY,2*i)]; pfW[X4(ixY,iyY,2*i)]=pfW[X4(ixY,iyY,2*iYX)]; pfW[X4(ixY,iyY,2*iYX)]=fTemp;
                    fTemp=pfW[X4(ixY,iyY,2*i+1)]; pfW[X4(ixY,iyY,2*i+1)]=pfW[X4(ixY,iyY,2*iYX+1)]; pfW[X4(ixY,iyY,2*iYX+1)]=fTemp;
                }
            }
        }
    }
    // END create symmetric W

	
	// Normalizing weights
    for (iy=0; iy<iNy; iy++)
        for(ix=0; ix<iNx; ix++)
    {
        fSum = 0.0;
        for (i=0; i<iNbNeigh; i++) fSum += pfW[X4(ix,iy,2*i)];
        if (fSum>1e-3)
            for (i=0; i<iNbNeigh; i++)
        {
            pfW[X4(ix,iy,2*i)] /= fSum;
            pfW[X4(ix,iy,2*i+1)] = SQRT(pfW[X4(ix,iy,2*i)]);
        }
    }
	
	free((float *) pfpx);
    free((float *) pfW2);
	free((float *) pfbuf);

 }