// SUBTRACT_MU.CPP
//
// X = SUBTRACT_MU( X, Mu ) subtracts bias term Mu from sparse data
// matrix X.

// Equivalent Matlab code:
//   M = spones(X);
//   X = X - repmat(Mu,1,size(X,2)).*M;

// This software is provided "as is", without warranty of any kind.
// Alexander Ilin, Tapani Raiko

#include <math.h>
#include <string.h>
#include "mex.h"
#include "matrix.h"

#define EPS 1e-15

#ifdef OLDMATLABAPI
typedef int mwIndex;
typedef int mwSize;
#endif

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, 
                 const mxArray *prhs[])
{
  const mxArray *mxX, *mxMu;
  mxArray       *mxXout;
  double        *X, *Mu, *Xout;
  mwIndex       r;
  mwSize        ndata;                  // Number of observed values
  mwSize        n1;                     // Dimensionalities of the
  mwSize        n2;                     //  the data matrix
  mwSize        ncomp;                  // Number of components
  int           numCPU;                 // Number of threads
  mwIndex       *ir, *jc;
  mwSize        nzmax;

  mxX = prhs[0];
  mxMu = prhs[1];

  n1 = mxGetM( mxX );
  n2 = mxGetN( mxX );

  if( mxGetM(mxMu) != n1 )
      mexErrMsgTxt(
      "SUTRACT_MU: Mu should contain the same number of columns as X" );

  X = (double *)mxGetPr( mxX );
  Mu = (double *)mxGetPr( mxMu );

  ir = mxGetIr( mxX );
  jc = mxGetJc( mxX );
  nzmax = mxGetNzmax(mxX);
  ndata = jc[n2];
  
  // Unsafe solution
  // mxXout = (mxArray *)prhs[0];

  // Copy the structure of matrix X to output matrix Xout
  mxXout = mxCreateSparse( n1, n2, nzmax, mxREAL );
  memcpy( mxGetIr( mxXout ), ir, nzmax*sizeof(mwIndex) );
  memcpy( mxGetJc( mxXout ), jc, (n2+1)*sizeof(mwIndex) );
  Xout = mxGetPr(mxXout);
  plhs[0] = mxXout;

  numCPU = 1;
  if( numCPU == 1 )
  {
      for( r=0; r < ndata; r++ )
      {
          Xout[r] = X[r] - Mu[ir[r]];
          if( !Xout[r] )
              Xout[r] = EPS;
      }
      return;
  }

}
