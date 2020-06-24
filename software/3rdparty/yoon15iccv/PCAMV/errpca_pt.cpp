// ERRPCA_PT.CPP
//
// errMx = ERRPCA_PT( X, A, S, numCPU ) computes a sparse matrix errMx
// of reconstruction errors (X - A*S) for the given sparse matrix X.
//
// numCPU specifies the number of CPUs used for parallel computing
// (default 1).
//
// See also COMPUTE_RMS.M, CF_PT.M

// Equivalent Matlab code:
//   M = spones(X);
//   errMx = (X - A*S).*M;

// This software is provided "as is", without warranty of any kind.
// Alexander Ilin, Tapani Raiko

#include <math.h>
#include <string.h>
#include "mex.h"
#include "matrix.h"
#ifndef NOTHREADS
#define USETHREADS
#include <pthread.h>
#endif

#ifdef OLDMATLABAPI
typedef int mwIndex;
typedef int mwSize;
#endif

typedef struct
{
  double     *A, *S, *X, *Err;
  mwIndex    *ir, *jc;
  mwSize     ndata;
  mwIndex    tfirst,jx;
  mwSize     ncomp;
  mwSize        n1;
}
TParams;

#ifdef USETHREADS
void *thread_function(void *);
void ThreadComputations(TParams*);
#endif

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, 
                 const mxArray *prhs[])
{
  const mxArray *mxX, *mxA, *mxS;
  mxArray       *mxErr;
  double        *X, *A, *S, *Err;
  mwIndex       r,jx,six,aix;
  mwSize        ndata;                  // Number of observed values
  mwSize        n1;                     // Dimensionalities of the
  mwSize        n2;                     //  the data matrix
  mwSize        ncomp;                  // Number of components
  double        res;
  int           numCPU;                 // Number of threads
  mwIndex       *ir, *jc, k;
  mwSize        nzmax;

  mxX = prhs[0];
  mxA = prhs[1];
  mxS = prhs[2];

  X = (double *)mxGetPr( mxX );
  A = (double *)mxGetPr( mxA );
  S = (double *)mxGetPr( mxS );

  n1 = mxGetM( mxX );
  n2 = mxGetN( mxX );
  ncomp = mxGetN( mxA );
  ir = mxGetIr( mxX );
  jc = mxGetJc( mxX );
  nzmax = mxGetNzmax(mxX);
  ndata = jc[n2];

  // printf( "ndata: %d, n1: %d, n2: %d, ncomp: %d\n",
  //    ndata, n1, n2, ncomp );

  // Copy the structure of matrix X to output matrix Err
  mxErr = mxCreateSparse( n1, n2, nzmax, mxREAL );
  memcpy( mxGetIr( mxErr ), ir, nzmax*sizeof(mwIndex) );
  memcpy( mxGetJc( mxErr ), jc, (n2+1)*sizeof(mwIndex) );
  Err = mxGetPr(mxErr);
  plhs[0] = mxErr;

#ifdef USETHREADS
  if( nrhs < 4 )
      numCPU = 1;
  else
      numCPU = (int)*(double *)mxGetPr( prhs[3] );
#else
  numCPU = 1;
#endif

  if( numCPU == 1 )
  {
      jx = 0;
      for( r=0; r < ndata; r++ )
      {
          res = 0;
          while( r == jc[jx+1] )
              jx++;
          // printf( "(%d %d)", ir[r]+1, jx+1 );
          aix = ir[r];
          six = jx*ncomp;
          for( k=0; k<ncomp; k++ )
          {
              res += A[ aix ] * S[ six ];
              six++;
              aix += n1;
          }
          Err[r] = X[r] - res;
      }
      return;
  }


#ifdef USETHREADS
  /*******************************************************************
                    Multi-thread implementation
  *******************************************************************/
  mwIndex          cfirst;                 // First column for a thread
  mwIndex          tmp;
  mwIndex          tlast;                  // Last value for a thread
  pthread_t        *mythread;
  TParams          *tp;
  int              i;

  mythread = (pthread_t *)malloc( numCPU*sizeof(pthread_t) );
  tp = (TParams *)malloc( numCPU*sizeof(TParams) );

  for( i=0; i < numCPU; i++ )
  {
      // Common thread arguments
      tp[i].A = A;
      tp[i].S = S;
      tp[i].ncomp = ncomp;
      tp[i].n1 = n1;
      tp[i].Err = Err;
      tp[i].X = X;
      tp[i].ir = ir;
      tp[i].jc = jc;

      // Thread specific arguments
      cfirst = i * (mwIndex)floor( (double)n2 / numCPU );
      tp[i].jx = cfirst;
      tp[i].tfirst = jc[cfirst];

      if( i == numCPU-1 )
          tlast = ndata;
      else
      {
          tmp = (i+1) * (mwIndex)floor( (double)n2 / numCPU );
          tlast = jc[tmp];
      }

      tp[i].ndata = tlast - tp[i].tfirst;

     if( i < numCPU-1 )
     {
         if( pthread_create( &(mythread[i]), NULL, thread_function,
                             (void*)(&tp[i]) ) )
         {
             mexErrMsgTxt("Error creating thread.");
         }
     }
     else
     {
         ThreadComputations( tp + numCPU-1 );
     }
  }

  for( i=0; i < numCPU-1; i++ )
  {
      if( pthread_join( mythread[i], NULL ) )
      {
          printf("Error joining thread\n");
      }
  }

  return;

#endif // USETHREADS

}


#ifdef USETHREADS

/*
**  Thread function
*/
void *thread_function(void *arg) 
{   
  TParams*   tp = (TParams*)arg;
  ThreadComputations(tp);
  return(NULL);
}

void ThreadComputations(TParams* tp)
{   
  double        *X, *A, *S, *Err;
  mwIndex       r,jx,six,aix;
  mwSize        ncomp,n1;
  mwSize        ndata;                  // Number of observed values
  double        res;
  mwIndex       *ir, *jc, k;
  mwIndex       tfirst;

  A = tp->A;
  S = tp->S;
  X = tp->X;
  ir = tp->ir;
  jc = tp->jc;
  Err = tp->Err;
  ndata = tp->ndata;
  ncomp = tp->ncomp;
  n1 = tp->n1;
  tfirst = tp->tfirst;

  jx = tp->jx;

  for( r=tfirst; r < tfirst+ndata; r++ )
  {
      res = 0;
      while( r == jc[jx+1] )
          jx++;
      // printf( "(%d %d)", ir[r]+1, jx+1 );
      aix = ir[r];
      six = jx*ncomp;
      for( k=0; k<ncomp; k++ )
      {
          res += A[ aix ] * S[ six ];
          six++;
          aix += n1;
      }
      Err[r] = X[r] - res;
  }
  return;

}

#endif
