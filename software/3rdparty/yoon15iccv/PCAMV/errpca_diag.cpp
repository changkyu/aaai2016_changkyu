// ERRPCA_DIAG.CPP
//
// [errMx,varcost] = ERRPCA_DIAG( X, A, S, numCPU ) computes a sparse
// matrix errMx of reconstruction errors (X - A*S) for the given
// sparse matrix X. Output parameter varcost is needed to compute the
// variance of the observation noise in PCA_DIAG. It is calculated
// using Sv and Av which is the posterior variances of S and A
// (assuming fully factorial posterior for both).
//
// numCPU specifies the number of CPUs used for parallel computing
// (default 1).
//
// See also CF_DIAG.M

// Equivalent Matlab code:
//   M = spones(X);
//   errMx = (X - A*S).*M;
//   if isempty(Muv)
//       varcost = ( (A.^2)*Sv ).*M;
//   else
//       varcost = ( Av*S.^2 + (A.^2)*Sv + Av*Sv ).*M;
//   end
//   varcost = full(sum(sum(varcost)));

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

#define METHOD_PT 0
#define METHOD_VBD 1
#define METHOD_PPCAD 2

#ifdef OLDMATLABAPI
typedef int mwIndex;
typedef int mwSize;
#endif

typedef struct
{
  double     *A, *S, *X, *ErrMx, *Sv, *Av, varcost;
  mwIndex    *ir, *jc;
  mwSize     ndata;
  mwIndex    tfirst,jx;
  mwSize     ncomp;
  mwSize     n1;
  int        method;
}
TParams;

#ifdef USETHREADS
void *thread_function(void *);
void ThreadComputations(TParams*);
#endif

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, 
                 const mxArray *prhs[])
{
  const mxArray *mxX, *mxA, *mxS, *mxAv, *mxSv;
  mxArray       *mxErrMx;
  double        *X, *A, *S, *Av, *Sv, *ErrMx;
  double        *pvarcost, varcost;
  mwIndex       r,jx,six,aix;
  mwSize        ndata;                  // Number of observed values
  mwSize        n1;                     // Dimensionalities of the
  mwSize        n2;                     //  the data matrix
  mwSize        ncomp;                  // Number of components
  double        res;
  int           numCPU;                 // Number of threads
  mwIndex       *ir, *jc, k;
  mwSize        nzmax;
  int           method;
  double        ak, sk, avk, svk;

  mxX = prhs[0];
  mxA = prhs[1];
  mxS = prhs[2];

  Sv = NULL; Av = NULL;
  if( nrhs > 3 )
  {
      mxSv = prhs[3];
      if( ! mxIsEmpty( mxSv ) )
          Sv = (double *)mxGetPr( mxSv );
  }
  if( nrhs > 4 )
  {
      mxAv = prhs[4];
      if( ! mxIsEmpty( mxAv ) )
          Av = (double *)mxGetPr( mxAv );
  }

  method = METHOD_PT;
  if( nlhs > 1 )
      if( Av != NULL && Sv != NULL )
          method = METHOD_VBD;
      else if( Av == NULL && Sv != NULL )
          method = METHOD_PPCAD;

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

  // Copy the structure of matrix X to output matrix ErrMx
  mxErrMx = mxCreateSparse( n1, n2, nzmax, mxREAL );
  memcpy( mxGetIr( mxErrMx ), ir, nzmax*sizeof(mwIndex) );
  memcpy( mxGetJc( mxErrMx ), jc, (n2+1)*sizeof(mwIndex) );
  ErrMx = mxGetPr(mxErrMx);
  plhs[0] = mxErrMx;

  if( nlhs > 1 )
  {
      plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
      pvarcost = mxGetPr(plhs[1]);
  }

#ifdef USETHREADS
  if( nrhs < 6 )
      numCPU = 1;
  else
      numCPU = (int)*(double *)mxGetPr( prhs[5] );
#else
  numCPU = 1;
#endif

  if( numCPU == 1 )
  {
      jx = 0;
      switch( method )
      {
      case METHOD_PT:
          for( r=0; r < ndata; r++ )
          {
              res = 0;
              while( r == jc[jx+1] ) jx++; // Move to next column
              aix = ir[r];
              six = jx*ncomp;
              
              for( k=0; k<ncomp; k++ )
              {
                  res += A[ aix ] * S[ six ];
                  six++; aix += n1;
              }
              ErrMx[r] = X[r] - res;
          }
          varcost = 0;
          break;

      case METHOD_VBD:
          varcost = 0;
          for( r=0; r < ndata; r++ )
          {
              res = 0;
              while( r == jc[jx+1] ) jx++; // Move to next column
              aix = ir[r];
              six = jx*ncomp;
              
              for( k=0; k<ncomp; k++ )
              {
                  ak = A[aix]; avk = Av[aix];
                  sk = S[six]; svk = Sv[six];
                  
                  res += ak *sk;
                  varcost += avk * sk*sk + ak*ak * svk + avk * svk;
                  six++; aix += n1;
              }
              ErrMx[r] = X[r] - res;
          }
          break;

      case METHOD_PPCAD:
          varcost = 0;
          for( r=0; r < ndata; r++ )
          {
              res = 0;
              while( r == jc[jx+1] ) jx++; // Move to next column
              aix = ir[r];
              six = jx*ncomp;
              
              for( k=0; k<ncomp; k++ )
              {
                  ak = A[aix];
                  sk = S[six]; svk = Sv[six];
                  
                  res += ak *sk;
                  varcost += ak*ak * svk;
                  six++; aix += n1;
              }
              ErrMx[r] = X[r] - res;
          }
          break;
      }

      if( nlhs > 1 )
          (*pvarcost) = varcost;
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
      tp[i].A = A; tp[i].Av = Av;
      tp[i].S = S; tp[i].Sv = Sv;
      tp[i].ncomp = ncomp;
      tp[i].n1 = n1;
      tp[i].ErrMx = ErrMx;
      tp[i].X = X;
      tp[i].ir = ir;
      tp[i].jc = jc;
      tp[i].method = method;

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

  if( nlhs > 1 )
  {
      (*pvarcost) = 0;
      for( i=0; i < numCPU; i++ )
          (*pvarcost) += tp[i].varcost;
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
  double        *X, *A, *S, *ErrMx, *Av, *Sv;
  mwIndex       r,jx,six,aix;
  mwSize        ncomp,n1;
  mwSize        ndata;                  // Number of observed values
  double        res;
  mwIndex       *ir, *jc, k;
  mwIndex       tfirst;
  double        ak, sk, avk, svk, varcost;

  A = tp->A; Av = tp->Av;
  S = tp->S; Sv = tp->Sv;
  X = tp->X;
  ir = tp->ir;
  jc = tp->jc;
  ErrMx = tp->ErrMx;
  ndata = tp->ndata;
  ncomp = tp->ncomp;
  n1 = tp->n1;
  tfirst = tp->tfirst;

  jx = tp->jx;
  switch( tp->method )
  {
  case METHOD_PT:
      for( r=tfirst; r < tfirst+ndata; r++ )
      {
          res = 0;
          while( r == jc[jx+1] ) jx++;
          aix = ir[r];
          six = jx*ncomp;
          for( k=0; k<ncomp; k++ )
          {
              res += A[ aix ] * S[ six ];
              six++; aix += n1;
          }
          ErrMx[r] = X[r] - res;
      }
      tp->varcost = 0;
      return;

  case METHOD_VBD:
      varcost = 0;
      for( r=tfirst; r < tfirst+ndata; r++ )
      {
          res = 0;
          while( r == jc[jx+1] ) jx++; // Move to next column
          aix = ir[r];
          six = jx*ncomp;
          
          for( k=0; k<ncomp; k++ )
          {
              ak = A[aix]; avk = Av[aix];
              sk = S[six]; svk = Sv[six];
              
              res += ak *sk;
              varcost += avk * sk*sk + ak*ak * svk + avk * svk;
              six++; aix += n1;
          }
          ErrMx[r] = X[r] - res;
      }
      tp->varcost = varcost;
      return;
      
  case METHOD_PPCAD:
      varcost = 0;
      for( r=tfirst; r < tfirst+ndata; r++ )
      {
          res = 0;
          while( r == jc[jx+1] ) jx++; // Move to next column
          aix = ir[r];
          six = jx*ncomp;
          
          for( k=0; k<ncomp; k++ )
          {
              ak = A[aix];
              sk = S[six]; svk = Sv[six];
              
              res += ak *sk;
              varcost += ak*ak * svk;
              six++; aix += n1;
          }
          ErrMx[r] = X[r] - res;
      }
      tp->varcost = varcost;
      return;
  }
}

#endif
