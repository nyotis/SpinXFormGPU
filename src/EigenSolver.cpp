// -----------------------------------------------------------------------------
//
// SpinXFormGPU -- EigenSolver.cpp
// 
// Code courtesy of Keenan Crane. 
//
// =============================================================================
// SpinXForm -- EigenSolver.cpp
// Keenan Crane
// August 16, 2011
//

#include "EigenSolver.h"
#include "LinearSolver.h"
#include <cmath>

void EigenSolver :: solve( QuaternionMatrix& A,
                           vector<Quaternion>& x )
// solves the eigenvalue problem Ax = cx for the
// eigenvector x with the smallest eigenvalue c
{
   // set the initial guess to the identity
   vector<Quaternion> b( x.size(), 1. );

   // perform a fixed number of inverse power iterations
   const int nIter = 3;
   for( int i = 0; i != nIter; i++ )
   {
      normalize( b );
      LinearSolver::solve( A, x, b, false );
      b = x;
   }

   // normalize the final solution
   normalize( x );
}

void EigenSolver :: normalize( vector<Quaternion>& x )
// rescales x to have unit length
{
   // compute length
   float norm = 0.;
   for( size_t i = 0; i != x.size(); i++ )
   {
      norm += x[i].norm2();
   }
   norm = sqrt( norm );

   // normalize
   for( size_t i = 0; i != x.size(); i++ )
   {
      x[i] /= norm;
   }
}

