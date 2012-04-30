// -----------------------------------------------------------------------------
//
// SpinXFormGPU -- EigenSolver.h
// 
// Code courtesy of Keenan Crane. 
//
// On the Power iteration method: 
// The method is based on the idea that if a given vector is repeatedly applied 
// to a matrix, and is properly normalized, then ultimately, it will lie in the
// direction of the eigenvector associated with the eigenvalues which are largest 
// in absolute value. The rate of convergence for the Power iteration depends on 
// the ratio of the second largest eigenvalue (in absolute value) to the largest 
// eigenvalue (in absolute value) and for many applications this leads to 
// unacceptably slow convergence. The method can be problematic if one wants to 
// compute a number of extremal eigenvalues. The Power iteration is still in use, 
// but most frequently as (implicit) part of more efficient techniques, e.g. 
// Krylov methods, inverse iteration, QR-method.
//
// Eigenvalue computation in the 20th century, G Golub (2000)
// Journal of Computational and Applied Mathematics, 123 (1-2), p. 35-65
//
// =============================================================================
// SpinXForm -- EigenSolver.h
// Keenan Crane
// August 16, 2011
//
// EigenSolver is used to solve the sparse eigenvalue problem
//
//    Ax = cx
//
// for the eigenvector x with smallest eigenvalue c, where A is a
// positive-definite matrix.  It implements the most basic inverse
// iteration scheme
//
//    x_{n+1} = A^-1 * x_n,
//
// which is equivalent to iteratively solving the linear system
//
//    A x_{n+1} = x_n, where x_0 is the starting estimate
//
// Currently the solver applies a fixed number of iterations and no
// other convergence criteria are considered.  Better results might be
// obtained by using a more sophisticated eigensolver such as SLEPc,
// PRIMME, or ARPACK (available as "eigs" in MATLAB).  Additionally, if
// a good initial guess x_0 for the eigenvector is available, faster
// convergence might be achieved by computing the eigenvalue estimate
//
//    c0 = (x'*A*x)/(x'*x), where x' == x transpose
//
// and applying the same inverse iteration scheme to the shifted
// matrix B = A-c0*I.
//


#ifndef SPINXFORM_EIGENSOLVER_H
#define SPINXFORM_EIGENSOLVER_H

#include "QuaternionMatrix.h"
#include <vector>

using namespace std;

class EigenSolver
{
   public:
      static void solve( QuaternionMatrix& A,
                         vector<Quaternion>& x );
      // solves the eigenvalue problem Ax = cx for the
      // eigenvector x with the smallest eigenvalue c

   protected:
      static void normalize( vector<Quaternion>& x );
      // rescales x to have unit length
};

#endif
