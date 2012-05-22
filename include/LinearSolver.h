// -----------------------------------------------------------------------------
//
// SpinXFormGPU -- LinearSolver.h
// 
// Code courtesy of Keenan Crane. For the original version, see:
// http://users.cms.caltech.edu/~keenan/src/spinxform_commandline.zip
// at http://users.cms.caltech.edu/~keenan/project_spinxform.html
//
// This version has been modified in several ways:
//    - A conjugate gradient solver from CUSP library was added.
//    - Routines related to solvers used by SpinXForm, i.e. a sparse Cholesky 
//      factorization and a simple conjugate gradient (CG) solver, were trimmed.
//
// ============================================================================
// SpinXForm -- LinearSolver.h
// Keenan Crane
// August 16, 2011
//
// LinearSolver is used to solve the sparse linear system
//
//    Ax = b
//
// where A is a positive-definite matrix.

#ifndef SPINXFORM_LINEAR_SOLVER_H
#define SPINXFORM_LINEAR_SOLVER_H

#include "QuaternionMatrix.h"
#include <vector>

class LinearSolver {
   public:
       
      // solves the linear system Ax = b where A is positive-semidefinite 
      // with a conjugate gradient solver from CUSP library 
      static void solve( QuaternionMatrix&        A,
			 std::vector<Quaternion>& x,
                         std::vector<Quaternion>& b,
                         bool precondition = true   );
      
      // converts vector from quaternion- to real-valued entries
      static void toReal( const std::vector<Quaternion>& uQuat,
			                 std::vector<float>& uReal );
      
      // converts vector from real- to quaternion-valued entries
      static void toQuat( std::vector<float>& uReal,
		                    std::vector<Quaternion>& uQuat );
};

#endif
