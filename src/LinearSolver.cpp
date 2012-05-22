// -----------------------------------------------------------------------------
//
// SpinXFormGPU -- LinearSolver.cpp
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
// =============================================================================
// SpinXForm -- LinearSolver.cpp
// Keenan Crane
// August 16, 2011
//

#include "LinearSolver.h"
#include <iostream>
#include <cassert>

#include "cusp_device.h"
#include <cusp/print.h>
#include <cusp/io/matrix_market.h>

#include <thrust/host_vector.h> //convert solution back to quaternions

using namespace std;

void LinearSolver :: solve( QuaternionMatrix&   A,
                            vector<Quaternion>& x,
                            vector<Quaternion>& b,
                            bool precondition     ) {
// solves the linear system Ax = b where A is positive-semidefinite with a 
// conjugate gradient solver from CUSP library       
      
   typedef cusp::coo_matrix<int, float, cusp::host_memory> coo_cusp;
        
   // initialize C (C holds the matrix of reals in COO format) 
   coo_cusp C(0,0,0);
        
   C = A.toRealCooFormat();
   size_t row_indices_size = C.row_indices.size();
      
   //cusp::print(C);
   //std::cout << "row_indices_size " << row_indices_size << "\n";
   
   // sanity check (cusp::coo format needs to be sorted by row/col)
   //std::cout << "is C sorted_by_row? " << C.is_sorted_by_row() << "\n";// TRUE
   
   // allocate space for result, rhs
   // result/rhs are sparse -> they contain as many rows as the num of rows in C
   // C * result = rhs -> dimensions of vectors should match the size of C elements 
   // rhs[n*4] till rhs[row_indices_size -1] = 0 (cannot get rid of this) 
   vector<float> result( row_indices_size );
   vector<float> rhs(    row_indices_size ); 
       
   // convert right-hand side(b) to real values (rhs)
   LinearSolver::toReal( b, rhs );
   
   //std::cout << "rhs.size() " << rhs.size() << "\n";
   
   // allocate array1d (CUSP's format for a dense matrix) on the host for rhs
   cusp::array1d<float, cusp::host_memory> rhs_host = rhs;
   //cusp::io::write_matrix_market_file(rhs_host, "rhs_host.mtx");
   //cusp::print(rhs_host);
   
   // sanity check (make sure rhs_host and rhs have the same size)
   thrust::host_vector<float> rhs_host_thrust( rhs_host.begin(), rhs_host.end() );
   assert( rhs_host_thrust.size() == rhs.size() );

   // allocate array1d on the host for result
   cusp::array1d<float, cusp::host_memory> result_host = result;
   //cusp::io::write_matrix_market_file( result_host, "result_host_before_cu.mtx" );
   
   // sanity check (make sure the sizes of result_host and result are the same)
   thrust::host_vector<float> result_host_thrust( result_host.begin(), 
                                                  result_host.end()   );
   assert( result_host_thrust.size() == result_host.size() );
   
   // calls cusp_device.cu and solves linear system on the device  
   solve_on_device( C, rhs_host, result_host );
   //cusp::io::write_matrix_market_file(result_host, "result_host_after_cu_no_views.mtx");
   
   //convert solution back to quaternions
   thrust::host_vector<float> thrust_result_to_quat( result_host.begin(), 
                                                     result_host.end()   );
   
   //vector<float> result_in_std_format( row_indices_size ); // same with below
   vector<float> result_in_std_format( thrust_result_to_quat.size() );
      
   thrust::copy( thrust_result_to_quat.begin(), thrust_result_to_quat.end(), 
                 result_in_std_format.begin() );
        
   //convert solution back to quaternions
   toQuat( result_in_std_format, x );
					
}

void LinearSolver :: toReal( const vector<Quaternion>& uQuat,
                             vector<float>& uReal )
// converts vector from quaternion- to real-valued entries
{
   for( size_t i = 0; i != uQuat.size(); i++ )
   {
      uReal[i*4+0] = uQuat[i].re();   // real
      uReal[i*4+1] = uQuat[i].im().x; // i
      uReal[i*4+2] = uQuat[i].im().y; // j
      uReal[i*4+3] = uQuat[i].im().z; // k
   }
}

//void LinearSolver :: toQuat( const vector<float>& uReal,
void LinearSolver :: toQuat( vector<float>& uReal,
                             vector<Quaternion>& uQuat )
// converts vector from real- to quaternion-valued entries
{
   for( size_t i = 0; i != uQuat.size(); i++ )
   {
      uQuat[i] = Quaternion( uReal[i*4+0],   // real
                             uReal[i*4+1],   // i
                             uReal[i*4+2],   // j
                             uReal[i*4+3] ); // k
   }
}