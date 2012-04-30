// -----------------------------------------------------------------------------
//
// SpinXFormGPU -- QuaternionMatrix.h
// 
// Code courtesy of Keenan Crane. For the original version, see:
// http://users.cms.caltech.edu/~keenan/src/spinxform_commandline.zip
// at http://users.cms.caltech.edu/~keenan/project_spinxform.html
//
// This version has been modified so that it returns a real matrix in CUSP's COO
// format.
//
// =============================================================================
// SpinXForm -- QuaternionMatrix.h
// Keenan Crane
// August 16, 2011
//
// QuaternionMatrix is a convenience class for representing sparse matrices with
// quaternion-valued entries.  Entries can be accessed via expressions like
//
//    A( i, j ) = Quaternion( 1., 2., 3., 4. );
//
// A QuaternionMatrix can be converted to a sparse matrix with real-valued
// entries by calling toRealCooFormat().
//

#ifndef SPINXFORM_QUATERNIONMATRIX_H
#define SPINXFORM_QUATERNIONMATRIX_H

#include <map>
#include <vector>
#include <iostream>
#include "Quaternion.h"
#include "sparse_matrix.h"

#include <cusp/coo_matrix.h>
#include <cusp/print.h>
#include <cusp/csr_matrix.h>
#include <cusp/io/matrix_market.h>

#include <thrust/device_ptr.h>
#include <thrust/device_malloc.h>
#include <thrust/device_free.h>
#include <thrust/device_vector.h>
#include <thrust/copy.h>

//#include <cusp/dia_matrix.h>
//#include <cusp/ell_matrix.h>
#include <cusp/hyb_matrix.h>


class QuaternionMatrix
{
   public:
      void resize( int m, int n );
      // allocates an mxn matrix of zeros
      
      int size( int dim ) const;
      // returns the size of the dimension specified by scalar dim

            Quaternion& operator()( int& row, int& col );
      const Quaternion& operator()( int row, int col ) const;
      // access element (row,col)
      // note: uses 0-based indexing
      
      // define type of CUSP's COO matrix --------------------------------------
      
      // which index type to use
      typedef size_t IndexType;

      // where to perform the computation
      typedef cusp::host_memory MemorySpace;

      // which floating point type to use
      typedef float ValueType;
   
      typedef cusp::coo_matrix<IndexType, ValueType, MemorySpace> coo_cusp; 
      // -----------------------------------------------------------------------
      
      coo_cusp toRealCooFormat( void );
      // returns real matrix in CUSP's COO format 
      // where each quaternion becomes a 4x4 block
      
   protected:
      typedef std::pair<int,int> EntryIndex; // NOTE: column THEN row! (makes it easier to build compressed format)
      typedef std::map<EntryIndex, Quaternion> EntryMap;

      EntryMap data;
      // non-zero entries

      int m, n;
      // rows, columns

      static Quaternion zero;
      // dummy value for const access of zeros
      
      SparseMatrixf A;
      // *NOTE* for doubles use SparseMatrixd A;
      
};

#endif
