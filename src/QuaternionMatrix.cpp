// -----------------------------------------------------------------------------
//
// SpinXFormGPU -- QuaternionMatrix.cpp
// 
// Code courtesy of Keenan Crane. For the original version, see:
// http://users.cms.caltech.edu/~keenan/src/spinxform_commandline.zip
// at http://users.cms.caltech.edu/~keenan/project_spinxform.html
//
// This version has been modified so that it returns a real matrix in CUSP's COO
// format.
//
// =============================================================================
// SpinXForm -- QuaternionMatrix.cpp
// Keenan Crane
// August 16, 2011
//

#include "QuaternionMatrix.h"
#include <iostream>
#include <fstream>

//#include <cusp/print.h>

using namespace std;

Quaternion QuaternionMatrix::zero( 0., 0., 0., 0. );
// dummy value for const access of zeros

void QuaternionMatrix :: resize( int _m, int _n )
// initialize an mxn matrix of zeros
{
   m = _m;
   n = _n;
   data.clear();
}

int QuaternionMatrix :: size( int dim ) const
// returns the size of the dimension specified by scalar dim
{
   if( dim == 1 ) return m;
   if( dim == 2 ) return n;
   return 0;
}

Quaternion& QuaternionMatrix :: operator()( int& row, int& col )
// return reference to element (row,col)
// note: uses 0-based indexing
{
   EntryIndex index( col, row ); // typedef std::pair<int,int> EntryIndex;
            
   EntryMap::const_iterator entry = data.find( index );
   if( entry == data.end())
   {
      data[ index ] = Quaternion( 0., 0., 0., 0. );
   }
 
   return data[ index ];
}

const Quaternion& QuaternionMatrix :: operator()( int row, int col ) const
// return const reference to element (row,col)
// note: uses 0-based indexing
{
   EntryIndex index( col, row );

   EntryMap::const_iterator entry = data.find( index );
   if( entry == data.end())
   {
      return zero;
   }

   return entry->second;
}

typedef cusp::coo_matrix<size_t, float, cusp::host_memory> coo_cusp;
coo_cusp QuaternionMatrix :: toRealCooFormat( void ) { // check const qualifier -> header
// return coo_cusp by value is *expensive*, but W is locally defined so cannot return it by reference    
    
   float Q[4][4];

   // convert quaternionic matrix to real matrix
   A.resize( n*4 );
  
   for( EntryMap::iterator e = data.begin(); e != data.end(); e++ )
   {
      int i = e->first.second; // row
      int j = e->first.first;  // column
      e->second.toMatrix( Q );

      for( int u = 0; u < 4; u++ )
      for( int v = 0; v < 4; v++ )
      {
         if( Q[u][v] != 0. )
         {
             A.set_element( i*4+u, j*4+v, Q[u][v] );
         }
      }
   }
   
   //std::cout <<  "dim of Sparse Matrix: " <<   A.n  << "\n";
  
   // save the index of each matrix's row for each pair of row and column indices in a nested 2D vector
   std::vector<int> rows;
   std::vector< std::vector<size_t> > cusp_rows = A.coo_format_rows( rows );
   //std::cout << cusp_rows.size() << "\n"; // 32016
   
   // hack to erase 0 content
   cusp_rows.erase( cusp_rows.begin(), cusp_rows.begin() + 16008 ); // hard-coded
   //cusp_rows.erase( cusp_rows.begin(), cusp_rows.begin() + n*4 ); // generic
   
   // uncomment the 3 lines below to get a feel how the cusp_rows index looks like
   // for ( size_t i = 0; i < cusp_rows.size(); ++i )
    //   for ( size_t int j =0; j < cusp_rows[i].size(); ++j )
     //      std::cout << "\n" << "cusp_rows[" << i << ", " << j <<  "]= " << cusp_rows[i][j] << " ";   
   
   // convert nested 2D cusp_rows to 1D vector
   std::vector<size_t> cusp_rows_1D; // *CAREFUL* not to use resize(), it adds a line of 0's in the beginning
   
   for ( size_t i=0; i < cusp_rows.size(); ++i )
      copy( cusp_rows[i].begin(), cusp_rows[i].end(), back_inserter(cusp_rows_1D) );
      
   cout << "\n" << "CUSP_COO_rows_1D_size: " << cusp_rows_1D.size() << "\n";
       
   // save the index of each matrix's column for each pair of row and column indices in a nested 2D vector
   std::vector<int> cols;
   std::vector< std::vector<size_t> > cusp_columns = A.coo_format_columns( cols );
   //std::cout << "cusp_columns size: " << cusp_columns.size() << "\n"; // 32016
   
   // hack to erase 0 content 
   cusp_columns.erase( cusp_columns.begin(), cusp_columns.begin() + 16008 ); // hard-coded
   //cusp_rows.erase( cusp_columns.begin(), cusp_columns.begin() + n*4 ); // generic
  
   /*
   // uncomment the 3 lines below to get a feel how the cusp_columns index looks like
   for ( size_t i = 0; i < cusp_columns.size(); ++i )
      for  ( size_t  int j =0; j < cusp_columns[i].size(); ++j )
         std::cout << "\n" << "cusp_columns[" << i << ", " << j <<  "]= " << cusp_columns[i][j] << " ";
   */
   //std::cout << "cusp_columns size: " << cusp_columns.size() << "\n";

   // convert nested 2D cusp columns to 1D vector
   std::vector<size_t> cusp_columns_1D; // *CAREFUL* not to use resize(), it adds a line of 0's in the beginning
   
   for ( size_t i=0; i < cusp_columns.size(); ++i )
      copy( cusp_columns[i].begin(), cusp_columns[i].end(), back_inserter( cusp_columns_1D ) );
      
   cout << "\n" << "CUSP_COO_columns_1D_size: " << cusp_columns_1D.size() << "\n";
   
   // CUSP's COO values --------------------------------------------------------
   
   // save each matrix's value for each pair of row 
   // and column indices in a nested 2D vector
   std::vector<float> vals;
   std::vector< std::vector<float> > cusp_values = A.coo_format_values( vals );   

   // hack to erase 0 content 
   cusp_values.erase( cusp_values.begin(), cusp_values.begin() + 16008 ); // hard-coded
   //cusp_values.erase( cusp_values.begin(), cusp_values.begin() + n*4 ); // generic
   
   // convert nested 2D cusp_values to 1D vector
   std::vector<float> cusp_values_1D; // careful with resize() adds a line of 0's in the beginning
     
   for ( size_t i=0; i < cusp_values.size(); ++i )
      copy( cusp_values[i].begin(), cusp_values[i].end(), back_inserter( cusp_values_1D ) );
	
   // uncomment the 3 lines to get a feel how the cusp_columns index looks like
   //for ( size_t i = 0; i < cusp_values.size(); ++i )    
     // for ( size_t j =0; j < cusp_values[i].size(); ++j )
       // std::cout << "\n" << "cusp_values[" << i << ", " << j <<  "]= " << cusp_values[i][j] << " ";

   cout << "\n" << "CUSP_COO_values_1D_size: " << cusp_values_1D.size() << "\n";
   
   // combine CUSP's COO rows, columns and values into a single COO matrix
   coo_cusp W( cusp_rows_1D.size(), cusp_columns_1D.size(), cusp_values_1D.size() ); 
   
   for( size_t i = 0; i < cusp_rows_1D.size(); ++i ) { // cusp_values_1D.size() gives same to cusp_rows_1D.size()   
        
      W.row_indices[i]    = cusp_rows_1D[i];
      W.column_indices[i] = cusp_columns_1D[i];
      W.values[i]         = cusp_values_1D[i];
   }
               
   //cusp::print(W);
   //cusp::io::write_matrix_market_file( W, "W.mtx" );
          
   return W;
   // return coo_cusp by value is *expensive*, but W is locally defined so cannot return it by reference  
} 
