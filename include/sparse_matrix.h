// =============================================================================
// sparse_matrix.h
//
// Code courtesy of Robert Bridson.  For the original version, see:
//    http://www.cs.ubc.ca/~rbridson/fluidsimulation/
//
// This version has been modified in several ways:
//    - Routines have been added for saving the sparse matrix in CUSP's COO format
//    - Fixed version of SparseMatrix and routines related to a different solver 
//      have been trimmed


#ifndef SPARSE_MATRIX_H
#define SPARSE_MATRIX_H

#include <iostream>
#include <vector>
#include <cassert>
#include <algorithm> // for std::max_element
#include "util.h"

//============================================================================
// Dynamic compressed sparse row matrix.

template<class T>
struct SparseMatrix
{
   size_t n; // dimension
   
   // for each row (with non-zero values), a list of all column indices (sorted)                          
   std::vector< std::vector<size_t> > index;
   std::vector< std::vector<T> >      value; // values corresponding to index
   
   explicit SparseMatrix(size_t n_=0, size_t expected_nonzeros_per_row=7)
      : n(n_), index(n_), value(n_)
   {
      for(size_t i=0; i<n; ++i){
         index[i].reserve(expected_nonzeros_per_row);
         value[i].reserve(expected_nonzeros_per_row);
      }
   }

   void clear(void)
   {
      n=0;
      index.clear();
      value.clear();
   }

   void zero(void)
   {
      for(size_t i=0; i<n; ++i){
         index[i].resize(0);
         value[i].resize(0);
      }
   }

   void resize(int n_)
   {
      n=n_;
      index.resize(n);
      value.resize(n);
   }

   T operator()(size_t i, size_t j) const
   {
      for(size_t k=0; k<index[i].size(); ++k){
         if(index[i][k]==j) return value[i][k];
         else if(index[i][k]>j) return 0;
      }
      return 0;
   }
   
   void set_element(size_t i, size_t j, T new_value)
   {
      size_t k=0;
      for(; k<index[i].size(); ++k){
         if(index[i][k]==j){
            value[i][k]=new_value;
            return;
         }else if(index[i][k]>j){
            insert(index[i], k, j);
            insert(value[i], k, new_value);
            return;
         }
      }
      index[i].push_back(j);
      value[i].push_back(new_value);
   }

   void add_to_element(size_t i, size_t j, T increment_value)
   {
      size_t k=0;
      for(; k<index[i].size(); ++k){
         if(index[i][k]==j){
            value[i][k]+=increment_value;
            return;
         }else if(index[i][k]>j){
            insert(index[i], k, j);
            insert(value[i], k, increment_value);
            return;
         }
      }
      index[i].push_back(j);
      value[i].push_back(increment_value);
   }

   // assumes indices is already sorted (handy for cusp::coo conversion)
   void add_sparse_row(size_t i, const std::vector<size_t> &indices, const std::vector<T> &values)
   {
      size_t j=0, k=0;
      while(j<indices.size() && k<index[i].size()){
         if(index[i][k]<indices[j]){
            ++k;
         }else if(index[i][k]>indices[j]){
            insert(index[i], k, indices[j]);
            insert(value[i], k, values[j]);
            ++j;
         }else{
            value[i][k]+=values[j];
            ++j;
            ++k;
         }
      }
      for(;j<indices.size(); ++j){
         index[i].push_back(indices[j]);
         value[i].push_back(values[j]);
      }
   }

   // assumes matrix has symmetric structure - so the indices in row i tell us which columns to delete i from
   void symmetric_remove_row_and_column(size_t i)
   {
      for(size_t a=0; a<index[i].size(); ++a){
         size_t j=index[i][a]; // 
         for(size_t b=0; b<index[j].size(); ++b){
            if(index[j][b]==i){
               erase(index[j], b);
               erase(value[j], b);
               break;
            }
         }
      }
      index[i].resize(0);
      value[i].resize(0);
   }
            
   // Sparse matrix in CUSP's COO format ---------------------------------------
   // (see [1] at the end of file) ---------------------------------------------
   
   // COO rows ---------------------------------------------------------------
   
   // store row indices with non-zero values (nnz) in a 2D vector
   // (each row index is saved as many times as the nnz values)
   std::vector< std::vector<size_t> > coo_format_rows( std::vector<int>& row_indices ) {
      
      ///*
      // compute size of each row (num of nnz elements) 
      // and then the one with max size
      std::vector<size_t> row_size(n);  
      for ( size_t k = 0; k < n; ++k ) {
           
         row_size[k] = index[k].size();  
         //std::cout << "size of row [" << k << "] " << row_size[k] << "\n";
      }
       
      // set the inner dimension of the nested vector coo_format_rows
      // equal to the row with max size
      size_t max_elements_per_row;
      max_elements_per_row = *std::max_element( row_size.begin(), row_size.end() );
      //*/ 
       
      std::vector< std::vector<size_t> > myvec( n, std::vector<size_t>( max_elements_per_row ) ); 
   
      for( size_t i = 0; i < n; ++i ){
         
         std::vector<size_t> row;
         for( size_t j = 0; j < index[i].size(); ++j ){ // condition _j<index[i].size();_ outputs the same
             
            row.push_back(i);
            //std::cout << row[j] << " ";             
         }
          myvec.push_back(row);
          //std::cout << row[i] << "\n"; 
      } return myvec; // myvec cannot be returned by reference (local variable)
   }
   
   // COO columns ------------------------------------------------------------
   
   // store column indices in a 2D vector (each column index corresponds 
   // to the row indices computed in the coo_format_rows)
   std::vector< std::vector<size_t> > coo_format_columns( std::vector<int>& column_indices ) {
        
      ///*  
      // compute size of each row (num of nnz elements) and then the one with max size
      std::vector<size_t> row_size(n);  
      for ( size_t k = 0; k < n; ++k ) {
           
         row_size[k] = index[k].size();  
         //std::cout << "size of row [" << k << "] " << row_size[k] << "\n";
      }
       
      // set the inner dim of the nested vector coo_format_rows equal to the row with max size
      size_t max_elements_per_row;
      max_elements_per_row = *std::max_element( row_size.begin(), row_size.end() );
      //*/
       
      std::vector< std::vector<size_t> > myvec( n, std::vector<size_t>( max_elements_per_row ) );   
      
      for( size_t i = 0; i < n; ++i ){
         
         std::vector<size_t> row;
         for( size_t j = 0; j < index[i].size(); ++j ){
             
            row.push_back( index[i][j] );
            //std::cout << row[j] << " ";             
         }
          myvec.push_back(row);
          //std::cout << row[i] << "\n"; 
      } return myvec; // myvec cannot be returned by reference (local variable)
       
   }
   
   // COO values -------------------------------------------------------------

   std::vector< std::vector<float> > coo_format_values( std::vector<float>& values ) {
         
      ///*  
      // compute size of each row (num of nnz elements) and then the one with max size 
      std::vector<size_t> row_size(n);  
      for ( size_t k = 0; k < n; ++k ) {
           
         row_size[k] = index[k].size();  
         //std::cout << "size of row [" << k << "] " << row_size[k] << "\n";
      }
       
      size_t max_elements_per_row;
      max_elements_per_row = *std::max_element( row_size.begin(), row_size.end() );
      //*/ 
       
      std::vector< std::vector< float > > myvec( n, std::vector<float>( max_elements_per_row ) );
         
      for( size_t i = 0; i < n; ++i ){

         std::vector< float > row; 
         for( size_t j=0; j<value[i].size(); ++j ){
             
            row.push_back( value[i][j] );
            //std::cout << row[j] << " ";             
         }
         myvec.push_back(row);
         //std::cout << row[i] << "\n"; 
      } return myvec; // myvec cannot be returned by reference (local variable) 
    
   }
    
};    
    
typedef SparseMatrix<float> SparseMatrixf;
//typedef SparseMatrix<double> SparseMatrixd;

#endif

   // [1] ----------------------------------------------------------------------------
   //
   // I repeat the exact code for computing max_elements_per_row in coo_format_rows(), 
   // coo_format_columns(), coo_format_values(), which is obviously a bad practice. 
   // when I substitute the repeated instances with the following lines I get errors
   
   // compute size of each row (num of nnz elements) and then the one with max size
   //std::vector<size_t> row_size(n);  
   //for ( size_t k = 0; k < n; ++k ) {
   
   //   row_size[k] = index[k].size();  
   //   //std::cout << "size of row [" << k << "] " << row_size[k] << "\n";
   //}
       
   // set the inner dim of the nested vector coo_format_rows equal to the row with max size
   //static size_t max_elements_per_row = 0;
   //max_elements_per_row = *std::max_element( row_size.begin(), row_size.end() );
   //
   // [1] ----------------------------------------------------------------------------
