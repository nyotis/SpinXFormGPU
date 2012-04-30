// -----------------------------------------------------------------------------
//
// SpinXFormGPU -- Vector.h
// 
// Code courtesy of Keenan Crane. 
//
// =============================================================================
// SpinXForm -- Vector.h
// Keenan Crane
// August 16, 2011
//
// Vector represents a three-dimensional point or vector.
//

#ifndef SPINXFORM_VECTOR_H
#define SPINXFORM_VECTOR_H

#include <ostream>

class Vector
{
   public:
      // CONSTRUCTORS ----------------------------------------------------------
      Vector( void );                        // initializes all components to zero
      Vector( float x, float y, float z); // initializes with specified components
      Vector( const Vector& v );             // initializes from existing vector

      // ACCESSORS -------------------------------------------------------------
            float& operator[] ( int index );       // returns reference to the specified component (0-based indexing: x, y, z )
      const float& operator[] ( int index ) const; // returns const reference to the specified component (0-based indexing: x, y, z )

      // VECTOR SPACE OPERATIONS -----------------------------------------------
      Vector operator+  ( const Vector& v ) const; // addition
      Vector operator-  ( const Vector& v ) const; // subtraction
      Vector operator-  ( void ) const;            // negation
      Vector operator*  ( const float& c ) const; // scalar multiplication
      Vector operator/  ( const float& c ) const; // scalar division
      void   operator+= ( const Vector& v );       // addition / assignment
      void   operator-= ( const Vector& v );       // subtraction / assignment
      void   operator*= ( const float& c );       // scalar multiplication / assignment
      void   operator/= ( const float& c );       // scalar division / assignment

      // ALGEBRAIC OPERATIONS --------------------------------------------------
      float operator*  ( const Vector& v ) const; // dot product
      Vector operator^  ( const Vector& v ) const; // cross product
      
      // NORMS -----------------------------------------------------------------
      float norm( void ) const;  // returns Euclidean length
      float norm2( void ) const; // returns Euclidean length squared
      Vector unit( void ) const;  // returns unit vector
      void normalize( void );     // divides by Euclidean length

      // STORAGE ---------------------------------------------------------------
      float x, y, z; // components
};

// VECTOR SPACE OPERATIONS -----------------------------------------------
Vector operator*( const float& c, const Vector& v ); // scalar multiplication

// I/O -------------------------------------------------------------------------
std::ostream& operator<<( std::ostream& os, const Vector& o ); // prints components

#endif

