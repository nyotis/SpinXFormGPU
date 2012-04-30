// -----------------------------------------------------------------------------
//
// SpinXFormGPU -- Quaternion.h
// 
// Code courtesy of Keenan Crane. 
//
// =============================================================================
// SpinXForm -- Quaternion.h
// Keenan Crane
// August 16, 2011
//
// Quaternion represents an element of the quaternions.  Operator overloading
// makes it possible to write expressions involving quaternions, vectors, and
// scalars.  For instance,
//
//    Quaternion q;
//    q = 1.23;
//
// sets the real part of q to 1.23 and the imaginary part to zero.  Similarly,
//
//    Quaternion q;
//    Vector v;
//    q = v;
//
// sets the imaginary part to v and the real part to zero.
//

#ifndef SPINXFORM_QUATERNION_H
#define SPINXFORM_QUATERNION_H

#include "Vector.h"
#include <ostream>

class Quaternion
{
   public:
      // CONSTRUCTORS ----------------------------------------------------------
      Quaternion( void );                                      // initializes all components to zero
      Quaternion( const Quaternion& q );                       // initializes from existing quaternion
      Quaternion( float s, float vi, float vj, float vk ); // initializes with specified real (s) and imaginary (v) components
      Quaternion( float s, const Vector& v );                 // initializes with specified real (s) and imaginary (v) components
      Quaternion( float s );                                  // initializes purely real quaternion with specified real (s) component
      Quaternion( const Vector& v );                           // initializes purely imaginary quaternion with specified imaginary (v) component
      
      // ASSIGNMENT OPERATORS --------------------------------------------------
      const Quaternion& operator=( float s ); // assigns a purely real quaternion with real value s
      const Quaternion& operator=( const Vector& v ); // assigns a purely real quaternion with imaginary value v

      // ACCESSORS -------------------------------------------------------------
            float& operator[]( int index );       // returns reference to the specified component (0-based indexing: r, i, j, k)
      const float& operator[]( int index ) const; // returns const reference to the specified component (0-based indexing: r, i, j, k)
      void toMatrix( float Q[4][4] ) const;       // builds 4x4 matrix Q representing (left) quaternion multiplication
            float& re( void );                    // returns reference to float part
      const float& re( void ) const;              // returns const reference to float part
            Vector& im( void );                    // returns reference to imaginary part
      const Vector& im( void ) const;              // returns const reference to imaginary part

      // VECTOR SPACE OPERATIONS -----------------------------------------------
      Quaternion operator+( const Quaternion& q ) const; // addition
      Quaternion operator-( const Quaternion& q ) const; // subtraction
      Quaternion operator-( void ) const;                // negation
      Quaternion operator*( float c ) const;            // scalar multiplication
      Quaternion operator/( float c ) const;            // scalar division
      void       operator+=( const Quaternion& q );      // addition / assignment
      void       operator+=( float c );                 // addition / assignment of pure real
      void       operator-=( const Quaternion& q );      // subtraction / assignment
      void       operator-=( float c );                 // subtraction / assignment of pure real
      void       operator*=( float c );                 // scalar multiplication / assignment
      void       operator/=( float c );                 // scalar division / assignment

      // ALGEBRAIC OPERATIONS --------------------------------------------------
      Quaternion operator*( const Quaternion& q ) const; // Hamilton product
      void       operator*=( const Quaternion& q );      // Hamilton product / assignment
      Quaternion operator~( void ) const;                // conjugation
      Quaternion inv( void ) const;                      // inverse
      
      // NORMS -----------------------------------------------------------------
      float norm( void ) const;     // returns Euclidean length
      float norm2( void ) const;    // returns Euclidean length squared
      Quaternion unit( void ) const; // returns unit quaternion
      void normalize( void );        // divides by Euclidean length

   protected:
      // STORAGE ---------------------------------------------------------------
      float s; // scalar (float) part
      Vector v; // vector (imaginary) part
};

// VECTOR SPACE OPERATIONS -----------------------------------------------
Quaternion operator*( float c, const Quaternion& q ); // scalar multiplication

// I/O -------------------------------------------------------------------------
std::ostream& operator<<( std::ostream& os, const Quaternion& q ); // prints components

#endif

