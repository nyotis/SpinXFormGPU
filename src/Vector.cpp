// -----------------------------------------------------------------------------
//
// SpinXFormGPU -- Vector.cpp
// 
// Code courtesy of Keenan Crane. 
//
// =============================================================================
// SpinXForm -- Vector.cpp
// Keenan Crane
// August 16, 2011
//

#include "Vector.h"
#include <math.h>

// CONSTRUCTORS ----------------------------------------------------------------

Vector :: Vector( void )
// initializes all components to zero
: x( 0. ),
  y( 0. ),
  z( 0. )
{}

Vector :: Vector( float x0,
                  float y0,
                  float z0 )
// initializes with specified components
: x( x0 ),
  y( y0 ),
  z( z0 )
{}

Vector :: Vector( const Vector& v )
// initializes from existing vector
: x( v.x ),
  y( v.y ),
  z( v.z )
{}


// ACCESSORS -------------------------------------------------------------------

float& Vector :: operator[]( int index )
  // returns reference to the specified component (0-based indexing: x, y, z )
{
   return ( &x )[ index ];
}

const float& Vector :: operator[]( int index ) const
  // returns const reference to the specified component (0-based indexing: x, y, z )
{
   return ( &x )[ index ];
}


// VECTOR SPACE OPERATIONS -----------------------------------------------------

Vector Vector :: operator+( const Vector& v ) const
  // addition
{
   return Vector( x + v.x,
                  y + v.y,
                  z + v.z );
}

Vector Vector :: operator-( const Vector& v ) const
  // subtraction
{
   return Vector( x - v.x,
                  y - v.y,
                  z - v.z );
}

Vector Vector :: operator-( void ) const
  // negation
{
   return Vector( -x,
                  -y,
                  -z );
}

Vector Vector :: operator*( const float& c ) const
  // scalar multiplication
{
   return Vector( x*c,
                  y*c,
                  z*c );
}

Vector operator*( const float& c, const Vector& v )
  // scalar multiplication
{
   return v*c;
}

Vector Vector :: operator/( const float& c ) const
  // scalar division
{
   return (*this) * ( 1./c );
}

void Vector :: operator+=( const Vector& v )
  // addition / assignment
{
   x += v.x;
   y += v.y;
   z += v.z;
}

void Vector :: operator-=( const Vector& v )
  // subtraction / assignment
{
   x -= v.x;
   y -= v.y;
   z -= v.z;
}

void Vector :: operator*=( const float& c )
  // scalar multiplication / assignment
{
   x *= c;
   y *= c;
   z *= c;
}

void Vector :: operator/=( const float& c )
  // scalar division / assignment
{
   (*this) *= ( 1./c );
}


// ALGEBRAIC OPERATIONS --------------------------------------------------------

float Vector :: operator*( const Vector& v ) const
  // dot product
{
   return x*v.x +
          y*v.y +
          z*v.z ;
}

Vector Vector :: operator^( const Vector& v ) const
  // cross product
{
   return Vector( y*v.z - z*v.y,
                  z*v.x - x*v.z,
                  x*v.y - y*v.x );
}


// NORMS -----------------------------------------------------------------------

float Vector :: norm( void ) const
  // returns Euclidean length
{
   return sqrt( norm2());
}

float Vector :: norm2( void ) const
  // returns Euclidean length squared
{
   return (*this) * (*this);
}

Vector Vector :: unit( void ) const
  // returns unit vector
{
   return (*this) / norm();
}

void Vector :: normalize( void )
  // divides by Euclidean length
{
   (*this) /= norm();
}

std::ostream& operator<<( std::ostream& os, const Vector& o )
  // scalar multiplication
{
   os << "[ "
      << o.x << " "
      << o.y << " "
      << o.z
      << " ]";

   return os;
}