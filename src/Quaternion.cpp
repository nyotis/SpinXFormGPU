// -----------------------------------------------------------------------------
//
// SpinXFormGPU -- Quaternion.cpp
// 
// Code courtesy of Keenan Crane. 
//
// =============================================================================
// SpinXForm -- Quaternion.cpp
// Keenan Crane
// August 16, 2011
//

#include "Quaternion.h"
#include <cmath>
#include <iostream>

using namespace std;

// CONSTRUCTORS ----------------------------------------------------------

Quaternion :: Quaternion( void )
// initializes all components to zero
: s( 0. ),
  v( 0., 0., 0. )
{}

Quaternion :: Quaternion( const Quaternion& q )
// initializes from existing quaternion
: s( q.s ),
  v( q.v )
{}

Quaternion :: Quaternion( float s_, float vi, float vj, float vk )
// initializes with specified float (s) and imaginary (v) components
: s( s_ ),
  v( vi, vj, vk )
{}

Quaternion :: Quaternion( float s_, const Vector& v_ )
// initializes with specified float(s) and imaginary (v) components
: s( s_ ),
  v( v_ )
{}

Quaternion :: Quaternion( float s_ )
: s( s_ )
{}

Quaternion :: Quaternion( const Vector& v_ )
: v( v_ )
{}


// ASSIGNMENT OPERATORS --------------------------------------------------

const Quaternion& Quaternion :: operator=( float _s )
// assigns a purely real quaternion with real value s
{
   s = _s;
   v = Vector( 0., 0., 0. );

   return *this;
}

const Quaternion& Quaternion :: operator=( const Vector& _v )
// assigns a purely real quaternion with imaginary value v
{
   s = 0.;
   v = _v;

   return *this;
}


// ACCESSORS -------------------------------------------------------------

float& Quaternion::operator[]( int index )
// returns reference to the specified component (0-based indexing: float, i, j, k)
{
   return ( &s )[ index ];
}

const float& Quaternion::operator[]( int index ) const
// returns const reference to the specified component (0-based indexing: float, i, j, k)
{
   return ( &s )[ index ];
}

void Quaternion::toMatrix( float Q[4][4] ) const
// returns 4x4 matrix representation
{
   Q[0][0] =   s; Q[0][1] = -v.x; Q[0][2] = -v.y; Q[0][3] = -v.z;
   Q[1][0] = v.x; Q[1][1] =    s; Q[1][2] = -v.z; Q[1][3] =  v.y;
   Q[2][0] = v.y; Q[2][1] =  v.z; Q[2][2] =    s; Q[2][3] = -v.x;
   Q[3][0] = v.z; Q[3][1] = -v.y; Q[3][2] =  v.x; Q[3][3] =    s;
}

float& Quaternion::re( void )
// returns reference to float part
{
   return s;
}

const float& Quaternion::re( void ) const
// returns const reference to float part
{
   return s;
}

Vector& Quaternion::im( void )
// returns reference to imaginary part
{
   return v;
}

const Vector& Quaternion::im( void ) const
// returns const reference to imaginary part
{
   return v;
}


// VECTOR SPACE OPERATIONS -----------------------------------------------

Quaternion Quaternion::operator+( const Quaternion& q ) const
// addition
{
   return Quaternion( s+q.s, v+q.v );
}

Quaternion Quaternion::operator-( const Quaternion& q ) const
// subtraction
{
   return Quaternion( s-q.s, v-q.v );
}

Quaternion Quaternion::operator-( void ) const
// negation
{
   return Quaternion( -s, -v );
}

Quaternion Quaternion::operator*( float c ) const
// scalar multiplication
{
   return Quaternion( s*c, v*c );
}

Quaternion operator*( float c, const Quaternion& q )
// scalar multiplication
{
   return q*c;
}

Quaternion Quaternion::operator/( float c ) const
// scalar division
{
   return Quaternion( s/c, v/c );
}

void Quaternion::operator+=( const Quaternion& q )
// addition / assignment
{
   s += q.s;
   v += q.v;
}

void Quaternion::operator+=( float c )
// addition / assignment of pure real
{
   s += c;
}

void Quaternion::operator-=( const Quaternion& q )
// subtraction / assignment
{
   s -= q.s;
   v -= q.v;
}

void Quaternion::operator-=( float c )
// subtraction / assignment of pure real
{
   s -= c;
}

void Quaternion::operator*=( float c )
// scalar multiplication / assignment
{
   s *= c;
   v *= c;
}

void Quaternion::operator/=( float c )
// scalar division / assignment
{
   s /= c;
   v /= c;
}


// ALGEBRAIC OPERATIONS --------------------------------------------------

Quaternion Quaternion::operator*( const Quaternion& q ) const
// Hamilton product
{
   const float& s1( s );
   const float& s2( q.s );
   const Vector& v1( v );
   const Vector& v2( q.v );

   return Quaternion( s1*s2 - v1*v2, s1*v2 + s2*v1 + (v1^v2) );
}

void Quaternion::operator*=( const Quaternion& q )
// Hamilton product / assignment
{
   *this = ( *this * q );
}

Quaternion Quaternion::operator~( void ) const
// conjugation
{
   return Quaternion( s, -v );
}

Quaternion Quaternion::inv( void ) const
{
   return ( ~( *this )) / this->norm2();
}


// NORMS -----------------------------------------------------------------

float Quaternion::norm( void ) const
// returns Euclidean length
{
   return sqrt( s*s + v.x*v.x + v.y*v.y + v.z*v.z );
}

float Quaternion::norm2( void ) const
// returns Euclidean length squared
{
   return s*s + v*v;
}

Quaternion Quaternion::unit( void ) const
// returns unit quaternion
{
   return *this / norm();
}

void Quaternion::normalize( void )
// divides by Euclidean length
{
   *this /= norm();
}


// GEOMETRIC OPERATIONS --------------------------------------------------

Quaternion slerp( const Quaternion& q0, const Quaternion& q1, float t )
// spherical-linear interpolation
{
   // interpolate length
   float m0 = q0.norm();
   float m1 = q1.norm();
   float m = (1-t)*m0 + t*m1;

   // interpolate direction
   Quaternion p0 = q0 / m0;
   Quaternion p1 = q1 / m1;
   float theta = acos(( (~p0)*p1 ).re() );
   Quaternion p = ( sin((1-t)*theta)*p0 + sin(t*theta)*p1 )/sin(theta);

   return m*p;
}


// I/O -------------------------------------------------------------------------

std::ostream& operator<<( std::ostream& os, const Quaternion& q )
// prints components
{
   os << "( " << q.re() << ", " << q.im() << " )";

   return os;
}

