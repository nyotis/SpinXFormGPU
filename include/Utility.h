// -----------------------------------------------------------------------------
//
// SpinXFormGPU -- Utility.h
// 
// Code courtesy of Keenan Crane. 
//
// =============================================================================
// SpinXForm -- Utility.h
// Keenan Crane
// August 16, 2011
//
// This file contains simple convenience functions used by multiple objects.
//

#ifndef SPINXFORM_UTILITY_H
#define SPINXFORM_UTILITY_H

#include <vector>

inline float sqr( float x )
{
   return x*x;
}

template <class T>
inline void removeMean( std::vector<T>& v )
{
   T mean = 0.;

   for( size_t i = 0; i < v.size(); i++ )
   {
      mean += v[i];
   }

   mean /= (float) v.size();

   for( size_t i = 0; i < v.size(); i++ )
   {
      v[i] -= mean;
   }
}

inline bool bigEndian( void )
{
   int n = 1;
   char *c = (char *) &n;

   if( c[0] == 1 )
   {
      return false;
   }
   return true;
}

inline void swapShort( short& x )
{
   x = (( x & 0x00FF ) << 8 ) +
       (( x & 0xFF00 ) >> 8 ) ;
}

#endif
