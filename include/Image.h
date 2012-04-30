// -----------------------------------------------------------------------------
//
// SpinXFormGPU -- Image.h
// 
// Code courtesy of Keenan Crane.
//
// For a C implementation with the 'ugly' fprintf instead of streams 
//    http://paulbourke.net/dataformats/tga/
// In-depth info about the TGA format 
//    http://people.sc.fsu.edu/~jburkardt/pdf/targa.pdf
//    http://people.sc.fsu.edu/~jburkardt/cpp_src/tga_io/tga_io.html
//
// =============================================================================
// SpinXForm -- Image.h
// Keenan Crane
// August 16, 2011
//
// Image represents a grayscale bitmapped image.  Standard usage might
// look something like
//
//    Image im;
//    im.load( "image.tga" );
//
//    // sample image at point p
//    float p[2] = { .5, 1.23 };
//    float value = im.sample( p[0], p[1] );
//

#ifndef SPINXFORM_IMAGE_H
#define SPINXFORM_IMAGE_H

#include <vector>
#include <string>

using namespace std;

class Image
{
   public:
            float& operator()( int x, int y );
      const float& operator()( int x, int y ) const;
      // accesses float (x,y)

      float sample( float x, float y ) const;
      // samples image at (x,y) using bilinear filtering

      int  width( void ) const;
      int height( void ) const;
      // returns image dimensions

      void read( const char* filename );
      // loads an image file in Truevision TGA format
      // (must be uncompressed RGB image with 24 or 32 bits per float)

      void reload( void );
      // updates image from disk

   protected:
      void clamp( int& x, int& y ) const;
      // clamps coordinates to range [0,w-1] x [0,h-1]

      string filename; // name of source file
      vector<float> pixels; // interleaved RGBA float data in range [0-1]
      int w; // width
      int h; // height
};

#endif
