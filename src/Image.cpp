// -----------------------------------------------------------------------------
//
// SpinXFormGPU -- EigenSolver.cpp
// 
// Code courtesy of Keenan Crane. 
//
// =============================================================================
// SpinXForm -- Image.cpp
// Keenan Crane
// August 16, 2011
//

#include <fstream>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include "Image.h"
#include "Utility.h"

float& Image :: operator()( int x, int y )
// accesses pixel (x,y)
{
   return pixels[ x + y*w ];
}

const float& Image :: operator()( int x, int y ) const
// accesses pixel (x,y)
{
   return pixels[ x + y*w ];
}

float Image :: sample( float x, float y ) const
// samples image at (x,y) using bilinear filtering
{
   const Image& I( *this );
   float ax = x - floor( x );
   float ay = y - floor( y );
   float bx = 1. - ax;
   float by = 1. - ay;
   int x0 = (int) floor( x );
   int y0 = (int) floor( y );
   int x1 = x0 + 1;
   int y1 = y0 + 1;

   clamp( x0, y0 );
   clamp( x1, y1 );

   return by * ( bx * I(x0,y0) + ax * I(x1,y0) ) +
          ay * ( bx * I(x0,y1) + ax * I(x1,y1) ) ;
}

int Image :: width( void ) const
// returns image width
{
   return w;
}

int Image :: height( void ) const
// returns image height
{
   return h;
}

class TGAHeader
// header format for Truevision TGA images
{
   public:
      char  idFieldSize;
      char  colorMapType;
      char  dataTypeCode;
      short colorMapOrigin;
      short colorMapLength;
      char  colorMapEntrySize;
      short xOrigin;
      short yOrigin;
      short width;
      short height;
      char  bitsPerPixel;
      char  imageSpecification;
};

void Image :: read( const char* _filename )
// loads an image file in Truevision TGA format
// (must be uncompressed RGB image with 24 or 32 bits per pixel)
{
   filename = string( _filename );
   ifstream in( filename.c_str(), ios_base::binary );

   if( !in.is_open() )
   {
      cerr << "Error: could not open file " << filename << " for input!" << endl;
      exit( 1 );
   }

   // read header
   TGAHeader header;
   in.read( (char*) &(header.idFieldSize),        1 );
   in.read( (char*) &(header.colorMapType),       1 );
   in.read( (char*) &(header.dataTypeCode),       1 );
   in.read( (char*) &(header.colorMapOrigin),     2 );
   in.read( (char*) &(header.colorMapLength),     2 );
   in.read( (char*) &(header.colorMapEntrySize),  1 );
   in.read( (char*) &(header.xOrigin),            2 );
   in.read( (char*) &(header.yOrigin),            2 );
   in.read( (char*) &(header.width),              2 );
   in.read( (char*) &(header.height),             2 );
   in.read( (char*) &(header.bitsPerPixel),       1 );
   in.read( (char*) &(header.imageSpecification), 1 );
   if( bigEndian() )
   {
      swapShort( header.colorMapOrigin );
      swapShort( header.colorMapLength );
      swapShort( header.xOrigin );
      swapShort( header.yOrigin );
      swapShort( header.width );
      swapShort( header.height );
   }

   w = header.width;
   h = header.height;

   // validate data type
   const char uncompressedGrayscale = 3;
   if( header.dataTypeCode != uncompressedGrayscale ||
       header.bitsPerPixel != 8 )
   {
      cerr << "Error: input must be uncompressed grayscale image with 8 bits per pixel." << endl;
      exit( 1 );
   }

   // read identification field (unused)
   vector<char> idField( header.idFieldSize );
   in.read( &idField[0], header.idFieldSize );

   // read color map data (unused)
   if( header.colorMapType == 1 )
   {
      int bytesPerEntry = header.colorMapEntrySize / 8;
      int colorMapSize = header.colorMapLength * bytesPerEntry;
      vector<char> colorMapData( colorMapSize );
      in.read( &colorMapData[0], colorMapSize );
   }

   // read pixel data
   vector<unsigned char> pixelData( w*h );
   in.read( (char*) &pixelData[0], w*h );

   // convert pixel data
   pixels.resize( w*h );
   for( int i = 0; i != w*h; ++i )
   {
      pixels[i] = (float) pixelData[i] / 255.;
   }
}

void Image :: reload( void )
// updates image from disk
{
   read( filename.c_str() );
}

void Image :: clamp( int& x, int& y ) const
// clamps coordinates to range [0,w-1] x [0,h-1]
{
   x = max( 0, min( w-1, x ));
   y = max( 0, min( h-1, y ));
}

