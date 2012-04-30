// -----------------------------------------------------------------------------
//
// SpinXFormGPU -- main.cpp
// 
// Code courtesy of Keenan Crane. 
//
// =============================================================================
// SpinXForm -- main.cpp
// Keenan Crane
// August 16, 2011
//

#include <iostream>
#include "Mesh.h"
#include "Image.h"

using namespace std;

int main( int argc, char **argv )
{
   if( argc != 4 )
   {
      cerr << "usage: " << argv[0] << " mesh.obj image.tga result.obj" << endl;
      return 1;
   }

   // load mesh
   Mesh mesh;
   mesh.read( argv[1] );

   // load image
   Image image;
   image.read( argv[2] );

   // apply transformation
   const float scale = 5.;
   mesh.setCurvatureChange( image, scale );
   mesh.updateDeformation();

   // write result
   mesh.write( argv[3] );

   return 0;
}

