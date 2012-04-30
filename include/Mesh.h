// -----------------------------------------------------------------------------
//
// SpinXFormGPU -- Mesh.h
// 
// Code courtesy of Keenan Crane. 
//
// ============================================================================
// SpinXForm -- Mesh.h
// Keenan Crane
// August 16, 2011
//
// Mesh is a very simple mesh data structure consisting of a list of vertices
// and a collection of triangles specified via indices into the vertex list.
// It also stores the data necessary to compute a spin transformation of the
// surface.
//
// To specify a deformation, the user should set a value of "rho" on each
// face corresponding to the desired change in curvature.  The deformation
// is computed by calling updateDeformation(), which puts the transformed
// vertices in the list "newVertices."
//

#ifndef SPINXFORM_MESH_H
#define SPINXFORM_MESH_H

#include <vector>
#include <string>
#include "Quaternion.h"
#include "QuaternionMatrix.h"
#include "Image.h"

using namespace std;

class Face
{
   public:
      int vertex[3]; // indices into vertex list
      Vector uv[3]; // texture coordinates (for visualization only)
};

class Mesh
{
   public:
      void read( const string& filename );
      // loads a triangle mesh in Wavefront OBJ format

      void write( const string& filename );
      // saves a triangle mesh in Wavefront OBJ format

      void setCurvatureChange( const Image& image, const float scale );
      // sets rho values by interpreting "image" as a square image
      // in the range [0,1] x [0,1] and mapping values to the
      // surface via vertex texture coordinates -- grayscale
      // values in the  range [0,1] get mapped (linearly) to values
      // in the range [-scale,scale]

      void updateDeformation( void );
      // computes a conformal deformation using the current rho

      void resetDeformation( void );
      // restores surface to its original configuration

      float area( int i );
      // returns area of triangle i in the original mesh

      vector<Face> faces;
      // list of triangles as indices into vertex list

      vector<Quaternion> vertices, newVertices;
      // original and deformed vertex coordinates

      vector<float> rho;
      // controls change in curvature (one value per face)

   protected:

      vector<Quaternion> lambda;
      // local similarity transformation (one value per vertex)

      vector<Quaternion> omega;
      // divergence of target edge vectors

      QuaternionMatrix L; // Laplace matrix
      QuaternionMatrix E; // matrix for eigenvalue problem

      void buildEigenvalueProblem( void );
      void buildPoissonProblem( void );
      void buildLaplacian( void );
      void buildOmega( void );
      void normalizeSolution( void );
};

#endif
