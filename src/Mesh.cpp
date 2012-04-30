// -----------------------------------------------------------------------------
//
// SpinXFormGPU -- Mesh.cpp
//
// Code courtesy of Keenan Crane. 
//
// ============================================================================
// SpinXForm -- Mesh.cpp
// Keenan Crane
// August 16, 2011
//

#include <fstream>
#include <sstream>
#include <cmath>
#include <ctime>
#include "Mesh.h"
#include "LinearSolver.h"
#include "EigenSolver.h"
#include "Utility.h"

#include <iostream>

void Mesh :: updateDeformation( void )
{
   int t0 = clock();

   // solve eigenvalue problem for local similarity transformation lambda
   buildEigenvalueProblem();
   EigenSolver::solve( E, lambda ); // E(4002 x 4002)

   // solve Poisson problem for new vertex positions
  buildPoissonProblem();
  LinearSolver::solve( L, newVertices, omega );
  normalizeSolution();

   int t1 = clock();
   cout << "time: " << (t1-t0)/(float) CLOCKS_PER_SEC << "s" << endl;
}

void Mesh :: resetDeformation( void )
{
   // copy original mesh vertices to current mesh
   for( size_t i = 0; i < vertices.size(); i++ )
   {
      newVertices[i] = vertices[i];
   }

   normalizeSolution();
}

void Mesh :: setCurvatureChange( const Image& image, const float scale )
// sets rho values by interpreting "image" as a square image
// in the range [0,1] x [0,1] and mapping values to the
// surface via vertex texture coordinates -- grayscale
// values in the  range [0,1] get mapped (linearly) to values
// in the range [-scale,scale]
{
   float w = (float) image.width();
   float h = (float) image.height();

   for( size_t i = 0; i < faces.size(); i++ )
   {
      // get handle to current value of rho
      rho[i] = 0.;

      // compute average value over the face
      for( int j = 0; j < 3; j++ )
      {
         Vector uv = faces[i].uv[j];
         rho[i] += image.sample( uv.x*w, uv.y*h ) / 3.;
      }

      // map value to range [-scale,scale]
      rho[i] = (2.*(rho[i]-.5)) * scale;
   }
}

float Mesh :: area( int i )
// returns area of triangle i in the original mesh
{
   Vector& p1 = vertices[ faces[i].vertex[0] ].im();
   Vector& p2 = vertices[ faces[i].vertex[1] ].im();
   Vector& p3 = vertices[ faces[i].vertex[2] ].im();

   return .5 * (( p2-p1 ) ^ ( p3-p1 )).norm();
}

void Mesh :: buildEigenvalueProblem( void )
{
   // allocate a sparse |V|x|V| matrix
   int nV = vertices.size();
   E.resize( nV, nV );

   // visit each face
   for( size_t k = 0; k < faces.size(); k++ )
   {
      float A = area(k);
      float a = -1. / (4.*A);
      float b = rho[k] / 6.;
      float c = A*rho[k]*rho[k] / 9.;

      // get vertex indices
      int I[3] =
      {
         faces[k].vertex[0],
         faces[k].vertex[1],
         faces[k].vertex[2]
      };

      // compute edges across from each vertex
      Quaternion e[3];
      for( int i = 0; i < 3; i++ )
      {
         e[i] = vertices[ I[ (i+2) % 3 ]] -
                vertices[ I[ (i+1) % 3 ]] ;
      }

      // increment matrix entry for each ordered pair of vertices
      for( int i = 0; i < 3; i++ )
      for( int j = 0; j < 3; j++ )
      {
         E(I[i],I[j]) += a*e[i]*e[j] + b*(e[j]-e[i]) + c;
      }
   }
  // std::cout << "first dim of Quatern matrix: " << E.size(1) << "second dim of Quatern matrix: " << E.size(2) << "\n";
}

void Mesh :: buildPoissonProblem( void )
{
   buildLaplacian();
   buildOmega();
}

void Mesh :: buildLaplacian( void )
// builds the cotan-Laplace operator
{
   // allocate a sparse |V|x|V| matrix
   int nV = vertices.size();
   L.resize( nV, nV );

   // visit each face
   for( size_t i = 0; i < faces.size(); i++ )
   {
      // visit each triangle corner
      for( int j = 0; j < 3; j++ )
      {
         // get vertex indices
         int k0 = faces[i].vertex[ (j+0) % 3 ];
         int k1 = faces[i].vertex[ (j+1) % 3 ];
         int k2 = faces[i].vertex[ (j+2) % 3 ];

         // get vertex positions
         Vector f0 = vertices[k0].im();
         Vector f1 = vertices[k1].im();
         Vector f2 = vertices[k2].im();

         // compute cotangent of the angle at the current vertex
         // (equal to cosine over sine, which equals the dot
         // product over the norm of the cross product)
         Vector u1 = f1 - f0;
         Vector u2 = f2 - f0;
         float cotAlpha = (u1*u2)/(u1^u2).norm();

         // add contribution of this cotangent to the matrix
         L( k1, k2 ) -= cotAlpha / 2.;
         L( k2, k1 ) -= cotAlpha / 2.;
         L( k1, k1 ) += cotAlpha / 2.;
         L( k2, k2 ) += cotAlpha / 2.;
      }
   }
}

void Mesh :: buildOmega( void )
{
   // clear omega
   for( size_t i = 0; i < omega.size(); i++ )
   {
      omega[i] = 0.;
   }

   // visit each face
   for( size_t i = 0; i < faces.size(); i++ )
   {
      // get indices of the vertices of this face
      int v[3] = { faces[i].vertex[0],
                   faces[i].vertex[1],
                   faces[i].vertex[2] };

      // visit each edge
      for( int j = 0; j < 3; j++ )
      {
         // get vertices
         Quaternion f0 = vertices[ v[ (j+0) % 3 ]];
         Quaternion f1 = vertices[ v[ (j+1) % 3 ]];
         Quaternion f2 = vertices[ v[ (j+2) % 3 ]];

         // determine orientation of this edge
         int a = v[ (j+1) % 3 ];
         int b = v[ (j+2) % 3 ];
         if( a > b )
         {
            swap( a, b );
         }

         // compute transformed edge vector
         Quaternion lambda1 = lambda[a];
         Quaternion lambda2 = lambda[b];
         Quaternion e = vertices[b] - vertices[a];
         Quaternion eTilde = (1./3.) * (~lambda1) * e * lambda1 +
                             (1./6.) * (~lambda1) * e * lambda2 +
                             (1./6.) * (~lambda2) * e * lambda1 +
                             (1./3.) * (~lambda2) * e * lambda2 ;

         // compute cotangent of the angle opposite the current edge
         Vector u1 = ( f1 - f0 ).im();
         Vector u2 = ( f2 - f0 ).im();
         float cotAlpha = (u1*u2)/(u1^u2).norm();

         // add contribution of this edge to the divergence at its vertices
         omega[a] -= cotAlpha * eTilde / 2.;
         omega[b] += cotAlpha * eTilde / 2.;
      }
   }

   removeMean( omega );
}

void Mesh :: normalizeSolution( void )
{
   // center vertices around the origin
   removeMean( newVertices );

   // find the vertex with the largest norm
   float r = 0.;
   for( size_t i = 0; i < vertices.size(); i++ )
   {
      r = max( r, newVertices[i].norm2() );
   }
   r = sqrt(r);

   // rescale so that vertices have norm at most one
   for( size_t i = 0; i < vertices.size(); i++ )
   {
      newVertices[i] /= r;
   }
}

// FILE I/O --------------------------------------------------------------------

void Mesh :: read( const string& filename )
// loads a triangle mesh in Wavefront OBJ format
{
   // open mesh file
   ifstream in( filename.c_str() );
   if( !in.is_open() )
   {
      cerr << "Error: couldn't open file ";
      cerr << filename;
      cerr << " for input!" << endl;
      exit( 1 );
   }

   // temporary list of vertex coordinates
   vector<Vector> uv;

   // parse mesh file
   string s;
   while( getline( in, s ))
   {
      stringstream line( s );
      string token;

      line >> token;

      if( token == "v" ) // vertex
      {
         float x, y, z;

         line >> x >> y >> z;

         vertices.push_back( Quaternion( 0., x, y, z ));
         newVertices.push_back( Quaternion( 0., x, y, z ));
      }
      if( token == "vt" ) // texture coordinate
      {
         float u, v;

         line >> u >> v;

         uv.push_back( Vector( u, v, 0. ));
      }
      else if( token == "f" ) // face
      {
         Face triangle;

         // iterate over vertices
         for( int i = 0; i < 3; i++ )
         {
            line >> s;
            stringstream item( s );

            int I[3] = { -1, -1, -1 };

            // iterate over v, vt, and vn indices
            for( int j = 0; getline( item, s, '/' ) && j < 3; j++ )
            {
               stringstream index( s );
               index >> I[j];
            }

            triangle.vertex[i] = I[0]-1;

            if( I[1] != -1 )
            {
               triangle.uv[i] = uv[ I[1]-1 ];
            }
         }

         faces.push_back( triangle );
      }
   }

   // allocate space for mesh attributes
   lambda.resize( vertices.size() );
   omega.resize( vertices.size() );
   rho.resize( faces.size() );
   normalizeSolution();
}

void Mesh :: write( const string& filename )
// saves a triangle mesh in Wavefront OBJ format
{
   ofstream out( filename.c_str() );

   if( !out.is_open() )
   {
      cerr << "Error: couldn't open file ";
      cerr << filename;
      cerr << " for output!" << endl;
      return;
   }

   for( size_t i = 0; i < vertices.size(); i++ )
   {
      out << "v " << newVertices[i].im().x << " "
                  << newVertices[i].im().y << " "
                  << newVertices[i].im().z << endl;
   }

   for( size_t i = 0; i < faces.size(); i++ )
   {
      out << "f " << 1+faces[i].vertex[0] << " "
                  << 1+faces[i].vertex[1] << " "
                  << 1+faces[i].vertex[2] << endl;
   }
}
