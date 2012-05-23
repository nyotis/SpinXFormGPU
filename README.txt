-------------------------------------------------------------------------------

                   SpinXFormGPU, Nikos Yiotis, May 2012

SpinXFormGPU is a GPU-aware version of SpinXForm [1] that uses NVidia's CUDA. 
SpinXForm [1] is a C++ implementation of an algorithm for computing conformal 
transformations of polygon meshes based on the paper

  K. Crane, U. Pinkall, and P. Schroeder, "Spin Transformations of
  Discrete Surfaces," ACM Transactions on Graphics (SIGGRAPH) 2011.

The basic functionality is to take a manifold, triangulated mesh (without 
boundary) plus a grayscale image map as input, and produce a new mesh with the 
same connectivity but modified vertex positions. Triangles in the new mesh will 
have nearly the same angles and aspect ratios as those from the original mesh, 
but with curvature that differs by the amount specified in the input image  
(see the paper above for more details).  

SpinXFormGPU solves the linear systems of SpinXForm [1] on the GPU using CUSP's 
[2] conjugate gradient solver but *without* any preconditioner. Depending on 
which distribution of SpinXFormGPU you have installed, SpinXFormGPU might use 
CUSP's views (see ./views.txt for an explanation) or not. In order to solve the 
linear systems using CUSP, I converted the sparse matrices of SpinXForm to the 
Coordinate (or triplet) format (COO) on the host (see COO_format.txt). Note that 
SpinXFormGPU and SpinXFormGPUViews only support single precision arithmetic 
since I developed SpinXFormGPU for a NVidia's card of Compute Compatibility 1.1,
which doesn't support double precision. I provide double precision support for 
CUDA devices of Compute Capability 1.3 and higher (completely untested, but 
likely to be much faster compared to the single precision), see [7], [8] for the 
equivalent to the SpinXFormGPU, SpinXFormGPU_Views versions respectively. 

Note that the code was written in ISO C++, however some strict ISO restrictions 
(in particular, compilation with ansi, pedantic flags) were relaxed due to a bug 
in Thrust, see the comments in the makefile for details. Apart from:

- CUDA v.4.1 or higher,
- Thrust v.1.5.1 (a parallel algorithms library resembling the C++ Standard 
  Template Library (STL), see [2]), or higher
- CUSP v.0.3.0 (a library for sparse linear algebra and graph computations on 
  CUDA, built on top of Thrust, see [3]) or higher.

there are no other external dependencies so the code should compile on any 
machine. The SpinXFormGPU executable takes two input arguments (i.e., one 
triangulated mesh and one grayscale image) and one optional output argument:

   spinxformgpu mesh.obj image.tga [result.obj]

I provide sphere.obj/sphere_original.obj and bumpy.tga files (courtesy of 
Keenan Crane) so one can compile and execute in command line with 

make clean
make
./spinxformgpu sphere.obj bumpy.tga result.obj

or 

make clean
make
./spinxformgpu sphere_original.obj bumpy.tga result.obj

Note that there is an issue in QuaternionMatrix.cpp. In the 

'cusp_columns.erase( cusp_columns.begin(), cusp_columns.begin() + 84648  );'

function - found in the 'CUSP's COO columns' part of the routine - the value 
84648, which represents 4-times the number of vertices of the input mesh 
(sphere_original.obj in this case), should be hard-coded. Note however, that the 

'cusp_rows.erase( cusp_rows.begin(), cusp_rows.begin() + n*4 );'      // generic 
and
'cusp_values.erase( cusp_values.begin(), cusp_values.begin() + n*4 );'// generic

functions found respectively in the 'CUSP's COO rows' and 'CUSP's COO values' 
parts of the routine work fine with generic values (n*4), that is there is no 
need to hard-code the value.


* Notes on Performance

For meshes that consist of < 5K vertices (like sphere.obj) timings on CPU/GPU 
(SpinXForm/SpinXFormGPU respectively) are approximately the same. For meshes 
that consist of ~20K vertices (like sphere_original.obj), SpinXFormGPU is 2-to-3 
times faster compared to SpinXForm [1] .

Another concern is that SpinXFormGPU was developed on a Mac OS X 10.6.8 (Snow 
Leopard) 2.66GHz Intel Core Duo, gcc version 4.2.1 with an NVidia GeForce 9400M 
card of Compute Compatibility 1.1. Since double precision is not supported for 
Compute Compatibility 1.1, I had to convert doubles to floats, which resulted in 
a severe loss in performance. Single precision accuracy also affects Conjugate 
Gradient's convergence on the GPU. For CUDA devices of Compute Capability 1.3 
and higher I provide double precision support for the SpinXFormGPU code at [7], 
[8] respectively.

Note that I didn't try to further accelerate in CUDA the "fast" version of 
SpinXForm (see [9]), which uses CHOLMOD [5], since sparse matrix factorization 
algorithms on the GPU are far from being mature and CHOLMOD is significantly 
more difficult to build/link on diverse platforms.

Performance of the software is memory-bound: its performance is limited by 
memory bandwidth, not floating point operations (the amount of computations 
per load/store is relatively small). This is due to the sparse-matrix-vector 
multiplications that occur when the conjugate solver is executing.

The iterative nature of the solution ('EigenSolver::solve' in EigenSolver.cpp 
applies three inverse power iterations) pushes data to the device and back to 
the host at the start and end of each iteration, which is an expensive operation. 
This makes it hard to optimize further (ideally transfer of data between the host 
and the device through the PCIe bus should be minimized. Data should reside on 
the GPU for as long as possible and ideally one should perform all the expensive 
computations on the device and only by the end of the computations should transfer 
the data back to the host). In particular, at the end of each iteration once the 
solver ('solve_on_device_views', which is called from the routine 
'LinearSolver::solve') on the device finishes, normalization of the vector 
and conversion from reals to quaternions takes place back to the host. Note that 
one could implement normalization on the device or skip it until completion of 
the three iterations, since the only concern is overflow after a large number of 
iterations (cheers Keenan Crane). Keeping data in real format on the device, 
however, so that we avoid the bottleneck of transferring data to the host is not 
straightforward.

Whether CUSP's views are used or not, SpinXFormGPU is not significantly faster 
compared to SpinXForm [1] mainly for two reasons:

1. In constrast to SpinXForm, no preconditioner is used in SpinXFormGPU. 
Preconditioners supported by CUSP do not work. For the linear systems to be 
solved, matrices of quaternions are converted into a system of reals, so it 
might be just that preconditioners for real matrices don't work for 
quaternionic matrices. 

2. SpinXFormGPU could have been significantly accelerated further by using 
prefactorization of the sparse matrix A in A x = b in the eigenvalue problem. 
The main advantage is that one can just factor once and then iteratively apply 
back-substitution, so essentially each power iteration is a back-substitution.  
This has already implemented for the CPU by Crane et al. in [4] by using 
CHOLMOD [5]. However, implementing sparse Cholesky factorization on the GPU 
tends to be a fairly complicated process, which does not always exploit the 
maximum potential for parallel computations. To the best of my knowledge, 
efficient and intuitive implementations of sparse Cholesky (or other)
factorization algorithms on the GPU do not actually exist.

Here are some pointers to the literature for the adventurous ones who would like
to implement a Cholesky factorization on the GPU (cheers to Vasily Volkov).

[-] Christen et al., 2007 General-Purpose Sparse Matrix Building Blocks using 
    the NVIDIA CUDA Technology Platform, 
    http://www.cs.jhu.edu/~misha/ReadingSeminar/Papers/Christen07.pdf
[-] Krawezik and Poole, 2009, Accelerating the ANSYS Direct Sparse Solver with 
    GPUs, 
    http://saahpc.ncsa.illinois.edu/09/papers/Krawezik_paper.pdf
[-] Yu et al., 2011, A CPU–GPU hybrid approach for the unsymmetric multifrontal 
    method, 
    http://www.sciencedirect.com/science/article/pii/S0167819111001293
[-] George et al., 2011, Multifrontal Factorization of Sparse SPD Matrices on 
    GPUs, http://goo.gl/NSuD2
[-] Lucas et al., 2012, Multifrontal Sparse Matrix Factorization on Graphics 
    Processing Units, ftp://ftp.isi.edu/isi-pubs/tr-677.pdf
  
* An alternative for recording time is to use cudaEvents and cudaEventElapsedTime()

* Note that the reference_solution_sphere.obj is computed from the original 
version of the code (SpinXForm)[1] by substituting double with float.

* ./identical_system.txt is a sanity check to make sure that the gpu version of 
the code (SpinXFormGPU) solves the same linear system with the original version 
of the code (SpinXForm [1]) when sphere.obj and bumpy.tga files are used as input.

[1] http://users.cms.caltech.edu/~keenan/src/spinxform_commandline.zip
    at http://users.cms.caltech.edu/~keenan/project_spinxform.html

[2] http://thrust.github.com/
    
[3] http://code.google.com/p/cusp-library/

[4] http://users.cms.caltech.edu/~keenan/src/spinxform_fast.zip
    at http://users.cms.caltech.edu/~keenan/project_spinxform.html

[5] CHOLMOD: routines for sparse Cholesky factorization
    http://www.cise.ufl.edu/research/sparse/cholmod/

[6] http://code.google.com/p/gpuocelot/

[7] https://github.com/nyiotis/SpinXFormGPU_double

[8] https://github.com/nyiotis/SpinXFormGPUViews_double

[9] http://users.cms.caltech.edu/~keenan/src/spinxform_fast.zip
    at http://users.cms.caltech.edu/~keenan/project_spinxform.html

Code courtesy of Keenan Crane unless otherwise stated. For the original version 
(SpinXForm), see [1]. Thanks to Keenan Crane for making his code publicly 
available and many useful discussions and suggestions. 


===============================================================================

                             SpinXForm v0.1

                              Keenan Crane
                             August 16, 2011

===============================================================================            


--------------
0. CONTENTS
--------------

  0. CONTENTS
  1. ABOUT
  2. DEPENDENCIES
  3. BUILDING
  4. TESTING
  5. USAGE
   A. INPUTS
   B. COMMAND LINE MODE
   C. INTERACTIVE MODE
  6. FILE FORMATS
   A. Wavefront OBJ
   B. Truevision TGA
  7. SOURCE CODE
  8. LICENSE


-----------
1. ABOUT
-----------

This archive contains a C++ implementation of an algorithm for computing
conformal transformations of polygon meshes.  It is based on the paper

  K. Crane, U. Pinkall, and P. Schroeder, "Spin Transformations of
  Discrete Surfaces," ACM Transactions on Graphics (SIGGRAPH) 2011.

The basic functionality is to take a manifold, triangulated mesh (without
boundary) plus a grayscale image map as input, and produce a new mesh
with the same connectivity but modified vertex positions.  Triangles in the new
mesh will have nearly the same angles and aspect ratios as those from the
original mesh, but with curvature that differs by the amount specified in the
input image.  (See the paper above for more details.)  In the GUI version of
SpinXForm the input image is reloaded from disk every time the deformation is
updated.  This way the user can "paint" a curvature change in some external
program (such as Adobe Photoshop) and immediately view the resulting effect
on the surface.

PLEASE NOTE that this is research code and has not been tested extensively.
Luckily, it was written by a friendly researcher who is happy to help you out
via email in your times of trouble: keenan@cs.caltech.edu

Thanks to Robert Bridson for providing the preconditioned conjugate gradient
code (found in the "pcg" subdirectory).


------------------
2. DEPENDENCIES
------------------

In order to make the code as easy as possible to build, several distributions
are available, each with an increasing number of dependencies:

  -spinxform_nodep  -- no dependencies; should build with any C++ compiler
  -spinxform_opengl -- depends only on OpenGL and GLUT
  -spinxform_fast   -- depends on OpenGL, GLUT, and SuiteSparse

The first two distributions are easiest to build, but also exhibit suboptimal
performance because they use a simple conjugate gradient solver (with only a
diagonal preconditioner).  The "fast" version depends on Tim Davis' CHOLMOD
cholesky factorization library:

http://www.cise.ufl.edu/research/sparse/cholmod/

which in turn depends on SuiteSparse and METIS:

http://www.cise.ufl.edu/research/sparse/SuiteSparse/
http://glaros.dtc.umn.edu/gkhome/views/metis

as well as some (hopefully optimized!) BLAS/LAPACK implementation.  On
UNIX-like systems you will probably end up needing the libraries

  bamd.a
  libcamd.a
  libcolamd.a
  libccolamd.a
  libcholmod.a

from SuiteSparse and

  libmetis.a

from METIS.  If you want to avoid compiling all of SuiteSparse, you can simply
type "make" in each of the appropriate Lib directories (e.g., AMD/Lib) after
setting up UFConfig and copying the resulting library (.a) files to
/usr/local/bin or some other appropriate place.  Further instructions on
building SuiteSparse and its dependencies can be found on the SuiteSparse home
page.

On Mac OS X, the easiest way to link to an efficient BLAS
implementation is by adding the framework

  -framework Accelerate

On other platforms, Kazushige Goto's GotoBLAS library is a popular choice:

  http://www.tacc.utexas.edu/tacc-projects/gotoblas2/


--------------
3. BUILDING
--------------

On UNIX-like systems, you should simply need to type

  make

at the command line in the root install directory.  This should build the
command-line utility "spinxform."  On Windows it may be simplest to install Cygwin:

  http://www.cygwin.com/

which includes an implementation of the GNUMake system which can be used to
execute the Makefiles for SpinXForm, SuiteSparse, and METIS.  If you need to build
code that does not depend on Cygwin DLLs, MINGW is an option:

  http://www.mingw.org/

Finally, a VisualStudio project has been included, but uses only the (slower)
conjugate gradient solver.  For compiling SuiteSparse, METIS, or BLAS
with VisualStudio, you may want to take a look at

  http://matrixprogramming.com/2008/05/umfpack-vc


-------------
4. TESTING
-------------

To test that SpinXForm is working properly, type

   make test

at the command line.  Results produced by SpinXForm will be checked against
reference solutions included in th distribution.  Depending on output
formatting and the libraries used for BLAS, linear solver, etc., your result
may not match the reference solution exactly, but it should be pretty close.


-----------
5. USAGE
-----------

SpinXForm can be executed from the command line (see Section 5B), by running
SpinXForm.app (on Mac OS X), or by running SpinXForm.bat (Windows).  On Mac OS
X the user will be prompted for the mesh file and then the corresponding image
file.  In Windows, inputs are specified by editing SpinXForm.bat.

5A. INPUTS
----------

Meshes should be encoded as Wavefront OBJ files; images should be encoded as
grayscale, 8 bits per pixel Truevision TGA files.  (See Section 6 for more
information on file formats.)  

The input image file represents the requested change in curvature as grayscale
values.  A value of 128 (~50%) represents no change; values greater than 128
indicates an increase in curvature, and values less than 128 indicate a
decrease.

** Note: the input mesh MUST have texture coordinates! **

Values are mapped from the image to the surface using texture coordinates in
the range [0,1] x [0,1].  This means that the image is treated as a SQUARE
image, independent of its actual dimensions.


5B. COMMAND LINE MODE
---------------------

The SpinXForm executable takes two input arguments and one optional output argument:

   spinxform mesh.obj image.tga [result.obj]

These arguments can be specified on the command line or by creating a .BAT file
in Windows.  Specifying an output argument runs SpinXForm in batch mode: rather
than displaying the GUI, results are simply written to the command line.


5C. INTERACTIVE MODE
--------------------

The graphical user interface (GUI) displays the current state of the mesh.  The
view is controlled by clicking and dragging anywhere in the main window.
Right-clicking brings up the main menu with the following options:

  [space] Transform  - reload the image from disk and update the deformation
      [r] Reset Mesh - restore the surface to its original configuration
      [w] Write Mesh - write the transformed mesh to "result.obj"
      [q] Exit       - quit SpinXForm

The "view" submenu contains the following options:

  [s] Smooth Shaded  - draws the surface using smooth shading
  [f] Wireframe      - draws the surface as flat-shaded triangles with outlined edges
  [e] QC Error       - displays the per-face deviation from original triangle shape
  [p] Rho            - displays the requested change in curvature
  [t] Textured       - displays a checkerboard pattern
  [-] Shrink Texture - decreases the size of the checkerboard pattern
  [+] Grow Texture   - increases the size of the checkerboard pattern

Each of these options has a keyboard equivalent, indicated in square brackets [].
Note that the checkerboard pattern will look good only if the initial texture
coordinates come from something like a conformal parameterization of the surface.


------------------
6. FILE FORMATS
------------------

6A. Wavefront OBJ
-----------------

A Wavefront OBJ file specifies vertex and face data for a polygonal surface
mesh, as well as vertex texture coordinates.  Vertex data is specified by lines
of the form

  v [x] [y] [z]

where [x], [y], and [z] are x, y and z vertex coordinates.  Texture coordinates
are similarly specified via

  vt [u] [v]

Face data is specified by lines of the form

  f [i1]/[j1] [i2]/[j2] ... [in]/[jn]

where [i1], [i2], etc., are 1-based indices into the vertex list and [j1],
j[2], etc., are indices into the list of texture coordinates.  For
instance, the following data specifies a simple mesh of a square
consisting of two triangles:

  v 0 0 0
  v 1 0 0
  v 1 1 0
  v 0 1 0
  vt 0 0
  vt 1 0
  vt 1 1
  vt 0 1
  f 1/1 3/3 4/4
  f 1/1 2/2 3/3

Note that SpinXForm supports only meshes where all faces are triangles.  Vertex
normals are ignored.

6B. Truevision TGA
------------------

TGA is a very simple bitmap image format consisting of a fixed amount of header
data followed by a chunk of color data.  A more detailed specification can be
found at

http://en.wikipedia.org/wiki/Truevision_TGA

Note that SpinXForm supports only uncompressed grayscale images at 8 bits per pixel.


-----------------
7. SOURCE CODE
-----------------

As an alternative to the command-line utility, you can of course interface with
the code directly.  Basic usage looks like:

  Mesh mesh;
  Image image;
  mesh.read( "mesh.obj" );                  // read a polygon mesh
  image.read( "image.tga" );                // read an image
  const double scale = 5.;                  // set the magnitude of rho
  mesh.setCurvatureChange( image, scale );  // set the values of rho
  mesh.updateDeformation();                 // compute the solution
  mesh.write( "result.obj" );               // write the transformed surface


-------------
8. LICENSE
-------------

SpinXForm is covered by the following free software license based on the 2-clause
BSD license.  This license is compatible with most other free software licenses,
including the GNU General Public License.

*
* Copyright 2011 Keenan Crane. All rights reserved.
* 
* Redistribution and use in source and binary forms, with or without modification,
* are permitted provided that the following conditions are met:
* 
* 1. Redistributions of source code must retain the above copyright notice, this
*    list of conditions and the following disclaimer.
* 2. Redistributions in binary form must reproduce the above copyright notice,
*    this list of conditions and the following disclaimer in the documentation
*    and/or other materials provided with the distribution.
* 
* THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR IMPLIED
* WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
* MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
* SHALL THE FREEBSD PROJECT OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
* INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
* LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
* PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
* LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
* OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
* ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
* 
* The views and conclusions contained in the software and documentation are those
* of the author and should not be interpreted as representing official policies,
* either expressed or implied, of any other person or institution.
*

