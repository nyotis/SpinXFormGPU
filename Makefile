# replaced "-pedantic -ansi" CFLAGS and LDFLAGS from original makefile with
# -Wno-long-long due to this bug in Thrust
# https://github.com/thrust/thrust/issues/11
# https://groups.google.com/forum/#!msg/thrust-users/jSujgymfbwM/MckpqfhvcFsJ

# when the -G debug flag is used for a pre-Fermi architecture too many warnings
# are triggered. there is nothing wrong, simply ignore them (this is a known 
# limitation of nvcc when targeting sm_13 and earlier hardware, 
# see http://goo.gl/3or39) 

#LLVM build removed -ansi from CFLAGS | LDFLAGS
# run ./spinxform ./examples/bumpy/sphere.obj ./examples/bumpy/bumpy.tga ./examples/bumpy/result.obj

TARGET = spinxformgpu

CFLAGS  = -Wall -Werror -Wno-long-long -O3 -Iinclude -I/usr/local/cuda/include/

#debug version
#CFLAGS  = -Wall -Werror -Wno-long-long -O0 -g -G -Iinclude -I/usr/local/cuda/include/
# in contrast to thrust, cusp not part of official CUDA release so need this -I/usr/local/cuda/include/

LDFLAGS = -Wall -Werror -O3
#LDFLAGS = -Wall -Werror -O0 -g -G
LIBS = -L/usr/local/cuda/lib -lcudart
OBJS = EigenSolver.o cusp_device.o Image.o LinearSolver.o Mesh.o Quaternion.o QuaternionMatrix.o Vector.o main.o

all: $(TARGET)

$(TARGET): $(OBJS)
	g++ $(OBJS) $(LDFLAGS) $(LIBS) -o $(TARGET)

EigenSolver.o: src/EigenSolver.cpp include/EigenSolver.h include/QuaternionMatrix.h include/Quaternion.h include/Vector.h include/sparse_matrix.h include/util.h include/LinearSolver.h
	g++ $(CFLAGS) -c src/EigenSolver.cpp

Image.o: src/Image.cpp include/Image.h
	g++ $(CFLAGS) -c src/Image.cpp

cusp_device.o: cusp_device.cu include/cusp_device.h
#	nvcc -m64 -G -c cusp_device.cu
	nvcc -m64 -c cusp_device.cu
	
LinearSolver.o: src/LinearSolver.cpp include/LinearSolver.h include/QuaternionMatrix.h include/Quaternion.h include/Vector.h include/sparse_matrix.h include/util.h
	g++ $(CFLAGS) -c src/LinearSolver.cpp
    	
Mesh.o: src/Mesh.cpp include/Mesh.h include/Quaternion.h include/Vector.h include/QuaternionMatrix.h include/sparse_matrix.h include/util.h include/Image.h include/LinearSolver.h include/EigenSolver.h include/Utility.h
	g++ $(CFLAGS) -c src/Mesh.cpp

Quaternion.o: src/Quaternion.cpp include/Quaternion.h include/Vector.h
	g++ $(CFLAGS) -c src/Quaternion.cpp

QuaternionMatrix.o: src/QuaternionMatrix.cpp include/QuaternionMatrix.h include/Quaternion.h include/Vector.h include/sparse_matrix.h include/util.h
	g++ $(CFLAGS) -c src/QuaternionMatrix.cpp

Vector.o: src/Vector.cpp include/Vector.h
	g++ $(CFLAGS) -c src/Vector.cpp

main.o: src/main.cpp
	g++ $(CFLAGS) -c src/main.cpp
	

clean:
	rm -f $(TARGET)
	rm -f *.o
	rm -f examples/bumpy/solution.obj
	rm -f examples/spacemonkey/solution.obj

