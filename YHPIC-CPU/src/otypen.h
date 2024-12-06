#ifndef _MYMATH_H
#define _MYMATH_H
#include <math.h>
#ifndef PI
#define PI 3.1415926
#define TWOPI (2.0*PI)
#endif
#define NSPECIES 2
#define DAMP_LENGTH 50
//#define LOWPRECISION
#ifdef PARALLEL
#include <mpi.h>
#endif 

#ifdef NOSCALE
#define SCALE 1
#else 
#define SCALE 1
#endif

#ifdef LOWPRECISION
#define Scalar float
#define FLOOR(x) floorf(x)
#define CEIL(x)  ceilf(x)
#define SQRT(x)  sqrtf(x)
#define POW(x,y) powf(x,y)
#define RINT(x)  rintf(x)
#define H5T_OOPIC_SCALAR H5T_NATIVE_FLOAT
#ifdef PARALLEL
#define MPI_SCALAR MPI_FLOAT
#endif
#else
#define Scalar double
#define FLOOR(x) floor(x)
#define CEIL(x)  ceil(x)
#define SQRT(x)  sqrt(x)
#define POW(x,y) pow(x,y)
#define RINT(x)  rint(x)
#define H5T_OOPIC_SCALAR H5T_NATIVE_DOUBLE
#ifdef PARALLEL
#define MPI_SCALAR MPI_DOUBLE
#endif
#endif

#define sqr(x) (x*x)
#define MIN(a,b)			((a<b) ? (a) : (b))
#define MAX(a,b)			((a>b) ? (a) : (b))
#endif
