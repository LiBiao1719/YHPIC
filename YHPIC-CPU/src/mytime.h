#ifndef MY_TIME_H
#define MY_TIME_H

#include <time.h>

#ifdef PARALLEL
#include <mpi.h>
#endif

double pic_time();

#endif





