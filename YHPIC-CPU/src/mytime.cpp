#include "mytime.h"

double pic_time()
{
  double x;
#ifdef PARALLEL
  return MPI_Wtime();
#else
  //  return cpu_time();
#endif  
}
