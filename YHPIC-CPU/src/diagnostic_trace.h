#ifndef DIAG_TRACE_H
#define DIAG_TRACE_H
#include "parparam.h"
#include "readfile.h"
#include "fieldinfo.h"

struct tracepositions
{
  int i;
  int j;
  int k;
};

class Diagnostic_trace
{
 public:
  readfile rf;
  ParParam* para;

  char name[100];
  int IMIN,IMAX,JMIN,JMAX,KMIN,KMAX;  //region size

  int Q;
  int t_start;
  int t_stop;
  int t_step;
  int traces;
  int ISTEP;
  struct tracepositions *tracepos;
  Scalar *ex;
  Scalar *ey;
  Scalar *ez;
  Scalar *jx;
  Scalar *jy;
  Scalar *jz;

  Diagnostic_trace(char* filename,ParParam* para);
  ~Diagnostic_trace();
  void store_traces(FieldInfo* field,Scalar t);
  void write_traces(FieldInfo* field,Scalar t);
};
#endif

