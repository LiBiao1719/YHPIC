#ifndef DIAGNOSTIC_H
#define DIAGNOSTIC_H
#include "parparam.h"
#include "cellinfo.h"

const int DIM = 400;

class Diagnostic
{
 public:
  char name[100];
  char name_err[100];
  char name_vtu[100];
  char name_vts[100];
  char name_pvts[100];
   char name1[100];
   char name_pvtu[100];
   char name2[100];
  int IMIN,IMAX,JMIN,JMAX,KMIN,KMAX;  //region size
  int x[DIM],y[DIM],z[DIM],a[DIM];
  Scalar vcut;
int iTime;
  int  diag_save;
  Scalar  diag_dt;

  ParParam* para;

  Diagnostic(ParParam* para);
  
  void diagnostic_title();
  void diagnostic_time();
  void diagnostic_velocity(int spec,CellInfo*,Scalar t);
  void diagnostic_field(FieldInfo* field,Scalar t);
  void diagnostic_write_pvtu();
  void diagnostic_write_pvd();
  void diagnostic_particle(CellInfo* cells,Scalar t);
  void diagnostic_energy(FieldInfo* field,CellInfo* cells,Scalar t);
  void diagnostic_dump(FieldInfo* fields,CellInfo* cells,Scalar t);
  void diagnostic_restore(FieldInfo* fields,CellInfo* cells);
  void diagnostic_dumpH5(FieldInfo* fields,CellInfo* cells,Scalar t);
};
#endif

