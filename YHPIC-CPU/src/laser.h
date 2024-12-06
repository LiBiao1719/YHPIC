#ifndef _LASER_H
#define _LASER_H
#include "parparam.h"
#include "fieldinfo.h"
#include "readfile.h"

class LaserParam
{
public:
    
  readfile rf;
  //laser pulse
  int Q;
  int Qtrain;
  int trains;
  Scalar *time_delay;

  Scalar spotwidth;  //half width of waist
  Scalar amplitude;  //
  int cells_per_wl;  //
  int focus_z;       //focus plane position
  int polarization;
  int shape;
  int raise;
  int duration;
  Scalar dt;
  Scalar dx;
  Scalar single_pulse(Scalar t);
  Scalar envelop(Scalar t);
  Scalar Ey(Scalar i,Scalar j,Scalar k,ParParam* p,Scalar t);
  Scalar Ex(Scalar i,Scalar j,Scalar k,ParParam* p,Scalar t);
    
  LaserParam(char *filename);
  ~LaserParam();
  void pulsecorrect(int TYPE,FieldInfo& f,ParParam *p,Scalar t); //
  void absorbcorrect(int TYPE,FieldInfo& f,ParParam *p,Scalar t); //
  void Murabsorbcorrect(int TYPE,FieldInfo& f,ParParam *p,Scalar t); //
  void boundarycorrect(int TYPE,FieldInfo& f,ParParam *p,Scalar t); 
};
#endif

