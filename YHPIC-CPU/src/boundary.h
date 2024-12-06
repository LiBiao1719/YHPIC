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

  Scalar spotwidth;  //half width of waist
  Scalar amplitude;  //
  int cells_per_wl;  //
  int focus_z;       //focus plane position
  int polarization;
  int shape;
  int raise;
  int duration;
  
  Scalar envelop(Scalar t);
  Scalar Ey(int i,int j,ParParam* p,Scalar t);
  Scalar Ex(int i,int j,ParParam* p,Scalar t);
    
  LaserParam(char *filename);
  void pulsecorrect(int TYPE,FieldInfo& f,ParParam *p,Scalar t); //
  void boundarycorrect(int TYPE,FieldInfo& f,ParParam *p,Scalar t); 
};
#endif

