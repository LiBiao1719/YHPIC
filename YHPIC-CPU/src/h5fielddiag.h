#ifndef _H5PARTDIAG_H
#define _H5PARTDIAG_H

#include "otypen.h"
#include "readfile.h"
#include "species.h"
#include "cellinfo.h"
#include "fieldinfo.h"

#include "H5Part.h"
#include "H5Block.h"

class H5FieldDiag
{
public:
    readfile rf;
    int Q;
    int t_start;               //time of start diag
    int t_stop;                //time of stop  diag
    int t_step;                //time step of  diag
    int t_count;
    int t_save;
    
    //box
    int cells_per_wl;
    int cells_axis;  //cells number of axis direction
    int cells_trans; //cells number of transverse direction

    int q_e;
    int q_b;
    int q_x;
    int q_y;

    int  XMIN;
    int  XMAX;
    int  YMIN;
    int  YMAX;
    int  ZMIN;
    int  ZMAX;

    int  imin;
    int  imax;
    int  jmin;
    int  jmax;
    int  kmin;
    int  kmax;

    char path[100];
    char h5filename[100];
    H5PartFile* h5partfile;

    int  h5index;

    FieldInfo*  field;
    ParParam* param;
   
    H5FieldDiag(char* filename,ParParam* p,FieldInfo* f);

    void dump_field(int step);
    void write_field_item(char *name,int qt);
    void write_density_item(char *name,int isp);

};

#endif
      
