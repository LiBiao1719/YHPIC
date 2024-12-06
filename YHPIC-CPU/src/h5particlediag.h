#ifndef _H5PARTICLEDIAG_H
#define _H5PARTICLEDIAG_H

#include "otypen.h"
#include "readfile.h"
#include "species.h"
#include "cellinfo.h"

#include "H5Part.h"
#include "H5Block.h"

class H5ParticleDiag
{
public:
    readfile rf;
    int Q;
    Scalar t_start;               //time of start diag
    Scalar t_stop;                //time of stop  diag
    Scalar t_step;                //time step of  diag
    int t_count;
    int t_save;
    
    //box
    int cells_per_wl;
    int cells_axis;  //cells number of axis direction
    int cells_trans; //cells number of transverse direction

    int  IMIN;
    int  IMAX;
    int  JMIN;
    int  JMAX;
    int  KMIN;
    int  KMAX;
    
    //species
    int species_num;
    struct species specy[NSPECIES];

    //file name and path
    char path[100];
    char h5filename[100];
    H5PartFile* h5partfile;

    CellInfo *cells;
    
    H5ParticleDiag(char* filename,ParParam* p, CellInfo* f);
    void dump_particle(int step);
    void write_particle_item(char *name,int qt);
};

#endif
      
