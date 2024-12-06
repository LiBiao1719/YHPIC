#ifndef _PARPARAM_H
#define _PARPARAM_H

#ifdef PARALLEL
#include "mpi.h"
#define ECOMM           0*6+1           // several communication IDs
#define PCOMM           1*6+1
#endif

#include "otypen.h"
#include "readfile.h"
#include "species.h"

class FieldInfo;
class LaserParam;

class SerParam
{
public:
    readfile rf;
    int XMIN,XMAX,YMIN,YMAX,ZMIN,ZMAX;

    //box
    int cells_per_wl;
    int cells_axis;  //cells number of axis direction
    int cells_trans; //cells number of transverse direction
    int cells_left;  //cells number of left vacuum
    int cells_gap;   //cells number of upper and down vacuum
    int cells_plasma;  //cells number of plasma region
    int cells_ramp;
    Scalar angle;    //

    Scalar n_ion_over_nc;  // maximum density/critical density    
    Scalar n_el_over_nc;
    Scalar n_ion_ramp;
    //prograte
    int prop_start;
    int prop_stop;
    int prop_save;
    int prop_restart;
    int flag_restart;

    //electron-0 ion-1
    struct species specy[4];
    //    int fix[2];
    //    int z[2];
    //    int m[2];
    int ppc[2];
    
    Scalar vtherm[2]; 
    //shift region
    int shift_dstep;
    int shift_step;
    int shift_flag;

    char *path;
    SerParam(char* filename);
};

class ParParam
{
public:

    int IMIN,JMIN,KMIN;
    int IMAX,JMAX,KMAX;

    int IDme;
    int ID_xr,ID_xl,ID_yl,ID_yr,ID_zl,ID_zr; //-1 means no neighbour
#ifdef PARALLEL
    MPI_Comm    MYCOMM;                   // MPI specific stuff
    MPI_Request srequest[6][2]; 
    MPI_Request rrequest[6][2];  
#endif

    Scalar *sbuf[6],*rbuf[6];
    int sbuflen[6],rbuflen[6];

    double *sbufpt,*rbufpt;
    int sbufptlth,rbufptlth;
    // input parameter
    SerParam *parameter;
    LaserParam* laser;
    
    ParParam(SerParam* spara,LaserParam* lpara);
    ~ParParam();
    void initialize(); 
    Scalar* ptrBuffer(int i){ return rbuf[i]; }
    Scalar* ptsBuffer(int i){ return sbuf[i]; }
    void commField(int dirction);
    void commParticle(int ID_s,int ID_r,int sendnum,int& recvnum);
    void correctField(FieldInfo &f,Scalar t);
};
#endif
      
