#ifndef _FIELD_H
#define _FIELD_H
#include <fstream>
#include <iostream>
#include "matrix.h"

#define ILOOP       for(i=IMIN;i<IMAX;i++)
#define JLOOP       for(j=JMIN;j<JMAX;j++)
#define KLOOP       for(k=KMIN;k<KMAX;k++)

#define ALLILOOP    for(i=IMIN-1;i<=IMAX;i++)
#define ALLJLOOP    for(j=JMIN-1;j<=JMAX;j++)
#define ALLKLOOP    for(k=KMIN-1;k<=KMAX;k++)

#define IJLOOP  ILOOP JLOOP
#define IKLOOP  ILOOP KLOOP
#define JKLOOP  JLOOP KLOOP
#define IJKLOOP ILOOP JLOOP KLOOP

#define ALLIJLOOP ALLILOOP ALLJLOOP
#define ALLIKLOOP ALLILOOP ALLKLOOP
#define ALLJKLOOP ALLJLOOP ALLKLOOP
#define ALLIJKLOOP ALLILOOP ALLJLOOP ALLKLOOP

using namespace std;

class ParParam;
class LaserParam;
class FieldInfo
{
public:
    ParParam *par;
    LaserParam *laser;
    int size[6];  //memory allocated size
    int IMIN,IMAX,JMIN,JMAX,KMIN,KMAX;  //region size
    int DAMPKMIN,DAMPKMAX;    //damp+simulation region size
    int DAMP_REGION_F;
    int DAMP_REGION_B;
    Scalar dt,dx,dtdx; //timestep and grid spacing
    Scalar poynting_left_p,poynting_right_p,poynting_left_n,poynting_right_n;
    Scalar laser_flux_left,laser_flux_right;
#ifdef DEBUGMOVE
    Vector3 EStatic;
    Vector3 BStatic;
#endif      
    Matrix<Vector3> E;
    Matrix<Vector3> B;
    Matrix<Vector3> Ic;
    Matrix<Vector3> ENode;
    Matrix<Vector3> BNode;
    Matrix<Scalar>  Dens[NSPECIES];

    Matrix<Vector3> E_older_min;
    Matrix<Vector3> E_oldest_min;
    Matrix<Vector3> E_older_max;
    Matrix<Vector3> E_oldest_max;

    Matrix<Vector3> Eave;
    Matrix<Vector3> Bave;

    FieldInfo(ParParam *p);
    Scalar get_dt(){ return dt; }

    void correctField();         // *****   message exchange member function ****
    void messagepack1st(int DI);  //     DI = 0: message sendrcv in x-direction
    void messageunpack1st(int DI);//     DI = 1: message sendrcv in y-direction 
    void messagepack2nd(int DI);  //     DI = 2: message sendrcv in z-direction
    void messageunpack2nd(int DI);//
    void exchange1st();          //     step = 1: exchange 1nd
    void exchange2nd();          //     step = 1: exchange 1st

    void exchangeDensity();          //     step = 1: exchange 1st
    void messagepackDensity(int DI);  //     DI = 0: message sendrcv in x-direction
    void messageunpackDensity(int DI);//     DI = 1: message sendrcv in y-direction 

    void advanceB();
    void advanceE();
    void advance(Scalar t);
    void clearQI();

    void front_flux_change(Scalar t);
    void back_flux_change(Scalar t);

    void shiftField(int direction);
    void average(Scalar t);
    Scalar damp_function(int j,int LENGTH,Scalar DAMP_R);
    Scalar energy();
    Scalar work_energy();
    Scalar memory(){ return (IMAX-IMIN+5)*(JMAX-JMIN+5)*(KMAX-KMIN+5)*18*8; }
    void dump(FILE* file);
    void restore(FILE* file);
    void gatherDataXz(Matrix<Vector3>& component,Scalar **data_out,int JM,int xmin,int xmax,int zmin,int zmax,int root,MPI_Comm oopic_comm);
    void gatherDataYz(Matrix<Vector3>& component,Scalar **data_out,int IM,int ymin,int ymax,int zmin,int zmax,int root,MPI_Comm oopic_comm);
    friend ofstream&
        operator<<(ofstream& os, const FieldInfo& f); 
};

// NEW TEST
#endif
