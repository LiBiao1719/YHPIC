#ifndef CELLINFO_H 
#define CELLINFO_H
#include "fieldinfo.h"
#include "particle.h"
#include "cell.h"
#include "parparam.h"
#include "ostack.h"
class Stack;

class CellInfo
{
 public:
  ParParam *para;
  FieldInfo *field;
  Vector3 enf[8];  //temper for accelerate particle
  Vector3 bnf[8];  //used for optimize
  Stack *stk;
  int size[6];
  int IMIN,IMAX,JMIN,JMAX,KMIN,KMAX;
  Scalar dx,dt;
  Matrix<cell> CellsM;
  int n_el;
  int n_ion;
  int n_part;

  Scalar cells_left_axis;   
  Scalar cells_plasma_axis; 
  Scalar cells_ramp_axis;   

#ifdef HIGHORDERCLOUD
  Plane<Scalar> S0;     //
  Plane<Scalar> S1;     //  used for highorder CIC cloud
  Matrix<Vector3> WW;   //
  Matrix<Vector3> Ipar; //
#endif

  CellInfo(FieldInfo* field);
  ~CellInfo();
  Scalar gauss_rand48();
  void InitCells();
  void ChainParticles();
  void InitParticles();
  void advance(Scalar t);
  void acceleratefield(int i,int j,int k);
  void accelerate(particle* part,int i,int j,int k);
  void move(particle* part,int i,int j,int k,volatile Scalar &pdx,volatile Scalar &pdy,volatile Scalar &pdz);
  void deposit_charge(particle* part,int i,int j,int k);
  void has_to_change_cell(particle* part,int i,int j,int k);
  void deposit_current(particle* part,int i,int j,int k,volatile Scalar pdx,volatile Scalar pdy,volatile Scalar pdz);
  void deposit_gamma(Scalar t);
  void do_change_cell();
  void reflect_particles();
  void periodic_particles();
  void exchange_particles();
  void shift_particles(int direction);
  void shift_chainparticles(int direction);
  void shift_initparticles(int direction);
  void ptsndrcv(int i0 ,int i1 ,int j0 ,int j1 ,int k0 ,int k1 ,int ID_s,
		int ii0,int ii1,int jj0,int jj1,int kk0,int kk1,int ID_r);
  void messagepack(int i0,int i1,int j0,int j1,int k0,int k1,int& send_num);
  void messageunpack(int i0,int i1,int j0,int j1,int k0,int k1,int& recv_num);
  Scalar particle_energy();
  Scalar energy();
  Scalar memory(){ return n_part*sizeof(particle); }

  void dump(FILE* file);
  void restore(FILE* file);

  friend ofstream&
    operator<<(ofstream& os,const CellInfo& cf);
};

#endif
