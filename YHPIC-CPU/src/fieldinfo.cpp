#include "fieldinfo.h"
#include "otypen.h"
#include "parparam.h"
#include "mytime.h"
#include "laser.h"

#ifdef PARALLEL
extern  MPI_Comm OOPIC_COMM;
extern int MY_RANK;
extern int NP_X;
extern int NP_Y;
extern int NP_Z;
MPI_Comm OOPIC_XZ_COMM;
MPI_Comm OOPIC_YZ_COMM;
#endif
#define DAMP 0

const Scalar DAMP_R = 1.0;

extern Scalar field_advance_time,field_exchange_time;

FieldInfo::FieldInfo(ParParam* para)
{
	par = para;
	laser = para->laser;
	dx = 1.0/par->parameter->cells_per_wl;
	dt = 0.5*dx;
	dtdx = 0.5;

	poynting_left_p = 0.0;
	poynting_left_n = 0.0;
	poynting_right_p = 0.0;
	poynting_right_n = 0.0;
	laser_flux_left = 0.0;
	laser_flux_right = 0.0;

    IMIN = par->IMIN;
    IMAX = par->IMAX;
    JMIN = par->JMIN;
    JMAX = par->JMAX;
    KMIN = par->KMIN;
    KMAX = par->KMAX;
    
    DAMPKMIN = KMIN;
    DAMPKMAX = KMAX;
    if(KMIN == par->parameter->ZMIN) DAMPKMIN = KMIN-DAMP_LENGTH;
    if(KMAX == par->parameter->ZMAX) DAMPKMAX = KMAX+DAMP_LENGTH;

    size[0] = IMIN-2;
    size[1] = IMAX+2;
    size[2] = JMIN-2;
    size[3] = JMAX+2;
    size[4] = DAMPKMIN-2;
    size[5] = DAMPKMAX+2;
  

    E.Init(size); ENode.Init(size);
    B.Init(size); BNode.Init(size);
    Eave.Init(size); Bave.Init(size);
    Ic.Init(size); 

    for(int isp=0;isp<NSPECIES;isp++)
      Dens[isp].Init(size);

    int size_t[6] = {IMIN-2,IMAX+2,JMIN-2,JMAX+2,0,1};
    if(KMIN==par->parameter->ZMIN){
        E_older_min.Init( size_t);
        E_oldest_min.Init(size_t);
    }
    if(KMAX==par->parameter->ZMAX){
        E_older_max.Init( size_t);
        E_oldest_max.Init(size_t);
    }

    for(int i=size[0];i<=size[1];i++)
      for(int j=size[2];j<=size[3];j++)
	for(int k=size[4];k<=size[5];k++){
	  E[i][j][k] = Vector3(0.0,0.0,0.0);
	  B[i][j][k] = Vector3(0.0,0.0,0.0);
	  ENode[i][j][k] = Vector3(0.0,0.0,0.0);
	  BNode[i][j][k] = Vector3(0.0,0.0,0.0);
	  Ic[i][j][k] = Vector3(0.0,0.0,0.0);
	}

    int aa;
}

void FieldInfo::advance(Scalar t)
{
  Scalar t0,t1,t2,t3,t4,t5,t6;
  //B(0)->B(0.5)
  advanceB();
  t1 = pic_time();
  exchange1st();
  //E(0)->E(1.0)
  t2 = pic_time();
  advanceE();
  par->correctField(*this,t);
    
  //B(0.5)->B(1.0)
  advanceB();
  t3 = pic_time();
  exchange2nd();
  t4 = pic_time();
  average(t);

  field_exchange_time += (Scalar)((t4-t3)+(t2-t1));

  //do flux computing
  front_flux_change(t);
  back_flux_change(t);

}

void FieldInfo::exchange1st()
{
  for(int i=0;i<3;i++){
    messagepack1st(i);
    par->commField(i);
    messageunpack1st(i);
  }
}

void FieldInfo::exchange2nd()
{
  for(int i=0;i<3;i++){
    messagepack2nd(i);
    par->commField(i);
    messageunpack2nd(i);
  }
}

void FieldInfo::exchangeDensity()
{
  for(int i=0;i<3;i++){
    messagepackDensity(i);
    par->commField(i);
    messageunpackDensity(i);
  }
}

void FieldInfo::messagepackDensity(int DIRECTION)
{
  int index;
  switch(DIRECTION){ 
  case 2: /*************************  z-dirction ******************************/
    {
      // front k=KMAX
      index=0;
      Scalar* buff0=par->ptsBuffer(0);
      for(int i=IMIN-2;i<=IMAX+2;i++)
	for(int j=JMIN-2;j<=JMAX+2;j++)
	  {
	    for(int isp=0;isp<NSPECIES;isp++){
	      buff0[index++]=Dens[isp][i][j][KMAX];
	      buff0[index++]=Dens[isp][i][j][KMAX+1];
	      buff0[index++]=Dens[isp][i][j][KMAX+2];
	    }
	  }
      par->sbuflen[0]=index++;
      
      //back k=KMIN
      index=0;
      Scalar* buff1=par->ptsBuffer(1);

      for(int i=IMIN-2;i<=IMAX+2;i++)
	for(int j=JMIN-2;j<=JMAX+2;j++)
	  {
	    for(int isp=0;isp<NSPECIES;isp++){
	      buff1[index++]=Dens[isp][i][j][KMIN];
	      buff1[index++]=Dens[isp][i][j][KMIN-1];
	      buff1[index++]=Dens[isp][i][j][KMIN-2];
	    }
	  }    
      par->sbuflen[1]=index++;
    }
    break;
  case 1: /*************************  y-dirction ******************************/
    {
      //right J=JMAX
      index=0;
      Scalar* buff4=par->ptsBuffer(4);

      for(int i=IMIN-2;i<=IMAX+2;i++)
	for(int k=DAMPKMIN-2;k<=DAMPKMAX+2;k++)
	  {
	    for(int isp=0;isp<NSPECIES;isp++){
	      buff4[index++]=Dens[isp][i][JMAX][k];
	      buff4[index++]=Dens[isp][i][JMAX+1][k];
	      buff4[index++]=Dens[isp][i][JMAX+2][k];
	    }
	  }
      par->sbuflen[4]=index++;
      //left J=JMIN
      index=0;
      Scalar* buff5=par->ptsBuffer(5);

      for(int i=IMIN-2;i<=IMAX+2;i++)
	for(int k=DAMPKMIN-2;k<=DAMPKMAX+2;k++)
	  {
	    for(int isp=0;isp<NSPECIES;isp++){
	      buff5[index++]=Dens[isp][i][JMIN][k];
	      buff5[index++]=Dens[isp][i][JMIN-1][k];
	      buff5[index++]=Dens[isp][i][JMIN-2][k];
	    }
	  }    
      par->sbuflen[5]=index++;
    }
    break;
  case 0: /*************************  x-dirction ******************************/
    {
      //right I=IMAX
      index=0;
      Scalar* buff2=par->ptsBuffer(2);

      for(int j=JMIN-2;j<=JMAX+2;j++)
	for(int k=DAMPKMIN-2;k<=DAMPKMAX+2;k++)
	  {
	    for(int isp=0;isp<NSPECIES;isp++){
	      buff2[index++]=Dens[isp][IMAX][j][k];
	      buff2[index++]=Dens[isp][IMAX+1][j][k];
	      buff2[index++]=Dens[isp][IMAX+2][j][k];
	    }
	  }
      par->sbuflen[2]=index++;
      //left I=IMIN
      index=0;
      Scalar* buff3=par->ptsBuffer(3);

      for(int j=JMIN-2;j<=JMAX+2;j++)
	for(int k=DAMPKMIN-2;k<=DAMPKMAX+2;k++)
	  {
	    for(int isp=0;isp<NSPECIES;isp++){
	      buff3[index++]=Dens[isp][IMIN][j][k];
	      buff3[index++]=Dens[isp][IMIN-1][j][k];
	      buff3[index++]=Dens[isp][IMIN-2][j][k];
	    }
	  }    
      par->sbuflen[3]=index++;
    }
    break;
  }
}

void FieldInfo::messageunpackDensity(int DIRECTION)
{
  int index;
  switch(DIRECTION){
  case 2:
    {
      // front k=KMIN
      index=0;
      
      Scalar* buff0=par->ptrBuffer(0);
      if( KMIN != par->parameter->ZMIN){

	for(int i=IMIN-2;i<=IMAX+2;i++)
	for(int j=JMIN-2;j<=JMAX+2;j++)  
	  {
	    for(int isp=0;isp<NSPECIES;isp++){
	      Dens[isp][i][j][KMIN]+=buff0[index++];
	      Dens[isp][i][j][KMIN+1]+=buff0[index++];
	      Dens[isp][i][j][KMIN+2]+=buff0[index++];
	    }
	  }
      }
      //back k=KMAX
      index=0;
      Scalar* buff1=par->ptrBuffer(1);
      if( KMAX != par->parameter->ZMAX){

	for(int i=IMIN-2;i<=IMAX+2;i++)
	for(int j=JMIN-2;j<=JMAX+2;j++)
	  {
	    for(int isp=0;isp<NSPECIES;isp++){
	      Dens[isp][i][j][KMAX]+=buff1[index++];
	      Dens[isp][i][j][KMAX-1]+=buff1[index++];
	      Dens[isp][i][j][KMAX-2]+=buff1[index++];
	    }
	  }
      }    
    }
    break;
  case 1:
    { 
      //left J=JMIN
      index=0;
      Scalar* buff4=par->ptrBuffer(4);

      for(int i=IMIN-2;i<=IMAX+2;i++)
        for(int k=DAMPKMIN-2;k<=DAMPKMAX+2;k++)
	  {
	    for(int isp=0;isp<NSPECIES;isp++){
	      Dens[isp][i][JMIN][k]+=buff4[index++];
	      Dens[isp][i][JMIN+1][k]+=buff4[index++];
	      Dens[isp][i][JMIN+2][k]+=buff4[index++];
	    }
	  }
      
      //right J=JMAX
      index=0;
      Scalar* buff5=par->ptrBuffer(5);

      for(int i=IMIN-2;i<=IMAX+2;i++)
      for(int k=DAMPKMIN-2;k<=DAMPKMAX+2;k++)
	{
	  for(int isp=0;isp<NSPECIES;isp++){
	    Dens[isp][i][JMAX][k]+=buff5[index++];
	    Dens[isp][i][JMAX-1][k]+=buff5[index++];
	    Dens[isp][i][JMAX-2][k]+=buff5[index++];
	  }
	} 
    }
    break;
  case 0:
    { 
      //left I=IMIN
      index=0;
      Scalar* buff2=par->ptrBuffer(2);

      for(int j=JMIN-2;j<=JMAX+2;j++)
        for(int k=DAMPKMIN-2;k<=DAMPKMAX+2;k++)
	  {
	    for(int isp=0;isp<NSPECIES;isp++){
	      Dens[isp][IMIN][j][k]+=buff2[index++];
	      Dens[isp][IMIN+1][j][k]+=buff2[index++];
	      Dens[isp][IMIN+2][j][k]+=buff2[index++];
	    }
	  }
      
      //right I=IMAX
      index=0;
      Scalar* buff3=par->ptrBuffer(3);

      for(int j=JMIN-2;j<=JMAX+2;j++)
      for(int k=DAMPKMIN-2;k<=DAMPKMAX+2;k++)
	{
	  for(int isp=0;isp<NSPECIES;isp++){
	    Dens[isp][IMAX][j][k]+=buff3[index++];
	    Dens[isp][IMAX+1][j][k]+=buff3[index++];
	    Dens[isp][IMAX+2][j][k]+=buff3[index++];
	  }
	} 
    }
    break;
  }                                                                                               
}

void FieldInfo::messagepack1st(int DIRECTION)
{
  int index;
  switch(DIRECTION){ 
  case 2: /*************************  z-dirction ******************************/
    {
      // front k=KMAX
      index=0;
      Scalar* buff0=par->ptsBuffer(0);
      for(int i=IMIN-2;i<=IMAX+2;i++)
	for(int j=JMIN-2;j<=JMAX+2;j++)
	  {
	    buff0[index++]=B[i][j][KMAX-1].e1();
	    buff0[index++]=B[i][j][KMAX-1].e2();
	    buff0[index++]=B[i][j][KMAX-1].e3();
	    buff0[index++]=Ic[i][j][KMAX].e1();
	    buff0[index++]=Ic[i][j][KMAX].e2();
	    buff0[index++]=Ic[i][j][KMAX].e3();
	    buff0[index++]=Ic[i][j][KMAX+1].e1();
	    buff0[index++]=Ic[i][j][KMAX+1].e2();
	    buff0[index++]=Ic[i][j][KMAX+1].e3();
#ifdef HIGHORDERCLOUD
	    buff0[index++]=Ic[i][j][KMAX+2].e1();
	    buff0[index++]=Ic[i][j][KMAX+2].e2();
	    buff0[index++]=Ic[i][j][KMAX+2].e3();
#endif
	  }
      par->sbuflen[0]=index++;
      
      //back k=KMIN
      index=0;
      Scalar* buff1=par->ptsBuffer(1);
      for(int i=IMIN-2;i<=IMAX+2;i++)
	for(int j=JMIN-2;j<=JMAX+2;j++)
	  {
	    buff1[index++]=B[i][j][KMIN].e1();
	    buff1[index++]=B[i][j][KMIN].e2();
	    buff1[index++]=B[i][j][KMIN].e3();
	    buff1[index++]=Ic[i][j][KMIN].e1();
	    buff1[index++]=Ic[i][j][KMIN].e2();
	    buff1[index++]=Ic[i][j][KMIN].e3();
	    buff1[index++]=Ic[i][j][KMIN-1].e1();
	    buff1[index++]=Ic[i][j][KMIN-1].e2();
	    buff1[index++]=Ic[i][j][KMIN-1].e3();
#ifdef HIGHORDERCLOUD
	    buff1[index++]=Ic[i][j][KMIN-2].e1();
	    buff1[index++]=Ic[i][j][KMIN-2].e2();
	    buff1[index++]=Ic[i][j][KMIN-2].e3();
#endif
	  }    
      par->sbuflen[1]=index++;
    }
    break;
  case 0: /*************************  x-dirction ******************************/
    {
      //up i=IMAX
      index=0;
      Scalar* buff2=par->ptsBuffer(2);
      for(int j=JMIN-2 ;j<=JMAX+2 ;j++)
        for(int k=DAMPKMIN-2;k<=DAMPKMAX+2;k++)
	  {
	    buff2[index++]=B[IMAX-1][j][k].e1();
	    buff2[index++]=B[IMAX-1][j][k].e2();
	    buff2[index++]=B[IMAX-1][j][k].e3();
	    buff2[index++]=Ic[IMAX][j][k].e1();
	    buff2[index++]=Ic[IMAX][j][k].e2();
	    buff2[index++]=Ic[IMAX][j][k].e3();
	    buff2[index++]=Ic[IMAX+1][j][k].e1();
	    buff2[index++]=Ic[IMAX+1][j][k].e2();
	    buff2[index++]=Ic[IMAX+1][j][k].e3();
#ifdef HIGHORDERCLOUD
	    buff2[index++]=Ic[IMAX+2][j][k].e1();
	    buff2[index++]=Ic[IMAX+2][j][k].e2();
	    buff2[index++]=Ic[IMAX+2][j][k].e3();
#endif
	  }     
      par->sbuflen[2]=index++;
  
      //down i=IMIN send to down region i=IMAX
      index=0;
      Scalar* buff3=par->ptsBuffer(3);
      for(int j=JMIN-2;j<=JMAX+2;j++)
        for(int k=DAMPKMIN-2;k<=DAMPKMAX+2;k++)
	  {
	    buff3[index++]=B[IMIN][j][k].e1();
	    buff3[index++]=B[IMIN][j][k].e2();
	    buff3[index++]=B[IMIN][j][k].e3();
	    buff3[index++]=Ic[IMIN][j][k].e1();
	    buff3[index++]=Ic[IMIN][j][k].e2();
	    buff3[index++]=Ic[IMIN][j][k].e3();
	    buff3[index++]=Ic[IMIN-1][j][k].e1();
	    buff3[index++]=Ic[IMIN-1][j][k].e2();
	    buff3[index++]=Ic[IMIN-1][j][k].e3();
#ifdef HIGHORDERCLOUD
	    buff3[index++]=Ic[IMIN-2][j][k].e1();
	    buff3[index++]=Ic[IMIN-2][j][k].e2();
	    buff3[index++]=Ic[IMIN-2][j][k].e3();
#endif
	  } 
      par->sbuflen[3]=index++;
    }
    break;
  case 1: /*************************  y-dirction ******************************/
    {
      //right J=JMAX
      index=0;
      Scalar* buff4=par->ptsBuffer(4);
      for(int i=IMIN-2;i<=IMAX+2;i++)
        for(int k=DAMPKMIN-2;k<=DAMPKMAX+2;k++)
	  {
	    buff4[index++]=B[i][JMAX-1][k].e1();
	    buff4[index++]=B[i][JMAX-1][k].e2();
	    buff4[index++]=B[i][JMAX-1][k].e3();
	    buff4[index++]=Ic[i][JMAX][k].e1();
	    buff4[index++]=Ic[i][JMAX][k].e2();
	    buff4[index++]=Ic[i][JMAX][k].e3();
	    buff4[index++]=Ic[i][JMAX+1][k].e1();
	    buff4[index++]=Ic[i][JMAX+1][k].e2();
	    buff4[index++]=Ic[i][JMAX+1][k].e3();
#ifdef HIGHORDERCLOUD
	    buff4[index++]=Ic[i][JMAX+2][k].e1();
	    buff4[index++]=Ic[i][JMAX+2][k].e2();
	    buff4[index++]=Ic[i][JMAX+2][k].e3();
#endif
	  }
      par->sbuflen[4]=index++;
      //left J=JMIN
      index=0;
      Scalar* buff5=par->ptsBuffer(5);
      for(int i=IMIN-2;i<=IMAX+2;i++)
        for(int k=DAMPKMIN-2;k<=DAMPKMAX+2;k++)
	  {
	    buff5[index++]=B[i][JMIN][k].e1();
	    buff5[index++]=B[i][JMIN][k].e2();
	    buff5[index++]=B[i][JMIN][k].e3();
	    buff5[index++]=Ic[i][JMIN][k].e1();
	    buff5[index++]=Ic[i][JMIN][k].e2();
	    buff5[index++]=Ic[i][JMIN][k].e3();
	    buff5[index++]=Ic[i][JMIN-1][k].e1();
	    buff5[index++]=Ic[i][JMIN-1][k].e2();
	    buff5[index++]=Ic[i][JMIN-1][k].e3();
#ifdef HIGHORDERCLOUD
	    buff5[index++]=Ic[i][JMIN-2][k].e1();
	    buff5[index++]=Ic[i][JMIN-2][k].e2();
	    buff5[index++]=Ic[i][JMIN-2][k].e3();
#endif
	  }    
      par->sbuflen[5]=index++;
    }
    break;
  }
}

void FieldInfo::messageunpack1st(int DIRECTION)
{
  int index;
  switch(DIRECTION){
  case 2:
    {
      // front k=KMIN
      index=0;
      
      Scalar* buff0=par->ptrBuffer(0);
      if( KMIN != par->parameter->ZMIN){
	for(int i=IMIN-2;i<=IMAX+2;i++)
	  for(int j=JMIN-2;j<=JMAX+2;j++)  
	    {
	      B[i][j][KMIN-1].set_e1(buff0[index++]);
	      B[i][j][KMIN-1].set_e2(buff0[index++]);
	      B[i][j][KMIN-1].set_e3(buff0[index++]);
	      Ic[i][j][KMIN]+=Vector3(buff0[index++],0,0);
	      Ic[i][j][KMIN]+=Vector3(0,buff0[index++],0);
	      Ic[i][j][KMIN]+=Vector3(0,0,buff0[index++]);
	      Ic[i][j][KMIN+1]+=Vector3(buff0[index++],0,0);
	      Ic[i][j][KMIN+1]+=Vector3(0,buff0[index++],0);
	      Ic[i][j][KMIN+1]+=Vector3(0,0,buff0[index++]);
#ifdef HIGHORDERCLOUD
	      Ic[i][j][KMIN+2]+=Vector3(buff0[index++],0,0);
	      Ic[i][j][KMIN+2]+=Vector3(0,buff0[index++],0);
	      Ic[i][j][KMIN+2]+=Vector3(0,0,buff0[index++]);
#endif
	    }
      }
      //back k=KMAX
      index=0;
      Scalar* buff1=par->ptrBuffer(1);
      if( KMAX != par->parameter->ZMAX){
	for(int i=IMIN-2;i<=IMAX+2;i++)
	  for(int j=JMIN-2;j<=JMAX+2;j++)
	    {
	      B[i][j][KMAX].set_e1(buff1[index++]);
	      B[i][j][KMAX].set_e2(buff1[index++]);
	      B[i][j][KMAX].set_e3(buff1[index++]);
	      Ic[i][j][KMAX]+=Vector3(buff1[index++],0,0);
	      Ic[i][j][KMAX]+=Vector3(0,buff1[index++],0);
	      Ic[i][j][KMAX]+=Vector3(0,0,buff1[index++]);
	      Ic[i][j][KMAX-1]+=Vector3(buff1[index++],0,0);
	      Ic[i][j][KMAX-1]+=Vector3(0,buff1[index++],0);
	      Ic[i][j][KMAX-1]+=Vector3(0,0,buff1[index++]);
#ifdef HIGHORDERCLOUD
	      Ic[i][j][KMAX-2]+=Vector3(buff1[index++],0,0);
	      Ic[i][j][KMAX-2]+=Vector3(0,buff1[index++],0);
	      Ic[i][j][KMAX-2]+=Vector3(0,0,buff1[index++]);
#endif
	    }
      }    
    }
    break;
  case 0:
    {
      //down i=IMIN
      index=0;
      Scalar* buff2=par->ptrBuffer(2);
      for(int j=JMIN-2;j<=JMAX+2;j++)
	for(int k=DAMPKMIN-2;k<=DAMPKMAX+2;k++)
	  {
	    B[IMIN-1][j][k].set_e1(buff2[index++]);
	    B[IMIN-1][j][k].set_e2(buff2[index++]);
	    B[IMIN-1][j][k].set_e3(buff2[index++]);
	    Ic[IMIN][j][k]+=Vector3(buff2[index++],0,0);
	    Ic[IMIN][j][k]+=Vector3(0,buff2[index++],0);
	    Ic[IMIN][j][k]+=Vector3(0,0,buff2[index++]);
	    Ic[IMIN+1][j][k]+=Vector3(buff2[index++],0,0);
	    Ic[IMIN+1][j][k]+=Vector3(0,buff2[index++],0);
	    Ic[IMIN+1][j][k]+=Vector3(0,0,buff2[index++]);
#ifdef HIGHORDERCLOUD
	    Ic[IMIN+2][j][k]+=Vector3(buff2[index++],0,0);
	    Ic[IMIN+2][j][k]+=Vector3(0,buff2[index++],0);
	    Ic[IMIN+2][j][k]+=Vector3(0,0,buff2[index++]);
#endif
	  }     
      
      //up i=IMAX; seceive from top region i=IMIN
      index=0;
      Scalar* buff3=par->ptrBuffer(3);
      for(int j=JMIN-2;j<=JMAX+2;j++)
        for(int k=DAMPKMIN-2;k<=DAMPKMAX+2;k++)
	  {
	    B[IMAX][j][k].set_e1(buff3[index++]);
	    B[IMAX][j][k].set_e2(buff3[index++]);
	    B[IMAX][j][k].set_e3(buff3[index++]);
	    Ic[IMAX][j][k]+=Vector3(buff3[index++],0,0);
	    Ic[IMAX][j][k]+=Vector3(0,buff3[index++],0);
	    Ic[IMAX][j][k]+=Vector3(0,0,buff3[index++]);
	    Ic[IMAX-1][j][k]+=Vector3(buff3[index++],0,0);
	    Ic[IMAX-1][j][k]+=Vector3(0,buff3[index++],0);
	    Ic[IMAX-1][j][k]+=Vector3(0,0,buff3[index++]);
#ifdef HIGHORDERCLOUD
	    Ic[IMAX-2][j][k]+=Vector3(buff3[index++],0,0);
	    Ic[IMAX-2][j][k]+=Vector3(0,buff3[index++],0);
	    Ic[IMAX-2][j][k]+=Vector3(0,0,buff3[index++]);
#endif
	  } 
    }
    break;
  case 1:
    { 
      //left J=JMIN
      index=0;
      Scalar* buff4=par->ptrBuffer(4);
      for(int i=IMIN-2;i<=IMAX+2;i++)
        for(int k=DAMPKMIN-2;k<=DAMPKMAX+2;k++)
	  {
	    B[i][JMIN-1][k].set_e1(buff4[index++]);
	    B[i][JMIN-1][k].set_e2(buff4[index++]);
	    B[i][JMIN-1][k].set_e3(buff4[index++]);
	    Ic[i][JMIN][k]+=Vector3(buff4[index++],0,0);
	    Ic[i][JMIN][k]+=Vector3(0,buff4[index++],0);
	    Ic[i][JMIN][k]+=Vector3(0,0,buff4[index++]);
	    Ic[i][JMIN+1][k]+=Vector3(buff4[index++],0,0);
	    Ic[i][JMIN+1][k]+=Vector3(0,buff4[index++],0);
	    Ic[i][JMIN+1][k]+=Vector3(0,0,buff4[index++]);
#ifdef HIGHORDERCLOUD
	    Ic[i][JMIN+2][k]+=Vector3(buff4[index++],0,0);
	    Ic[i][JMIN+2][k]+=Vector3(0,buff4[index++],0);
	    Ic[i][JMIN+2][k]+=Vector3(0,0,buff4[index++]);
#endif
	  }
      
      //right J=JMAX
      index=0;
      Scalar* buff5=par->ptrBuffer(5);
      for(int i=IMIN-2;i<=IMAX+2;i++)
        for(int k=DAMPKMIN-2;k<=DAMPKMAX+2;k++)
	  {
	    B[i][JMAX][k].set_e1(buff5[index++]);
	    B[i][JMAX][k].set_e2(buff5[index++]);
	    B[i][JMAX][k].set_e3(buff5[index++]);
	    Ic[i][JMAX][k]+=Vector3(buff5[index++],0,0);
	    Ic[i][JMAX][k]+=Vector3(0,buff5[index++],0);
	    Ic[i][JMAX][k]+=Vector3(0,0,buff5[index++]);
	    Ic[i][JMAX-1][k]+=Vector3(buff5[index++],0,0);
	    Ic[i][JMAX-1][k]+=Vector3(0,buff5[index++],0);
	    Ic[i][JMAX-1][k]+=Vector3(0,0,buff5[index++]);
#ifdef HIGHORDERCLOUD
	    Ic[i][JMAX-2][k]+=Vector3(buff5[index++],0,0);
	    Ic[i][JMAX-2][k]+=Vector3(0,buff5[index++],0);
	    Ic[i][JMAX-2][k]+=Vector3(0,0,buff5[index++]);
#endif
	  } 
    }
    break;
  }                                                                                               
}


void FieldInfo::messagepack2nd(int DIRECTION)
{
  int index;
  switch(DIRECTION){
  case 2: 
    {
      // front k=KMAX
      index=0;
      Scalar* buff0=par->ptsBuffer(0);
      for(int i=IMIN-2;i<=IMAX+2;i++)
        for(int j=JMIN-2;j<=JMAX+2;j++)
	  {
	    buff0[index++]=B[i][j][KMAX-1].e1();
	    buff0[index++]=B[i][j][KMAX-1].e2();
	    buff0[index++]=B[i][j][KMAX-1].e3();
	    buff0[index++]=E[i][j][KMAX-1].e1();
	    buff0[index++]=E[i][j][KMAX-1].e2();
	    buff0[index++]=E[i][j][KMAX-1].e3();

	    buff0[index++]=B[i][j][KMAX-2].e1();
	    buff0[index++]=B[i][j][KMAX-2].e2();
	    buff0[index++]=B[i][j][KMAX-2].e3();
	    buff0[index++]=E[i][j][KMAX-2].e1();
	    buff0[index++]=E[i][j][KMAX-2].e2();
	    buff0[index++]=E[i][j][KMAX-2].e3();
	  }
      par->sbuflen[0]=index++;
      
      //back k=KMIN
      index=0;
      Scalar* buff1=par->ptsBuffer(1);
      for(int i=IMIN-2;i<=IMAX+2;i++)
        for(int j=JMIN-2;j<=JMAX+2;j++)
	  {
	    buff1[index++]=B[i][j][KMIN].e1();
	    buff1[index++]=B[i][j][KMIN].e2();
	    buff1[index++]=B[i][j][KMIN].e3();
	    buff1[index++]=E[i][j][KMIN].e1();
	    buff1[index++]=E[i][j][KMIN].e2();
	    buff1[index++]=E[i][j][KMIN].e3();

	    buff1[index++]=B[i][j][KMIN+1].e1();
	    buff1[index++]=B[i][j][KMIN+1].e2();
	    buff1[index++]=B[i][j][KMIN+1].e3();
	    buff1[index++]=E[i][j][KMIN+1].e1();
	    buff1[index++]=E[i][j][KMIN+1].e2();
	    buff1[index++]=E[i][j][KMIN+1].e3();
	  }    
      par->sbuflen[1]=index++;
    }
    break;
  case 0:
    {
      //up i=IMAX
      index=0;
      Scalar* buff2=par->ptsBuffer(2);
      for(int j=JMIN-2 ;j<=JMAX+2 ;j++)
        for(int k=DAMPKMIN-2;k<=DAMPKMAX+2;k++)
	  {
	    buff2[index++]=B[IMAX-1][j][k].e1();
	    buff2[index++]=B[IMAX-1][j][k].e2();
	    buff2[index++]=B[IMAX-1][j][k].e3();
	    buff2[index++]=E[IMAX-1][j][k].e1();
	    buff2[index++]=E[IMAX-1][j][k].e2();
	    buff2[index++]=E[IMAX-1][j][k].e3();

	    buff2[index++]=B[IMAX-2][j][k].e1();
	    buff2[index++]=B[IMAX-2][j][k].e2();
	    buff2[index++]=B[IMAX-2][j][k].e3();
	    buff2[index++]=E[IMAX-2][j][k].e1();
	    buff2[index++]=E[IMAX-2][j][k].e2();
	    buff2[index++]=E[IMAX-2][j][k].e3();
	  }     
      par->sbuflen[2]=index++;
    
      //down i=IMIN send to down region i=IMAX
      index=0;
      Scalar* buff3=par->ptsBuffer(3);
      for(int j=JMIN-2;j<=JMAX+2;j++)
        for(int k=DAMPKMIN-2;k<=DAMPKMAX+2;k++)
	  {
	    buff3[index++]=B[IMIN][j][k].e1();
	    buff3[index++]=B[IMIN][j][k].e2();
	    buff3[index++]=B[IMIN][j][k].e3();
	    buff3[index++]=E[IMIN][j][k].e1();
	    buff3[index++]=E[IMIN][j][k].e2();
	    buff3[index++]=E[IMIN][j][k].e3();

	    buff3[index++]=B[IMIN+1][j][k].e1();
	    buff3[index++]=B[IMIN+1][j][k].e2();
	    buff3[index++]=B[IMIN+1][j][k].e3();
	    buff3[index++]=E[IMIN+1][j][k].e1();
	    buff3[index++]=E[IMIN+1][j][k].e2();
	    buff3[index++]=E[IMIN+1][j][k].e3();
	  } 
      par->sbuflen[3]=index++;
    }
    break;
  case 1:
    {   
      //right J=JMAX
      index=0;
      Scalar* buff4=par->ptsBuffer(4);
      for(int i=IMIN-2;i<=IMAX+2;i++)
        for(int k=DAMPKMIN-2;k<=DAMPKMAX+2;k++)
	  {
	    buff4[index++]=B[i][JMAX-1][k].e1();
	    buff4[index++]=B[i][JMAX-1][k].e2();
	    buff4[index++]=B[i][JMAX-1][k].e3();
	    buff4[index++]=E[i][JMAX-1][k].e1();
	    buff4[index++]=E[i][JMAX-1][k].e2();
	    buff4[index++]=E[i][JMAX-1][k].e3();

	    buff4[index++]=B[i][JMAX-2][k].e1();
	    buff4[index++]=B[i][JMAX-2][k].e2();
	    buff4[index++]=B[i][JMAX-2][k].e3();
	    buff4[index++]=E[i][JMAX-2][k].e1();
	    buff4[index++]=E[i][JMAX-2][k].e2();
	    buff4[index++]=E[i][JMAX-2][k].e3();
	  }
      par->sbuflen[4]=index++;
      //left J=JMIN
      index=0;
      Scalar* buff5=par->ptsBuffer(5);
      for(int i=IMIN-2;i<=IMAX+2;i++)
        for(int k=DAMPKMIN-2;k<=DAMPKMAX+2;k++)
	  {
	    buff5[index++]=B[i][JMIN][k].e1();
	    buff5[index++]=B[i][JMIN][k].e2();
	    buff5[index++]=B[i][JMIN][k].e3();
	    buff5[index++]=E[i][JMIN][k].e1();
	    buff5[index++]=E[i][JMIN][k].e2();
	    buff5[index++]=E[i][JMIN][k].e3();

	    buff5[index++]=B[i][JMIN+1][k].e1();
	    buff5[index++]=B[i][JMIN+1][k].e2();
	    buff5[index++]=B[i][JMIN+1][k].e3();
	    buff5[index++]=E[i][JMIN+1][k].e1();
	    buff5[index++]=E[i][JMIN+1][k].e2();
	    buff5[index++]=E[i][JMIN+1][k].e3();
	  }         
    
      par->sbuflen[5]=index++;
    }
    break;
  }
}

void FieldInfo::messageunpack2nd(int DIRECTION)
{
  int index;
  switch(DIRECTION){
  case 2:
    {
      // front k=KMIN
      index=0;
      Scalar* buff0=par->ptrBuffer(0);
      if( KMIN != par->parameter->ZMIN){
	for(int i=IMIN-2;i<=IMAX+2;i++)
	  for(int j=JMIN-2;j<=JMAX+2;j++)
	    {
	      B[i][j][KMIN-1].set_e1(buff0[index++]);
	      B[i][j][KMIN-1].set_e2(buff0[index++]);
	      B[i][j][KMIN-1].set_e3(buff0[index++]);
	      E[i][j][KMIN-1].set_e1(buff0[index++]);
	      E[i][j][KMIN-1].set_e2(buff0[index++]);
	      E[i][j][KMIN-1].set_e3(buff0[index++]);

	      B[i][j][KMIN-2].set_e1(buff0[index++]);
	      B[i][j][KMIN-2].set_e2(buff0[index++]);
	      B[i][j][KMIN-2].set_e3(buff0[index++]);
	      E[i][j][KMIN-2].set_e1(buff0[index++]);
	      E[i][j][KMIN-2].set_e2(buff0[index++]);
	      E[i][j][KMIN-2].set_e3(buff0[index++]);
	    }
      }
      
      //back k=KMAX
      index=0;
      Scalar* buff1=par->ptrBuffer(1);
      if( KMAX != par->parameter->ZMAX){
	for(int i=IMIN-2;i<=IMAX+2;i++)
	  for(int j=JMIN-2;j<=JMAX+2;j++)
	    {
	      B[i][j][KMAX].set_e1(buff1[index++]);
	      B[i][j][KMAX].set_e2(buff1[index++]);
	      B[i][j][KMAX].set_e3(buff1[index++]);
	      E[i][j][KMAX].set_e1(buff1[index++]);
	      E[i][j][KMAX].set_e2(buff1[index++]);
	      E[i][j][KMAX].set_e3(buff1[index++]);

	      B[i][j][KMAX+1].set_e1(buff1[index++]);
	      B[i][j][KMAX+1].set_e2(buff1[index++]);
	      B[i][j][KMAX+1].set_e3(buff1[index++]);
	      E[i][j][KMAX+1].set_e1(buff1[index++]);
	      E[i][j][KMAX+1].set_e2(buff1[index++]);
	      E[i][j][KMAX+1].set_e3(buff1[index++]);
	    }
      }    
    }
    break;
  case 0:
    {
      //down i=IMIN
      index=0;
      Scalar* buff2=par->ptrBuffer(2);
      for(int j=JMIN-2;j<=JMAX+2;j++)
        for(int k=DAMPKMIN-2;k<=DAMPKMAX+2;k++)
	  {
	    B[IMIN-1][j][k].set_e1(buff2[index++]);
	    B[IMIN-1][j][k].set_e2(buff2[index++]);
	    B[IMIN-1][j][k].set_e3(buff2[index++]);
	    E[IMIN-1][j][k].set_e1(buff2[index++]);
	    E[IMIN-1][j][k].set_e2(buff2[index++]);
	    E[IMIN-1][j][k].set_e3(buff2[index++]);

	    B[IMIN-2][j][k].set_e1(buff2[index++]);
	    B[IMIN-2][j][k].set_e2(buff2[index++]);
	    B[IMIN-2][j][k].set_e3(buff2[index++]);
	    E[IMIN-2][j][k].set_e1(buff2[index++]);
	    E[IMIN-2][j][k].set_e2(buff2[index++]);
	    E[IMIN-2][j][k].set_e3(buff2[index++]);
	  }     
      
      //up i=IMAX; seceive from top region i=IMIN
      index=0;
      Scalar* buff3=par->ptrBuffer(3);
      for(int j=JMIN-2;j<=JMAX+2;j++)
        for(int k=DAMPKMIN-2;k<=DAMPKMAX+2;k++)
	  {
	    B[IMAX][j][k].set_e1(buff3[index++]);
	    B[IMAX][j][k].set_e2(buff3[index++]);
	    B[IMAX][j][k].set_e3(buff3[index++]);
	    E[IMAX][j][k].set_e1(buff3[index++]);
	    E[IMAX][j][k].set_e2(buff3[index++]);
	    E[IMAX][j][k].set_e3(buff3[index++]);

	    B[IMAX+1][j][k].set_e1(buff3[index++]);
	    B[IMAX+1][j][k].set_e2(buff3[index++]);
	    B[IMAX+1][j][k].set_e3(buff3[index++]);
	    E[IMAX+1][j][k].set_e1(buff3[index++]);
	    E[IMAX+1][j][k].set_e2(buff3[index++]);
	    E[IMAX+1][j][k].set_e3(buff3[index++]);
	  } 
    }
    break;
  case 1:
    {   
      //left J=JMIN
      index=0;
      Scalar* buff4=par->ptrBuffer(4);
      for(int i=IMIN-2;i<=IMAX+2;i++)
        for(int k=DAMPKMIN-2;k<=DAMPKMAX+2;k++)
	  {
	    B[i][JMIN-1][k].set_e1(buff4[index++]);
	    B[i][JMIN-1][k].set_e2(buff4[index++]);
	    B[i][JMIN-1][k].set_e3(buff4[index++]);
	    E[i][JMIN-1][k].set_e1(buff4[index++]);
	    E[i][JMIN-1][k].set_e2(buff4[index++]);
	    E[i][JMIN-1][k].set_e3(buff4[index++]);

	    B[i][JMIN-2][k].set_e1(buff4[index++]);
	    B[i][JMIN-2][k].set_e2(buff4[index++]);
	    B[i][JMIN-2][k].set_e3(buff4[index++]);
	    E[i][JMIN-2][k].set_e1(buff4[index++]);
	    E[i][JMIN-2][k].set_e2(buff4[index++]);
	    E[i][JMIN-2][k].set_e3(buff4[index++]);
	  }
      
      //right J=JMAX
      index=0;
      Scalar* buff5=par->ptrBuffer(5);
      for(int i=IMIN-2;i<=IMAX+2;i++)
        for(int k=DAMPKMIN-2;k<=DAMPKMAX+2;k++)
	  {
	    B[i][JMAX][k].set_e1(buff5[index++]);
	    B[i][JMAX][k].set_e2(buff5[index++]);
	    B[i][JMAX][k].set_e3(buff5[index++]);
	    E[i][JMAX][k].set_e1(buff5[index++]);
	    E[i][JMAX][k].set_e2(buff5[index++]);
	    E[i][JMAX][k].set_e3(buff5[index++]);

	    B[i][JMAX+1][k].set_e1(buff5[index++]);
	    B[i][JMAX+1][k].set_e2(buff5[index++]);
	    B[i][JMAX+1][k].set_e3(buff5[index++]);
	    E[i][JMAX+1][k].set_e1(buff5[index++]);
	    E[i][JMAX+1][k].set_e2(buff5[index++]);
	    E[i][JMAX+1][k].set_e3(buff5[index++]);
	  } 
    }
    break;
  }                                
}

Scalar FieldInfo::damp_function(int k,int DAMPLENGTH, Scalar DAMP_R)
{
  /*
  Scalar ZMAX = par->parameter->ZMAX;
  Scalar ZMIN = par->parameter->ZMIN;
  Scalar offset = fabs(k-0.5*ZMAX);
  if(offset <= 0.5*ZMAX)
    return 1.0;
  else if(offset <= 0.5*ZMAX+DAMPLENGTH)
    return 1.0-DAMP_R*DAMP_R*(offset-0.5*ZMAX)*(offset-0.5*ZMAX)/(DAMPLENGTH*DAMPLENGTH);
  else
    return 0.0;
  */
  Scalar ZMAX = par->parameter->ZMAX;
  Scalar ZMIN = par->parameter->ZMIN;
  return (
	  k<ZMIN-DAMPLENGTH    ? 0.0:
	  k<ZMIN               ? (1.0-(DAMP_R*DAMP_R*(k-ZMIN)*(k-ZMIN)/(DAMPLENGTH*DAMPLENGTH)) ):
	  k<ZMAX               ? 1.0:
	  k<ZMAX+DAMPLENGTH    ? (1.0-(DAMP_R*DAMP_R*(k-ZMAX)*(k-ZMAX)/(DAMPLENGTH*DAMPLENGTH)) ):
	  0.0
	  );
}

void FieldInfo::advanceB()
{
    int i,j,k;
    int I0,I1,J0,J1,K0,K1;
    I0 = IMIN; I1=IMAX;
    J0 = JMIN; J1=JMAX;
    K0 = DAMPKMIN; K1=DAMPKMAX;

    Scalar dtdx2=0.5*dtdx;
    //intenal point
    IJKLOOP
      //      B[i][j][k] = B[i][j][k]-dtdx2*ROATE_p(E,i,j,k);
      B[i][j][k] = ( B[i][j][k]-dtdx2*ROATE_p(E,i,j,k)*damp_function(k,DAMP_LENGTH,0.0) )*damp_function(k,DAMP_LENGTH,1.0);

    if(KMIN == par->parameter->ZMIN){
      for(i=I0-2; i<=I1+1; i++)
	for(j=J0-2; j<=J1+1; j++)
	  for(k=K0-2; k<K0; k++){
	    B[i][j][k].set_e1(0.0);
	    B[i][j][k].set_e2(0.0);
	    B[i][j][k].set_e3(0.0);
	  }
    }

    if(KMAX == par->parameter->ZMAX){
      for(i=I0-2; i<=I1+1; i++)
	for(j=J0-2; j<=J1+1; j++)
	  for(k=K1; k<K1+2; k++){
	    B[i][j][k].set_e1(0.0);
	    B[i][j][k].set_e2(0.0);
	    B[i][j][k].set_e3(0.0);
	  }
    }
}

void FieldInfo::advanceE()
{
    int i,j,k;
    int I0,I1,J0,J1,K0,K1;
    I0 = IMIN; I1=IMAX;
    J0 = JMIN; J1=JMAX;
    K0 = DAMPKMIN; K1=DAMPKMAX;
//    if(KMIN == par->parameter->ZMIN) K0=KMIN+1;
//    if(KMAX == par->parameter->ZMAX) K1=KMAX-1;
    if(DAMP == 1 && KMIN == par->parameter->ZMIN){
      for(i=I0; i<=I1; i++)
	for(j=J0; j<=J1; j++)
	  for(k=K0-2; k<K0; k++)
	    E[i][j][k] = (E[i][j][k]+dtdx*ROATE_n(B,i,j,k)-dt*TWOPI*Ic[i][j][k])*(1.0-DAMP_R*sqr(k*1.0/DAMP_REGION_F));
    }

    for(i=I0; i<=I1; i++)
      for(j=J0; j<=J1; j++)
	for(k=K0; k<=K1; k++)
	  //    	  E[i][j][k] = E[i][j][k]+dtdx*ROATE_n(B,i,j,k)-dt*TWOPI*Ic[i][j][k];
	  E[i][j][k] = ( E[i][j][k]+(dtdx*ROATE_n(B,i,j,k)-dt*TWOPI*Ic[i][j][k])*damp_function(k,DAMP_LENGTH,0.0) )*damp_function(k,DAMP_LENGTH,1.0);


    if(DAMP == 1 && KMAX == par->parameter->ZMAX){
      for(i=I0; i<=I1; i++)
	for(j=J0; j<=J1; j++)
	  for(k=K1+1; k<=KMAX+DAMP_REGION_B; k++)
	    E[i][j][k] = (E[i][j][k]+dtdx*ROATE_n(B,i,j,k)-dt*TWOPI*Ic[i][j][k])*(1-DAMP_R*sqr((k*1.0-KMAX)/DAMP_REGION_B));
    }
}

void FieldInfo::average(Scalar t)
{
  int i,j,k;

  if(KMIN == par->parameter->ZMIN){
    ALLIJLOOP{
      //absorbwave in KMIN information save
      E_oldest_min[i][j][0] = E_older_min[i][j][0];
      E_older_min[i][j][0]  = E[i][j][0];
      E_oldest_min[i][j][1] = E_older_min[i][j][1];
      E_older_min[i][j][1]  = E[i][j][1];
    }
  }
  if(KMAX == par->parameter->ZMAX){
    ALLIJLOOP{
      //absorbwave in KMAX information save
      E_oldest_max[i][j][0]=E_older_max[i][j][0];
      E_older_max[i][j][0] =E[i][j][KMAX];
      E_oldest_max[i][j][1]=E_older_max[i][j][1];
      E_older_max[i][j][1] =E[i][j][KMAX-1];
    }
  }

  for(i=IMIN-1;i<=IMAX+1;i++)
    for(j=JMIN-1;j<=JMAX+1;j++)
      for(k=KMIN-1;k<=KMAX+1;k++){
	ENode[i][j][k].set_e1(0.5*(E[i][j][k].e1()+E[i-1][j][k].e1()));
	ENode[i][j][k].set_e2(0.5*(E[i][j][k].e2()+E[i][j-1][k].e2()));
	ENode[i][j][k].set_e3(0.5*(E[i][j][k].e3()+E[i][j][k-1].e3()));

	BNode[i][j][k].set_e1(0.25*(B[i][j][k].e1()+B[i][j][k-1].e1() +
			            B[i][j-1][k].e1()+B[i][j-1][k-1].e1()));
	BNode[i][j][k].set_e2(0.25*(B[i][j][k].e2()+B[i][j][k-1].e2() +
			            B[i-1][j][k].e2()+B[i-1][j][k-1].e2()));
	BNode[i][j][k].set_e3(0.25*(B[i][j][k].e3()+B[i][j-1][k].e3() +
			            B[i-1][j][k].e3()+B[i-1][j-1][k].e3()));
      }
#ifdef DEBUGMOVE
  EStatic = Vector3(0,0,0);
  BStatic = Vector3(0.0,0.0,0.0);
  for(i=IMIN;i<=IMAX;i++)
    for(j=JMIN;j<=JMAX;j++)
      for(k=KMIN;k<=KMAX;k++){
	ENode[i][j][k] += EStatic;
	BNode[i][j][k] += BStatic;
      }
#endif
  int istep = int(t/dt);

  if(fmod(istep,2*par->parameter->cells_per_wl) == 1){
    ALLIJKLOOP{
      Eave[i][j][k] = ENode[i][j][k]/(2.0*par->parameter->cells_per_wl) ;
      Bave[i][j][k] = BNode[i][j][k]/(2.0*par->parameter->cells_per_wl) ;
    }
  }
  else{
    ALLIJKLOOP{
      Eave[i][j][k] += ENode[i][j][k]/(2.0*par->parameter->cells_per_wl) ;
      Bave[i][j][k] += BNode[i][j][k]/(2.0*par->parameter->cells_per_wl) ;
    }
  }

}

void FieldInfo::clearQI()
{
  for(int i=IMIN-2;i<=IMAX+2;i++)
    for(int j=JMIN-2;j<=JMAX+2;j++)
      for(int k=KMIN-2;k<=KMAX+2;k++){
          Ic[i][j][k] = Vector3(0.0,0.0,0.0);
	  for(int isp=0;isp<NSPECIES;isp++)
	    Dens[isp][i][j][k] = 0.0;
      }
}

void FieldInfo::shiftField(int direction)
{
    int i,j,k;
    int I0,I1,J0,J1,K0,K1;

    I0 = IMIN-2; I1=IMAX+2;
    J0 = JMIN-2; J1=JMAX+2;
    K0 = KMIN-2; K1=KMAX+1;
    //copy k+1 cell infomation to k cell

    if(KMAX == par->parameter->ZMAX){
      for(i=I0; i<=I1; i++)
        for(j=J0; j<=J1; j++){
  	  E[i][j][KMAX].set_e1(0.0);
	  E[i][j][KMAX].set_e2(0.0);
	  E[i][j][KMAX].set_e3(0.0);
  	}
    }

    for(i=I0; i<=I1; i++)
      for(j=J0; j<=J1; j++)
        for(k=K0; k<=K1; k++){
	  E[i][j][k].set_e1(E[i][j][k+1].e1());
	  E[i][j][k].set_e2(E[i][j][k+1].e2());
	  E[i][j][k].set_e3(E[i][j][k+1].e3());
	  B[i][j][k].set_e1(B[i][j][k+1].e1());
	  B[i][j][k].set_e2(B[i][j][k+1].e2());
	  B[i][j][k].set_e3(B[i][j][k+1].e3());
	}

  if(KMAX == par->parameter->ZMAX){
    for(i=I0; i<=I1; i++)
      for(j=J0; j<=J1; j++){
	  E[i][j][KMAX].set_e1(0.0);
	  E[i][j][KMAX].set_e2(0.0);
	  E[i][j][KMAX].set_e3(0.0);
	}
  }
}

Scalar FieldInfo::energy()
{
  int i,j,k;
  int k0 = KMIN;
  int k1 = KMAX;
  Scalar energy0 = 0.0;
  if(KMIN == par->parameter->ZMIN){
    k0 = KMIN+4;
  }
  if(KMAX == par->parameter->ZMAX){
    k1 = KMAX-4;
  }

  IJKLOOP{
    energy0 += 0.5*(E[i][j][k].energy()+B[i][j][k].energy());
  }
  return energy0;
}

Scalar FieldInfo::work_energy()
{
  int i,j,k;
  Scalar energy0 = 0.0;
  IJKLOOP{
    energy0 = energy0 + E[i][j][k]*Ic[i][j][k]*dt;
  }
  return energy0;
}

void FieldInfo::front_flux_change(Scalar t)
{
  int i,j,k;
  //  flux0 = 0.0;
  if(KMIN == par->parameter->ZMIN){
    IJLOOP{
      poynting_left_p += fmax(-0.5*(ENode[i][j][KMIN+5].e2()*BNode[i][j][KMIN+5].e1()) + 0.5*(ENode[i][j][KMIN+5].e1()*BNode[i][j][KMIN+5].e2()), 0.0);
      poynting_left_n +=fmin( -0.5*(ENode[i][j][KMIN+5].e2()*BNode[i][j][KMIN+5].e1()) + 0.5*(ENode[i][j][KMIN+5].e1()*BNode[i][j][KMIN+5].e2()) ,0.0);
      laser_flux_left += 0.5*(sqr(laser->Ey(i,j,KMIN+4,par,t))+sqr(laser->Ex(i,j,KMIN+4,par,t)));
    }
  }
  else{
    laser_flux_left = 0.0;
    poynting_left_p = 0.0;
    poynting_left_n = 0.0;
  }
}

void FieldInfo::back_flux_change(Scalar t)
{
  int i,j,k;
  if(KMAX == par->parameter->ZMAX){
    IJLOOP{
      poynting_right_p +=fmax( 0.5*(ENode[i][j][KMAX-5].e2()*BNode[i][j][KMAX-5].e1()) - 0.5*(ENode[i][j][KMAX-5].e1()*BNode[i][j][KMAX-5].e2()) ,0.0);
      poynting_right_n +=fmin( 0.5*(ENode[i][j][KMAX-5].e2()*BNode[i][j][KMAX-5].e1()) - 0.5*(ENode[i][j][KMAX-5].e1()*BNode[i][j][KMAX-5].e2()) ,0.0);
      laser_flux_right += 0.5*(sqr(laser->Ey(i,j,KMAX-4,par,t))+sqr(laser->Ex(i,j,KMAX-4,par,t)));
    }
  }
  else{
    laser_flux_right = 0.0;
    poynting_right_p = 0.0;
    poynting_right_n = 0.0;
  }
}

void FieldInfo::dump(FILE* file)
{
  Scalar ee1,ee2,ee3,eb1,eb2,eb3,ei1,ei2,ei3;

  Scalar temp1 = 1.1111; Scalar temp2 = 2.2222;
  fwrite(&temp1,sizeof(Scalar),1,file);
  for(int k=KMIN-2;k<=KMAX+1;k++)
    for(int j=JMIN-2;j<=JMAX+1;j++)
      for(int i=IMIN-2;i<=IMAX+1;i++){
	ee1 = E[i][j][k].e1();
	ee2 = E[i][j][k].e2();
	ee3 = E[i][j][k].e3();
	eb1 = B[i][j][k].e1();
	eb2 = B[i][j][k].e2();
	eb3 = B[i][j][k].e3();
	ei1 = Ic[i][j][k].e1();
	ei2 = Ic[i][j][k].e2();
	ei3 = Ic[i][j][k].e3();
	fwrite(&ee1,sizeof(Scalar),1,file);
	fwrite(&ee2,sizeof(Scalar),1,file);
	fwrite(&ee3,sizeof(Scalar),1,file);
	fwrite(&eb1,sizeof(Scalar),1,file);
	fwrite(&eb2,sizeof(Scalar),1,file);
	fwrite(&eb3,sizeof(Scalar),1,file);
	fwrite(&ei1,sizeof(Scalar),1,file);
	fwrite(&ei2,sizeof(Scalar),1,file);
	fwrite(&ei3,sizeof(Scalar),1,file);
      }

  //dump mur boundary data
  if(KMIN==par->parameter->ZMIN){
    for(int k=0;k<=1;k++)
      for(int j=JMIN-2;j<=JMAX+2;j++)
	for(int i=IMIN-2;i<=IMAX+2;i++){
	  ee1 = E_older_min[i][j][k].e1();
	  ee2 = E_older_min[i][j][k].e2();
	  ee3 = E_older_min[i][j][k].e3();
	  eb1 = E_oldest_min[i][j][k].e1();
	  eb2 = E_oldest_min[i][j][k].e2();
	  eb3 = E_oldest_min[i][j][k].e3();
	  fwrite(&ee1,sizeof(Scalar),1,file);
	  fwrite(&ee2,sizeof(Scalar),1,file);
	  fwrite(&ee3,sizeof(Scalar),1,file);
	  fwrite(&eb1,sizeof(Scalar),1,file);
	  fwrite(&eb2,sizeof(Scalar),1,file);
	  fwrite(&eb3,sizeof(Scalar),1,file);
	}
  }
  if(KMAX==par->parameter->ZMAX){
    for(int k=0;k<=1;k++)
      for(int j=JMIN-2;j<=JMAX+2;j++)
	for(int i=IMIN-2;i<=IMAX+2;i++){
	  ee1 = E_older_max[i][j][k].e1();
	  ee2 = E_older_max[i][j][k].e2();
	  ee3 = E_older_max[i][j][k].e3();
	  eb1 = E_oldest_max[i][j][k].e1();
	  eb2 = E_oldest_max[i][j][k].e2();
	  eb3 = E_oldest_max[i][j][k].e3();
	  fwrite(&ee1,sizeof(Scalar),1,file);
	  fwrite(&ee2,sizeof(Scalar),1,file);
	  fwrite(&ee3,sizeof(Scalar),1,file);
	  fwrite(&eb1,sizeof(Scalar),1,file);
	  fwrite(&eb2,sizeof(Scalar),1,file);
	  fwrite(&eb3,sizeof(Scalar),1,file);
	}
  }

  fwrite(&temp2,sizeof(Scalar),1,file);
}

void FieldInfo::restore(FILE* file)
{
  Scalar ee1,ee2,ee3,eb1,eb2,eb3,ei1,ei2,ei3;
  Scalar temp1,temp2;
  fread(&temp1,sizeof(Scalar),1,file);
  if(temp1 != 1.1111){
    printf("read E file start wrong\n");
    exit(0);
  }
  for(int k=KMIN-2;k<=KMAX+1;k++)
    for(int j=JMIN-2;j<=JMAX+1;j++)
      for(int i=IMIN-2;i<=IMAX+1;i++){
	fread(&ee1,sizeof(Scalar),1,file);
	fread(&ee2,sizeof(Scalar),1,file);
	fread(&ee3,sizeof(Scalar),1,file);
	fread(&eb1,sizeof(Scalar),1,file);
	fread(&eb2,sizeof(Scalar),1,file);
	fread(&eb3,sizeof(Scalar),1,file);
	fread(&ei1,sizeof(Scalar),1,file);
	fread(&ei2,sizeof(Scalar),1,file);
	fread(&ei3,sizeof(Scalar),1,file);
	E[i][j][k] = Vector3(ee1,ee2,ee3);
	B[i][j][k] = Vector3(eb1,eb2,eb3);
	Ic[i][j][k] = Vector3(ei1,ei2,ei3);
      }
  //restore mur boundary 
  //dump mur boundary data
  if(KMIN==par->parameter->ZMIN){
    for(int k=0;k<=1;k++)
      for(int j=JMIN-2;j<=JMAX+2;j++)
	for(int i=IMIN-2;i<=IMAX+2;i++){
	  fread(&ee1,sizeof(Scalar),1,file);
	  fread(&ee2,sizeof(Scalar),1,file);
	  fread(&ee3,sizeof(Scalar),1,file);
	  fread(&eb1,sizeof(Scalar),1,file);
	  fread(&eb2,sizeof(Scalar),1,file);
	  fread(&eb3,sizeof(Scalar),1,file);
	  E_older_min[i][j][k].set_e1(ee1);
	  E_older_min[i][j][k].set_e2(ee2);
	  E_older_min[i][j][k].set_e3(ee3);
	  E_oldest_min[i][j][k].set_e1(eb1);
	  E_oldest_min[i][j][k].set_e2(eb2);
	  E_oldest_min[i][j][k].set_e3(eb3);
	}
  }
  if(KMAX==par->parameter->ZMAX){
    for(int k=0;k<=1;k++)
      for(int j=JMIN-2;j<=JMAX+2;j++)
	for(int i=IMIN-2;i<=IMAX+2;i++){
	  fread(&ee1,sizeof(Scalar),1,file);
	  fread(&ee2,sizeof(Scalar),1,file);
	  fread(&ee3,sizeof(Scalar),1,file);
	  fread(&eb1,sizeof(Scalar),1,file);
	  fread(&eb2,sizeof(Scalar),1,file);
	  fread(&eb3,sizeof(Scalar),1,file);
	  E_older_max[i][j][k].set_e1(ee1);
	  E_older_max[i][j][k].set_e2(ee2);
	  E_older_max[i][j][k].set_e3(ee3);
	  E_oldest_max[i][j][k].set_e1(eb1);
	  E_oldest_max[i][j][k].set_e2(eb2);
	  E_oldest_max[i][j][k].set_e3(eb3);
	}
  }


  fread(&temp2,sizeof(Scalar),1,file);
  if(temp2 != 2.2222){
    printf("read E file end wrong\n");
    exit(0);
  }

}

ofstream&
operator<<(ofstream& os, const FieldInfo& f) 
{

  os<< "<?xml version=\"1.0\"?>" <<"\n";
os<< "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\"  compressor=\"vtkZLibDataCompressor\"> "<<"\n";
os<<"<StructuredGrid WholeExtent=\" "<< f.IMIN << " "<< f.IMAX<<" "<< f.JMIN<<" " <<f.JMAX<<" " <<f.KMIN<<" "<< f.KMAX<<" \"  > \n";
os<< " <Piece Extent=\" "<< f.IMIN << " "<< f.IMAX<<" "<< f.JMIN<<" " <<f.JMAX<<" " <<f.KMIN<<" "<< f.KMAX <<" \"  > \n";
 os<<" <PointData> \n";
os<< " <DataArray type=\"Float64\" Name=\"energy\" format=\"ascii\"> "<<" \n";
  for(int k=f.KMIN;k<=f.KMAX;k++)
    for(int j=f.JMIN;j<=f.JMAX;j++)
      for(int i=f.IMIN;i<=f.IMAX;i++){
	os<<f.ENode[i][j][k].e1()<<" ";
	 
      }
os<<"\n";
os<< "</DataArray>"<<"\n";
os<< " <DataArray type=\"Float64\" Name=\"Dens0\" format=\"ascii\"> "<<" \n";
  for(int k=f.KMIN;k<=f.KMAX;k++)
    for(int j=f.JMIN;j<=f.JMAX;j++)
      for(int i=f.IMIN;i<=f.IMAX;i++){
	os<<f.ENode[i][j][k].e3()<<" ";
	 
      }
os<<"\n";
os<< "</DataArray>"<<"\n";
os<< " <DataArray type=\"Float64\" Name=\"Dens1\" format=\"ascii\"> "<<" \n";
  for(int k=f.KMIN;k<=f.KMAX;k++)
    for(int j=f.JMIN;j<=f.JMAX;j++)
      for(int i=f.IMIN;i<=f.IMAX;i++){
//	os<<f.Dens[0][i][j][k]<<" ";
	 
      }
os<<"\n";
os<< "</DataArray>"<<"\n";
os<< "</PointData>"<<"\n"<<" <CellData>"<<"\n"<<"  </CellData>"<<"\n";
      os<<" <Points> "<<"\n";
      os<<"  <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\" > "<<"\n";
for(int k=f.KMIN;k<=f.KMAX;k++)
    for(int j=f.JMIN;j<=f.JMAX;j++)
      for(int i=f.IMIN;i<=f.IMAX;i++){
	os<<i<<" "<<j<<" "<<k<<" 	";
	}
    os<< " </DataArray>"<<"\n";
   os<< " </Points> "<<"\n";
  

os<< " </Piece> "<<"\n";
os<< "</StructuredGrid>"<<"\n" << "</VTKFile>"<<"\n";
 		
    return os;
     
}

