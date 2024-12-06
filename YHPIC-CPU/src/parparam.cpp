#include "parparam.h"
#include "otypen.h"
#include "laser.h"
#include "fieldinfo.h"
#include <string.h>
#include <stdio.h>

#ifdef PARALLEL
extern  int MY_RANK;
extern  MPI_Comm OOPIC_COMM;
extern int NPROC;
extern int NP_X;
extern int NP_Y;
extern int NP_Z;
#endif

SerParam::SerParam(char* input_file_name)
{
  path = new char[100];
  rf.openinput(input_file_name);

  flag_restart = 0;
  // propagate -----------------------------------------------
  prop_start = atoi( rf.setget("&propagate", "prop_start") );
  prop_stop  = atoi( rf.setget("&propagate", "prop_stop") );
  prop_save  = atoi( rf.setget("&propagate", "prop_save") );
  prop_restart  = atoi( rf.setget("&propagate", "prop_restart") );
  if(prop_start > 0) flag_restart=1; 

  // parallel ------------------------------------------------
#ifdef PARALLEL
  NPROC = atoi( rf.setget("&parallel", "nproc") );
  NP_Z  = atoi( rf.setget("&parallel", "np_z" ) );
  NP_X  = atoi( rf.setget("&parallel", "np_x" ) );
  NP_Y  = atoi( rf.setget("&parallel", "np_y" ) );
#endif
  // box ------------------------------------------------------
  cells_per_wl = atoi( rf.setget("&box", "cells_per_wl") );
  cells_axis   = int( atof( rf.setget("&box", "cells_axis"  ) ) );
  cells_trans  = int( atof( rf.setget("&box", "cells_trans" ) ) );  
  cells_left   = int( atof( rf.setget("&box", "cells_left" ) )  ); 
  cells_gap    = int( atof( rf.setget("&box", "cells_gap" ) )   ); 
  cells_ramp   = int( atof( rf.setget("&box", "cells_ramp" ) )  );  
  cells_plasma = int( atof( rf.setget("&box", "cells_plasma" ) ));  
  n_ion_over_nc  = atof( rf.setget("&box", "n_ion_over_nc" ) );           
  n_ion_ramp     = atof( rf.setget("&box", "n_ion_ramp" ) );           
  angle        = atof( rf.setget("&pulse", "angle") );

  // ions ----------------------------------------------------------
  ppc[1]         = atoi( rf.setget( "&box", "i_ppc" ) );
  vtherm[1]      = atof( rf.setget( "&box", "i_vtherm" ) );

  specy[1].fix         = atoi( rf.setget( "&box", "i_fix" ) );
  specy[1].ze          = atoi( rf.setget( "&box", "i_z" ) );
  specy[1].m           = atoi( rf.setget( "&box", "i_m" ) );
  specy[1].zm = 1.0*specy[1].ze/specy[1].m;
  specy[1].n  = SCALE*SCALE*n_ion_over_nc / ppc[1];
  specy[1].zn = specy[1].n * specy[1].ze;

  // electrons -----------------------------------------------------
  ppc[0]         = atoi( rf.setget( "&box", "e_ppc" ) );
  vtherm[0]      = atof( rf.setget( "&box", "e_vtherm" ) );
  
  specy[0].fix         = atoi( rf.setget( "&box", "e_fix" ) );
  specy[0].ze          = -1;        // DEFAULT, SHOULD NOT BE CHANGED
  specy[0].m           = +1;        // DEFAULT, SHOULD NOT BE CHANGED
  specy[0].zm = 1.0*specy[0].ze/specy[0].m;
  specy[0].n  = SCALE*SCALE*fabs( 1.0 * specy[1].ze / specy[0].ze ) * n_ion_over_nc / ppc[0];
  specy[0].zn = specy[0].n * specy[0].ze;

  // shift  -------- -----------------------------------------------
  shift_step   = 2*atoi( rf.setget("&shift", "shift_strt") );
  shift_dstep  = 2*atoi( rf.setget("&shift", "shift_dx") );

  //output ---------------------------------------------------------
  strcpy( path, rf.setget( "&output", "path" ) ); 

  if ( specy[1].ze == 0 && ppc[0] > 0 ) { 
    cout << " WARNING     : You selected neutral atoms and free electrons" << endl;
    cout << "               Setting # MacroElectrons := 0" << endl;
    ppc[0] = 0;
  }

  n_el_over_nc = 1.0 * specy[1].ze * n_ion_over_nc;
  angle = angle*PI/180;

  prop_start   = atoi( rf.setget("&propagate", "prop_start"  ) );
  prop_stop    = atoi( rf.setget("&propagate", "prop_stop"  ) );
  if(prop_start*cells_per_wl > shift_step/2.0) 
    shift_flag   = 1;
  else
    shift_flag   = 0;

  XMIN=0; YMIN=0; ZMIN=0;
  XMAX=cells_trans;
  YMAX=cells_trans;
#ifdef DIM_TWO
  YMAX = 2;
#endif
  ZMAX=cells_axis;

  //// set in-out file name //////////////////////////////////////////////////////////////
  char outname[100];
  sprintf( outname, "%s/output.lpi", path );
#ifdef PARALLEL
  if(MY_RANK == 0){
#endif
      ofstream outfile(outname);
      outfile<<"&specy"<<"\n";
      outfile<<"     num = "<<2<<"\n";
      for(int i=0;i<2;i++){
          outfile<<"     zn"<<i<<" = "<<specy[i].zn<<"\n";
          outfile<<"      n"<<i<<" = "<<specy[i].n<<"\n";
      }
#ifdef PARALLEL
  }
#endif


  rf.closeinput();    
}

ParParam::ParParam(SerParam* spara,LaserParam* lpara)
{
	parameter=spara;
	laser=lpara;
}

void ParParam::initialize()
{	
  IMIN=parameter->XMIN;
  IMAX=parameter->XMAX;
  JMIN=parameter->YMIN;
  JMAX=parameter->YMAX;
  KMIN=parameter->ZMIN;
  KMAX=parameter->ZMAX;
  ID_xl=ID_xr=ID_yl=ID_yr=0;ID_zl=ID_zr=-1;
#ifdef PARALLEL
  ID_xl=ID_xr=ID_yl=ID_yr=ID_zl=ID_zr=-1;
  int nPROC,ndims[3],periods[3],mycoords[3];
  MPI_Comm_size(MPI_COMM_WORLD,&nPROC);
  if(nPROC != NPROC){
    printf("number of processes wrong");
    exit(0);
  }
  // creat cart coordinate
  ndims[0]=NP_X;
  ndims[1]=NP_Y;
  ndims[2]=NP_Z;
  periods[0]=0;
  periods[1]=0;
  periods[2]=0;

  MPI_Cart_create(MPI_COMM_WORLD,3,ndims,periods,true,&OOPIC_COMM);
  MPI_Comm_rank(OOPIC_COMM,&MY_RANK);
  MPI_Cart_coords(OOPIC_COMM,MY_RANK,3,mycoords);
  int mepx = mycoords[0];
  int mepy = mycoords[1];
  int mepz = mycoords[2];
  // get IMIN,IMAX,JMIN,JMAX,KMIN,KMAX
  int DXX = int((parameter->XMAX-parameter->XMIN)/NP_X);
  int DYY = int((parameter->YMAX-parameter->YMIN)/NP_Y);
  int DZZ = int((parameter->ZMAX-parameter->ZMIN)/NP_Z);

  IMIN = parameter->XMIN + mepx*DXX;
  IMAX = parameter->XMIN + (mepx+1)*DXX;
  JMIN = parameter->YMIN + mepy*DYY;
  JMAX = parameter->YMIN + (mepy+1)*DYY;
  KMIN = parameter->ZMIN + mepz*DZZ;
  KMAX = parameter->ZMIN + (mepz+1)*DZZ;

  if(mepx == (NP_X-1) && IMAX != parameter->XMAX)
    IMAX = parameter->XMAX;
  if(mepy == (NP_Y-1) && JMAX != parameter->YMAX)
    JMAX = parameter->YMAX;
  if(mepz == (NP_Z-1) && KMAX != parameter->ZMAX)
    KMAX = parameter->ZMAX;

  int DAMPKMIN = KMIN;
  int DAMPKMAX = KMAX;
  if(KMIN == parameter->ZMIN) DAMPKMIN = KMIN-DAMP_LENGTH;
  if(KMAX == parameter->ZMAX) DAMPKMAX = KMAX+DAMP_LENGTH;

  //get neiightbor rank
  int neib1[3];

  neib1[1] = mepy;
  neib1[2] = mepz;
  neib1[0] = mepx - 1;
  if( neib1[0] == -1) neib1[0] = NP_X-1;
  MPI_Cart_rank( OOPIC_COMM,neib1,&ID_xl);
  neib1[0] = mepx + 1;
  if( neib1[0] == NP_X) neib1[0] = 0;
  MPI_Cart_rank( OOPIC_COMM,neib1,&ID_xr);
  
  neib1[0] = mepx;
  neib1[2] = mepz;
  neib1[1] = mepy - 1;
  if( neib1[1] == -1) neib1[1] = NP_Y-1;
  MPI_Cart_rank( OOPIC_COMM,neib1,&ID_yl);
  neib1[1] = mepy + 1;
  if( neib1[1] == NP_Y) neib1[1] = 0;
  MPI_Cart_rank( OOPIC_COMM,neib1,&ID_yr);
  
  neib1[0] = mepx;
  neib1[1] = mepy;
  neib1[2] = mepz - 1;
  if( neib1[2] == -1) neib1[2] = NP_Z-1;
  MPI_Cart_rank( OOPIC_COMM,neib1,&ID_zl);
  neib1[2] = mepz + 1;
  if( neib1[2] == NP_Z) neib1[2] = 0;
  MPI_Cart_rank( OOPIC_COMM,neib1,&ID_zr);
  if(mepz == 0)      ID_zl = -1;
  if(mepz == NP_Z-1) ID_zr = -1;
//  if(mepx == 0)      ID_xl = -1; //add for test
//  if(mepx == NP_X-1) ID_xr = -1; //add for test
//  if(mepy == 0)      ID_yl = -1; //add for test
//  if(mepy == NP_Y-1) ID_yr = -1; //add for test

#endif
  //  sbufptlth = (IMAX-IMIN+2)*(JMAX-JMIN+2)*(KMAX-KMIN+2)*10*8;
  //  rbufptlth = (IMAX-IMIN+2)*(JMAX-JMIN+2)*(KMAX-KMIN+2)*10*8;
  //  sbufptlth = 2e8;
  //  rbufptlth = 2e8;
  //  sbufpt  = new Scalar[sbufptlth];
  //  rbufpt  = new Scalar[rbufptlth];
  //  memset(sbufpt,0,sbufptlth*sizeof(Scalar));
  //  memset(rbufpt,0,rbufptlth*sizeof(Scalar));
  sbufpt = NULL;
  rbufpt = NULL;
  rbuf[0] = new Scalar[2*12*(IMAX-IMIN+10)*(JMAX-JMIN+10)];
  rbuf[1] = new Scalar[2*12*(IMAX-IMIN+10)*(JMAX-JMIN+10)];
  rbuf[2] = new Scalar[2*12*(JMAX-JMIN+10)*(DAMPKMAX-DAMPKMIN+10)];
  rbuf[3] = new Scalar[2*12*(JMAX-JMIN+10)*(DAMPKMAX-DAMPKMIN+10)];
  rbuf[4] = new Scalar[2*12*(IMAX-IMIN+10)*(DAMPKMAX-DAMPKMIN+10)];
  rbuf[5] = new Scalar[2*12*(IMAX-IMIN+10)*(DAMPKMAX-DAMPKMIN+10)];
  
  sbuf[0] = new Scalar[2*12*(IMAX-IMIN+10)*(JMAX-JMIN+10)];
  sbuf[1] = new Scalar[2*12*(IMAX-IMIN+10)*(JMAX-JMIN+10)];
  sbuf[2] = new Scalar[2*12*(JMAX-JMIN+10)*(DAMPKMAX-DAMPKMIN+10)];
  sbuf[3] = new Scalar[2*12*(JMAX-JMIN+10)*(DAMPKMAX-DAMPKMIN+10)];
  sbuf[4] = new Scalar[2*12*(IMAX-IMIN+10)*(DAMPKMAX-DAMPKMIN+10)];
  sbuf[5] = new Scalar[2*12*(IMAX-IMIN+10)*(DAMPKMAX-DAMPKMIN+10)];

  for(int im=0;im<6;im++)
    sbuflen[im] = 0;
}

ParParam::~ParParam()
{
  for(int i=0;i<6;i++){
    delete[] sbuf[i];
    delete[] rbuf[i];
  }
  delete[] sbufpt;
  delete[] rbufpt;
  delete[] parameter->path;
  delete   parameter;
}

#ifdef PARALLEL

void ParParam::commField(int direction)
{
    int         i, j, k, z, gc, row;
    MPI_Status  status;
    switch(direction){
    case 2:{
      if(ID_zr!=-1) {
	MPI_Isend(sbuf[0],sbuflen[0],MPI_SCALAR,ID_zr,ECOMM+0,OOPIC_COMM,&srequest[0][0]); 
	MPI_Irecv(rbuf[1],sbuflen[0],MPI_SCALAR,ID_zr,ECOMM+1,OOPIC_COMM,&rrequest[1][0]);
      }
      if(ID_zl!=-1) {
	MPI_Isend(sbuf[1],sbuflen[1],MPI_SCALAR,ID_zl,ECOMM+1,OOPIC_COMM,&srequest[1][0]);
	MPI_Irecv(rbuf[0],sbuflen[1],MPI_SCALAR,ID_zl,ECOMM+0,OOPIC_COMM,&rrequest[0][0]);
      }
      if (ID_zr!=-1) {
	MPI_Wait(&srequest[0][0],&status);
	MPI_Wait(&rrequest[1][0],&status);
      }
      if(ID_zl!=-1) {
	MPI_Wait(&srequest[1][0],&status);
	MPI_Wait(&rrequest[0][0],&status);  
      }
    }
      break;
    case 0:{
      if(ID_xr!=-1) {
	MPI_Isend(sbuf[2],sbuflen[2],MPI_SCALAR,ID_xr,ECOMM+2,OOPIC_COMM,&srequest[2][0]);
	MPI_Irecv(rbuf[3],sbuflen[2],MPI_SCALAR,ID_xr,ECOMM+3,OOPIC_COMM,&rrequest[3][0]);
      }
      if(ID_xl!=-1) {
	MPI_Isend(sbuf[3],sbuflen[3],MPI_SCALAR,ID_xl,ECOMM+3,OOPIC_COMM,&srequest[3][0]);
	MPI_Irecv(rbuf[2],sbuflen[3],MPI_SCALAR,ID_xl,ECOMM+2,OOPIC_COMM,&rrequest[2][0]);
      }
      if(ID_xr!=-1) {
	MPI_Wait(&srequest[2][0],&status);
	MPI_Wait(&rrequest[3][0],&status);  
      }
      if(ID_xl!=-1) {
	MPI_Wait(&srequest[3][0],&status);
	MPI_Wait(&rrequest[2][0],&status);  
      }
    }
      break;
    case 1:{
      if(ID_yr!=-1) {
	MPI_Isend(sbuf[4],sbuflen[4],MPI_SCALAR,ID_yr,ECOMM+4,OOPIC_COMM,&srequest[4][0]);
	MPI_Irecv(rbuf[5],sbuflen[4],MPI_SCALAR,ID_yr,ECOMM+5,OOPIC_COMM,&rrequest[5][0]);
      }
      if(ID_yl!=-1) {
	MPI_Isend(sbuf[5],sbuflen[5],MPI_SCALAR,ID_yl,ECOMM+5,OOPIC_COMM,&srequest[5][0]);
	MPI_Irecv(rbuf[4],sbuflen[5],MPI_SCALAR,ID_yl,ECOMM+4,OOPIC_COMM,&rrequest[4][0]);
      }
      if(ID_yr!=-1) {
	MPI_Wait(&srequest[4][0],&status);
	MPI_Wait(&rrequest[5][0],&status);  
      }
      if(ID_yl!=-1) {
	MPI_Wait(&srequest[5][0],&status);
	MPI_Wait(&rrequest[4][0],&status);
      }
    }
      break;
    }
}

void ParParam::commParticle(int ID_s,int ID_r,int sendnum,int& recvnum)
{
  int rnum;
  int snum = sendnum;
  MPI_Status  status;
  int ID_sendto  = ID_s;
  int ID_rcvfrom = ID_r;
  if(ID_sendto  == -1) ID_sendto  = MPI_PROC_NULL;
  if(ID_rcvfrom == -1) ID_rcvfrom = MPI_PROC_NULL;

  MPI_Isend(&snum,1,MPI_INT,ID_sendto, PCOMM+100,OOPIC_COMM,&srequest[0][0]); 
  MPI_Irecv(&rnum,1,MPI_INT,ID_rcvfrom,PCOMM+100,OOPIC_COMM,&rrequest[1][0]);
  MPI_Wait(&rrequest[1][0],&status);
  MPI_Wait(&srequest[0][0],&status);
  recvnum = rnum;

  if(ID_rcvfrom == MPI_PROC_NULL) recvnum = 0;
  //  cout<<"num sendrecv finish"<<sendnum<<" "<<recvnum<<endl;

  sbufptlth = sendnum*8+20;
  rbufptlth = recvnum*8+20;
  
//  sbufpt = new double[sbufptlth];  //only for test! rm it for ok
  rbufpt = new double[rbufptlth];
  MPI_Isend(sbufpt,sbufptlth,MPI_DOUBLE,ID_sendto, PCOMM+101,OOPIC_COMM,&srequest[2][0]); 
  MPI_Irecv(rbufpt,rbufptlth,MPI_DOUBLE,ID_rcvfrom,PCOMM+101,OOPIC_COMM,&rrequest[3][0]);
  MPI_Wait(&rrequest[3][0],&status);
  MPI_Wait(&srequest[2][0],&status);
}
#else

void ParParam::commParticle(int ID_s,int ID_r,int sendnum,int& recvnum)
{
  if(sendnum==0) return;
  recvnum = sendnum;
  sbufptlth = sendnum*8+2;
  rbufptlth = recvnum*8+2;
  rbufpt = new double[sbufptlth];
  memcpy(rbufpt,sbufpt,sizeof(double)*sbufptlth);
}

void ParParam::commField(int direction)
{
  switch(direction){
  case 2:{
    memset(rbuf[0],0    ,sizeof(Scalar)*sbuflen[0]);
    memset(rbuf[1],0    ,sizeof(Scalar)*sbuflen[1]);
  }
    break;
  case 0:{
    memcpy(rbuf[3],sbuf[3],sizeof(Scalar)*sbuflen[3]);
    memcpy(rbuf[2],sbuf[2],sizeof(Scalar)*sbuflen[2]);
  }
    break;
  case 1:{
    memcpy(rbuf[5],sbuf[5],sizeof(Scalar)*sbuflen[5]);
    memcpy(rbuf[4],sbuf[4],sizeof(Scalar)*sbuflen[4]);
  }
    break;
  }
}

#endif

//apply laser pulse in boundary condition
void ParParam::correctField(FieldInfo& field,Scalar t)
{
    // laser input correct
    if(KMIN==parameter->ZMIN && parameter->shift_flag == 0){
      laser->Murabsorbcorrect(0,field,this,t);
      laser->pulsecorrect(0,field,this,t);
    }
    //laser output correct
    if(KMAX==parameter->ZMAX && parameter->shift_flag == 0){
      laser->Murabsorbcorrect(1,field,this,t);
//      laser->pulsecorrect(1,field,this,t);
    }
    //IMAX boundary correct
    //    laser->boundarycorrect(0,field,this,t);

    //JMAX boundary correct
    //    laser->boundarycorrect(1,field,this,t);
}


