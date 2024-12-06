#include "h5particlediag.h"
#include "otypen.h"
#include "laser.h"
#include "cellinfo.h"
#include <string.h>
#include <stdio.h>

#ifdef PARALLEL
#include "mpi.h"
extern MPI_Comm OOPIC_COMM;
extern int MY_RANK;
extern int NP_X;
extern int NP_Y;
extern int NP_Z;
#endif

H5ParticleDiag::H5ParticleDiag(char* input_file_name,ParParam* para, CellInfo* _cell)
{
  char var_name[100];
  rf.openinput(input_file_name);

  // specie --------------------------------------------------
  for(int i=0;i<NSPECIES;i++) {
    specy[i].ze = para->parameter->specy[i].ze;
    specy[i].fix  = para->parameter->specy[i].fix;
    specy[i].m = para->parameter->specy[i].m;
  }

  
  // diagnostic region --------------------------------------
  Q           = atoi( rf.setget("&particle", "Q") );
  if(Q!=0){    
      t_start     = atof( rf.setget("&particle", "t_start"  ) );
      t_stop      = atof( rf.setget("&particle", "t_stop"  ) );
      t_step      = atof( rf.setget("&particle", "t_step"  ) );
      
      // box ------------------------------------------------------
      cells_per_wl = atoi( rf.setget("&box", "cells_per_wl") );
      cells_axis   = atoi( rf.setget("&box", "cells_axis"  ) );
      cells_trans  = atoi( rf.setget("&box", "cells_trans" ) );  
  }
  
  //output ---------------------------------------------------------
  strcpy( path, rf.setget( "&output", "path" ) ); 

  t_save      = t_start*cells_per_wl*2;
  t_count     = 0;
  rf.closeinput();    

  //initialize
  IMIN        = para->IMIN;
  IMAX        = para->IMAX;
  JMIN        = para->JMIN;
  JMAX        = para->JMAX;
  KMIN        = para->KMIN;
  KMAX        = para->KMAX;
  cells=_cell;

  //clear lpi-particle.h5
  sprintf(h5filename,"%s/lpi-particle.h5",path);
  h5partfile = H5PartOpenFileParallel(h5filename, H5PART_WRITE, OOPIC_COMM);
  H5PartCloseFile(h5partfile);

}

void H5ParticleDiag::dump_particle(int step)
{
    char pname[100];
    if(Q == 0) return;   

    if(t_count >= t_start*cells_per_wl*2 &&
       t_count <=  t_stop*cells_per_wl*2    ){
      if(t_count >= t_save){
	h5partfile = H5PartOpenFileParallel(h5filename, H5PART_APPEND, OOPIC_COMM);
	H5PartSetStep(h5partfile, int(t_save/(cells_per_wl*2.0)));
	write_particle_item("ex",0);
	H5PartCloseFile(h5partfile);
	t_save += t_step*cells_per_wl*2;
      }
    }
    t_count++;
}


void H5ParticleDiag::write_particle_item(char* name,int step)
{
  int i,j,k;
  particle* part;
  int num=0;
  int index=0;

  //tem array for data
  int *specy;
  float *x,*y,*z,*ux,*uy,*uz,*gamma;
 
  // comput particle number
  IJKLOOP{
      if(cells->CellsM[i][j][k].first!=NULL){
          for(part=cells->CellsM[i][j][k].first; part!=NULL;part=part->next){
              num++;
          }
      }
  }

  //allocate memory
  x  = new float[num];
  y  = new float[num];
  z  = new float[num];
  ux = new float[num];
  uy = new float[num];
  uz = new float[num];
  gamma  = new float[num];
  specy = new int[num];
  
  // fill particle array data
  IJKLOOP{
      if(cells->CellsM[i][j][k].first!=NULL){
          for(part=cells->CellsM[i][j][k].first; part!=NULL;part=part->next){
              x[index]   =(h5part_float32_t)part->x;
              y[index]   =(h5part_float32_t)part->y;
              z[index]   =(h5part_float32_t)part->z;
              ux[index]  =(h5part_float32_t)part->ux;
              uy[index]  =(h5part_float32_t)part->uy;
              uz[index]  =(h5part_float32_t)part->uz;
	      specy[index] = (h5part_int32_t)part->species;
              gamma[index]  =(h5part_float32_t)1.0/part->igamma;
              index++;
          }
      }
  }

  H5PartSetNumParticles(h5partfile, num);
  H5PartWriteDataFloat32(h5partfile,"x",x);
  H5PartWriteDataFloat32(h5partfile,"y",y);
  H5PartWriteDataFloat32(h5partfile,"z",z);      
  H5PartWriteDataFloat32(h5partfile,"ux",ux);
  H5PartWriteDataFloat32(h5partfile,"uy",uy);
  H5PartWriteDataFloat32(h5partfile,"uz",uz);      
  H5PartWriteDataFloat32(h5partfile,"gamma",gamma);
  H5PartWriteDataInt32(h5partfile,"specy",specy);

  delete[] x;
  delete[] y;
  delete[] z;
  delete[] ux;
  delete[] uy;
  delete[] uz;
  delete[] gamma;
  delete[] specy;
}

