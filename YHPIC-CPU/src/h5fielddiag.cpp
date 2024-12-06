#include "h5fielddiag.h"
#include "otypen.h"
#include "laser.h"
#include "fieldinfo.h"
#include "cellinfo.h"
#include <string.h>
#include <stdio.h>
#include "hdf5.h"

#ifdef PARALLEL
#include "mpi.h"
extern MPI_Comm OOPIC_COMM;
extern int MY_RANK;
extern int NP_X;
extern int NP_Y;
extern int NP_Z;
#endif

#ifdef PARALLEL_IO
#include "H5Part.h"
#include "H5Block.h"
#endif

H5FieldDiag::H5FieldDiag(char* input_file_name,ParParam* para,FieldInfo* f)
{
  char var_name[100];
  rf.openinput(input_file_name);

  // diagnostic region --------------------------------------
  Q           = atoi( rf.setget("&field", "Q") );
  if(Q!=0){   
      t_start     = atoi( rf.setget("&field", "t_start"  ) );
      t_stop      = atoi( rf.setget("&field", "t_stop"  ) );
      t_step      = atoi( rf.setget("&field", "t_step"  ) );
      
      // box ------------------------------------------------------
      cells_per_wl = atoi( rf.setget("&box", "cells_per_wl") );
      cells_axis   = atoi( rf.setget("&box", "cells_axis"  ) );
      cells_trans  = atoi( rf.setget("&box", "cells_trans" ) );  

      q_e = atoi( rf.setget("&field", "q_e") );
      q_b = atoi( rf.setget("&field", "q_b") );
      q_x = atoi( rf.setget("&field", "q_x") );
      q_y = atoi( rf.setget("&field", "q_y") );
  }
  
  //output ---------------------------------------------------------
  strcpy( path, rf.setget( "&output", "path" ) ); 

  t_save      = t_start*cells_per_wl*2;
  t_count     = 0;

  h5index     = 0;

  field       = f;

  imin        = para->IMIN;
  imax        = para->IMAX;
  jmin        = para->JMIN;
  jmax        = para->JMAX;
  kmin        = para->KMIN;
  kmax        = para->KMAX;
  XMIN        = para->parameter->XMIN;
  XMAX        = para->parameter->XMAX;
  YMIN        = para->parameter->YMIN;
  YMAX        = para->parameter->YMAX;
  ZMIN        = para->parameter->ZMIN;
  ZMAX        = para->parameter->ZMAX;

  rf.closeinput();    

  sprintf(h5filename,"%s/lpi-field.h5",path);

  h5partfile = H5PartOpenFileParallel(h5filename, H5PART_WRITE, OOPIC_COMM);

  H5PartCloseFile(h5partfile);


}

void H5FieldDiag::dump_field(int step)
{
    char pname[100];
    if(Q==0) return;   

    if(t_count >= t_start*cells_per_wl*2 &&
       t_count <=  t_stop*cells_per_wl*2    ){
        if(t_count == t_save){

	  h5partfile = H5PartOpenFileParallel(h5filename, H5PART_APPEND, OOPIC_COMM);
	  H5PartSetStep(h5partfile, int(t_save/(cells_per_wl*2.0)));
          cout<<"save at"<<t_count<<endl;
	  if(q_e==1) write_field_item("ex",0);
	  if(q_e==1) write_field_item("ey",1);
	  if(q_e==1) write_field_item("ez",2);
	  if(q_b==1) write_field_item("bx",3);
	  if(q_b==1) write_field_item("by",4);
	  if(q_b==1) write_field_item("bz",5);
	  if(q_x==1) write_field_item("xx",6);
	  if(q_x==1) write_field_item("xy",7);
	  if(q_x==1) write_field_item("xz",8);
	  if(q_y==1) write_field_item("yx",9);
	  if(q_y==1) write_field_item("yy",10);
	  if(q_y==1) write_field_item("yz",11);

	  write_density_item("d0",0);
          H5PartCloseFile(h5partfile);

	  t_save += t_step*cells_per_wl*2;
        }
    }
    t_count++;
}

void H5FieldDiag::write_field_item(char *name,int qt)
{
  int index,i,j,k;

  int IMIN = field->IMIN;
  int IMAX = field->IMAX;
  int JMIN = field->JMIN;
  int JMAX = field->JMAX;
  int KMIN = field->KMIN;
  int KMAX = field->KMAX;
  int XMIN = field->par->parameter->XMIN;
  int XMAX = field->par->parameter->XMAX;
  int YMIN = field->par->parameter->YMIN;
  int YMAX = field->par->parameter->YMAX;
  int ZMIN = field->par->parameter->ZMIN;
  int ZMAX = field->par->parameter->ZMAX;

  h5part_float32_t *data = new h5part_float32_t[(IMAX-IMIN)*(JMAX-JMIN)*(KMAX-KMIN)];
  
  switch(qt){
    case 0:
      index=0;    
      IJKLOOP
	{data[index] = (h5part_float32_t) field->ENode[i][j][k].e1(); index++;}
      break;
    case 1:
      index=0;    
      IJKLOOP
	{data[index] =  (h5part_float32_t) field->ENode[i][j][k].e2(); index++;}
      break;
    case 2:
      index=0;    
      IJKLOOP
	{data[index] =  (h5part_float32_t) field->ENode[i][j][k].e3(); index++;}
      break;
    case 3:
      index=0;    
      IJKLOOP
	{data[index] =  (h5part_float32_t) field->BNode[i][j][k].e1(); index++;}
      break;
    case 4:
      index=0;    
      IJKLOOP
	{data[index] =  (h5part_float32_t) field->BNode[i][j][k].e2(); index++;}
      break;
    case 5:
      index=0;    
      IJKLOOP
	{data[index] =  (h5part_float32_t) field->BNode[i][j][k].e3(); index++;}
      break;
    case 6:
      index=0;    
      IJKLOOP
	{data[index] =  (h5part_float32_t) field->Eave[i][j][k].e1(); index++;}
      break;
    case 7:
      index=0;    
      IJKLOOP
	{data[index] =  (h5part_float32_t) field->Eave[i][j][k].e2(); index++;}
      break;
    case 8:
      index=0;    
      IJKLOOP
	{data[index] =  (h5part_float32_t) field->Eave[i][j][k].e3(); index++;}
      break;
    case 9:
      index=0;    
      IJKLOOP
	{data[index] =  (h5part_float32_t) field->Bave[i][j][k].e1(); index++;}
      break;
    case 10:
      index=0;    
      IJKLOOP
	{data[index] =  (h5part_float32_t) field->Bave[i][j][k].e2(); index++;}
      break;
    case 11:
      index=0;    
      IJKLOOP
	{data[index] =  (h5part_float32_t) field->Bave[i][j][k].e3(); index++;}
      break;
  }
  
  H5BlockDefine3DFieldLayout(h5partfile, KMIN, KMAX-1, JMIN, JMAX-1, IMIN, IMAX-1);
  H5Block3dWriteScalarFieldFloat32(h5partfile, name, data);
  free(data);
}

void H5FieldDiag::write_density_item(char *name,int ispe)
{
  int index,i,j,k;

  int IMIN = field->IMIN;
  int IMAX = field->IMAX;
  int JMIN = field->JMIN;
  int JMAX = field->JMAX;
  int KMIN = field->KMIN;
  int KMAX = field->KMAX;

  h5part_float32_t *data = new h5part_float32_t[(IMAX-IMIN)*(JMAX-JMIN)*(KMAX-KMIN)];
  
  index=0;    
  IJKLOOP
      {data[index] = (h5part_float32_t) fabs(field->Dens[ispe][i][j][k]); index++;}
  
  H5BlockDefine3DFieldLayout(h5partfile, KMIN, KMAX-1, JMIN, JMAX-1, IMIN, IMAX-1);
  H5Block3dWriteScalarFieldFloat32(h5partfile, name, data);
  free(data);
}

