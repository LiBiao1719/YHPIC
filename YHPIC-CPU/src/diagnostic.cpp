#include "diagnostic.h"
#include "particle.h"
#include <stdlib.h>
#include <time.h>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "hdf5.h"

#ifdef PARALLEL
#include "mpi.h"
extern MPI_Comm OOPIC_COMM;
extern int MY_RANK;
extern int NP_X;
extern int NP_Y;
extern int NP_Z;
#endif
extern Scalar init_time,end_time,dump_time;
extern Scalar cell_advance_time,cell_exchange_time;
extern Scalar field_advance_time,field_exchange_time;

using namespace std;

void Diagnostic::diagnostic_title()
{
  printf("|====================================================================================|\n"); fflush(stdout);
  printf("|                                                                                    |\n"); fflush(stdout);
  printf("|            THREE DIMENSIONAL PARTICLE_IN_CELL CODE BENCHMARK                       |\n"); fflush(stdout);
  printf("|                                                                                    |\n"); fflush(stdout);
  printf("|====================================================================================|\n"); fflush(stdout);
  printf("|                                                                                    |\n"); fflush(stdout);
  printf("|%9s|%20s|%18s|%18s|%15s|\n","time","particle energy","field energy","flux_l energy","flux_r energy"); fflush(stdout);
  printf("|                                                                                    |\n"); fflush(stdout);
}

Diagnostic::Diagnostic(ParParam* param)
{
  para = param;
  IMIN = para->IMIN;
  IMAX = para->IMAX;
  JMIN = para->JMIN;
  JMAX = para->JMAX;
  KMIN = para->KMIN;
  KMAX = para->KMAX;

  sprintf(name_err,"%s/lpi-error.log",para->parameter->path);
  ofstream errfile(name_err,ios::out);

  errfile.precision( 3 );
  errfile.setf( ios::showpoint | ios::scientific );
    
  errfile << "#" << setw(11) << "time"
       << setw(12) << "particle_en" 
       << setw(12) << "fields_en"
       << setw(12) << "laser_left"
       << setw(12) << "kinetic" << endl;

  errfile.close();
}

/*
void Diagnostic::diagnostic_field(FieldInfo* field,Scalar time)
{
#ifdef PARALLEL
  sprintf(name,"%s/field-%d-%.3f", "./", MY_RANK, time);
#else
  sprintf(name,"%s/field-%.3f", "./", time);
#endif
  ofstream outfile(name);
  outfile<<(*field);
}
*/

void Diagnostic::diagnostic_field(FieldInfo* field,Scalar time)
{
int prop_stop = field->par->parameter->prop_stop;
int prop_save = field->par->parameter->prop_save;
 diag_save = prop_save;   
 diag_dt =   field->dt;
int n = 0;
n = prop_stop/prop_save;

#ifdef PARALLEL
  sprintf(name,"%s/result/piece-data/field-%d-%.3f.vts", ".", MY_RANK, time);
   sprintf( name_vts,"%.3f.vts",   time);
  if(MY_RANK==0){
  iTime = iTime + 1;
  }
#else
  sprintf(name,"%s/field-%.3f.vts", "./", time);
#endif
 //cout<<"write vts file begin!";ield-3.681.pvts
  ofstream outfile(name);
  if( !outfile )
  {
     cerr<<"error:Unable to open outfile!"<<"\n";
  }
  outfile<<(*field);
 // cout<<"write vts file OK!";
 // ofstream tfile("xgy");
  //tfile<<"This is a field diagnostic program!";
  /*
diagnostic_write_pvts();
return ;
*/

#ifdef PARALLEL
int  MY_SIZE;
MPI_Comm my_comm;
//my_comm = field->par->MYCOMM;
MPI_Comm_size(OOPIC_COMM,&MY_SIZE);


 if(MY_RANK==0 /* MY_RANK id is 0*/){
 
    
     /*write .pvts fille*/
     sprintf(name1,"%s/result/field-%.3f.pvts", ".", time);
     sprintf(name_pvts,"field-%.3f.pvts",  time);
     //  cout<<"Hello "<<name_pvts<<" hello";
     ofstream out_pvts_file(name1);
     if( !out_pvts_file )
     {
          cerr<<"error:Unable to open out_pvts_file!"<<"\n";
     }
     out_pvts_file<<"<?xml version=\"1.0\"?>"<<"\n"
     <<"<VTKFile type=\"PStructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">"<<"\n"
     <<" <PStructuredGrid WholeExtent=\" "<< field->par->parameter->XMIN << " "<< field->par->parameter->XMAX
     <<" "<<  field->par->parameter->YMIN << " "<< field->par->parameter->YMAX
     <<" "<< field->par->parameter->ZMIN << " "<< field->par->parameter->ZMAX <<" \" GhostLevel=\"0\">"<<"\n"
     <<"  <PPointData>"<<"\n"
     <<"  <PDataArray type=\"Float64\" Name=\"energy\"/>"<<"\n"
     <<"  <PDataArray type=\"Float64\" Name=\"Dens1\"/>"<<"\n"
     <<"  <PDataArray type=\"Float64\" Name=\"Dens0\"/>"<<"\n"
     <<" </PPointData>"<<"\n"
     <<" <PCellData>"<<"\n"
  	 <<" </PCellData>"<<"\n"
     <<"  <PPoints>"<<"\n"
     <<"   <PDataArray type=\"Float32\" NumberOfComponents=\"3\"/>"<<"\n"
     <<"  </PPoints>"<<"\n";
     
     /*
     MPI size ? 
     */
     /*create three arrays*/
     int	*XRANGE, *YRANGE, *ZRANGE;
     int	m,n,l,i,j,k;

     //m = (IMAX - IMIN)/NP_X;		//ÔöÁ¿  [0, m, 2m]
     //n = (JMAX - JMIN)/NP_Y;
     //l = (KMAX - KMIN)/NP_Z;
     XRANGE = new int[NP_X + 1];	//Õâ¸öÊý×é°üº¬ÁË¿ÉÄÜµÄÈ¡Öµ,Èç0,50,100,150µÈ
     YRANGE = new int[NP_Y + 1];
     ZRANGE = new int[NP_Z + 1];
     /*intialize the values of three arrays*/
     for(i = 0; i < NP_X; i++)
     	XRANGE[i] = 0 + i * (IMAX - IMIN);
     XRANGE[NP_X] = field->par->parameter->XMAX;
     for(j = 0; j < NP_Y; j++)
     	YRANGE[j] = 0 + j * (JMAX - JMIN);
     YRANGE[NP_Y] = field->par->parameter->YMAX;
     for(k = 0; k < NP_Z; k++)
     	ZRANGE[k] = 0 + k * (KMAX - KMIN);
     ZRANGE[NP_Z] = field->par->parameter->ZMAX;
     int tt=0;
      for(i = 0; i < NP_X; i++)
      	for(j = 0; j < NP_Y; j++)
      		for(k = 0; k < NP_Z; k++)
      			out_pvts_file<<" <Piece Extent=\""  
      			<<XRANGE[i]<<" "<<XRANGE[i+1]<<" "
      			<<YRANGE[j]<<" "<<YRANGE[j+1]<<" "
      			<<ZRANGE[k]<<" "<<ZRANGE[k+1]<<"\""<<"  "
      			<<"Source="<<"\""<<  "/piece-data/field-"<< tt++ <<"-"<< name_vts<<"\"/>"<<"\n";
      			
  /*    
     for(int i=0; i < MY_SIZE; i++)
     {
       out_pvts_file<<" <Piece Extent="<<0 5 0 10 0 10"Source="<<"\""<<  "field-"<< i <<"-"<< name_vts<<"\"/>"<<"\n";
      //
     		
     }*/
     out_pvts_file<<"</PStructuredGrid>"<<"\n"
     <<"</VTKFile>"<<"\n";
  
}
   #else
     ;
     #endif
     
    return ;
}


/*void Diagnostic::diagnostic_write_pvts()
{

}*/

void Diagnostic::diagnostic_write_pvd()
{
 
     /*write .pvd fille*/
     /*
     in the last to write .pvd file and
     How to diagnostic last time?
     */
        float  save_dt = diag_dt*diag_save;
       
   	  sprintf(name2,"%s/field.pvd", "./");
    	  ofstream ofile(name2);
    	  if( !ofile )
    	  {
    	       cerr<<"error:Unable to open ofile!"<<"\n";
    	  }
    	 ofile<<"<?xml version=\"1.0\"?>"<<"\n"
    	 <<"<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">"<<"\n"
    	 <<" <Collection>"<<"\n";
   	  for (int i = 0; i < iTime; i++){
	    sprintf(name_pvts,"field-%.3f.pvts",save_dt*(i+1));
       	   ofile<<" <DataSet timestep=\""<< i <<"\" group=\"LegaVTK0\" part=\"0\" file=\""<< name_pvts <<"\"/>"<<"\n";
      	}
     	ofile<<"</Collection>"<<"\n"
   	    <<"</VTKFile>"<<"\n";
   
}


void Diagnostic::diagnostic_particle(CellInfo* cells,Scalar time)
{
#ifdef PARALLEL
  sprintf(name,"%s/particle-%d-%.3f", "./", MY_RANK, time);
#else
  sprintf(name,"%s/particle-%.3f", "./", time);
#endif
  ofstream outfile(name);
  outfile<<(*cells);
}

void Diagnostic::diagnostic_velocity(int species,CellInfo* cells,Scalar time)
{
  int i,j,k;
  Scalar vx, vy, vz, v, absolut;
  int bvx, bvy, bvz, bv;
  FILE *file;
  vcut = 1.0;
  struct particle *part;
  
  for(int num=0;num<DIM;num++)
    x[num]=y[num]=z[num]=a[num]=0;

  IJKLOOP{
    for(part = cells->CellsM[i][j][k].first;part!=NULL;part=part->next){
      if (part->species == species) {
	vx = part->ux * part->igamma;
	vy = part->uy * part->igamma;
	vz = part->uz * part->igamma;
	
	absolut = sqrt( sqr(vx) + sqr(vy) + sqr(vz) );
	
	bvx = (int) FLOOR( 0.5 * DIM * (1.0 + vx/vcut) + 0.5 );
	bvy = (int) FLOOR( 0.5 * DIM * (1.0 + vy/vcut) + 0.5 );
	bvz = (int) FLOOR( 0.5 * DIM * (1.0 + vz/vcut) + 0.5 );
	bv  = (int) FLOOR( 0.5 * DIM * (1.0 + absolut/vcut) + 0.5 );
	
	if (bvx>=0 && bvx<=DIM) x[bvx]++;
	else { printf("velocity bin out of range\n"); exit(1); }
	if (bvy>=0 && bvy<=DIM) y[bvy]++;
	else { printf( "velocity bin out of range\n" ); exit(1); }
	if (bvz>=0 && bvz<=DIM) z[bvz]++;
	else { printf( "velocity bin out of range\n" ); exit(1); }
	if (bv>=0 && bv<=DIM)   a[bv]++;
	else  { printf( "velocity bin out of range\n" ); exit(1); }
	
      }
    }
  }
#ifdef PARALLEL
  sprintf(name,"%s/field-%d-%.3f", "./", MY_RANK, time);
#else
  sprintf(name,"%s/velocity-sp-%d-%.3f", "./", species, time); 
#endif
  file = fopen( name, "w" );
  
  for( i=0; i<=DIM; i++ ) {
    v = -vcut + 2.0*vcut*i/DIM;
    fprintf( file, "\n %.4e  %d %d %d %d", v, x[i], y[i], z[i], a[i] ); 
  }
  
  fclose( file );
}

void Diagnostic::diagnostic_energy(FieldInfo* fields,CellInfo* cells,Scalar t)
{
  Scalar energy[8];
  energy[0] = cells->particle_energy();
  energy[1] = fields->energy();
  energy[2] = fields->poynting_left_p;
  energy[3] = fields->poynting_left_n;
  energy[4] = fields->laser_flux_left;
  energy[5] = fields->poynting_right_p;
  energy[6] = fields->poynting_right_n;
  energy[7] = fields->laser_flux_right;

#ifdef PARALLEL
  Scalar temp_energy[8];
  for(int i=0;i<8;i++) temp_energy[i] = energy[i];
  MPI_Allreduce(temp_energy,energy,8,MPI_SCALAR,MPI_SUM,OOPIC_COMM);
  if(MY_RANK == 0){
#endif
    Scalar err_energy = ((energy[2]-energy[3] + energy[5]-energy[4]) - energy[0] - energy[1]);
    printf("%12.4f%17.4e%17.4e%17.4e%17.4e\n",t,energy[0],energy[1],energy[2]+energy[3],energy[5]+energy[6]);
#ifdef PARALLEL  
  }
  if(MY_RANK == 0){
#endif
  ofstream errfile(name_err,ios::app);
  errfile.precision( 3 );
  errfile.setf( ios::showpoint | ios::scientific );

  errfile << setw(11) << t 
              << setw(12) << energy[0] 
              << setw(12) << energy[1]
	      << setw(12) << energy[2]
	      << setw(12) << energy[3]
	      << setw(12) << energy[4]
	      << setw(12) << energy[5] << endl;

  errfile.close();
#ifdef PARALLEL  
  }
#endif
}

void Diagnostic::diagnostic_time()
{
  time_t time_now;
  time(&time_now);

  printf("                                                                         \n"); fflush(stdout);
  printf("|===================<< performance results >>===========================|\n"); fflush(stdout);
  printf("                                                                         \n"); fflush(stdout);
#ifdef PARALLEL
  printf("%30s%3d%3s%3d%3s%3d\n","nproc               =",NP_X,"*",NP_Y,"*",NP_Z); fflush(stdout);
#endif
  printf("%30s%15.2f\n","loop time           =",end_time);  fflush(stdout);
  printf("%30s%15.2f\n","dump time           =",dump_time); fflush(stdout);
  printf("%30s%15.2f\n","cell advance  time  =",cell_advance_time); fflush(stdout);
  printf("%30s%15.2f\n","field advance  time =",field_advance_time); fflush(stdout);
  printf("Time now is :%s\n",ctime(&time_now)); fflush(stdout);
}

void Diagnostic::diagnostic_dump(FieldInfo* fields,CellInfo* cells,Scalar t)
{
  char *path = para->parameter->path;
  char name[200];
#ifdef PARALLEL
  sprintf(name,"%s/lpi-restart-%.3f%s-%d", path, t, ".dmp", MY_RANK);
#else
  sprintf(name,"%s/lpi-restart-%.3f%s", path, t, ".dmp");
#endif
  FILE* file = fopen(name,"w");
  fields->dump(file);
  cells->dump(file);
  fclose(file);
}

void Diagnostic::diagnostic_dumpH5(FieldInfo* fields,CellInfo* cells,Scalar t)
{
  char *path = para->parameter->path;
  char fname[200];
  char pname[200];
#ifdef PARALLEL
  sprintf(fname,"%s/lpi-fields-%.3f-%d", path , t,MY_RANK);
  sprintf(pname,"%s/lpi-particles-%.3f-%d", path, t,MY_RANK);
#else
  sprintf(fname,"%s/lpi-fields-%.3f%s", "./", t, ".h5");
  sprintf(pname,"%s/lpi-particles-%.3f%s", "./", t, ".h5");
#endif
}

void Diagnostic::diagnostic_restore(FieldInfo* fields,CellInfo* cells)
{
  char *path = para->parameter->path;
  Scalar dt = fields->get_dt();
  Scalar t  = para->parameter->prop_start;

  if(t==0) return;
  char name[200];
#ifdef PARALLEL
  sprintf(name,"%s/lpi-restart-%.3f%s-%d", path, t, ".dmp", MY_RANK);
#else
  sprintf(name,"%s/lpi-restart-%.3f%s", path, t, ".dmp");
#endif
  FILE* file = fopen(name,"r");
  fields->restore(file);
  cells->restore(file);
  fclose(file);
}
