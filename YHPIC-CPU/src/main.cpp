#include <stdlib.h>
#include "laser.h"
#include "parparam.h"
#include "cellinfo.h"
#include "diagnostic.h"
#include "diagnostic_trace.h"
#include "mytime.h"
#ifdef PARALLEL
#include <mpi.h>

int MY_RANK;
int NPROC;
int NP_X;
int NP_Y;
int NP_Z;
MPI_Comm OOPIC_COMM;
MPI_Group MPI_GROUP_WORLD;
#endif

//#include "h5fielddiag.h"
//#include "h5particlediag.h"

Scalar init_time,end_time,dump_time;
Scalar cell_advance_time,cell_exchange_time,field_advance_time,field_exchange_time;

int main(int argc, char *argv[])
{
  Scalar start_t,init_end_t,run_end_t;
  Scalar ts,t0,t1,t2,t3,t4,t5,t6,t7,t8,t9;;
  cell_advance_time=cell_exchange_time=field_advance_time=field_exchange_time=0.0;

  char filename[100];
  if(argc==1)
    sprintf(filename,"%s","../input/lpi.inp");
  else
    sprintf(filename,"%s",argv[1]);

#ifdef PARALLEL
  MPI_Init(&argc,&argv);
#endif
  start_t = MPI_Wtime();         //get start time

  LaserParam* laser=new LaserParam(filename);
  SerParam*   spara=new SerParam(filename);
  ParParam*   ppara=new ParParam(spara,laser);

  ppara->initialize();

  FieldInfo* field=new FieldInfo(ppara);
  CellInfo*  cells=new CellInfo(field);

//  H5FieldDiag*  fdiag   = new H5FieldDiag(filename,ppara,field);
//  H5ParticleDiag*  pdiag   = new H5ParticleDiag(filename,ppara,cells);

  Diagnostic* diags=new Diagnostic(ppara);
  Diagnostic_trace* diag_trace=new Diagnostic_trace(filename,ppara);

  diags->diagnostic_restore(field,cells);

#ifdef PARALLEL
  if(MY_RANK == 0)
#endif
    diags->diagnostic_title();
  
  Scalar t=0.0; Scalar dt = field->get_dt();

  int dstep = int((ppara->parameter->prop_save+0.0001)/dt);
  int rstep = int((ppara->parameter->prop_restart+0.0001)/dt);
  int shift_step =  ppara->parameter->shift_step;
  int shift_dstep = ppara->parameter->shift_dstep;

  for(t=ppara->parameter->prop_start; t<=ppara->parameter->prop_stop+dt; t+=dt)
    {
      if(fmod(int(t/dt),10) == 0 && t!=0 ){
	diags->diagnostic_energy(field,cells,t);
      }

      int istep = int((t+0.0001)/dt); 

      if(fmod(istep,rstep) == 0 && t!=0 ){
	diags->diagnostic_dump(field,cells,t);
      }
      /*  advance field */
      t0 = MPI_Wtime();
      field->advance(t);
      t1 = MPI_Wtime();

      /*  parallel IO             */
      /*  anotation for benchmark */
      //fdiag->dump_field(0); 
      //pdiag->dump_particle(0);
      t2 = MPI_Wtime();
      
      /*  shift fields and particles */
      /*  anotaion for benchmark
      if(fmod(istep,shift_dstep) == 0 && istep>shift_step){
        ppara->parameter->shift_flag = 1;
        field->shiftField(1);
		cells->shift_particles(1);
      }
      */

      /*  advance particles  */
      cells->advance(t);
      cells->reflect_particles();
      cells->periodic_particles();
      field->exchangeDensity();
      t3 = MPI_Wtime();

      field_advance_time += (Scalar)(t1-t0);
      dump_time          += (Scalar)(t2-t1);
      cell_advance_time  += (Scalar)(t3-t2);

      diag_trace->store_traces(field,t);
      diag_trace->write_traces(field,t);
    }

    end_time = MPI_Wtime()-start_t;

#ifdef PARALLEL
  if(MY_RANK == 0)
#endif
    diags->diagnostic_time();

  delete laser;
  delete ppara;
  delete field;
  delete cells;
  delete diags;
//  delete fdiag;
//  delete pdiag;
  
#ifdef PARALLEL
  MPI_Finalize();
#endif
  return 0;
}
