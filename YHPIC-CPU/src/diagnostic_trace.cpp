#include "diagnostic_trace.h"
#ifdef PARALLEL
extern MPI_Comm OOPIC_COMM;
extern int MY_RANK;
#endif

Diagnostic_trace::Diagnostic_trace(char* input_file_name,ParParam* para)
{
  IMIN = para->IMIN;
  IMAX = para->IMAX;
  JMIN = para->JMIN;
  JMAX = para->JMAX;
  KMIN = para->KMIN;
  KMAX = para->KMAX;

  rf.openinput( input_file_name );

  Q         = atoi( rf.setget( "&traces", "Q" ) );
  t_start   = atoi( rf.setget( "&traces", "t_start" ) );
  t_stop    = atoi( rf.setget( "&traces", "t_stop" ) );
  t_step    = 2*para->parameter->cells_per_wl;
  ISTEP     = 0;
  traces            = atoi( rf.setget( "&traces", "traces" ) );
  if(Q==0) return;

  tracepos = new struct tracepositions[traces];
  for(int i=0;i<traces;i++) {
    sprintf(name,"t%d_x",i);        // trace positions expected in variables t0_x, t1_x, t2_x, ...
    tracepos[i].i = atoi( rf.setget( "&traces", name ) );
    sprintf(name,"t%d_y",i);        // trace positions expected in variables t0_y, t1_y, t2_y, ...
    tracepos[i].j = atoi( rf.setget( "&traces", name ) );
    sprintf(name,"t%d_z",i);        // trace positions expected in variables t0_z, t1_z, t2_z, ...
    tracepos[i].k = atoi( rf.setget( "&traces", name ) );
  }

  rf.closeinput();

  ex = new Scalar[traces*t_step];
  ey = new Scalar[traces*t_step];
  ez = new Scalar[traces*t_step];
  jx = new Scalar[traces*t_step];
  jy = new Scalar[traces*t_step];
  jz = new Scalar[traces*t_step];
}

void Diagnostic_trace::store_traces(FieldInfo* field,Scalar t)
{
  if(Q==0) return;
  if(t>t_start){
    for(int index=0;index<traces;index++){
      int i=tracepos[index].i;
      int j=tracepos[index].j;
      int k=tracepos[index].k;
      if(i>=IMIN&&i<IMAX&&j>=JMIN&&j<JMAX&&k>=KMIN&&k<KMAX)
	{
	  ex[index+ISTEP*traces]=field->ENode[i][j][k].e1();
	  ey[index+ISTEP*traces]=field->ENode[i][j][k].e2();
	  ez[index+ISTEP*traces]=field->ENode[i][j][k].e3();
	  jx[index+ISTEP*traces]=field->Ic[i][j][k].e1();
	  jy[index+ISTEP*traces]=field->Ic[i][j][k].e2();
	  jz[index+ISTEP*traces]=field->Ic[i][j][k].e3();
	}else
	{
	  ex[index+ISTEP*traces]=0.0;
	  ey[index+ISTEP*traces]=0.0;
	  ez[index+ISTEP*traces]=0.0;
	  jx[index+ISTEP*traces]=0.0;
	  jy[index+ISTEP*traces]=0.0;
	  jz[index+ISTEP*traces]=0.0;
	}
    }
	  ISTEP++;
  }

}

void Diagnostic_trace::write_traces(FieldInfo* field,Scalar t)
{
  
  FILE* file;
  if(Q==0) return;
  if(ISTEP>=t_step){
#ifdef PARALLEL        
    Scalar *ex_t = new Scalar[traces*t_step];
    Scalar *ey_t = new Scalar[traces*t_step];
    Scalar *ez_t = new Scalar[traces*t_step];
    Scalar *jx_t = new Scalar[traces*t_step];
    Scalar *jy_t = new Scalar[traces*t_step];
    Scalar *jz_t = new Scalar[traces*t_step];
    MPI_Allreduce(ex,ex_t,traces*t_step,MPI_SCALAR,MPI_SUM,OOPIC_COMM);
    MPI_Allreduce(ey,ey_t,traces*t_step,MPI_SCALAR,MPI_SUM,OOPIC_COMM);
    MPI_Allreduce(ez,ez_t,traces*t_step,MPI_SCALAR,MPI_SUM,OOPIC_COMM);
    MPI_Allreduce(jx,jx_t,traces*t_step,MPI_SCALAR,MPI_SUM,OOPIC_COMM);
    MPI_Allreduce(jy,jy_t,traces*t_step,MPI_SCALAR,MPI_SUM,OOPIC_COMM);
    MPI_Allreduce(jz,jz_t,traces*t_step,MPI_SCALAR,MPI_SUM,OOPIC_COMM);
    if(MY_RANK == 0){
#endif
      sprintf(name,"%s/trace-%.3f", ".", t);
      file = fopen(name,"w");
#ifdef PARALLEL
      fwrite(ex_t,sizeof(Scalar),traces*t_step,file);
      fwrite(ey_t,sizeof(Scalar),traces*t_step,file);
      fwrite(ez_t,sizeof(Scalar),traces*t_step,file);
      fwrite(jx_t,sizeof(Scalar),traces*t_step,file);
      fwrite(jy_t,sizeof(Scalar),traces*t_step,file);
      fwrite(jz_t,sizeof(Scalar),traces*t_step,file);
#else
      fwrite(ex,sizeof(Scalar),traces*t_step,file);
      fwrite(ey,sizeof(Scalar),traces*t_step,file);
      fwrite(ez,sizeof(Scalar),traces*t_step,file);
      fwrite(jx,sizeof(Scalar),traces*t_step,file);
      fwrite(jy,sizeof(Scalar),traces*t_step,file);
      fwrite(jz,sizeof(Scalar),traces*t_step,file);
#endif
      fclose(file);
#ifdef PARALLEL
    }
#endif
    ISTEP=0;
#ifdef PARALLEL
    delete []ex_t;
    delete []ey_t;
    delete []ez_t;
    delete []jx_t;
    delete []jy_t;
    delete []jz_t;
#endif
  }
}

Diagnostic_trace::~Diagnostic_trace()
{
  if(Q==0) return;
    delete []ex;
    delete []ey;
    delete []ez;
    delete []jx;
    delete []jy;
    delete []jz;
}
