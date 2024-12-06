#include "laser.h"

LaserParam::LaserParam(char* input_file_name)
{
    char name[100];
    rf.openinput(input_file_name);
    Q         = atoi( rf.setget("&pulse", "Q"          ) );
    spotwidth = atof( rf.setget("&pulse", "spotwidth"  ) );
    amplitude = atof( rf.setget("&pulse", "amplitude"  ) );
    focus_z   = atoi( rf.setget("&pulse", "focus_z" ) ); 
    polarization    = atoi( rf.setget("&pulse", "polarization"   ) );
    shape           = atoi( rf.setget("&pulse", "shape"          ) );
    raise           = atoi( rf.setget("&pulse", "raise"          ) );
    duration        = atoi( rf.setget("&pulse", "duration"       ) );
    int cells_per_wl    = atoi( rf.setget("&box", "cells_per_wl") );

    Qtrain        = atoi( rf.setget("&pulse", "Qtrain"       ) );
    if(Qtrain==1){
	trains    = atoi( rf.setget("&pulse", "trains"       ) );
	if(trains>0){
	    time_delay = new Scalar[trains];
	    for(int i=1;i<=trains;i++) {
		sprintf(name,"t%d",i);        // train delay expected in variables t0, t1, t2, ...
		time_delay[i] = atoi( rf.setget( "&pulse", name ) );
	    }
	}
    }else{
	trains = 0;
    }
    
    rf.closeinput();    
    
    dx = 1.0/cells_per_wl;
    dt = 0.5*dx;
}

LaserParam::~LaserParam()
{
  if(Qtrain==1)
      if(trains>0)
	  delete []time_delay;
}

Scalar LaserParam::single_pulse(Scalar t)
{
  return (
	  shape == 1 ?
	  (
           t < 0? 0.0:
	   t < raise ?            ( t / raise ) :
	   t < duration - raise ? 1.0 :
	   t < duration ?         ( duration - t ) / raise :
	   0
	   ) :
	   shape == 2 ?
	   (
           t < 0? 0.0:
	    t < raise ? 	   sin( 0.5 * PI * t / raise ) :
	    t < duration - raise ? 1.0 :
	    t < duration ?         sin( 0.5 * PI * (duration - t) / raise ) :
	    0
	    ) :
       shape == 3 ?
       (
           t < 0? 0.0:
	     t < raise ?         sqr( sin( 0.5 * PI * t / raise ) ) :
	     t < duration - raise ? 1.0 :
	     t < duration ?      sqr( sin( 0.5 * PI * (duration - t) / raise ) ) :
	     0
	    ):
	   // assume shape == 4
        exp(  -(t-2.5*raise)*(t-2.5*raise)/(2*raise*raise)  ) 

     );
}

Scalar LaserParam::envelop(Scalar t)
{
    //pay attension : lpi.inp trains delay in periods
    Scalar value = single_pulse(t);
    for(int i=1;i<=trains;i++)
	value+=single_pulse(t-time_delay[i]);
    return value;
}

Scalar LaserParam::Ex(Scalar i,Scalar j,Scalar k,ParParam* p,Scalar t)
{
    int XMIN = p->parameter->XMIN;
    int XMAX = p->parameter->XMAX;
    int YMIN = p->parameter->YMIN;
    int YMAX = p->parameter->YMAX;
    
    Scalar XMIDDLE = XMIN+0.5*(XMAX-XMIN);
    Scalar YMIDDLE = YMIN+0.5*(YMAX-YMIN);

//    Scalar dx = 1.0/p->parameter->cells_per_wl;
    Scalar PHI0 = 0.0;
    Scalar w0 = spotwidth;
    Scalar x0 = XMIDDLE*dx;
    Scalar y0 = YMIDDLE*dx;
    Scalar z0 = focus_z;
    Scalar ZR = PI*w0*w0;
    Scalar zz = z0;
    Scalar xx = i*dx-x0;
    Scalar yy = j*dx-y0;
    Scalar RC = zz*(1+sqr(ZR/zz));
    Scalar wb = w0*sqrt(1+zz*zz/(ZR*ZR));
    Scalar phase_G = TWOPI*(xx*xx+yy*yy)/(2.0*RC)-atan(zz/ZR)-PHI0;
#ifdef DIM_TWO    
    yy = 0.0; 
#endif   
    if(Q == 1 || Q == 2){
	    if(polarization == 3)
	      return 0.7071*amplitude*envelop(t)*
                cos(TWOPI*t+phase_G)*
                exp( -(xx*xx)/(wb*wb))*
                exp( -(yy*yy)/(wb*wb));
	    else if(polarization == 1) 
	      return amplitude*envelop(t)*
                sin(TWOPI*t+phase_G)*
                exp( -(xx*xx)/(wb*wb))*
                exp( -(yy*yy)/(wb*wb));
            else return 0.0;
    }    
    else
      return 0  ;
}

Scalar LaserParam::Ey(Scalar i,Scalar j,Scalar k,ParParam* p,Scalar t)
{
    int XMIN = p->parameter->XMIN;
    int XMAX = p->parameter->XMAX;
    int YMIN = p->parameter->YMIN;
    int YMAX = p->parameter->YMAX;
    
    Scalar XMIDDLE = XMIN+0.5*(XMAX-XMIN);
    Scalar YMIDDLE = YMIN+0.5*(YMAX-YMIN);

//    Scalar dx = 1.0/p->parameter->cells_per_wl;
    Scalar PHI0 = 0.0;
    Scalar w0 = spotwidth;
    Scalar x0 = XMIDDLE*dx;
    Scalar y0 = YMIDDLE*dx;
    Scalar z0 = focus_z;
    Scalar ZR = PI*w0*w0;
    Scalar zz = z0;
    Scalar xx = i*dx-x0;
    Scalar yy = j*dx-y0;
    Scalar RC = zz*(1+sqr(ZR/zz));
    Scalar wb = w0*sqrt(1+zz*zz/(ZR*ZR));
    Scalar phase_G = TWOPI*(xx*xx+yy*yy)/(2.0*RC)-atan(zz/ZR)-PHI0;

#ifdef DIM_TWO    
    yy = 0.0; 
#endif   
    if(Q == 1 || Q == 2){
	    if(polarization == 3)
	      return 0.7071*amplitude*envelop(t)*
                sin(TWOPI*t+phase_G)*
                exp( -(xx*xx)/(wb*wb))*
                exp( -(yy*yy)/(wb*wb));
	    else if(polarization == 2) 
	      return amplitude*envelop(t)*
                sin(TWOPI*t+phase_G)*
                exp( -(xx*xx)/(wb*wb))*
                exp( -(yy*yy)/(wb*wb));
            else return 0.0;
    }    
    else
      return 0  ;
}

void LaserParam::pulsecorrect(int TYPE,FieldInfo& fields,ParParam* p,Scalar t)
{
  int i,j,k;
  int IMIN = p->IMIN;
  int IMAX = p->IMAX;
  int JMIN = p->JMIN;
  int JMAX = p->JMAX;
  
  int ZMIN = p->parameter->ZMIN;
  int ZMAX = p->parameter->ZMAX;
  
  Scalar dtdx=fields.dtdx;
  switch(TYPE){
  case 0:
    if(Q == 1){
      for(i=IMIN;i<=IMAX;i++)
	for(j=JMIN;j<=JMAX;j++){
	  Scalar exi      = Ex(i+0.5,j,ZMIN,p,t+0.5*dt);
	  Scalar eyi      = Ey(i,j+0.5,ZMIN,p,t+0.5*dt);
	  Scalar byi      = Ex(i+0.5,j,ZMIN-0.5,p,t+0.5*dt);
	  Scalar bxi      = -Ey(i,j+0.5,ZMIN-0.5,p,t+0.5*dt);
	  fields.E[i][j][ZMIN] +=  0.5*Vector3(byi,-bxi,0.0);
	  fields.B[i][j][ZMIN-1] += 0.5*Vector3(-eyi,exi,0.0);
	}
    }
    break;
  case 1:
    if(Q == 2){
      for(i=IMIN;i<=IMAX;i++)
	for(j=JMIN;j<=JMAX;j++){
	  Scalar exi      = Ex(i+0.5,j,ZMIN,p,t);
	  Scalar eyi      = Ey(i,j+0.5,ZMIN,p,t);
	  Scalar byi      = -Ex(i+0.5,j,ZMIN-0.5,p,t);
	  Scalar bxi      = +Ey(i,j+0.5,ZMIN-0.5,p,t);
	  fields.E[i][j][ZMAX] += 0.5* Vector3(-byi,+bxi,0.0);
	  fields.B[i][j][ZMAX-1] += 0.5* Vector3(eyi,-exi,0.0);
	}
    }

    break;
  }
}

void LaserParam::absorbcorrect(int TYPE,FieldInfo& fields,ParParam* p,Scalar t)
{
  int i,j,k;
  int IMIN = p->IMIN;
  int IMAX = p->IMAX;
  int JMIN = p->JMIN;
  int JMAX = p->JMAX;
  
  int ZMIN = p->parameter->ZMIN;
  int ZMAX = p->parameter->ZMAX;
  
  Scalar dtdx=fields.dtdx;
  Scalar a[14];

  switch(TYPE){
  case 0:
    for(i=IMIN;i<IMAX;i++)
      for(j=JMIN;j<JMAX;j++){
	/* set Ex value */
	a[1] = fields.E_oldest_min[i][j][1].e1();
	a[2] = fields.E[i][j][1].e1();
	a[3] = fields.E_oldest_min[i][j][0].e1();
	a[4] = fields.E_older_min[i][j][1].e1();
	a[5] = fields.E_older_min[i][j][0].e1();

	a[6] = fields.E_older_min[i+1][j][0].e1();
	a[7] = fields.E_older_min[i-1][j][0].e1();
	a[8] = fields.E_older_min[i][j+1][0].e1();
	a[9] = fields.E_older_min[i][j-1][0].e1();

	a[10] = fields.E_older_min[i+1][j][1].e1();
	a[11] = fields.E_older_min[i-1][j][1].e1();
	a[12] = fields.E_older_min[i][j+1][1].e1();
	a[13] = fields.E_older_min[i][j-1][1].e1();

#ifdef DIM_TWO
	a[0] = -a[1]+(0.5-1.0)/(0.5+1.0)*(a[2]+a[3])
	  +2.0/(0.5+1.0)*(a[4]+a[5])
	  +0.25/3.0*(a[6]-2*a[5]+a[7])
          +0.25/3.0*(a[10]-2*a[4]+a[11]);
#else
	a[0] = -a[1]+(0.5-1.0)/(0.5+1.0)*(a[2]+a[3])
	  -2.0*(0.5-1.0)*(a[4]+a[5])
	  +0.25/3.0*(a[6]+a[7]+a[8]+a[9]+a[10]+a[11]+a[12]+a[13]);
#endif
	fields.E[i][j][ZMIN].set_e1(a[0]);

	/* set Ey value */
	a[1] = fields.E_oldest_min[i][j][1].e2();
	a[2] = fields.E[i][j][1].e2();
	a[3] = fields.E_oldest_min[i][j][0].e2();
	a[4] = fields.E_older_min[i][j][1].e2();
	a[5] = fields.E_older_min[i][j][0].e2();

	a[6] = fields.E_older_min[i+1][j][0].e2();
	a[7] = fields.E_older_min[i-1][j][0].e2();
	a[8] = fields.E_older_min[i][j+1][0].e2();
	a[9] = fields.E_older_min[i][j-1][0].e2();

	a[10] = fields.E_older_min[i+1][j][1].e2();
	a[11] = fields.E_older_min[i-1][j][1].e2();
	a[12] = fields.E_older_min[i][j+1][1].e2();
	a[13] = fields.E_older_min[i][j-1][1].e2();

#ifdef DIM_TWO
	a[0] = -a[1]+(0.5-1.0)/(0.5+1.0)*(a[2]+a[3])
	  +2.0/(0.5+1.0)*(a[4]+a[5])
	  +0.25/3.0*(a[6]-2*a[5]+a[7])
          +0.25/3.0*(a[10]-2*a[4]+a[11]);
#else
	a[0] = -a[1]+(0.5-1.0)/(0.5+1.0)*(a[2]+a[3])
	  -2.0*(0.5-1.0)*(a[4]+a[5])
	  +0.25/3.0*(a[6]+a[7]+a[8]+a[9]+a[10]+a[11]+a[12]+a[13]);
#endif
	fields.E[i][j][ZMIN].set_e2(a[0]);

      }
    /*move to fieldinfo::average
    for(i=IMIN;i<=IMAX;i++)
      for(j=JMIN;j<=JMAX;j++)
	for(int k=0;k<=1;k++){
	  fields.E_oldest_min[i][j][k] = fields.E_older_min[i][j][k];
	  fields.E_older_min[i][j][k]  = fields.E[i][j][k];
	}
    */
    break;

  case 1:
    for(i=IMIN;i<IMAX;i++)
      for(j=JMIN;j<JMAX;j++){
	/* set Ex value */
	a[1] = fields.E_oldest_max[i][j][0].e1();
	a[2] = fields.E[i][j][ZMAX-1].e1();
	a[3] = fields.E_oldest_max[i][j][1].e1();
	a[4] = fields.E_older_max[i][j][0].e1();
	a[5] = fields.E_older_max[i][j][1].e1();

	a[6] = fields.E_older_max[i+1][j][1].e1();
	a[7] = fields.E_older_max[i-1][j][1].e1();
	a[8] = fields.E_older_max[i][j+1][1].e1();
	a[9] = fields.E_older_max[i][j-1][1].e1();

	a[10] = fields.E_older_max[i+1][j][0].e1();
	a[11] = fields.E_older_max[i-1][j][0].e1();
	a[12] = fields.E_older_max[i][j+1][0].e1();
	a[13] = fields.E_older_max[i][j-1][0].e1();

#ifdef DIM_TWO
	a[0] = -a[1]+(0.5-1.0)/(0.5+1.0)*(a[2]+a[3])
	  +2.0/(0.5+1.0)*(a[4]+a[5])
	  +0.25/3.0*(a[6]-2*a[5]+a[7])
          +0.25/3.0*(a[10]-2*a[4]+a[11]);
#else
	a[0] = -a[1]+(0.5-1.0)/(0.5+1.0)*(a[2]+a[3])
	  -2.0*(0.5-1.0)*(a[4]+a[5])
	  +0.25/3.0*(a[6]+a[7]+a[8]+a[9]+a[10]+a[11]+a[12]+a[13]);

//		a[0] = a[4]+(0.5-1)/(0.5+1)*(a[2]-a[5]);
#endif
	fields.E[i][j][ZMAX].set_e1(a[0]);

	
        /* set Ey value */
	a[1] = fields.E_oldest_max[i][j][0].e2();
	a[2] = fields.E[i][j][ZMAX-1].e2();
	a[3] = fields.E_oldest_max[i][j][1].e2();
	a[4] = fields.E_older_max[i][j][0].e2();
	a[5] = fields.E_older_max[i][j][1].e2();

	a[6] = fields.E_older_max[i+1][j][1].e2();
	a[7] = fields.E_older_max[i-1][j][1].e2();
	a[8] = fields.E_older_max[i][j+1][1].e2();
	a[9] = fields.E_older_max[i][j-1][1].e2();

	a[10] = fields.E_older_max[i+1][j][0].e2();
	a[11] = fields.E_older_max[i-1][j][0].e2();
	a[12] = fields.E_older_max[i][j+1][0].e2();
	a[13] = fields.E_older_max[i][j-1][0].e2();

#ifdef DIM_TWO
	a[0] = -a[1]+(0.5-1.0)/(0.5+1.0)*(a[2]+a[3])
	  +2.0/(0.5+1.0)*(a[4]+a[5])
	  +0.25/3.0*(a[6]-2*a[5]+a[7])
          +0.25/3.0*(a[10]-2*a[4]+a[11]);
#else
	a[0] = -a[1]+(0.5-1.0)/(0.5+1.0)*(a[2]+a[3])
	  -2.0*(0.5-1.0)*(a[4]+a[5])
	  +0.25/3.0*(a[6]+a[7]+a[8]+a[9]+a[10]+a[11]+a[12]+a[13]);
#endif
	fields.E[i][j][ZMAX].set_e2(a[0]);

      }
    /* move to fieldinfo::average
    for(i=IMIN-2;i<=IMAX+2;i++)
      for(j=JMIN-2;j<=JMAX+2;j++){
	  fields.E_oldest_max[i][j][0]=fields.E_older_max[i][j][0];
	  fields.E_older_max[i][j][0] =fields.E[i][j][ZMAX-1];
	  fields.E_oldest_max[i][j][1]=fields.E_older_max[i][j][1];
	  fields.E_older_max[i][j][1] =fields.E[i][j][ZMAX];
      }
    */
    break;
  }
}

void LaserParam::Murabsorbcorrect(int TYPE,FieldInfo& fields,ParParam* p,Scalar t)
{
  int i,j,k;
  int IMIN = p->IMIN;
  int IMAX = p->IMAX;
  int JMIN = p->JMIN;
  int JMAX = p->JMAX;
  
  int ZMIN = p->parameter->ZMIN;
  int ZMAX = p->parameter->ZMAX;
  
  Scalar dtdx=fields.dtdx;
  Scalar a[14];

  switch(TYPE){
  case 0:
    for(i=IMIN;i<IMAX;i++)
      for(j=JMIN;j<JMAX;j++){
	/* set Ex value */
	a[1] = fields.E_oldest_min[i][j][1].e1();

	a[2] = fields.E[i][j][ZMIN+1].e1();
	a[3] = fields.E_oldest_min[i][j][0].e1();

	a[4] = fields.E_older_min[i][j][0].e1();
	a[5] = fields.E_older_min[i][j][1].e1();

	a[6] = fields.E_older_min[i+1][j][0].e1();
	a[7] = fields.E_older_min[i-1][j][0].e1();
	a[8] = fields.E_older_min[i][j+1][0].e1();
	a[9] = fields.E_older_min[i][j-1][0].e1();

	a[10] = fields.E_older_min[i][j][0].e1();

	a[0] = -a[1]+(0.5-1.0)/(0.5+1.0)*(a[2]+a[3])
	  +2.0/(0.5+1.0)*(a[4]+a[5])
	  +0.25/(0.5+1.0)*(a[6]+a[7]+a[8]+a[9]-4.0*a[10]);

	fields.E[i][j][ZMIN].set_e1(a[0]);

	/* set Ey value */
	a[1] = fields.E_oldest_min[i][j][1].e2();

	a[2] = fields.E[i][j][ZMIN+1].e2();
	a[3] = fields.E_oldest_min[i][j][0].e2();

	a[4] = fields.E_older_min[i][j][0].e2();
	a[5] = fields.E_older_min[i][j][1].e2();

	a[6] = fields.E_older_min[i+1][j][0].e2();
	a[7] = fields.E_older_min[i-1][j][0].e2();
	a[8] = fields.E_older_min[i][j+1][0].e2();
	a[9] = fields.E_older_min[i][j-1][0].e2();

	a[10] = fields.E_older_min[i][j][0].e2();

	a[0] = -a[1]+(0.5-1.0)/(0.5+1.0)*(a[2]+a[3])
	  +2.0/(0.5+1.0)*(a[4]+a[5])
	  +0.25/(0.5+1.0)*(a[6]+a[7]+a[8]+a[9]-4.0*a[10]);

	fields.E[i][j][ZMIN].set_e2(a[0]);

      }
    /*move to fieldinfo::average
    for(i=IMIN;i<=IMAX;i++)
      for(j=JMIN;j<=JMAX;j++)
	for(int k=0;k<=1;k++){
	  fields.E_oldest_min[i][j][k] = fields.E_older_min[i][j][k];
	  fields.E_older_min[i][j][k]  = fields.E[i][j][k];
	}
    */
    break;

  case 1:
    for(i=IMIN;i<IMAX;i++)
      for(j=JMIN;j<JMAX;j++){
	/* set Ex value */
	a[1] = fields.E_oldest_max[i][j][1].e1();

	a[2] = fields.E[i][j][ZMAX-1].e1();
	a[3] = fields.E_oldest_max[i][j][0].e1();

	a[4] = fields.E_older_max[i][j][0].e1();
	a[5] = fields.E_older_max[i][j][1].e1();

	a[6] = fields.E_older_max[i+1][j][0].e1();
	a[7] = fields.E_older_max[i-1][j][0].e1();
	a[8] = fields.E_older_max[i][j+1][0].e1();
	a[9] = fields.E_older_max[i][j-1][0].e1();

	a[10] = fields.E_older_max[i][j][0].e1();

	a[0] = -a[1]+(0.5-1.0)/(0.5+1.0)*(a[2]+a[3])
	  +2.0/(0.5+1.0)*(a[4]+a[5])
	  +0.25/(0.5+1.0)*(a[6]+a[7]+a[8]+a[9]-4.0*a[10]);

	fields.E[i][j][ZMAX].set_e1(a[0]);

	
        /* set Ey value */
	a[1] = fields.E_oldest_max[i][j][1].e2();

	a[2] = fields.E[i][j][ZMAX-1].e2();
	a[3] = fields.E_oldest_max[i][j][0].e2();

	a[4] = fields.E_older_max[i][j][0].e2();
	a[5] = fields.E_older_max[i][j][1].e2();

	a[6] = fields.E_older_max[i+1][j][0].e2();
	a[7] = fields.E_older_max[i-1][j][0].e2();
	a[8] = fields.E_older_max[i][j+1][0].e2();
	a[9] = fields.E_older_max[i][j-1][0].e2();

	a[10] = fields.E_older_max[i][j][0].e2();

	a[0] = -a[1]+(0.5-1.0)/(0.5+1.0)*(a[2]+a[3])
	  +2.0/(0.5+1.0)*(a[4]+a[5])
	  +0.25/(0.5+1.0)*(a[6]+a[7]+a[8]+a[9]-4.0*a[10]);


	fields.E[i][j][ZMAX].set_e2(a[0]);

      }
    /* move to fieldinfo::average
    for(i=IMIN-2;i<=IMAX+2;i++)
      for(j=JMIN-2;j<=JMAX+2;j++){
	  fields.E_oldest_max[i][j][0]=fields.E_older_max[i][j][0];
	  fields.E_older_max[i][j][0] =fields.E[i][j][ZMAX-1];
	  fields.E_oldest_max[i][j][1]=fields.E_older_max[i][j][1];
	  fields.E_older_max[i][j][1] =fields.E[i][j][ZMAX];
      }
    */
    break;
  }
}



void LaserParam::boundarycorrect(int TYPE,FieldInfo& fields,ParParam *p,Scalar t)
{
    int i,j,k;
    Scalar dtdx = fields.dtdx;
    int IMIN = p->IMIN;
    int IMAX = p->IMAX;
    int JMIN = p->JMIN;
    int JMAX = p->JMAX;
    int KMIN = p->KMIN;
    int KMAX = p->KMAX;
    
    int XMIN = p->parameter->XMIN;    
    int XMAX = p->parameter->XMAX;    
    int YMIN = p->parameter->YMIN;    
    int YMAX = p->parameter->YMAX; 
    int ZMIN = p->parameter->ZMIN;
    int ZMAX = p->parameter->ZMAX;   
    
    switch(TYPE)
    {
    case 0:  // IMAX plane E correct
        i=IMAX;
	for(j=JMIN;j<JMAX;j++)
	  for(k=KMIN;k<KMAX;k++)
	    fields.E[i][j][k] = fields.E[i][j][k]+\
	      dtdx*ROATE_n(fields.B,i,j,k)-\
	      dtdx*TWOPI*fields.Ic[i][j][k];
	break;

    case 1:  // JMAX plane E correct
      j=JMAX;
      for(i=IMIN;i<IMAX;i++)
	for(k=KMIN;k<KMAX;k++)
	  fields.E[i][j][k] = fields.E[i][j][k]+\
	    dtdx*ROATE_n(fields.B,i,j,k)-\
	    dtdx*TWOPI*fields.Ic[i][j][k];
      break;
    }    
}


 
