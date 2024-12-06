#include "laser.h"

LaserParam::LaserParam(char* input_file_name)
{
  rf.openinput(input_file_name);
  Q         = atoi( rf.setget("&pulse", "Q"          ) );
  spotwidth = atof( rf.setget("&pulse", "spotwidth"  ) );
  amplitude = atof( rf.setget("&pulse", "amplitude"  ) );
  focus_z   = atoi( rf.setget("&pulse", "focus_z" ) ); 
  polarization    = atoi( rf.setget("&pulse", "polarization"   ) );
  shape           = atoi( rf.setget("&pulse", "shape"          ) );
  raise           = atoi( rf.setget("&pulse", "raise"          ) );
  duration        = atoi( rf.setget("&pulse", "duration"       ) );
  rf.closeinput();    
}

Scalar LaserParam::envelop(Scalar t)
{
  return (
	  shape == 1 ?
	  (
	   t < raise ?            ( t / raise ) :
	   t < duration - raise ? 1.0 :
	   t < duration ?         ( duration - t ) / raise :
	   0
	   ) :
	   shape == 2 ?
	   (
	    t < raise ? 	   sin( 0.5 * PI * t / raise ) :
	    t < duration - raise ? 1.0 :
	    t < duration ?         sin( 0.5 * PI * (duration - t) / raise ) :
	    0
	    ) :
       shape == 3 ?
       (
	     t < raise ?         sqr( sin( 0.5 * PI * t / raise ) ) :
	     t < duration - raise ? 1.0 :
	     t < duration ?      sqr( sin( 0.5 * PI * (duration - t) / raise ) ) :
	     0
	    ):
	   // assume shape == 4
        exp(  -sqr(t-2.5*raise)/(2*sqr(raise))  ) 

     );
}

Scalar LaserParam::Ex(int i,int j,ParParam* p,Scalar t)
{
    int XMIN = p->parameter->XMIN;
    int XMAX = p->parameter->XMAX;
    int YMIN = p->parameter->YMIN;
    int YMAX = p->parameter->YMAX;
    
    Scalar XMIDDLE = XMIN+0.5*(XMAX-XMIN);
    Scalar YMIDDLE = YMIN+0.5*(YMAX-YMIN);

    Scalar dx = 1.0/p->parameter->cells_per_wl;
    Scalar PHI0 = 0.0;
    Scalar w0 = spotwidth;
    Scalar x0 = XMIDDLE*dx;
    Scalar y0 = YMIDDLE*dx;
    Scalar z0 = focus_z;
    Scalar ZR = PI*w0*w0*0.5;
    Scalar zz = 0.0 -z0;
    Scalar xx = i*dx-x0;
    Scalar yy = j*dx-y0;
    Scalar RZ = zz*(1+sqr(ZR/zz));
    Scalar ww = w0*sqrt(1+zz*zz/(ZR*ZR));
    Scalar phase = TWOPI*(t-zz);
    Scalar phase_G = phase-zz*(xx*xx+yy*yy)/(ww*ww*ZR)+atan(zz/ZR)-PHI0;
    Scalar phase_G1 = phase_G+atan(zz/ZR);
#ifdef DIM_TWO    
    yy = 0.0; 
#endif   
    if(Q == 1){
	    if(polarization == 3)
	      return amplitude*envelop(t)*(w0/ww)*
                cos(phase_G)*
                exp( -(xx*xx)/(ww*ww))*
                exp( -(yy*yy)/(ww*ww));
	    else 
	      return sqrt(2)*amplitude*envelop(t)*(w0/ww)*
                cos(phase_G)*
                exp( -(xx*xx)/(ww*ww))*
                exp( -(yy*yy)/(ww*ww));
    }    
    else
      return 0  ;
}

Scalar LaserParam::Ey(int i,int j,ParParam* p,Scalar t)
{
    int XMIN = p->parameter->XMIN;
    int XMAX = p->parameter->XMAX;
    int YMIN = p->parameter->YMIN;
    int YMAX = p->parameter->YMAX;
    
    Scalar XMIDDLE = XMIN+0.5*(XMAX-XMIN);
    Scalar YMIDDLE = YMIN+0.5*(YMAX-YMIN);

    Scalar dx = 1.0/p->parameter->cells_per_wl;
    Scalar PHI0 = 0.0;
    Scalar w0 = spotwidth;
    Scalar x0 = XMIDDLE*dx;
    Scalar y0 = YMIDDLE*dx;
    Scalar z0 = focus_z;
    Scalar ZR = PI*w0*w0*0.5;
    Scalar zz = 0.0 -z0;
    Scalar xx = i*dx-x0;
    Scalar yy = j*dx-y0;
    Scalar RZ = zz*(1+sqr(ZR/zz));
    Scalar ww = w0*sqrt(1+zz*zz/(ZR*ZR));
    Scalar phase = TWOPI*(t-zz);
    Scalar phase_G = phase-zz*(xx*xx+yy*yy)/(ww*ww*ZR)+atan(zz/ZR)-PHI0;
    Scalar phase_G1 = phase_G+atan(zz/ZR);
#ifdef DIM_TWO    
    yy = 0.0; 
#endif   
    if(Q == 1){
	    if(polarization == 3)
	      return amplitude*envelop(t)*(w0/ww)*
                sin(phase_G)*
                exp( -(xx*xx)/(ww*ww))*
                exp( -(yy*yy)/(ww*ww));
	    else 
	      return sqrt(2)*amplitude*envelop(t)*(w0/ww)*
                sin(phase_G)*
                exp( -(xx*xx)/(ww*ww))*
                exp( -(yy*yy)/(ww*ww));
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
    for(i=IMIN;i<=IMAX;i++)
      for(j=JMIN;j<=JMAX;j++){
	//	fields.E[i][j][ZMIN] = fields.E[i][j][ZMIN]+\
	//	            dtdx*ROATE_n(fields.B,i,j,ZMIN)-\
	//                    dtdx*TWOPI*fields.Ic[i][j][ZMIN];
	Scalar g0      = 2*Ex(i,j,p,t);		
	Scalar Ex_old	= fields.E[i][j][ZMIN].e1();			
	Scalar By_old  = fields.E[i][j][ZMIN].e2();
	Scalar By_new  = (g0-2*Ex_old-(1-dtdx)*By_old)/(1+dtdx);
	fields.E[i][j][ZMIN].set_e1(0.5*g0);
	fields.E[i][j][ZMIN].set_e2(0.0);
	fields.E[i][j][ZMIN].set_e3(0.0);
      }
    break;
  case 1:
    for(i=IMIN;i<=IMAX;i++)
      for(j=JMIN;j<=JMAX;j++){
	//	fields.E[i][j][ZMAX] = fields.E[i][j][ZMAX]+\
	//	            dtdx*ROATE_n(fields.B,i,j,ZMAX)-\
	//                    dtdx*TWOPI*fields.Ic[i][j][ZMAX];

	Scalar g0      = 0.0;		
	Scalar Ex_old	= fields.E[i][j][ZMAX].e1();			
	Scalar By_old  = fields.B[i][j][ZMAX].e2();
	Scalar By_new  = (g0+2*Ex_old-(1-dtdx)*By_old)/(1+dtdx);
	fields.E[i][j][ZMAX].set_e1(0.0);
	fields.E[i][j][ZMAX].set_e2(0.0);
	fields.E[i][j][ZMAX].set_e3(0.0);
      }
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


 
