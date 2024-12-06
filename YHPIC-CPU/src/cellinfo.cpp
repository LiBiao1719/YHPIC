#include "cellinfo.h"
#include "fieldinfo.h"
#include <algorithm>
#ifdef PARALLEL
#include "mpi.h"
extern  MPI_Comm OOPIC_COMM;
extern  int MY_RANK;
#endif

using namespace std;

CellInfo::CellInfo(FieldInfo* _field)
{
  field = _field;
  dx = field->dx;
  dt = field->dt;
  stk = new Stack();
  para = field->par;
  IMIN = para->IMIN;
  IMAX = para->IMAX;
  JMIN = para->JMIN;
  JMAX = para->JMAX;
  KMIN = para->KMIN;
  KMAX = para->KMAX;
  
  size[0] = IMIN-2;
  size[1] = IMAX+2;
  size[2] = JMIN-2;
  size[3] = JMAX+2;
  size[4] = KMIN-2;
  size[5] = KMAX+2;

  CellsM.Init(size);

  InitCells();
  if(para->parameter->flag_restart == 0){
    ChainParticles();
    InitParticles();
  }
#ifdef HIGHORDERCLOUD
  int size_s0[4]={-2,2,1,3};
  int size_w1[6]={-2,2,-2,2,-2,2};
  int size_w2[6]={-3,2,-3,2,-3,2};
  S0.Init(size_s0);
  S1.Init(size_s0);
  WW.Init(size_w1);
  Ipar.Init(size_w2);
#endif
}

CellInfo::~CellInfo()
{
  int i,j,k,number;
  struct particle *pi;
  IJKLOOP{
    pi=CellsM[i][j][k].first;
    while(pi!=NULL){
      CellsM[i][j][k].first=pi->next;
      delete pi;
      pi=CellsM[i][j][k].first;
    }
  }
  delete stk;
}

void CellInfo::InitCells()
{
  int i,j,k;
  int cells_left = para->parameter->cells_left;
  int cells_ramp = para->parameter->cells_ramp;
  int cells_gap  = para->parameter->cells_gap;
  int cells_plasma = para->parameter->cells_plasma;
  int flag_restart = para->parameter->flag_restart;
  Scalar n_ion_ramp = para->parameter->n_ion_ramp/para->parameter->n_ion_over_nc;
  n_part = 0;


  Scalar XMIDDLE = para->parameter->XMIN + 0.5*(para->parameter->XMAX-para->parameter->XMIN);
  Scalar YMIDDLE = para->parameter->YMIN + 0.5*(para->parameter->YMAX-para->parameter->YMIN);
  Scalar ZMIDDLE = para->parameter->ZMIN + 0.5*(para->parameter->ZMAX-para->parameter->ZMIN);
/*
  int cells_r0 = 10;
  int cells_r1 = 5;
  int cells_tk = 2;
  IJKLOOP{
        
    if(k<cells_left){
      CellsM[i][j][k].dens[0]=0.0;
      CellsM[i][j][k].dens[1]=0.0;
    }
    else if(k<cells_left+cells_ramp ){
      CellsM[i][j][k].dens[0]=(k-cells_left)*n_ion_ramp/cells_ramp;
      CellsM[i][j][k].dens[1]=CellsM[i][j][k].dens[0];
    }
    else if(k<cells_left+cells_plasma ){
      Scalar rs = (cells_left+cells_plasma-k)*(cells_r0-cells_r1)/cells_plasma+cells_r1;
      Scalar ri = sqrt((i-XMIDDLE)*(i-XMIDDLE)+(j-YMIDDLE)*(j-YMIDDLE));
      if(ri>rs && ri<rs+5){
	CellsM[i][j][k].dens[0]=1.0;
	CellsM[i][j][k].dens[1]=CellsM[i][j][k].dens[0];
      }
    }
    else{
      CellsM[i][j][k].dens[0]=0.0;
      CellsM[i][j][k].dens[1]=CellsM[i][j][k].dens[0];
    }
  */


#ifdef PLASMABALL
/* plasma ball init
  Scalar XMIDDLE = para->parameter->XMIN + 0.5*(para->parameter->XMAX-para->parameter->XMIN);
  Scalar YMIDDLE = para->parameter->YMIN + 0.5*(para->parameter->YMAX-para->parameter->YMIN);
  Scalar ZMIDDLE = para->parameter->ZMIN + 0.5*(para->parameter->ZMAX-para->parameter->ZMIN);
  Scalar R0 = 0.5*cells_plasma;
  IJKLOOP{
    if((i-XMIDDLE)*(i-XMIDDLE)+(j-YMIDDLE)*(j-YMIDDLE)+(k-ZMIDDLE)*(k-ZMIDDLE)<R0*R0){
      CellsM[i][j][k].dens[0]=1.0;
      CellsM[i][j][k].dens[1]=1.0;
    }
    else{
      CellsM[i][j][k].dens[0]=0.0;
      CellsM[i][j][k].dens[1]=0.0;
    }
  }
*/ 
/*  IJKLOOP{
    if(k<(para->parameter->ZMAX-cells_plasma-0.02*(i-0.5*para->parameter->XMAX)*(i-0.5*para->parameter->XMAX))){
      CellsM[i][j][k].dens[0]=0.0;
      CellsM[i][j][k].dens[1]=0.0;
    }
    else{
      CellsM[i][j][k].dens[0]=1.0;
      CellsM[i][j][k].dens[1]=1.0;
    }
  }
*/
#else
    IJKLOOP{
        
    if(k<cells_left){
      CellsM[i][j][k].dens[0]=0.0;
      CellsM[i][j][k].dens[1]=0.0;
    }
    else if(k<cells_left+cells_ramp ){
      CellsM[i][j][k].dens[0]=(k-cells_left)*n_ion_ramp/cells_ramp;
      CellsM[i][j][k].dens[1]=CellsM[i][j][k].dens[0];
    }
    else if(k<cells_left+cells_plasma ){
      CellsM[i][j][k].dens[0]=1.0;
      CellsM[i][j][k].dens[1]=CellsM[i][j][k].dens[0];
    }
    else{
      CellsM[i][j][k].dens[0]=0.0;
      CellsM[i][j][k].dens[1]=CellsM[i][j][k].dens[0];
    }
    
    /*
    Scalar z0 = cells_left+0.5*cells_plasma;
    Scalar x0 = 0.5*(para->parameter->XMAX - para->parameter->XMIN);
    Scalar y0 = 0.5*(para->parameter->YMAX - para->parameter->YMIN);

    Scalar z_left = z0 - 0.5*cells_plasma*cos(para->parameter->angle);
    Scalar x_left = x0 + 0.5*cells_plasma*sin(para->parameter->angle);
    Scalar z_right = z0 + 0.5*cells_plasma*cos(para->parameter->angle);
    Scalar x_right = x0 - 0.5*cells_plasma*sin(para->parameter->angle);

    Scalar z_down = z0 - (x0-cells_gap)*sin(para->parameter->angle);
    Scalar x_down = x0 - (x0-cells_gap)*cos(para->parameter->angle);
    Scalar z_top  = z0 + (x0-cells_gap)*sin(para->parameter->angle);
    Scalar x_top  = x0 + (x0-cells_gap)*cos(para->parameter->angle);
    Scalar y_down = y0 - (y0-cells_gap);
    Scalar y_top  = y0 + (y0-cells_gap);

    Scalar z_plasma_left = z_left+(i-x_left)*tan(para->parameter->angle);
    Scalar z_plasma_right= z_right +(i-x_right)*tan(para->parameter->angle);

    Scalar x_plasma_down = x_down+(z_down-k)*(tan(para->parameter->angle));
    Scalar x_plasma_top  = x_top +(z_top -k)*(tan(para->parameter->angle));

    if(k>z_plasma_left && k<z_plasma_right && i>x_plasma_down && i<x_plasma_top  && flag_restart==0){
#ifndef DIM_TWO
if(j>y_down&&j<y_top)
{
#endif
      CellsM[i][j][k].dens[0] = 1.0;
      CellsM[i][j][k].dens[1] = CellsM[i][j][k].dens[0];
#ifndef DIM_TWO
}
#endif
    }	
    else{
      CellsM[i][j][k].dens[0]=0.0;
      CellsM[i][j][k].dens[1]=0.0;
    }
    */ 
   }
#endif
#ifdef DEBUGMOVE
  CellsM[IMAX/2][JMAX/2][10].dens[0]=1.0;
#endif                                                                             
}

void CellInfo::ChainParticles()
{
  int i,j,k,ppc[2];
  int ITemp1,ITemp2,ITemp3;
  Scalar delta,delta2,sqrt3,sqrt32;
  struct particle *pn,*po;

  ppc[0] = para->parameter->ppc[0];
  ppc[1] = para->parameter->ppc[1];

  IJKLOOP{
    if(i==4)
    {
	    int xxx=0;
    }
    CellsM[i][j][k].first = NULL;
    CellsM[i][j][k].last  = NULL;

    po = CellsM[i][j][k].first;

    for(int ispe=0;ispe<2;ispe++){
      int number = (int) FLOOR(CellsM[i][j][k].dens[ispe]*ppc[ispe]+0.5);
      if(para->parameter->specy[ispe].fix==1) number = 0;
      CellsM[i][j][k].np[ispe] = number;
      if(ispe==0) n_el      += number;
      else        n_ion     += number;
      CellsM[i][j][k].npart += number;
      n_part                += number;
      if(number!=0){
	sqrt3 = int(CEIL(POW(number,1.0/3.0)));
	sqrt32= sqrt3*sqrt3;
	delta = 1.0 / sqrt3;
	
	for(int index=1;index<=number;index++){
	  pn          = new ( struct particle );
	  if (!pn) {printf("allocation error\n"); exit(1);}
	  
	  pn->prev    = po;
	  if (po==NULL) CellsM[i][j][k].first = pn;
	  else          po->next    = pn;
	  pn->species = ispe;
	  ITemp1      = (int)FLOOR((index-1)/sqrt32);
	  ITemp2      = (int)FLOOR((index-1)/sqrt3)-ITemp1*sqrt3;
	  ITemp3      = MAX(0,int(index -1 -ITemp1*sqrt32 - ITemp2*sqrt3));
	  pn->x       = Scalar(i) + ((Scalar)ITemp1+0.49999999) * delta;
	  pn->y       = Scalar(j) + ((Scalar)ITemp2+0.49999999) * delta;
	  pn->z       = Scalar(k) + ((Scalar)ITemp3+0.49999999) * delta;
#ifdef DEBUG
	  if((int)FLOOR(pn->x)!=i || (int)FLOOR(pn->y)!=j || (int)FLOOR(pn->z)!=k){
	    printf("initial particle out cell\n");
	    exit(1);
	  }
#endif 
	  po          = pn;
	}
      }
    }
    if (po!=NULL) {
      po->next = NULL;
      CellsM[i][j][k].last = po;
    }
  }
#ifdef PARALLEL
      int temp;
      MPI_Allreduce(&n_part,&temp,1,MPI_INT,MPI_SUM,OOPIC_COMM);
      n_part = temp;
#endif

}

void CellInfo::InitParticles()
{
  int i,j,k,index;
  int number[2];
  Scalar ux,uy,uz;
  struct particle *part;
  SerParam *parameter = para->parameter;

  IJKLOOP{
    if(CellsM[i][j][k].first!=NULL){
      for(part=CellsM[i][j][k].first; part!=NULL;part=part->next){
#ifdef DEBUG
	if(!part) { printf("segmentation\n");exit(1); }
	if((int)FLOOR(part->x)!=i || (int)FLOOR(part->y)!=j || (int)FLOOR(part->z)!=k ){
 	  printf("InitParticle:particle not int cell \n");
	  cout<<i<<" "<<j<<" "<<k<<" "<<part->x<<" "<<part->y<<" "<<part->z<<"\n";
	  exit(1);
	}
#endif
	int species= part->species;
	int    fix = parameter->specy[species].fix;
	Scalar ze  = parameter->specy[species].ze;
	Scalar m   = parameter->specy[species].m;
	Scalar zm  = parameter->specy[species].zm;
	// part->zn = charge state * density / critical density
	// ---------------------------------------------------------
	// Lorentz-Transformation: part->zn is scaled up with Gamma^3!
	// L-Contraction in y-direction leads to n_M = Gamma n_L
	// and Doppler shift leads to n_c_M = 1/Gamma^2 * n_c_L
	// ---------------------------------------------------------
	// part->zn is designed such that the sum of MacroParticle charges
	// (electrons and ions) is zero in each cell initially
	
	
	//	if (ze == 0) {     // neutral atoms
	//	  parameter->specy[part->species].n  = parameter->n_ion_over_nc / parameter->ppc[part->species];
	//	  parameter->specy[part->species].zn = 0;
	//	}
	//	else {                  // electrons or ions
	//	  parameter->specy[part->species].n  = fabs( 1.0 * parameter->specy[1].ze / ze ) 
	//	    * parameter->n_ion_over_nc / parameter->ppc[species];
	//	  parameter->specy[part->species].zn = ze * parameter->specy[part->species].n;
	//	}
		
	ux = parameter->vtherm[part->species]*gauss_rand48();
	uy = parameter->vtherm[part->species]*gauss_rand48();
	uz = parameter->vtherm[part->species]*gauss_rand48();
	/*	
	if(part->species == 0){
	  vx = 0.1;
	  vy = 0.0;
	  vz = 0.1;
	}else{
	  vx = 0.0;
	  vy = 0.0;
	  vz = 0.0;
	}
	*/
	// determine gamma*v
//	part->igamma  = sqrt( 1.0 - vx*vx - vy*vy - vz*vz );
	part->igamma = (1.0 / sqrt(1.0 + (ux*ux + uy*uy + uz*uz)));
	part->ux      = ux ;
	part->uy      = uy ;
	part->uz      = uz ;
	
      }
    }
  }
}
//////////////////////////////////////////////////////////////////////////////////////////


Scalar CellInfo::gauss_rand48( void )
{
/*
  Scalar r1, r2;

  r1 = drand48();
  r2 = drand48();

  return sqrt( -2.0 * log( 1.0 - r1 ) ) * sin( 2*PI*r2 );
*/

    static double V1, V2, S;
    static int phase = 0;
    double X;
    
    if ( phase == 0 ) {
        do {
            double U1 = (double)rand() / RAND_MAX;
            double U2 = (double)rand() / RAND_MAX;
            
            V1 = 2 * U1 - 1;
            V2 = 2 * U2 - 1;
            S = V1 * V1 + V2 * V2;
        } while(S >= 1 || S == 0);
        
        X = V1 * sqrt(-2 * log(S) / S);
    } else
        X = V2 * sqrt(-2 * log(S) / S);
    
    phase = 1 - phase;
    
    return X;
}

void CellInfo::acceleratefield(int i,int j,int k)
{
  enf[0] = field->ENode[i][j][k];
  enf[1] = field->ENode[i][j][k+1];
  enf[2] = field->ENode[i][j+1][k];
  enf[3] = field->ENode[i][j+1][k+1];;
  enf[4] = field->ENode[i+1][j][k];
  enf[5] = field->ENode[i+1][j][k+1];
  enf[6] = field->ENode[i+1][j+1][k];
  enf[7] = field->ENode[i+1][j+1][k+1];

  bnf[0] = field->BNode[i][j][k];
  bnf[1] = field->BNode[i][j][k+1];
  bnf[2] = field->BNode[i][j+1][k];
  bnf[3] = field->BNode[i][j+1][k+1];;
  bnf[4] = field->BNode[i+1][j][k];
  bnf[5] = field->BNode[i+1][j][k+1];
  bnf[6] = field->BNode[i+1][j+1][k];
  bnf[7] = field->BNode[i+1][j+1][k+1];
}

void CellInfo::accelerate(particle* part,int i,int j,int k)
{
#ifdef DEBUG
  if((int)FLOOR(part->x)!=i||(int)FLOOR(part->y)!=j||(int)FLOOR(part->z)!=k){
    printf("accelerate particle is out cell\n");
    exit(1);
  }
#endif
  //
  // this function is currently not used
  //
  // full acceleration according Boris (in Birdsall, Langdon) 

  register Scalar zmpidt = para->parameter->specy[part->species].zm * PI * dt;
#ifdef HIGHORDERCLOUD

  //clear S0 S1
  for(int mm=-2;mm<=2;mm++)
    for(int nn=1;nn<=3;nn++){
      S0[mm][nn]=0.0;
      S1[mm][nn]=0.0;
    }

  int ii = RINT(part->x);
  int jj = RINT(part->y);
  int kk = RINT(part->z);

  //x direction
  Scalar delta_x = ii - part->x;
  S0[0][1]  = 0.75-delta_x*delta_x;
  S0[1][1]  = 0.5*(0.5-delta_x)*(0.5-delta_x);
  S0[-1][1] = 0.5*(0.5+delta_x)*(0.5+delta_x);

  //y direction
  Scalar delta_y = jj - part->y;
  S0[0][2]  = 0.75-delta_y*delta_y;
  S0[1][2]  = 0.5*(0.5-delta_y)*(0.5-delta_y);
  S0[-1][2] = 0.5*(0.5+delta_y)*(0.5+delta_y);
  
  //z direction
  Scalar delta_z = kk - part->z;
  S0[0][3]  = 0.75-delta_z*delta_z;
  S0[1][3]  = 0.5*(0.5-delta_z)*(0.5-delta_z);
  S0[-1][3] = 0.5*(0.5+delta_z)*(0.5+delta_z);

  Scalar ex,ey,ez;
  Scalar bx,by,bz;
  ex=ey=ez=bx=by=bz=0.0;

  for(int im=-1;im<=1;im++)
    for(int jm=-1;jm<=1;jm++)
      for(int km=-1;km<=1;km++)
	{
	  ex += S0[im][1]*S0[jm][2]*S0[km][3]*field->ENode[ii+im][jj+jm][kk+km].e1();
	  ey += S0[im][1]*S0[jm][2]*S0[km][3]*field->ENode[ii+im][jj+jm][kk+km].e2();
	  ez += S0[im][1]*S0[jm][2]*S0[km][3]*field->ENode[ii+im][jj+jm][kk+km].e3();
	  bx += S0[im][1]*S0[jm][2]*S0[km][3]*field->BNode[ii+im][jj+jm][kk+km].e1();
	  by += S0[im][1]*S0[jm][2]*S0[km][3]*field->BNode[ii+im][jj+jm][kk+km].e2();
	  bz += S0[im][1]*S0[jm][2]*S0[km][3]*field->BNode[ii+im][jj+jm][kk+km].e3();
	} 
#else
  register Scalar w1    = part->x-i;
  register Scalar w2    = part->y-j;
  register Scalar w3    = part->z-k;
  
  register Scalar ex = (1-w1)*(1-w2)*(1-w3) * enf[0].e1() 
                     + (1-w1)*(1-w2)*w3     * enf[1].e1()
                     + (1-w1)*w2*(1-w3)     * enf[2].e1()
                     + (1-w1)*w2*w3         * enf[3].e1()
                     + w1*(1-w2)*(1-w3)     * enf[4].e1() 
                     + w1*(1-w2)*w3         * enf[5].e1()
                     + w1*w2*(1-w3)         * enf[6].e1()
                     + w1*w2*w3             * enf[7].e1();
       // interpolate fields
  register Scalar ey = (1-w1)*(1-w2)*(1-w3) * enf[0].e2() 
                     + (1-w1)*(1-w2)*w3     * enf[1].e2()
                     + (1-w1)*w2*(1-w3)     * enf[2].e2()
                     + (1-w1)*w2*w3         * enf[3].e2()
                     + w1*(1-w2)*(1-w3)     * enf[4].e2() 
                     + w1*(1-w2)*w3         * enf[5].e2()
                     + w1*w2*(1-w3)         * enf[6].e2()
                     + w1*w2*w3             * enf[7].e2();

  register Scalar ez = (1-w1)*(1-w2)*(1-w3) * enf[0].e3() 
                     + (1-w1)*(1-w2)*w3     * enf[1].e3()
                     + (1-w1)*w2*(1-w3)     * enf[2].e3()
                     + (1-w1)*w2*w3         * enf[3].e3()
                     + w1*(1-w2)*(1-w3)     * enf[4].e3() 
                     + w1*(1-w2)*w3         * enf[5].e3()
                     + w1*w2*(1-w3)         * enf[6].e3()
                     + w1*w2*w3             * enf[7].e3();

  register Scalar bx = (1-w1)*(1-w2)*(1-w3) * bnf[0].e1() 
                     + (1-w1)*(1-w2)*w3     * bnf[1].e1()
                     + (1-w1)*w2*(1-w3)     * bnf[2].e1()
                     + (1-w1)*w2*w3         * bnf[3].e1()
                     + w1*(1-w2)*(1-w3)     * bnf[4].e1() 
                     + w1*(1-w2)*w3         * bnf[5].e1()
                     + w1*w2*(1-w3)         * bnf[6].e1()
                     + w1*w2*w3             * bnf[7].e1();

  register Scalar by = (1-w1)*(1-w2)*(1-w3) * bnf[0].e2() 
                     + (1-w1)*(1-w2)*w3     * bnf[1].e2()
                     + (1-w1)*w2*(1-w3)     * bnf[2].e2()
                     + (1-w1)*w2*w3         * bnf[3].e2()
                     + w1*(1-w2)*(1-w3)     * bnf[4].e2() 
                     + w1*(1-w2)*w3         * bnf[5].e2()
                     + w1*w2*(1-w3)         * bnf[6].e2()
                     + w1*w2*w3             * bnf[7].e2();

  register Scalar bz = (1-w1)*(1-w2)*(1-w3) * bnf[0].e3() 
                     + (1-w1)*(1-w2)*w3     * bnf[1].e3()
                     + (1-w1)*w2*(1-w3)     * bnf[2].e3()
                     + (1-w1)*w2*w3         * bnf[3].e3()
                     + w1*(1-w2)*(1-w3)     * bnf[4].e3() 
                     + w1*(1-w2)*w3         * bnf[5].e3()
                     + w1*w2*(1-w3)         * bnf[6].e3()
                     + w1*w2*w3             * bnf[7].e3();

#endif  

  register Scalar ux = part->ux + ex * zmpidt;                     // half acceleration
  register Scalar uy = part->uy + ey * zmpidt;
  register Scalar uz = part->uz + ez * zmpidt;
  
  register Scalar igamma = (1.0 / sqrt(1.0 + (ux*ux + uy*uy + uz*uz)));

  register Scalar tx = bx * zmpidt * igamma;
  register Scalar ty = by * zmpidt * igamma;                       // rotation 
  register Scalar tz = bz * zmpidt * igamma;
  register Scalar t_1 = 1.0/(1.0+tx*tx + ty*ty + tz*tz);
  register Scalar tx2 = tx*tx;
  register Scalar txy = tx*ty;
  register Scalar txz = tx*tz;
  register Scalar ty2 = ty*ty;
  register Scalar tyz = ty*tz;
  register Scalar tz2 = tz*tz;
  //  register Scalar sx = 2.0 * tz / ( 1.0 + t2 );
  //  register Scalar sy = 2.0 * ty / ( 1.0 + t2 );
  //  register Scalar sz = 2.0 * tz / ( 1.0 + t2 );
  
  register Scalar ux2 = t_1*( ux*(1.0+tx2-ty2-tz2) + 2.0*uy*(txy+tz) + 2.0*uz*(txz-ty) );
  register Scalar uy2 = t_1*( 2.0*ux*(txy-tz) + uy*(1.0-tx2+ty2-tz2) + 2.0*uz*(tyz+tx) );
  register Scalar uz2 = t_1*( 2.0*ux*(txz+ty) + 2.0*uy*(tyz-tx) + uz*(1.0-tx2-ty2+tz2) );

  part->ux = ux = ux2 + ex * zmpidt;           // second half acceleration
  part->uy = uy = uy2 + ey * zmpidt;
  part->uz = uz = uz2 + ez * zmpidt;
  
  part->igamma = (1.0 / sqrt(1.0 + (ux*ux + uy*uy + uz*uz)));

  // we take the square roots for all particles per cell in a separate loop !
  
}

void CellInfo::move(struct particle* part,int i,int j,int k,volatile Scalar &pdx,volatile Scalar &pdy,volatile Scalar &pdz)
{
  if ( para->parameter->specy[part->species].fix==1 ) {  }
  else {

    pdx = part->ux * part->igamma * 0.5; //change dx to normalize( * dt/dx=0.5) 
    pdy = 1.0/SCALE*part->uy * part->igamma * 0.5; // 
    pdz = 1.0/SCALE*part->uz * part->igamma * 0.5; // 

    part->x += pdx;
    part->y += pdy;
    part->z += pdz;

    if((int)FLOOR(part->x!=i) ||(int)FLOOR(part->y)!=j ||(int)FLOOR(part->z)!=k)
	    stk->put_on_stack(part,i,j,k);

#ifdef DEBUG
    if ( fabs(pdx) > 1.0 || fabs(pdy) > 1.0 || fabs(pdz) > 1.0){
      printf( "particle displacement larger than grid spacing!\n" );
      exit(1);
    }
#endif

  }
}

void CellInfo::has_to_change_cell(particle* part,int i,int j,int k)
{
}

void CellInfo::deposit_charge(particle* part,int i,int j,int k)
{
#ifdef DEBUG
  if( (int)floor(part->x)!=i || (int)floor(part->y)!=j || (int)floor(part->z)!=k ){
    printf("deposit_charge:particle out cell\n");
    exit(1);
  }
#endif
  register double w1    = part->x-i;
  register double w2    = part->y-j;
  register double w3    = part->z-k;
  register double dn     = para->parameter->specy[part->species].n;

  field->Dens[part->species][i][j][k]      += (1-w1)*(1-w2)*(1-w3) * dn; 
  field->Dens[part->species][i][j][k+1]    += (1-w1)*(1-w2)*w3     * dn;
  field->Dens[part->species][i][j+1][k]    += (1-w1)*w2*(1-w3)     * dn;
  field->Dens[part->species][i][j+1][k+1]  += (1-w1)*w2*w3         * dn;
  field->Dens[part->species][i+1][j][k]    += w1*(1-w2)*(1-w3)     * dn;
  field->Dens[part->species][i+1][j][k+1]  += w1*(1-w2)*w3         * dn;
  field->Dens[part->species][i+1][j+1][k]  += w1*w2*(1-w3)         * dn;
  field->Dens[part->species][i+1][j+1][k+1]+= w1*w2*w3             * dn;
}

void CellInfo::deposit_current(particle* part,int i1,int j1,int k1,volatile Scalar pdx
,volatile Scalar pdy,volatile Scalar pdz)
{
  Scalar iscale = 1.0/SCALE;
#ifdef HIGHORDERCLOUD
  Scalar xr,yr,zr;
  Scalar x1=part->x-pdx;
  Scalar y1=part->y-pdy;
  Scalar z1=part->z-pdz;

  Scalar x2=part->x;
  Scalar y2=part->y;
  Scalar z2=part->z;
  Scalar zn=para->parameter->specy[part->species].zn;
  Scalar zn2=2.0*para->parameter->specy[part->species].zn;
  Scalar qvx = zn*part->ux*part->igamma ; 

// dx/dt=2.0
  int ii0 = RINT(x1);
  int jj0 = RINT(y1);
  int kk0 = RINT(z1);

  int ii = RINT(x2);
  int jj = RINT(y2);
  int kk = RINT(z2);

  if(ii0==ii)
    xr = 0.5*(x1+x2);
  else
    xr = 0.5*(ii0+ii);
  if(jj0==jj)
    yr = 0.5*(y1+y2);
  else
    yr = 0.5*(jj0+jj);
  if(kk0==kk)
    zr = 0.5*(z1+z2);
  else
    zr = 0.5*(kk0+kk);

  Scalar Fx1 = iscale*iscale*zn2*(xr-x1);
  Scalar Fy1 = iscale*zn2*(yr-y1);
  Scalar Fz1 = iscale*zn2*(zr-z1);
  Scalar Fx2 = iscale*iscale*zn2*(x2-xr);
  Scalar Fy2 = iscale*zn2*(y2-yr);
  Scalar Fz2 = iscale*zn2*(z2-zr);

  //x direction for 1st
  Scalar delta_x1   = ii0-0.5*(xr+x1);
  Scalar wi1_       = 0.5*(0.5+delta_x1)*(0.5+delta_x1);
  Scalar wi1        = 0.75-delta_x1*delta_x1;
  Scalar wi_1       = 0.5*(0.5-delta_x1)*(0.5-delta_x1);

  //y direction for 1st
  Scalar delta_y1   = jj0-0.5*(yr+y1);
  Scalar wj1_       = 0.5*(0.5+delta_y1)*(0.5+delta_y1);
  Scalar wj1        = 0.75-delta_y1*delta_y1;
  Scalar wj_1       = 0.5*(0.5-delta_y1)*(0.5-delta_y1);

  //z direction for 1st
  Scalar delta_z1   = kk0-0.5*(zr+z1);
  Scalar wk1_       = 0.5*(0.5+delta_z1)*(0.5+delta_z1);
  Scalar wk1        = 0.75-delta_z1*delta_z1;
  Scalar wk_1       = 0.5*(0.5-delta_z1)*(0.5-delta_z1);

  //x direction for 2nd
  Scalar delta_x2   = ii-0.5*(xr+x2);
  Scalar wi2_       = 0.5*(0.5+delta_x2)*(0.5+delta_x2);
  Scalar wi2        = 0.75-delta_x2*delta_x2;
  Scalar wi_2       = 0.5*(0.5-delta_x2)*(0.5-delta_x2);

  //y direction for 2nd
  Scalar delta_y2   = jj-0.5*(yr+y2);
  Scalar wj2_       = 0.5*(0.5+delta_y2)*(0.5+delta_y2);
  Scalar wj2        = 0.75-delta_y2*delta_y2;
  Scalar wj_2       = 0.5*(0.5-delta_y2)*(0.5-delta_y2);

  //z direction for 2nd
  Scalar delta_z2   = kk-0.5*(zr+z2);
  Scalar wk2_       = 0.5*(0.5+delta_z2)*(0.5+delta_z2);
  Scalar wk2        = 0.75-delta_z2*delta_z2;
  Scalar wk_2       = 0.5*(0.5-delta_z2)*(0.5-delta_z2);

  // ist current *****************************************************************
  // along y
  field->Ic[ii0-1][jj0-1][kk0-1] += Vector3(0.0,Fy1*(0.5+delta_y1)*wi1_*wk1_,0.0);
  field->Ic[ii0-1][jj0-1][kk0  ] += Vector3(0.0,Fy1*(0.5+delta_y1)*wi1_*wk1,0.0);
  field->Ic[ii0-1][jj0-1][kk0+1] += Vector3(0.0,Fy1*(0.5+delta_y1)*wi1_*wk_1,0.0);
  field->Ic[ii0-1][jj0  ][kk0-1] += Vector3(0.0,Fy1*(0.5-delta_y1)*wi1_*wk1_,0.0);
  field->Ic[ii0-1][jj0  ][kk0  ] += Vector3(0.0,Fy1*(0.5-delta_y1)*wi1_*wk1,0.0);
  field->Ic[ii0-1][jj0  ][kk0+1] += Vector3(0.0,Fy1*(0.5-delta_y1)*wi1_*wk_1,0.0);

  field->Ic[ii0][jj0-1][kk0-1] += Vector3(0.0,Fy1*(0.5+delta_y1)*wi1*wk1_,0.0);
  field->Ic[ii0][jj0-1][kk0  ] += Vector3(0.0,Fy1*(0.5+delta_y1)*wi1*wk1,0.0);
  field->Ic[ii0][jj0-1][kk0+1] += Vector3(0.0,Fy1*(0.5+delta_y1)*wi1*wk_1,0.0);
  field->Ic[ii0][jj0  ][kk0-1] += Vector3(0.0,Fy1*(0.5-delta_y1)*wi1*wk1_,0.0);
  field->Ic[ii0][jj0  ][kk0  ] += Vector3(0.0,Fy1*(0.5-delta_y1)*wi1*wk1,0.0);
  field->Ic[ii0][jj0  ][kk0+1] += Vector3(0.0,Fy1*(0.5-delta_y1)*wi1*wk_1,0.0);

  field->Ic[ii0+1][jj0-1][kk0-1] += Vector3(0.0,Fy1*(0.5+delta_y1)*wi_1*wk1_,0.0);
  field->Ic[ii0+1][jj0-1][kk0  ] += Vector3(0.0,Fy1*(0.5+delta_y1)*wi_1*wk1,0.0);
  field->Ic[ii0+1][jj0-1][kk0+1] += Vector3(0.0,Fy1*(0.5+delta_y1)*wi_1*wk_1,0.0);
  field->Ic[ii0+1][jj0  ][kk0-1] += Vector3(0.0,Fy1*(0.5-delta_y1)*wi_1*wk1_,0.0);
  field->Ic[ii0+1][jj0  ][kk0  ] += Vector3(0.0,Fy1*(0.5-delta_y1)*wi_1*wk1,0.0);
  field->Ic[ii0+1][jj0  ][kk0+1] += Vector3(0.0,Fy1*(0.5-delta_y1)*wi_1*wk_1,0.0);
  
  // along z
  field->Ic[ii0-1][jj0-1][kk0-1] += Vector3(0.0,0.0,Fz1*(0.5+delta_z1)*wi1_*wj1_);
  field->Ic[ii0-1][jj0  ][kk0-1] += Vector3(0.0,0.0,Fz1*(0.5+delta_z1)*wi1_*wj1 );
  field->Ic[ii0-1][jj0+1][kk0-1] += Vector3(0.0,0.0,Fz1*(0.5+delta_z1)*wi1_*wj_1);
  field->Ic[ii0-1][jj0-1][kk0  ] += Vector3(0.0,0.0,Fz1*(0.5-delta_z1)*wi1_*wj1_);
  field->Ic[ii0-1][jj0  ][kk0  ] += Vector3(0.0,0.0,Fz1*(0.5-delta_z1)*wi1_*wj1 );
  field->Ic[ii0-1][jj0+1][kk0  ] += Vector3(0.0,0.0,Fz1*(0.5-delta_z1)*wi1_*wj_1);
  
  field->Ic[ii0][jj0-1][kk0-1] += Vector3(0.0,0.0,Fz1*(0.5+delta_z1)*wi1*wj1_);
  field->Ic[ii0][jj0  ][kk0-1] += Vector3(0.0,0.0,Fz1*(0.5+delta_z1)*wi1*wj1 );
  field->Ic[ii0][jj0+1][kk0-1] += Vector3(0.0,0.0,Fz1*(0.5+delta_z1)*wi1*wj_1);
  field->Ic[ii0][jj0-1][kk0  ] += Vector3(0.0,0.0,Fz1*(0.5-delta_z1)*wi1*wj1_);
  field->Ic[ii0][jj0  ][kk0  ] += Vector3(0.0,0.0,Fz1*(0.5-delta_z1)*wi1*wj1 );
  field->Ic[ii0][jj0+1][kk0  ] += Vector3(0.0,0.0,Fz1*(0.5-delta_z1)*wi1*wj_1);
  
  field->Ic[ii0+1][jj0-1][kk0-1] += Vector3(0.0,0.0,Fz1*(0.5+delta_z1)*wi_1*wj1_);
  field->Ic[ii0+1][jj0  ][kk0-1] += Vector3(0.0,0.0,Fz1*(0.5+delta_z1)*wi_1*wj1 );
  field->Ic[ii0+1][jj0+1][kk0-1] += Vector3(0.0,0.0,Fz1*(0.5+delta_z1)*wi_1*wj_1);
  field->Ic[ii0+1][jj0-1][kk0  ] += Vector3(0.0,0.0,Fz1*(0.5-delta_z1)*wi_1*wj1_);
  field->Ic[ii0+1][jj0  ][kk0  ] += Vector3(0.0,0.0,Fz1*(0.5-delta_z1)*wi_1*wj1 );
  field->Ic[ii0+1][jj0+1][kk0  ] += Vector3(0.0,0.0,Fz1*(0.5-delta_z1)*wi_1*wj_1);
  
  //2nd current ******************************************************************
  // along y
  field->Ic[ii-1][jj-1][kk-1] += Vector3(0.0,Fy2*(0.5+delta_y2)*wi2_*wk2_,0.0);
  field->Ic[ii-1][jj-1][kk  ] += Vector3(0.0,Fy2*(0.5+delta_y2)*wi2_*wk2,0.0);
  field->Ic[ii-1][jj-1][kk+1] += Vector3(0.0,Fy2*(0.5+delta_y2)*wi2_*wk_2,0.0);
  field->Ic[ii-1][jj  ][kk-1] += Vector3(0.0,Fy2*(0.5-delta_y2)*wi2_*wk2_,0.0);
  field->Ic[ii-1][jj  ][kk  ] += Vector3(0.0,Fy2*(0.5-delta_y2)*wi2_*wk2,0.0);
  field->Ic[ii-1][jj  ][kk+1] += Vector3(0.0,Fy2*(0.5-delta_y2)*wi2_*wk_2,0.0);

  field->Ic[ii][jj-1][kk-1] += Vector3(0.0,Fy2*(0.5+delta_y2)*wi2*wk2_,0.0);
  field->Ic[ii][jj-1][kk  ] += Vector3(0.0,Fy2*(0.5+delta_y2)*wi2*wk2,0.0);
  field->Ic[ii][jj-1][kk+1] += Vector3(0.0,Fy2*(0.5+delta_y2)*wi2*wk_2,0.0);
  field->Ic[ii][jj  ][kk-1] += Vector3(0.0,Fy2*(0.5-delta_y2)*wi2*wk2_,0.0);
  field->Ic[ii][jj  ][kk  ] += Vector3(0.0,Fy2*(0.5-delta_y2)*wi2*wk2,0.0);
  field->Ic[ii][jj  ][kk+1] += Vector3(0.0,Fy2*(0.5-delta_y2)*wi2*wk_2,0.0);

  field->Ic[ii+1][jj-1][kk-1] += Vector3(0.0,Fy2*(0.5+delta_y2)*wi_2*wk2_,0.0);
  field->Ic[ii+1][jj-1][kk  ] += Vector3(0.0,Fy2*(0.5+delta_y2)*wi_2*wk2,0.0);
  field->Ic[ii+1][jj-1][kk+1] += Vector3(0.0,Fy2*(0.5+delta_y2)*wi_2*wk_2,0.0);
  field->Ic[ii+1][jj  ][kk-1] += Vector3(0.0,Fy2*(0.5-delta_y2)*wi_2*wk2_,0.0);
  field->Ic[ii+1][jj  ][kk  ] += Vector3(0.0,Fy2*(0.5-delta_y2)*wi_2*wk2,0.0);
  field->Ic[ii+1][jj  ][kk+1] += Vector3(0.0,Fy2*(0.5-delta_y2)*wi_2*wk_2,0.0);

  // along z
  field->Ic[ii-1][jj-1][kk-1] += Vector3(0.0,0.0,Fz2*(0.5+delta_z2)*wi2_*wj2_);
  field->Ic[ii-1][jj  ][kk-1] += Vector3(0.0,0.0,Fz2*(0.5+delta_z2)*wi2_*wj2 );
  field->Ic[ii-1][jj+1][kk-1] += Vector3(0.0,0.0,Fz2*(0.5+delta_z2)*wi2_*wj_2);
  field->Ic[ii-1][jj-1][kk  ] += Vector3(0.0,0.0,Fz2*(0.5-delta_z2)*wi2_*wj2_);
  field->Ic[ii-1][jj  ][kk  ] += Vector3(0.0,0.0,Fz2*(0.5-delta_z2)*wi2_*wj2 );
  field->Ic[ii-1][jj+1][kk  ] += Vector3(0.0,0.0,Fz2*(0.5-delta_z2)*wi2_*wj_2);

  field->Ic[ii][jj-1][kk-1] += Vector3(0.0,0.0,Fz2*(0.5+delta_z2)*wi2*wj2_);
  field->Ic[ii][jj  ][kk-1] += Vector3(0.0,0.0,Fz2*(0.5+delta_z2)*wi2*wj2 );
  field->Ic[ii][jj+1][kk-1] += Vector3(0.0,0.0,Fz2*(0.5+delta_z2)*wi2*wj_2);
  field->Ic[ii][jj-1][kk  ] += Vector3(0.0,0.0,Fz2*(0.5-delta_z2)*wi2*wj2_);
  field->Ic[ii][jj  ][kk  ] += Vector3(0.0,0.0,Fz2*(0.5-delta_z2)*wi2*wj2 );
  field->Ic[ii][jj+1][kk  ] += Vector3(0.0,0.0,Fz2*(0.5-delta_z2)*wi2*wj_2);

  field->Ic[ii+1][jj-1][kk-1] += Vector3(0.0,0.0,Fz2*(0.5+delta_z2)*wi_2*wj2_);
  field->Ic[ii+1][jj  ][kk-1] += Vector3(0.0,0.0,Fz2*(0.5+delta_z2)*wi_2*wj2 );
  field->Ic[ii+1][jj+1][kk-1] += Vector3(0.0,0.0,Fz2*(0.5+delta_z2)*wi_2*wj_2);
  field->Ic[ii+1][jj-1][kk  ] += Vector3(0.0,0.0,Fz2*(0.5-delta_z2)*wi_2*wj2_);
  field->Ic[ii+1][jj  ][kk  ] += Vector3(0.0,0.0,Fz2*(0.5-delta_z2)*wi_2*wj2 );
  field->Ic[ii+1][jj+1][kk  ] += Vector3(0.0,0.0,Fz2*(0.5-delta_z2)*wi_2*wj_2);

  // along x ******************************************************************
  // 1st 
  field->Ic[ii0-1][jj0-1][kk0-1]  += Vector3(Fx1*(0.5+delta_x1)*wj1_*wk1_,0.0,0.0);
  field->Ic[ii0-1][jj0-1][kk0]    += Vector3(Fx1*(0.5+delta_x1)*wj1_*wk1 ,0.0,0.0);
  field->Ic[ii0-1][jj0-1][kk0+1]  += Vector3(Fx1*(0.5+delta_x1)*wj1_*wk_1,0.0,0.0);
  field->Ic[ii0-1][jj0  ][kk0-1]  += Vector3(Fx1*(0.5+delta_x1)*wj1 *wk1_,0.0,0.0);
  field->Ic[ii0-1][jj0  ][kk0]    += Vector3(Fx1*(0.5+delta_x1)*wj1 *wk1 ,0.0,0.0);
  field->Ic[ii0-1][jj0  ][kk0+1]  += Vector3(Fx1*(0.5+delta_x1)*wj1 *wk_1,0.0,0.0);
  field->Ic[ii0-1][jj0+1][kk0-1]  += Vector3(Fx1*(0.5+delta_x1)*wj_1*wk1_,0.0,0.0);
  field->Ic[ii0-1][jj0+1][kk0]    += Vector3(Fx1*(0.5+delta_x1)*wj_1*wk1 ,0.0,0.0);
  field->Ic[ii0-1][jj0+1][kk0+1]  += Vector3(Fx1*(0.5+delta_x1)*wj_1*wk_1,0.0,0.0);

  field->Ic[ii0][jj0-1][kk0-1]  += Vector3(Fx1*(0.5-delta_x1)*wj1_*wk1_,0.0,0.0);
  field->Ic[ii0][jj0-1][kk0]    += Vector3(Fx1*(0.5-delta_x1)*wj1_*wk1 ,0.0,0.0);
  field->Ic[ii0][jj0-1][kk0+1]  += Vector3(Fx1*(0.5-delta_x1)*wj1_*wk_1,0.0,0.0);
  field->Ic[ii0][jj0  ][kk0-1]  += Vector3(Fx1*(0.5-delta_x1)*wj1 *wk1_,0.0,0.0);
  field->Ic[ii0][jj0  ][kk0]    += Vector3(Fx1*(0.5-delta_x1)*wj1 *wk1 ,0.0,0.0);
  field->Ic[ii0][jj0  ][kk0+1]  += Vector3(Fx1*(0.5-delta_x1)*wj1 *wk_1,0.0,0.0);
  field->Ic[ii0][jj0+1][kk0-1]  += Vector3(Fx1*(0.5-delta_x1)*wj_1*wk1_,0.0,0.0);
  field->Ic[ii0][jj0+1][kk0]    += Vector3(Fx1*(0.5-delta_x1)*wj_1*wk1 ,0.0,0.0);
  field->Ic[ii0][jj0+1][kk0+1]  += Vector3(Fx1*(0.5-delta_x1)*wj_1*wk_1,0.0,0.0);

  // 2nd
  field->Ic[ii-1][jj-1][kk-1]  += Vector3(Fx2*(0.5+delta_x2)*wj2_*wk2_,0.0,0.0);
  field->Ic[ii-1][jj-1][kk]    += Vector3(Fx2*(0.5+delta_x2)*wj2_*wk2 ,0.0,0.0);
  field->Ic[ii-1][jj-1][kk+1]  += Vector3(Fx2*(0.5+delta_x2)*wj2_*wk_2,0.0,0.0);
  field->Ic[ii-1][jj  ][kk-1]  += Vector3(Fx2*(0.5+delta_x2)*wj2 *wk2_,0.0,0.0);
  field->Ic[ii-1][jj  ][kk]    += Vector3(Fx2*(0.5+delta_x2)*wj2 *wk2 ,0.0,0.0);
  field->Ic[ii-1][jj  ][kk+1]  += Vector3(Fx2*(0.5+delta_x2)*wj2 *wk_2,0.0,0.0);
  field->Ic[ii-1][jj+1][kk-1]  += Vector3(Fx2*(0.5+delta_x2)*wj_2*wk2_,0.0,0.0);
  field->Ic[ii-1][jj+1][kk]    += Vector3(Fx2*(0.5+delta_x2)*wj_2*wk2 ,0.0,0.0);
  field->Ic[ii-1][jj+1][kk+1]  += Vector3(Fx2*(0.5+delta_x2)*wj_2*wk_2,0.0,0.0);

  field->Ic[ii][jj-1][kk-1]  += Vector3(Fx2*(0.5-delta_x2)*wj2_*wk2_,0.0,0.0);
  field->Ic[ii][jj-1][kk]    += Vector3(Fx2*(0.5-delta_x2)*wj2_*wk2 ,0.0,0.0);
  field->Ic[ii][jj-1][kk+1]  += Vector3(Fx2*(0.5-delta_x2)*wj2_*wk_2,0.0,0.0);
  field->Ic[ii][jj  ][kk-1]  += Vector3(Fx2*(0.5-delta_x2)*wj2 *wk2_,0.0,0.0);
  field->Ic[ii][jj  ][kk]    += Vector3(Fx2*(0.5-delta_x2)*wj2 *wk2 ,0.0,0.0);
  field->Ic[ii][jj  ][kk+1]  += Vector3(Fx2*(0.5-delta_x2)*wj2 *wk_2,0.0,0.0);
  field->Ic[ii][jj+1][kk-1]  += Vector3(Fx2*(0.5-delta_x2)*wj_2*wk2_,0.0,0.0);
  field->Ic[ii][jj+1][kk]    += Vector3(Fx2*(0.5-delta_x2)*wj_2*wk2 ,0.0,0.0);
  field->Ic[ii][jj+1][kk+1]  += Vector3(Fx2*(0.5-delta_x2)*wj_2*wk_2,0.0,0.0);

#else
  Scalar xr,yr,zr;
  Scalar x1=part->x-pdx;
  Scalar y1=part->y-pdy;
  Scalar z1=part->z-pdz;
  Scalar x2=part->x;
  Scalar y2=part->y;
  Scalar z2=part->z;
  Scalar zn=para->parameter->specy[part->species].zn;

//  i1=(int)FLOOR(x1);
//  j1=(int)FLOOR(y1);
//  k1=(int)FLOOR(z1);

  int i2=(int)FLOOR(x2);
  int j2=(int)FLOOR(y2);
  int k2=(int)FLOOR(z2);
  
  
//  Scalar xr = MIN( MIN(i1,i2)+1 , MAX( MAX(i1,i2),0.5*(x1+x2) ) );
//  Scalar yr = MIN( MIN(j1,j2)+1 , MAX( MAX(j1,j2),0.5*(y1+y2) ) );
//  Scalar zr = MIN( MIN(k1,k2)+1 , MAX( MAX(k1,k2),0.5*(z1+z2) ) );
  
  if(i1==i2)
    xr = 0.5*(x1+x2);
  else
    xr = max(i1,i2);
  if(j1==j2)
    yr = 0.5*(y1+y2);
  else
    yr = max(j1,j2);
  if(k1==k2)
    zr = 0.5*(z1+z2);
  else
    zr = max(k1,k2);
  
  register Scalar wx1 = 0.5*(x1+xr)-i1;
  register Scalar wy1 = 0.5*(y1+yr)-j1;
  register Scalar wz1 = 0.5*(z1+zr)-k1;
  register Scalar wx1_ = 1.0-wx1;
  register Scalar wy1_ = 1.0-wy1;
  register Scalar wz1_ = 1.0-wz1;

  register Scalar wx2 = 0.5*(xr+x2)-i2;
  register Scalar wy2 = 0.5*(yr+y2)-j2;
  register Scalar wz2 = 0.5*(zr+z2)-k2;
  register Scalar wx2_ = 1.0-wx2;
  register Scalar wy2_ = 1.0-wy2;
  register Scalar wz2_ = 1.0-wz2;

  register Scalar Fx1 = iscale*iscale*zn*(xr-x1)*2.0; // dx/dt=2.0
  register Scalar Fy1 = iscale*zn*(yr-y1)*2.0;
  register Scalar Fz1 = iscale*zn*(zr-z1)*2.0;

  register Scalar Fx2 = iscale*iscale*zn*(x2-xr)*2.0;
  register Scalar Fy2 = iscale*zn*(y2-yr)*2.0;
  register Scalar Fz2 = iscale*zn*(z2-zr)*2.0;
  
  int i1st=(int)FLOOR(0.5*(x1+xr));
  int j1st=(int)FLOOR(0.5*(y1+yr));
  int k1st=(int)FLOOR(0.5*(z1+zr));

  int i2nd=(int)FLOOR(0.5*(x2+xr));
  int j2nd=(int)FLOOR(0.5*(y2+yr));
  int k2nd=(int)FLOOR(0.5*(z2+zr));
  
  field->Ic[i1st][j1st][k1st]      += Vector3(Fx1*wy1_*wz1_,Fy1*wx1_*wz1_,Fz1*wx1_*wy1_);
  field->Ic[i1st][j1st][k1st+1]    += Vector3(Fx1*wy1_*wz1,Fy1*wx1_*wz1,0.0);
  field->Ic[i1st][j1st+1][k1st]    += Vector3(Fx1*wy1*wz1_,0.0,Fz1*wx1_*wy1);
  field->Ic[i1st][j1st+1][k1st+1]  += Vector3(Fx1*wy1*wz1,0.0,0.0);
  field->Ic[i1st+1][j1st][k1st]    += Vector3(0.0,Fy1*wx1*wz1_,Fz1*wx1*wy1_);
  field->Ic[i1st+1][j1st][k1st+1]  += Vector3(0.0,Fy1*wx1*wz1,0.0);
  field->Ic[i1st+1][j1st+1][k1st]  += Vector3(0.0,0.0,Fz1*wx1*wy1);

  
  field->Ic[i2nd][j2nd][k2nd]      += Vector3(Fx2*wy2_* wz2_,Fy2*wx2_*wz2_,Fz2*wx2_*wy2_);
  field->Ic[i2nd][j2nd][k2nd+1]    += Vector3(Fx2*wy2_*wz2,Fy2*wx2_*wz2,0.0);
  field->Ic[i2nd][j2nd+1][k2nd]    += Vector3(Fx2*wy2*wz2_,0.0,Fz2*wx2_*wy2);
  field->Ic[i2nd][j2nd+1][k2nd+1]  += Vector3(Fx2*wy2*wz2,0.0,0.0);
  field->Ic[i2nd+1][j2nd][k2nd]    += Vector3(0.0,Fy2*wx1*wz2_,Fz2*wx2*wy2_);
  field->Ic[i2nd+1][j2nd][k2nd+1]  += Vector3(0.0,Fy2*wx1*wz2,0.0);
  field->Ic[i2nd+1][j2nd+1][k2nd]  += Vector3(0.0,0.0,Fz2*wx2*wy2);
#endif
}

void CellInfo::do_change_cell()
{
  stack_member *stack, *stack_delete;
  struct particle *part;
  int i,j,k;
  for(i=IMIN-1;i<=IMAX;i++)
    for(j=JMIN-1;j<=JMAX;j++)
      for(k=KMIN-1;k<=KMAX;k++){
	CellsM[i][j][k].insert = CellsM[i][j][k].first;
      }  

  stack = stk->zero->next;
  while( stack != stk->hole ){
    part = stack->part;
    int i0 = stack->i;
    int j0 = stack->j;
    int k0 = stack->k;
    stk->insert_particle(*this,part,i0,j0,k0);
    stack_delete = stack;
    stk->remove_from_stack(stack_delete);
    stack = stk->zero->next;
  }
}

void CellInfo::advance(Scalar t)
{
  int i,j,k,num=0;
  struct particle *part;

  volatile Scalar pdx,pdy,pdz;

  field->clearQI();
  
  IJKLOOP{
    if(CellsM[i][j][k].first!=NULL){
#ifndef HIGHORDERCLOUD
      acceleratefield(i,j,k);
#endif

      for(part=CellsM[i][j][k].first; part!=NULL;part=part->next){
#ifdef DEBUG
	if(!part) { printf("segmentation\n");exit(1); }
	if((int)FLOOR(part->x)!=i || (int)FLOOR(part->y)!=j || (int)FLOOR(part->z)!=k ){
 	  printf("advance:particle not int cell \n");
	  cout<<i<<" "<<j<<" "<<k<<" "<<part->x<<" "<<part->y<<" "<<part->z<<"\n";
	  exit(1);
	}
#endif
	deposit_charge(part,i,j,k);   //not necessary for the local algorithm
	                                      //charge distribution of the preceeding
	                                      //half time step
	accelerate(part,i,j,k);
	
	move( part,i,j,k,pdx,pdy,pdz);
		
	has_to_change_cell(part,i,j,k);
	
	deposit_current(part,i,j,k,pdx,pdy,pdz);
	
	num++;
      }
    }
  }
/*
#ifdef DEBUG
#ifdef PARALLEL
  int temp=num;
  MPI_Allreduce(&temp,&num,1,MPI_INT,MPI_SUM,OOPIC_COMM);
  if(MY_RANK ==0) 
#endif
    cout<<"num="<<num<<" ";
  if(n_part!=num){ printf(" advance: num of particle is run away"); exit(1);}
#endif
*/
  do_change_cell();
}

void CellInfo::shift_particles(int direction)
{
    int i,j,k;
    int I0,I1,J0,J1,K0,K1;
    struct particle *part;

    I0 = IMIN-1; I1=IMAX;
    J0 = JMIN-1; J1=JMAX;
    K0 = KMIN-1; K1=KMAX-2;
    if(KMIN == para->parameter->ZMIN) K0=KMIN;
    //copy k+1 cell infomation to k cell
    for(i=I0; i<=I1; i++)
      for(j=J0; j<=J1; j++)
	for(k=K0; k<=K1; k++)
    	  CellsM[i][j][k] = CellsM[i][j][k+1];

    //change part->z
    for(i=I0; i<=I1; i++)
    for(j=J0; j<=J1; j++)
    for(k=K0; k<=K1; k++)
    {
    if(CellsM[i][j][k].first!=NULL){
      for(part=CellsM[i][j][k].first; part!=NULL;part=part->next){
	part->z = part->z-1.0;
#ifdef DEBUG
	if((int)FLOOR(part->x)!=i || (int)FLOOR(part->y)!=j || (int)FLOOR(part->z)!=k ){
 	  printf("shift particle :particle not int cell \n");
	  cout<<i<<" "<<j<<" "<<k<<" "<<part->x<<" "<<part->y<<" "<<part->z<<"\n";
	  exit(1);
	}
#endif
      }
    }
    }

    //load particle in cell of ZMAX-1
    if(KMAX == para->parameter->ZMAX) shift_chainparticles(1);
    if(KMAX == para->parameter->ZMAX) shift_initparticles(1);

    int ID_zl = para->ID_zl; 
    int ID_zr = para->ID_zr;
    // exchange particle along z-dirction
    ptsndrcv(IMIN-1,  IMAX,  JMIN-1,  JMAX,  KMIN-1,  KMIN-1,  ID_zl,
	     IMIN-1,  IMAX,  JMIN-1,  JMAX,  KMAX-1,  KMAX-1,  ID_zr);

}

void CellInfo::shift_chainparticles(int di)
{
  int i,j,k,ppc[2];
  int ITemp1,ITemp2,ITemp3;
  Scalar delta,delta2,sqrt3,sqrt32;
  struct particle *pn,*po;

  ppc[0] = para->parameter->ppc[0];
  ppc[1] = para->parameter->ppc[1];
  
  k=KMAX-1;
  for(i=IMIN;i<IMAX;i++)
  for(j=JMIN;j<JMAX;j++)
  {
#ifdef DEBUG
    if(CellsM[i][j][k].first !=NULL || CellsM[i][j][k].last !=NULL){
         printf("shift particle, cell of KMAX-1 is not empty");
         exit(1);
    }
#endif
    CellsM[i][j][k].first = NULL;
    CellsM[i][j][k].last  = NULL;

    po = CellsM[i][j][k].first;

    for(int ispe=0;ispe<2;ispe++){
      int number = (int) FLOOR(CellsM[i][j][k].dens[ispe]*ppc[ispe]+0.5);
      CellsM[i][j][k].np[ispe] = number;
      if(ispe==0) n_el      += number;
      else        n_ion     += number;
      CellsM[i][j][k].npart += number;
      n_part                += number;
      if(number!=0){
	sqrt3 = int(CEIL(POW(number,1.0/3.0)));
	sqrt32= sqrt3*sqrt3;
	delta = 1.0 / sqrt3;
	
	for(int index=1;index<=number;index++){
	  pn          = new ( struct particle );
	  if (!pn) {printf("allocation error\n"); exit(1);}
	  
	  pn->prev    = po;
	  if (po==NULL) CellsM[i][j][k].first = pn;
	  else          po->next    = pn;
	  pn->species = ispe;
	  ITemp1      = (int)FLOOR((index-1)/sqrt32);
	  ITemp2      = (int)FLOOR((index-1)/sqrt3)-ITemp1*sqrt3;
	  ITemp3      = MAX(0,int(index -1 -ITemp1*sqrt32 - ITemp2*sqrt3));
	  pn->x       = Scalar(i) + ((Scalar)ITemp1+0.49999999) * delta;
	  pn->y       = Scalar(j) + ((Scalar)ITemp2+0.49999999) * delta;
	  pn->z       = Scalar(k) + ((Scalar)ITemp3+0.49999999) * delta;
#ifdef DEBUG
	  if((int)FLOOR(pn->x)!=i || (int)FLOOR(pn->y)!=j || (int)FLOOR(pn->z)!=k){
	    printf("initial particle out cell\n");
	    exit(1);
	  }
#endif 
	  po          = pn;
	}
      }
    }
    if (po!=NULL) {
      po->next = NULL;
      CellsM[i][j][k].last = po;
    }
  }
  //#ifdef PARALLEL
  //      int temp;
  //      MPI_Allreduce(&n_part,&temp,1,MPI_INT,MPI_SUM,OOPIC_COMM);
  //      n_part = temp;
  //#endif
}

void CellInfo::shift_initparticles(int di)
{
  int i,j,k,index;
  int number[2];
  Scalar vx,vy,vz;
  struct particle *part;
  SerParam *parameter = para->parameter;

  k=KMAX-1;
  for(i=IMIN;i<IMAX;i++)
  for(j=JMIN;j<JMAX;j++)
  {
    if(CellsM[i][j][k].first!=NULL){
      for(part=CellsM[i][j][k].first; part!=NULL;part=part->next){
#ifdef DEBUG
	if(!part) { printf("segmentation\n");exit(1); }
	if((int)FLOOR(part->x)!=i || (int)FLOOR(part->y)!=j || (int)FLOOR(part->z)!=k ){
 	  printf("InitParticle:particle not int cell \n");
	  cout<<i<<" "<<j<<" "<<k<<" "<<part->x<<" "<<part->y<<" "<<part->z<<"\n";
	  exit(1);
	}
#endif
		
	vx = parameter->vtherm[part->species]*gauss_rand48();
	vy = parameter->vtherm[part->species]*gauss_rand48();
	vz = parameter->vtherm[part->species]*gauss_rand48();
	if(fabs(vx)>=1 || fabs(vy)>=1 || fabs(vz)>=1){
	  printf("initial velocity is wrong!\n");
	  exit(1);
	}
	// determine gamma*v
	
	part->igamma  = sqrt( 1.0 - vx*vx - vy*vy - vz*vz );
	part->ux      = vx / part->igamma;
	part->uy      = vy / part->igamma;
	part->uz      = vz / part->igamma;
	
      }
    }
  }
}

void CellInfo::periodic_particles()
{
  particle *part;

  //************************periodic along x-direction at XMIN
  if(IMIN == para->parameter->XMIN){
    for(int j=JMIN-1;j<=JMAX;j++)
      for(int k=KMIN-1;k<=KMAX;k++)
	for(part=CellsM[IMIN-1][j][k].first; part!=NULL; part=part->next){
	  part->x  = part->x+para->parameter->XMAX ;
	}
  }
  //periodic along x-direction at XMAX
  if(IMAX == para->parameter->XMAX){
    for(int j=JMIN-1;j<=JMAX;j++)
      for(int k=KMIN-1;k<=KMAX;k++)
	for(part=CellsM[IMAX][j][k].first; part!=NULL; part=part->next){
	  part->x  = part->x-para->parameter->XMAX;
	}
  }

  int ID_xl = para->ID_xl; 
  int ID_xr = para->ID_xr;

  // exchange particle along x-direction
  ptsndrcv(IMIN-1,  IMIN-1,JMIN-1,  JMAX,  KMIN-1,  KMAX,  ID_xl,
	   IMAX-1,  IMAX-1,JMIN-1,  JMAX,  KMIN-1,  KMAX,  ID_xr);

  ptsndrcv(IMAX ,  IMAX,   JMIN-1,  JMAX,  KMIN-1,  KMAX,  ID_xr,
	   IMIN ,  IMIN,   JMIN-1,  JMAX,  KMIN-1,  KMAX,  ID_xl);


  //*******************periodic along y-direction at YMIN
  if(JMIN == para->parameter->YMIN){
    for(int i=IMIN-1;i<=IMAX;i++)
      for(int k=KMIN-1;k<=KMAX;k++)
	for(part=CellsM[i][JMIN-1][k].first; part!=NULL; part=part->next){
	  part->y  = part->y+(para->parameter->YMAX - para->parameter->YMIN);
	}
  }
  //periodic along y-direction at YMAX
  if(JMAX == para->parameter->YMAX){
    for(int i=IMIN-1;i<=IMAX;i++)
      for(int k=KMIN-1;k<=KMAX;k++)
	for(part=CellsM[i][JMAX][k].first; part!=NULL; part=part->next){
	  part->y  = part->y-(para->parameter->YMAX - para->parameter->YMIN);
	}
  }

  int ID_yl = para->ID_yl; 
  int ID_yr = para->ID_yr;

  // exchange particle along y-dirction
  ptsndrcv(IMIN-1,  IMAX,  JMIN-1,  JMIN-1,KMIN-1,  KMAX,  ID_yl,
  	   IMIN-1,  IMAX,  JMAX-1,  JMAX-1,KMIN-1,  KMAX,  ID_yr);
  
  ptsndrcv(IMIN-1,  IMAX,  JMAX,  JMAX,  KMIN-1,  KMAX,  ID_yr,
  	   IMIN-1,  IMAX,  JMIN,  JMIN,  KMIN-1,  KMAX,  ID_yl);
}

void CellInfo::reflect_particles()
{
  struct particle *part;
  //left face reflect
  if(KMIN == para->parameter->ZMIN){
    for(int i=IMIN-2;i<=IMAX+1;i++)
      for(int j=JMIN-2;j<=JMAX+1;j++)
          CellsM[i][j][KMIN-1].delete_particles();
/*
	for(part=CellsM[i][j][KMIN-1].first; part!=NULL; part=part->next){
	  part->uz = -part->uz;
	  part->dx = 0.0;
	  part->dy = 0.0;
	  part->dz = Scalar(KMIN)+1.0e-8 - part->z;
	  part->z  = Scalar(KMIN)+1.0e-8;
	  stk->put_on_stack(part);
	}
*/
  }
  //right face reflect
  if(KMAX == para->parameter->ZMAX){
    for(int i=IMIN-2;i<=IMAX+1;i++)
      for(int j=JMIN-2;j<=JMAX+1;j++)
          CellsM[i][j][KMAX].delete_particles();
/*
	for(part=CellsM[i][j][KMAX].first; part!=NULL; part=part->next){
	  part->uz = -part->uz;
	  part->dx = 0.0;
	  part->dy = 0.0;
	  part->dz = Scalar(KMAX)-(1.0e-8) - part->z;
	  part->z  = Scalar(KMAX)-(1.0e-8);
	  stk->put_on_stack(part);
	}
*/
  }


  do_change_cell();

  int ID_zl = para->ID_zl; 
  int ID_zr = para->ID_zr;
  // exchange particle along z-dirction
  ptsndrcv(IMIN-1,  IMAX,  JMIN-1,  JMAX,  KMIN-1,  KMIN-1,  ID_zl,
      	   IMIN-1,  IMAX,  JMIN-1,  JMAX,  KMAX-1,  KMAX-1,  ID_zr);
  
  ptsndrcv(IMIN-1,  IMAX,  JMIN-1,  JMAX,  KMAX,  KMAX,  ID_zr,
      	   IMIN-1,  IMAX,  JMIN-1,  JMAX,  KMIN,  KMIN,  ID_zl);
  

}

void CellInfo::exchange_particles()
{
  int ID_xl = para->ID_xl; 
  int ID_xr = para->ID_xr;
  int ID_yl = para->ID_yl; 
  int ID_yr = para->ID_yr;
  int ID_zl = para->ID_zl; 
  int ID_zr = para->ID_zr;

  // exchange particle along x-direction
  ptsndrcv(IMIN-1,  IMIN-1,JMIN-1,  JMAX,  KMIN-1,  KMAX,  ID_xl,
	   IMAX-1,  IMAX-1,JMIN-1,  JMAX,  KMIN-1,  KMAX,  ID_xr);

  ptsndrcv(IMAX ,  IMAX,   JMIN-1,  JMAX,  KMIN-1,  KMAX,  ID_xr,
	   IMIN ,  IMIN,   JMIN-1,  JMAX,  KMIN-1,  KMAX,  ID_xl);

  // exchange particle along y-dirction
  ptsndrcv(IMIN-1,  IMAX,  JMIN-1,  JMIN-1,KMIN-1,  KMAX,  ID_yl,
  	   IMIN-1,  IMAX,  JMAX-1,  JMAX-1,KMIN-1,  KMAX,  ID_yr);
  
  ptsndrcv(IMIN-1,  IMAX,  JMAX,  JMAX,  KMIN-1,  KMAX,  ID_yr,
  	   IMIN-1,  IMAX,  JMIN,  JMIN,  KMIN-1,  KMAX,  ID_yl);

  
}

void CellInfo::ptsndrcv(int i0 ,int i1 ,int j0 ,int j1 ,int k0 ,int k1 ,int ID_s,
			int ii0,int ii1,int jj0,int jj1,int kk0,int kk1,int ID_r)
{
  int sendnum=0; int recvnum=0;
//  if(ID_s != -1)
    messagepack(i0,i1,j0,j1,k0,k1,sendnum);
  para->commParticle(ID_s,ID_r,sendnum,recvnum);
//  if(ID_r != -1)
    messageunpack(ii0,ii1,jj0,jj1,kk0,kk1,recvnum);

  if(para->sbufpt != NULL)  {delete[] para->sbufpt; para->sbufpt = NULL;}
  if(para->rbufpt != NULL)  {delete[] para->rbufpt; para->rbufpt = NULL;}
}

void CellInfo::messagepack(int i0,int i1,int j0,int j1,int k0,int k1,int& send_num)
{
  int i,j,k;
  int index,cellnum;
  send_num = 0;

  for(i=i0;i<=i1;i++)
    for(j=j0;j<=j1;j++)
      for(k=k0;k<=k1;k++){
	for(particle* part=CellsM[i][j][k].first; part!=NULL; part=part->next){
	  send_num = send_num+1;
	}
      }
  
  int buflength = send_num*8+20;
  //  if(buflength>para->sbufptlth){
  //    printf("particle exchange message buff length is small < %d\n",buflength);
  //    exit(1);
  //  }

  para->sbufpt = new double[buflength];
  double* ptsbuf = para->sbufpt;
  index = 0;
  

  for(i=i0;i<=i1;i++)
    for(j=j0;j<=j1;j++)
      for(k=k0;k<=k1;k++){
	for(particle* part=CellsM[i][j][k].first; part!=NULL; part=part->next){
	  ptsbuf[index++] = part->species;	  
	  ptsbuf[index++] = part->x;
	  ptsbuf[index++] = part->y;
	  ptsbuf[index++] = part->z;
	  ptsbuf[index++] = part->ux;
	  ptsbuf[index++] = part->uy;
	  ptsbuf[index++] = part->uz;
	  ptsbuf[index++] = part->igamma;
	}
      }
  
  /*  
  int position=0;
  MPI_Pack(&send_num,1,MPI_INT,ptsbuf,para->sbufptlth,&position,OOPIC_COMM);
  for(i=i0;i<=i1;i++)
    for(j=j0;j<=j1;j++)
      for(k=k0;k<=k1;k++){
	for(particle* part=CellsM[i][j][k].first; part!=NULL; part=part->next){
	  MPI_Pack(&part->species,1,MPI_INT,ptsbuf,para->sbufptlth,&position,OOPIC_COMM);	  
	  MPI_Pack(&part->x,1,MPI_SCALAR,ptsbuf,para->sbufptlth,&position,OOPIC_COMM);
	  MPI_Pack(&part->y,1,MPI_SCALAR,ptsbuf,para->sbufptlth,&position,OOPIC_COMM);
	  MPI_Pack(&part->z,1,MPI_SCALAR,ptsbuf,para->sbufptlth,&position,OOPIC_COMM);
	  MPI_Pack(&part->ux,1,MPI_SCALAR,ptsbuf,para->sbufptlth,&position,OOPIC_COMM);
	  MPI_Pack(&part->uy,1,MPI_SCALAR,ptsbuf,para->sbufptlth,&position,OOPIC_COMM);
	  MPI_Pack(&part->uz,1,MPI_SCALAR,ptsbuf,para->sbufptlth,&position,OOPIC_COMM);
	  MPI_Pack(&part->igamma,1,MPI_SCALAR,ptsbuf,para->sbufptlth,&position,OOPIC_COMM);
	}
      }
  */
  // delete particles
  for(i=i0;i<=i1;i++)
    for(j=j0;j<=j1;j++)
      for(k=k0;k<=k1;k++){
	CellsM[i][j][k].delete_particles();
      }
}

void CellInfo::messageunpack(int i0,int i1,int j0,int j1,int k0,int k1,int& recv_num)
{
  int i,j,k;
  int index,cellnum,num=0;

  double* ptrbuf = para->rbufpt;
  
  index = 0; cellnum = 0;
  
  //  recv_num = ptrbuf[index++];
  for(int inum=0; inum<recv_num; inum++){
    particle* part = new particle;
    part->species = ptrbuf[index++];	  
    part->x = ptrbuf[index++];
    part->y = ptrbuf[index++];
    part->z = ptrbuf[index++];
    part->ux = ptrbuf[index++];
    part->uy = ptrbuf[index++];
    part->uz = ptrbuf[index++];
    part->igamma = ptrbuf[index++];

    i = (int)FLOOR(part->x);
    j = (int)FLOOR(part->y);
    k = (int)FLOOR(part->z);
#ifdef DEBUG
    if(i<i0 || i>i1 || j<j0 || j>j1 || k<k0 || k>k1 ){
      printf("unpack:particle not int cell \n");
      cout<<i0<<" "<<i1<<")("<<j0<<" "<<j1<<")("<<k0<<" "<<k1<<")"<<part->x<<" "<<part->y<<" "<<part->z<<"\n";
      exit(1);
    }
#endif
    CellsM[i][j][k].insert_particles(part);
  }
}

Scalar CellInfo::energy()
{
  Scalar result=particle_energy()+ field->energy();
  return result;
}

Scalar CellInfo::particle_energy()
{
  int i,j,k;
  particle* part;
  SerParam* parameter=para->parameter; 
  Scalar kinetic = 0.0;
  IJKLOOP{
    if(CellsM[i][j][k].first!=NULL){
      for(part=CellsM[i][j][k].first; part!=NULL;part=part->next){
	kinetic +=parameter->specy[part->species].n * parameter->specy[part->species].m*(1.0/part->igamma-1.0);
      }
    }
  }
  return kinetic;
}

void CellInfo::dump(FILE* file)
{
  char title_start[10]="00000";
  char title_end[10]  ="11111";

  int i,j,k;
  particle* part;
  int num = 0;
  IJKLOOP{
    num += CellsM[i][j][k].npart;
  }

  fwrite(title_start,sizeof(char),10,file);	
  fwrite(&num,sizeof(int),1,file);	

  IJKLOOP{
    if(CellsM[i][j][k].first!=NULL){
      for(part=CellsM[i][j][k].first; part!=NULL;part=part->next){
	int species = part->species;
	Scalar x = part->x;
	Scalar y = part->y;
	Scalar z = part->z;
	Scalar igamma = part->igamma;
	Scalar ux = part->ux;
	Scalar uy = part->uy;
	Scalar uz = part->uz;

	fwrite(&species,sizeof(int),1,file);	
	fwrite(&x  ,sizeof(Scalar),1,file);	
	fwrite(&y  ,sizeof(Scalar),1,file);	
	fwrite(&z  ,sizeof(Scalar),1,file);	
	fwrite(&igamma,sizeof(Scalar),1,file);	
	fwrite(&ux ,sizeof(Scalar),1,file);	
	fwrite(&uy ,sizeof(Scalar),1,file);	
	fwrite(&uz,sizeof(Scalar),1,file);	
      }
    }
  }
  fwrite(title_end,sizeof(char),10,file);	
}

void CellInfo::restore(FILE* file)
{
  char title_start[10];
  char title_end[10];
  int species,fix;
  Scalar ze,m,zm;
  Scalar igamma,x,y,z,ux,uy,uz,n,zn;

  int i,j,k;
  int number=0;

  fread(title_start,sizeof(char),10,file);
  if(strcmp(title_start,"00000")){
    printf("read lpi-particle file title wrong!!");
    exit(0);
  }
  fread(&number,sizeof(int),1,file);	
  
  for(int inum=0; inum<number; inum++){
    fread(&species,sizeof(int),1,file);	
    fread(&x  ,sizeof(Scalar),1,file);	
    fread(&y  ,sizeof(Scalar),1,file);	
    fread(&z  ,sizeof(Scalar),1,file);	
    fread(&igamma,sizeof(Scalar),1,file);	
    fread(&ux ,sizeof(Scalar),1,file);	
    fread(&uy ,sizeof(Scalar),1,file);	
    fread(&uz,sizeof(Scalar),1,file);	
    
    particle* part = new particle;
    part->species = species;	  
    part->x = x;
    part->y = y;
    part->z = z;
    part->ux = ux;
    part->uy = uy;
    part->uz = uz;
    part->igamma = igamma;
    i = (int)FLOOR(part->x);
    j = (int)FLOOR(part->y);
    k = (int)FLOOR(part->z);
#ifdef DEBUG
    if(i<IMIN || i>IMAX-1 || j<JMIN || j>JMAX-1 || k<KMIN || k>KMAX-1 ){
      printf("restore:particle not int cell \n");
      cout<<IMIN<<" "<<IMAX<<")("<<JMIN<<" "<<JMAX<<")("<<KMIN<<" "<<KMAX<<")"<<part->x<<" "<<part->y<<" "<<part->z<<"\n";
      exit(1);
    }
#endif
    CellsM[i][j][k].insert_particles(part);
  }

  fread(title_end,sizeof(char),10,file);
  if(strcmp(title_end,"11111")){
    printf("read lpi-particle file end wrong!!");
    exit(0);
  }

  n_part = number;
#ifdef PARALLEL
      int temp;
      MPI_Allreduce(&n_part,&temp,1,MPI_INT,MPI_SUM,OOPIC_COMM);
      n_part = temp;
#endif
}

ofstream&
operator<<(ofstream& os,const CellInfo& f)
{
  int num=0;
  struct particle *ptpar;
  //  os << "ZONE F=POINT ,I="<<f.IMAX-f.IMIN<<",J="<<f.JMAX-f.JMIN<<",K="<<f.KMAX-f.KMIN<<"\n";
  for(int k=f.KMIN;k<f.KMAX;k++)
    for(int j=f.JMIN;j<f.JMAX;j++)
      for(int i=f.IMIN;i<f.IMAX;i++){
	ptpar = f.CellsM[i][j][k].first;
	while(ptpar!=NULL){
	  os<<ptpar->x<<" "<<ptpar->y<<" "<<ptpar->z<<" "<<ptpar->ux<<" "<<ptpar->uy<<" "<<ptpar->uz<<" "<<"\n";
	  ptpar = ptpar->next;
	  num++;
	}
      }
  //#ifdef DEBUG
  //  if(num !=f.n_part) { printf("num is wrong!\n"),exit(1);}
  //#endif
    return os;
}

