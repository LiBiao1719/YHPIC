#include <hthread_device.h>
#include "kernel.h"
#include "m_ijk.h"

/* This part of the code runs on the MT-3000 accelerator array. 
 * The function performs particle advection and higher-order interpolation of electromagnetic field data.
 */
__global__ void dsp_mesh_accelerate(unsigned long coreNum, unsigned long extern_iter, unsigned long extern_num, unsigned long inner_num,unsigned long particle_each_iter,unsigned long *size, unsigned long *region_max_index, double *mt_EB1, double *mt_EB2,double *mt_EB3,double *mt_x, double *mt_y, double *mt_z, double *mt_ux, double *mt_uy, double *mt_uz, double *mt_ig)
{
	unsigned long core_id  = get_thread_id();	
	unsigned long index    = 0;
	int i_cell,j_cell,k_cell;

	long in[16];
	unsigned long offset,length;
	int it,i,j,k,im,jm,km;
	int ch_no1,ch_no2,ch_no3,ch_no4,ch_no5,ch_no6,ch_no7;

	lvector double zmpidt, ux,uy,uz,igamma,tx,ty,tz;
	lvector double tx2,ty2,tz2,txy,tyz,txz,t_1,ux2,uy2,uz2,t_1_ux,t_1_uy,t_1_uz;
	lvector double delta_x,delta_y,delta_z;
	lvector double ex,ey,ez,bx,by,bz;
	lvector double m_ii,m_jj,m_kk;
	lvector long   n_in,ii,jj,kk,IMIN,IMAX,JMIN,JMAX,KMIN,KMAX;
	lvector double dst,tmp,tmp1,tmp2,tmp3,tmp4,comp_1,comp_2,comp_3;
	lvector double j_dist,k_dist, d_in,d_IMIN,d_IMAX,d_JMIN,d_JMAX,d_KMIN,d_KMAX;

	/* 0x400000000 is the starting address of the AM cache, and the size of the AM space is 768KB.*/

	/* Allocate particle data in the AM space, including three components 
	 * each for spatial coordinates and velocity coordinates.
	 */
	// 970 vector = 0x1EF00
	lvector double * part_x  = (lvector double *)0x400000000;
	lvector double * part_y  = (lvector double *)0x40001E500;
	lvector double * part_z  = (lvector double *)0x40003CA00;
	lvector double * part_ux = (lvector double *)0x40005AF00;
	lvector double * part_uy = (lvector double *)0x400079400;
	lvector double * part_uz = (lvector double *)0x400097900; //tail 0x4000B5E00
	lvector double * S0      = (lvector double *)0x4000B5E00;  // 9 vector
	lvector double * BE_node = (lvector double *)0x4000B6280;  // 6	vector 

	/* Considering higher-order interpolation methods, the interpolation calculation 
	 * for a 4x4x4 grid region requires additional ghost layers, and the actual required 
	 * grid size is a 7x7x7 grid region. 
	 */

  	// Allocate electromagnetic field data storage on the AM, with the storage layout as 
	// shown at the top of Figure 6 in the paper, where the electric field E and magnetic 
	// field B are stored alternately and continuously.
	lvector double * EB1     = (lvector double *)0x4000BB680;  
	lvector double * EB2     = (lvector double *)0x4000BCBF0;
	lvector double * EB3     = (lvector double *)0x4000BE160; // tail 0x4000BF6D0

	//27*6 vectors.
	// The vector storage layout of the electromagnetic fields, as shown in the second 
	// part of Figure 6 in the paper, where these vectors can be mapped one-to-one to vector 
	//registers, facilitating subsequent vector interpolation calculations.
	lvector double * vec 	 = (lvector double *)0x4000B6580;  // tail0x4000BF280

	// 3 vector
	lvector double *mijk     = (lvector double *)0x4000BF6D0;  //tail 0x4000BF850
	long node_id[16];
	lvector double PI 	 	 = vec_movi(3.1415926);
	lvector double DT 		 = vec_movi(0.062500);
	lvector double zm		 = vec_movi(-1.0);

	i_cell = (size[1]-size[0])/4;
	j_cell = (size[3]-size[2])/4;
	k_cell = (size[5]-size[4])/4;

	IMIN   = vec_svbcast(size[0]);
	IMAX   = vec_svbcast(size[1]);
	JMIN   = vec_svbcast(size[2]);
	JMAX   = vec_svbcast(size[3]);
	KMIN   = vec_svbcast(size[4]);
	KMAX   = vec_svbcast(size[5]);
	
	d_IMIN = vec_itf(IMIN);
	d_IMAX = vec_itf(IMAX);
	d_JMIN = vec_itf(JMIN);
	d_JMAX = vec_itf(JMAX);
	d_KMIN = vec_itf(KMIN);
	d_KMAX = vec_itf(KMAX);

	long am_len = 686;  // 7x7x7x2
	unsigned long dsp_iter,dsp_res;
	long start_in;
	lvector long vec_in_reg;
	int e_it,i_it,in_reg;

	
	for(e_it=0;e_it<extern_iter;e_it++)	 // extern iteration
	{
		in_reg = core_id + e_it*coreNum;  // sub-region corresponding to each accelerator array
		long ii = in_reg/(j_cell*k_cell);
		long jj = (in_reg-ii*i_cell*k_cell)/k_cell;
		long kk = in_reg - ii*i_cell*k_cell - jj*k_cell;
		long i_min = 4*ii-1;
		long j_min = 4*jj-1;
		long k_min = 4*kk-1;

		lvector long vec_i_min = vec_svbcast(i_min);
		lvector long vec_j_min = vec_svbcast(j_min);
		lvector long vec_k_min = vec_svbcast(k_min);
		lvector double d_vec_i_min = vec_itf(vec_i_min);
		lvector double d_vec_j_min = vec_itf(vec_j_min);
		lvector double d_vec_k_min = vec_itf(vec_k_min);

		long extern_offset = e_it*extern_num;
		long inner_offset  = inner_num *core_id;

		dsp_iter = region_max_index[core_id]/particle_each_iter;
		dsp_res  = region_max_index[core_id]%particle_each_iter;

		start_in = in_reg * am_len;

		//Load the electromagetic field data, stored alternately and continuously in CPU memory, into the AM space via DMA
		// dma mt_EB1/2/3 in sub-region of in_reg;
		ch_no1 = dma_p2p(&mt_EB1[start_in],1,am_len*sizeof(double),0,EB1,1,am_len*sizeof(double),0,0,0);
		ch_no2 = dma_p2p(&mt_EB2[start_in],1,am_len*sizeof(double),0,EB2,1,am_len*sizeof(double),0,0,0);
		ch_no3 = dma_p2p(&mt_EB3[start_in],1,am_len*sizeof(double),0,EB3,1,am_len*sizeof(double),0,0,0);
 		dma_wait(ch_no1);
		dma_wait(ch_no2);
		dma_wait(ch_no3);
		
		for(i_it=0;i_it<=dsp_iter;i_it++)  // inner iteration
		{
			if(i_it < dsp_iter)  // divisible parts
			{
				length = particle_each_iter;
			}
			else   // the remainder
			{
				if(dsp_res%16 == 0)
					length = dsp_res;
				else
					length = (dsp_res/16+1)*16;
			}
			offset = extern_offset + inner_offset + particle_each_iter * i_it;

			//Load the particle data, stored contiguously in CPU memory, into the AM space via DMA.
			ch_no1 = raw_dma_p2p(&mt_x [offset],1,length*sizeof(double),0,part_x ,1,length*sizeof(double),0,false,0,0,0) ;
			ch_no2 = raw_dma_p2p(&mt_y [offset],1,length*sizeof(double),0,part_y ,1,length*sizeof(double),0,false,0,0,1) ;
			ch_no3 = raw_dma_p2p(&mt_z [offset],1,length*sizeof(double),0,part_z ,1,length*sizeof(double),0,false,0,0,2) ;
			ch_no4 = raw_dma_p2p(&mt_ux[offset],1,length*sizeof(double),0,part_ux,1,length*sizeof(double),0,false,0,0,3) ;
			ch_no5 = raw_dma_p2p(&mt_uy[offset],1,length*sizeof(double),0,part_uy,1,length*sizeof(double),0,false,0,0,4) ;
			ch_no6 = raw_dma_p2p(&mt_uz[offset],1,length*sizeof(double),0,part_uz,1,length*sizeof(double),0,false,0,0,5) ;
 			dma_wait(ch_no1);
 			dma_wait(ch_no2);
 			dma_wait(ch_no3);
 			dma_wait(ch_no4);
 			dma_wait(ch_no5);
 			dma_wait(ch_no6);

			lvector double dd_1,dd_2,dd_3;

			for(index=0;index < length/16;index++)
			{
				tmp1   = vec_muli(PI,DT);
				zmpidt = vec_muli(zm, tmp1);

				m_ijk(&part_x[index],&part_y[index],&part_z[index],mijk); //vm_rintd16
				
				// mijk.h
				// Implement the vector rint() operation using assembly instructions of the MT-3000 architecture
				m_ii = mijk[0];    //  m_ii = vector_rintd16(part_x[index]);
				m_jj = mijk[1];    //  m_jj = vector_rintd16(part_y[index]);
				m_kk = mijk[2];    //  m_kk = vector_rintd16(part_z[index]);

				//d_in = (m_ii - d_vec_i_min)*7*7 + (m_jj - d_vec_j_min)*7 + m_kk - d_vec_k_min
				dd_1 = m_ii - d_vec_i_min;
				d_in = dd_1*49.0;
				dd_2 = m_jj - d_vec_j_min;
				dd_2 = dd_2*7.0;
				d_in = d_in + dd_2;
				dd_3 = m_kk -d_KMIN;
				d_in = d_in + dd_3;
				d_in = 2*d_in;         // E and B ,need multiply 2
				n_in = vec_fti(d_in);
   
				mov_to_svr_v16di(n_in);  // the index of E or B
				node_id[0]  = mov_from_svr0();
				node_id[1]  = mov_from_svr1();
				node_id[2]  = mov_from_svr2();
				node_id[3]  = mov_from_svr3();
				node_id[4]  = mov_from_svr4();
				node_id[5]  = mov_from_svr5();
				node_id[6]  = mov_from_svr6();
				node_id[7]  = mov_from_svr7();
				node_id[8]  = mov_from_svr8();
				node_id[9]  = mov_from_svr9();
				node_id[10] = mov_from_svr10();
				node_id[11] = mov_from_svr11();
				node_id[12] = mov_from_svr12();
				node_id[13] = mov_from_svr13();
				node_id[14] = mov_from_svr14();
				node_id[15] = mov_from_svr15();

				delta_x = m_ii - part_x[index];
				S0[1] = 0.75 -   delta_x   *   delta_x; 
				S0[2] = 0.5*(0.5-delta_x)*(0.5-delta_x); 
				S0[0] = 0.5*(0.5+delta_x)*(0.5+delta_x); 
				delta_y = m_jj - part_y[index];
				S0[4] = 0.75-    delta_y *     delta_y; 
				S0[5] = 0.5*(0.5-delta_y)*(0.5-delta_y); 
				S0[3] = 0.5*(0.5+delta_y)*(0.5+delta_y); 
				delta_z = m_kk - part_z[index];
				S0[7] = 0.75-    delta_z *     delta_z; 
				S0[8] = 0.5*(0.5-delta_z)*(0.5-delta_z); 
				S0[6] = 0.5*(0.5+delta_z)*(0.5+delta_z); 

				// Implement the operation shown in Figure 6 of the paper using assembly instructions of the MT-3000 architecture.
				matrix_trans(EB1,EB2,EB3,vec,node_id);

				// Higher-order interpolation vector computation based on assembly instructions of the MT-3000 architecture.
				S0_multiply_eb(vec,BE_node,S0);

				ex     = BE_node[0];
				ey     = BE_node[2];
				ez     = BE_node[4];
				bx     = BE_node[1];
				by     = BE_node[3];
				bz     = BE_node[5];


				ux     = part_ux[index] + ex * zmpidt;
				uy     = part_uy[index] + ey * zmpidt;
				uz     = part_uz[index] + ez * zmpidt;
				igamma = vec_frsq(1.0 + (ux*ux+uy*uy+uz*uz));

				tx     = bx *zmpidt*igamma;
				ty     = by *zmpidt*igamma;
				tz     = bz *zmpidt*igamma;
				tx2    = tx*tx;
				ty2    = ty*ty;
				tz2    = tz*tz;

				t_1    = vec_frcp(1.0+ tx2+ty2+tz2);
				txy    = tx*ty;
				txz    = tx*tz;
				tyz    = ty*tz;
				ux2    = t_1*( ux*(1.0+tx2-ty2-tz2) + 2.0*uy*(txy+tz) + 2.0*uz*(txz-ty) );
				uy2    = t_1*( 2.0*ux*(txy-tz) + uy*(1.0-tx2+ty2-tz2) + 2.0*uz*(tyz+tx) );
				uz2    = t_1*( 2.0*ux*(txz+ty) + 2.0*uy*(tyz-tx) + uz*(1.0-tx2-ty2+tz2) );
				ux     = ux2 + ex * zmpidt;
				uy     = uy2 + ey * zmpidt;
				uz     = uz2 + ez * zmpidt;
				
//				igamma         = vm_sqrtd16(1.0+ux*ux+uy*uy+uz*uz);
//				part_x[index]  = vm_frecd16(igamma);
				part_x[index]  = vec_frsq(1.0+ux*ux+uy*uy+uz*uz); // mt_igamma[index]
				part_ux[index] = ux;
				part_uy[index] = uy;
				part_uz[index] = uz;
			} //for length
			vector_store(part_ux, &mt_ux[offset], length*sizeof(double));
	  		vector_store(part_uy, &mt_uy[offset], length*sizeof(double));
			vector_store(part_uz, &mt_uz[offset], length*sizeof(double));
			vector_store(part_x,  &mt_ig[offset], length*sizeof(double));
		} //for i_it
	} //for e_it
}


// two-buffer version
#define compute(ptr,l_index,EB1,EB2,EB3) \
{\
		for(index=0;index < half_len/16;index++){\
			tmp1   = vec_muli(PI,DT);\
			zmpidt = vec_muli((part_zm+ptr)[index], tmp1);\
			m_ii   = vm_rintd16((part_x+ptr)[index]);\
			m_jj   = vm_rintd16((part_y+ptr)[index]);\
			m_kk   = vm_rintd16((part_z+ptr)[index]);\
\
			dd_1 = m_ii - d_IMIN;\
			dd_1 = dd_1 +1.0;\
			d_in = dd_1*225.0;\
			\
			dd_2 = m_jj - d_JMIN;\
			dd_2 = dd_2 +1.0;\
			dd_2 = dd_2*15.0;\
			d_in = d_in + dd_2;\
			\
			dd_3 = m_kk -d_KMIN;\
			dd_3 = dd_3 + 1.0;\
			d_in = d_in + dd_3;\
			d_in = 2*d_in-l_index; \
			n_in = vec_fti(d_in);\
   \
			mov_to_svr_v16di(n_in);\
			node_id[0]  = mov_from_svr0();\
			node_id[1]  = mov_from_svr1();\
			node_id[2]  = mov_from_svr2();\
			node_id[3]  = mov_from_svr3();\
			node_id[4]  = mov_from_svr4();\
			node_id[5]  = mov_from_svr5();\
			node_id[6]  = mov_from_svr6();\
			node_id[7]  = mov_from_svr7();\
			node_id[8]  = mov_from_svr8();\
			node_id[9]  = mov_from_svr9();\
			node_id[10] = mov_from_svr10();\
			node_id[11] = mov_from_svr11();\
			node_id[12] = mov_from_svr12();\
			node_id[13] = mov_from_svr13();\
			node_id[14] = mov_from_svr14();\
			node_id[15] = mov_from_svr15();\
\
			delta_x = m_ii - (part_x+ptr)[index];\
			S0[1] = 0.75 -   delta_x   *   delta_x;\
			S0[2] = 0.5*(0.5-delta_x)*(0.5-delta_x);\
			S0[0] = 0.5*(0.5+delta_x)*(0.5+delta_x);\
			delta_y = m_jj - (part_y+ptr)[index];\
			S0[4] = 0.75-    delta_y *     delta_y;\
			S0[5] = 0.5*(0.5-delta_y)*(0.5-delta_y);\
			S0[3] = 0.5*(0.5+delta_y)*(0.5+delta_y);\
			delta_z = m_kk - (part_z+ptr)[index];\
			S0[7] = 0.75-    delta_z *     delta_z;\
			S0[8] = 0.5*(0.5-delta_z)*(0.5-delta_z);\
			S0[6] = 0.5*(0.5+delta_z)*(0.5+delta_z);\
			matrix_trans(EB1,EB2,EB3,vec,node_id);\
			S0_multiply_eb(vec,BE_node,S0);\
\
			ex     = BE_node[0];\
			ey     = BE_node[2];\
			ez     = BE_node[4];\
			bx     = BE_node[1];\
			by     = BE_node[3];\
			bz     = BE_node[5];\
    \
			ux     = (part_ux+ptr)[index] + ex * zmpidt;\
			uy     = (part_uy+ptr)[index] + ey * zmpidt;\
			uz     = (part_uz+ptr)[index] + ez * zmpidt;\
			igamma = vec_frsq(1.0 + (ux*ux+uy*uy+uz*uz));\
			tx     = bx *zmpidt*igamma;\
			ty     = by *zmpidt*igamma;\
			tz     = bz *zmpidt*igamma;\
			tx2    = tx*tx;\
			ty2    = ty*ty;\
			tz2    = tz*tz;\
			t_1    = vec_frcp(1.0+ tx2+ty2+tz2);\
			txy    = tx*ty;\
			txz    = tx*tz;\
			tyz    = ty*tz;\
			ux2    = t_1*( ux*(1.0+tx2-ty2-tz2) + 2.0*uy*(txy+tz) + 2.0*uz*(txz-ty) );\
			uy2    = t_1*( 2.0*ux*(txy-tz) + uy*(1.0-tx2+ty2-tz2) + 2.0*uz*(tyz+tx) );\
			uz2    = t_1*( 2.0*ux*(txz+ty) + 2.0*uy*(tyz-tx) + uz*(1.0-tx2-ty2+tz2) );\
			ux     = ux2 + ex * zmpidt;\
			uy     = uy2 + ey * zmpidt;\
			uz     = uz2 + ez * zmpidt;\
			\
			(part_x +ptr)[index]  = vec_frsq(1.0+ux*ux+uy*uy+uz*uz);\
			(part_ux+ptr)[index] = ux;\
			(part_uy+ptr)[index] = uy;\
		  (part_uz+ptr)[index] = uz;\
		}\
}

__global__ void mt_stream_node(unsigned long coreNum, unsigned long dsp_iter_num, unsigned long dsp_iter, unsigned long *size, double *mt_EB1, double *mt_EB2,double *mt_EB3,double *mt_x, double *mt_y, double *mt_z, double *mt_ux, double *mt_uy, double *mt_uz, double *mt_ig, double *mt_zm)
{
	unsigned long core_id  = get_thread_id();	
	unsigned long length   = dsp_iter_num/coreNum;
	unsigned long offset = core_id * length;
	unsigned long index    = 0;
	unsigned long offset_0,offset_1;
	unsigned long row_step1,row_step2,row_step3;
	unsigned long half_len = length/2;

	long in[16];
	long len0,len1,s_ii,s_jj,s_kk,s0_index,start_0;
	long t_ii,t_jj,t_kk,s1_index,start_1;

	int it,i,j,k,im,jm,km;
	int In0_no0,In0_no1,In0_no2,In0_no3,In0_no4,In0_no5,In0_no6,In0_no7;
	int Out0_no0,Out0_no1,Out0_no2;
	int In1_no0,In1_no1,In1_no2,In1_no3,In1_no4,In1_no5,In1_no6,In1_no7;
	int Out1_no0,Out1_no1,Out1_no2;

	double l_start_0,l_start_1;

	lvector double l_index_0,l_index_1,dd_1,dd_2,dd_3;
	lvector double zmpidt, ux,uy,uz,igamma,tx,ty,tz;
	lvector double tx2,ty2,tz2,txy,tyz,txz,t_1,ux2,uy2,uz2,t_1_ux,t_1_uy,t_1_uz;
	lvector double delta_x,delta_y,delta_z;
	lvector double ex,ey,ez,bx,by,bz;
	lvector double m_ii,m_jj,m_kk;
	lvector long   n_in,ii,jj,kk,IMIN,IMAX,JMIN,JMAX,KMIN,KMAX;
	lvector double dst,tmp,tmp1,tmp2,tmp3,tmp4,comp_1,comp_2,comp_3;
	lvector double j_dist,k_dist, d_in,d_IMIN,d_IMAX,d_JMIN,d_JMAX,d_KMIN,d_KMAX;

	// 780 vector = 0x18600
	unsigned long  half_par  = 2730; 
	unsigned long  par_0 		 = 0;
  // half particles buffer = 0x555000,  390 vectors = 0xC300
	// particles 0x400000000~0x4000AAA00;
	lvector double * part_x  = (lvector double *)0x400000000;
	lvector double * part_y  = (lvector double *)0x40000C300;
	lvector double * part_z  = (lvector double *)0x400018600;
	lvector double * part_ux = (lvector double *)0x400024900;
	lvector double * part_uy = (lvector double *)0x400030C00;
	lvector double * part_uz = (lvector double *)0x40003CF00;
	lvector double * part_zm = (lvector double *)0x400049200; // tail 0x400055500


	lvector double * S0      = (lvector double *)0x4000AAA00;  // 9 vector
	lvector double * BE_node = (lvector double *)0x4000AAE80;  // 6	vector 

	//0x5100 bytes = 27*6 vectors
	lvector double * vec 		 = (lvector double *)0x4000AB180;    

  	//  half = 0x2330  
	lvector double *		EB1  = (lvector double *)0x4000B0280;  
	lvector double *		EB2  = (lvector double *)0x4000B48E0;  
	lvector double *		EB3  = (lvector double *)0x4000B8F40;  

	lvector double *	h_EB1  = (lvector double *)0x4000B25B0;  
	lvector double *	h_EB2  = (lvector double *)0x4000B6C10;  
	lvector double *	h_EB3  = (lvector double *)0x4000BB270;  

	long node_id[16];
	lvector double PI 			 = vec_movi(3.1415926);
	lvector double DT 			 = vec_movi(0.062500);

	IMIN   = vec_svbcast(size[0]);
	IMAX   = vec_svbcast(size[1]);
	JMIN   = vec_svbcast(size[2]);
	JMAX   = vec_svbcast(size[3]);
	KMIN   = vec_svbcast(size[4]);
	KMAX   = vec_svbcast(size[5]);
	
	d_IMIN = vec_itf(IMIN);
	d_IMAX = vec_itf(IMAX);
	d_JMIN = vec_itf(JMIN);
	d_JMAX = vec_itf(JMAX);
	d_KMIN = vec_itf(KMIN);
	d_KMAX = vec_itf(KMAX);


//stream 
	offset_0 =  offset;
	offset_1 =  offset  + half_len;
	s_ii = (long)trunc(mt_x[offset_0]);
	s_jj = (long)trunc(mt_y[offset_0]);
	s_kk = (long)trunc(mt_z[offset_0]);

	t_ii = (long)trunc(mt_x[offset_1]);
	t_jj = (long)trunc(mt_y[offset_1]);
	t_kk = (long)trunc(mt_z[offset_1]);

	// E1,B1,E2,B2,E3,B3 edge of the mesh
	if(s_ii == size[1] -1)
	{
		// 4x15x15
		len0 = 2*450;   
	}
	else
	{
		// 5x15x15
		len0 = 2*563; 
	}

	if(t_ii == size[1] -1)
	{
		// 4x15x15
		len1 = 2*450;   
	}
	else
	{
		// 5x15x15
		len1 = 2*563; 
	}

	// (s_ii-1, JMIN-1, KMIN-1)
	s0_index = (s_ii-1-size[0]+1)*225;
	s1_index = (t_ii-1-size[0]+1)*225;
	start_0 = s0_index*2;
	start_1 = s1_index*2;
	l_start_0 = (double)start_0;
	l_start_1 = (double)start_1;
	l_index_0 = vec_svbcast(l_start_0);
	l_index_1 = vec_svbcast(l_start_1);

//stream 0
	row_step1 = &mt_y [offset_0] - &mt_x [offset_0+length];
	row_step2 = &mt_ux[offset_0] - &mt_z [offset_0+length];
	row_step3 = &mt_uz[offset_0] - &mt_uy[offset_0+length];


	In0_no0 = raw_dma_p2p(&mt_x [offset_0],2,half_len*sizeof(double),row_step1,part_x  ,1,length*sizeof(double)  ,0,false,0,STREAM0_INPUT_PRIORITY,0) ;
	In0_no1 = raw_dma_p2p(&mt_z [offset_0],2,half_len*sizeof(double),row_step2,part_z  ,1,length*sizeof(double)  ,0,false,0,STREAM0_INPUT_PRIORITY,1) ;
	In0_no2 = raw_dma_p2p(&mt_uy[offset_0],2,half_len*sizeof(double),row_step3,part_uy ,1,length*sizeof(double)  ,0,false,0,STREAM0_INPUT_PRIORITY,2) ;
	In0_no3 = raw_dma_p2p(&mt_zm[offset_0],1,half_len*sizeof(double),        0,part_zm ,1,half_len*sizeof(double),0,false,0,STREAM0_INPUT_PRIORITY,3) ;

	In0_no4 = raw_dma_p2p(&mt_EB1[start_0],1,len0*sizeof(double),0,EB1,1,len0*sizeof(double),0,false,0,STREAM0_INPUT_PRIORITY,4);
	In0_no5 = raw_dma_p2p(&mt_EB2[start_0],1,len0*sizeof(double),0,EB2,1,len0*sizeof(double),0,false,0,STREAM0_INPUT_PRIORITY,5);
	In0_no6 = raw_dma_p2p(&mt_EB3[start_0],1,len0*sizeof(double),0,EB3,1,len0*sizeof(double),0,false,0,STREAM0_INPUT_PRIORITY,6);

//	unsigned long s_row_step1 = &mt_y [offset_1] - &mt_x [offset_1+length];
//	unsigned long s_row_step2 = &mt_ux[offset_1] - &mt_z [offset_1+length];
//	unsigned long s_row_step3 = &mt_uz[offset_1] - &mt_uy[offset_1+length];
//	hthread_printf("(%lu,%lu,%lu), (%lu,%lu,%lu)\n",row_step1,row_step2,row_step3,s_row_step1,s_row_step2,s_row_step3);

// stream 1
	In1_no0 = raw_dma_p2p(&mt_x [offset_1],2,half_len*sizeof(double),row_step1,part_x +half_par,1,length*sizeof(double)  ,0,false,0,STREAM1_INPUT_PRIORITY,10) ;
	In1_no1 = raw_dma_p2p(&mt_z [offset_1],2,half_len*sizeof(double),row_step2,part_z +half_par,1,length*sizeof(double)  ,0,false,0,STREAM1_INPUT_PRIORITY,11) ;
	In1_no2 = raw_dma_p2p(&mt_uy[offset_1],2,half_len*sizeof(double),row_step3,part_uy+half_par,1,length*sizeof(double)  ,0,false,0,STREAM1_INPUT_PRIORITY,12) ;
	In1_no3 = raw_dma_p2p(&mt_zm[offset_1],1,half_len*sizeof(double),        0,part_zm+half_par,1,half_len*sizeof(double),0,false,0,STREAM1_INPUT_PRIORITY,13) ;

	In1_no4 = raw_dma_p2p(&mt_EB1[start_1],1,len1*sizeof(double),0,h_EB1,1,len1*sizeof(double),0,false,0,STREAM1_INPUT_PRIORITY,14);
	In1_no5 = raw_dma_p2p(&mt_EB2[start_1],1,len1*sizeof(double),0,h_EB2,1,len1*sizeof(double),0,false,0,STREAM1_INPUT_PRIORITY,15);
	In1_no6 = raw_dma_p2p(&mt_EB3[start_1],1,len1*sizeof(double),0,h_EB3,1,len1*sizeof(double),0,false,0,STREAM1_INPUT_PRIORITY,16);

	for(i=0;i<dsp_iter-1;i++)	
	{
		//stream0 wait
		raw_dma_wait((int)In0_no0);
		raw_dma_wait((int)In0_no1);
		raw_dma_wait((int)In0_no2);
		raw_dma_wait((int)In0_no3);
		raw_dma_wait((int)In0_no4);
		raw_dma_wait((int)In0_no5);
		raw_dma_wait((int)In0_no6);

		//stream0 compute
		compute(par_0,l_index_0,EB1,EB2,EB3);
		//stream0 out
		Out0_no0 = raw_dma_p2p(part_x, 1,half_len*sizeof(double),0,&mt_ig[offset_0],1,half_len*sizeof(double),0        ,false,0,STREAM0_OUTPUT_PRIORITY,7);
		Out0_no1 = raw_dma_p2p(part_ux,1,half_len*sizeof(double),0,&mt_ux[offset_0],1,half_len*sizeof(double),0        ,false,0,STREAM0_OUTPUT_PRIORITY,8);
		Out0_no2 = raw_dma_p2p(part_uy,1,length  *sizeof(double),0,&mt_uy[offset_0],2,half_len*sizeof(double),row_step3,false,0,STREAM0_OUTPUT_PRIORITY,9);

		//start the next flow stream0
		offset_0 =  (i+1)*dsp_iter_num + offset;
		s_ii = (long)trunc(mt_x[offset_0]);
		s_jj = (long)trunc(mt_y[offset_0]);
		s_kk = (long)trunc(mt_z[offset_0]);

		// E1,B1,E2,B2,E3,B3 edge of the mesh
		if(s_ii == size[1] -1)
		{
			// 4x15x15
			len0 = 2*450;   
		}
		else
		{
			// 5x15x15
			len0 = 2*563; 
		}
  
		// (s_ii-1, JMIN-1, KMIN-1)
		s0_index = (s_ii-1-size[0]+1)*225;
		start_0 = s0_index*2;
		l_start_0 = (double)start_0;
		l_index_0 = vec_svbcast(l_start_0);

		In0_no0 = raw_dma_p2p(&mt_x [offset_0],2,half_len*sizeof(double),row_step1,part_x  ,1,length*sizeof(double)  ,0,false,0,STREAM0_INPUT_PRIORITY,0) ;
		In0_no1 = raw_dma_p2p(&mt_z [offset_0],2,half_len*sizeof(double),row_step2,part_z  ,1,length*sizeof(double)  ,0,false,0,STREAM0_INPUT_PRIORITY,1) ;
		In0_no2 = raw_dma_p2p(&mt_uy[offset_0],2,half_len*sizeof(double),row_step3,part_uy ,1,length*sizeof(double)  ,0,false,0,STREAM0_INPUT_PRIORITY,2) ;
		In0_no3 = raw_dma_p2p(&mt_zm[offset_0],1,half_len*sizeof(double),        0,part_zm ,1,half_len*sizeof(double),0,false,0,STREAM0_INPUT_PRIORITY,3) ;
  
		In0_no4 = raw_dma_p2p(&mt_EB1[start_0],1,len0*sizeof(double),0,EB1,1,len0*sizeof(double),0,false,0,STREAM0_INPUT_PRIORITY,4);
		In0_no5 = raw_dma_p2p(&mt_EB2[start_0],1,len0*sizeof(double),0,EB2,1,len0*sizeof(double),0,false,0,STREAM0_INPUT_PRIORITY,5);
		In0_no6 = raw_dma_p2p(&mt_EB3[start_0],1,len0*sizeof(double),0,EB3,1,len0*sizeof(double),0,false,0,STREAM0_INPUT_PRIORITY,6);

		//stream1 wait
		raw_dma_wait((int)In1_no0);
		raw_dma_wait((int)In1_no1);
		raw_dma_wait((int)In1_no2);
		raw_dma_wait((int)In1_no3);
		raw_dma_wait((int)In1_no4);
		raw_dma_wait((int)In1_no5);
		raw_dma_wait((int)In1_no6);

		//stream1 compute
		compute(half_par,l_index_1,h_EB1,h_EB2,h_EB3);

		//stream1 out
		Out1_no0 = raw_dma_p2p(part_x +half_par,1,half_len*sizeof(double),0,&mt_ig[offset_1],1,half_len*sizeof(double),0        ,false,0,STREAM1_OUTPUT_PRIORITY,17);
		Out1_no1 = raw_dma_p2p(part_ux+half_par,1,half_len*sizeof(double),0,&mt_ux[offset_1],1,half_len*sizeof(double),0        ,false,0,STREAM1_OUTPUT_PRIORITY,18);
		Out1_no2 = raw_dma_p2p(part_uy+half_par,1,length  *sizeof(double),0,&mt_uy[offset_1],2,half_len*sizeof(double),row_step3,false,0,STREAM1_OUTPUT_PRIORITY,19);

		//start the next flow stream1
		offset_1 =  (i+1)*dsp_iter_num + offset + half_len;
		t_ii = (long)trunc(mt_x[offset_1]);
		t_jj = (long)trunc(mt_y[offset_1]);
		t_kk = (long)trunc(mt_z[offset_1]);
  
		// E1,B1,E2,B2,E3,B3 edge of the mesh
		if(t_ii == size[1] -1)
		{
			// 4x15x15
			len1 = 2*450;   
		}
		else
		{
			// 5x15x15
			len1 = 2*563; 
		}
  
		// (s_ii-1, JMIN-1, KMIN-1)
		s1_index = (t_ii-1-size[0]+1)*225;
		start_1 = s1_index*2;
		l_start_1 = (double)start_1;
		l_index_1 = vec_svbcast(l_start_1);

		// stream 1
		In1_no0 = raw_dma_p2p(&mt_x [offset_1],2,half_len*sizeof(double),row_step1,part_x +half_par,1,length*sizeof(double)  ,0,false,0,STREAM1_INPUT_PRIORITY,10) ;
		In1_no1 = raw_dma_p2p(&mt_z [offset_1],2,half_len*sizeof(double),row_step2,part_z +half_par,1,length*sizeof(double)  ,0,false,0,STREAM1_INPUT_PRIORITY,11) ;
		In1_no2 = raw_dma_p2p(&mt_uy[offset_1],2,half_len*sizeof(double),row_step3,part_uy+half_par,1,length*sizeof(double)  ,0,false,0,STREAM1_INPUT_PRIORITY,12) ;
		In1_no3 = raw_dma_p2p(&mt_zm[offset_1],1,half_len*sizeof(double),        0,part_zm+half_par,1,half_len*sizeof(double),0,false,0,STREAM1_INPUT_PRIORITY,13) ;
  
		In1_no4 = raw_dma_p2p(&mt_EB1[start_1],1,len1*sizeof(double),0,h_EB1,1,len1*sizeof(double),0,false,0,STREAM1_INPUT_PRIORITY,14);
		In1_no5 = raw_dma_p2p(&mt_EB2[start_1],1,len1*sizeof(double),0,h_EB2,1,len1*sizeof(double),0,false,0,STREAM1_INPUT_PRIORITY,15);
		In1_no6 = raw_dma_p2p(&mt_EB3[start_1],1,len1*sizeof(double),0,h_EB3,1,len1*sizeof(double),0,false,0,STREAM1_INPUT_PRIORITY,16);

		//free resource
		raw_dma_wait((int)Out0_no0);
		raw_dma_wait((int)Out0_no1);
		raw_dma_wait((int)Out0_no2);
		raw_dma_wait((int)Out1_no0);
		raw_dma_wait((int)Out1_no1);
		raw_dma_wait((int)Out1_no2);
	}

	// the last iteration of  stream0 and stream1
		//stream0 wait
		raw_dma_wait((int)In0_no0);
		raw_dma_wait((int)In0_no1);
		raw_dma_wait((int)In0_no2);
		raw_dma_wait((int)In0_no3);
		raw_dma_wait((int)In0_no4);
		raw_dma_wait((int)In0_no5);
		raw_dma_wait((int)In0_no6);

		//stream0 compute
		compute(par_0,l_index_0,EB1,EB2,EB3)

		//stream0 out
		Out0_no0 = raw_dma_p2p(part_x, 1,half_len*sizeof(double),0,&mt_ig[offset_0],1,half_len*sizeof(double),0        ,false,0,STREAM0_OUTPUT_PRIORITY,7);
		Out0_no1 = raw_dma_p2p(part_ux,1,half_len*sizeof(double),0,&mt_ux[offset_0],1,half_len*sizeof(double),0        ,false,0,STREAM0_OUTPUT_PRIORITY,8);
		Out0_no2 = raw_dma_p2p(part_uy,1,length  *sizeof(double),0,&mt_uy[offset_0],2,half_len*sizeof(double),row_step3,false,0,STREAM0_OUTPUT_PRIORITY,9);

		//stream1 wait
		raw_dma_wait((int)In1_no0);
		raw_dma_wait((int)In1_no1);
		raw_dma_wait((int)In1_no2);
		raw_dma_wait((int)In1_no3);
		raw_dma_wait((int)In1_no4);
		raw_dma_wait((int)In1_no5);
		raw_dma_wait((int)In1_no6);

		//stream1 compute
		compute(half_par,l_index_1,h_EB1,h_EB2,h_EB3)

		//stream1 out
		Out1_no0 = raw_dma_p2p(part_x +half_par,1,half_len*sizeof(double),0,&mt_ig[offset_1],1,half_len*sizeof(double),0        ,false,0,STREAM0_OUTPUT_PRIORITY,17);
		Out1_no1 = raw_dma_p2p(part_ux+half_par,1,half_len*sizeof(double),0,&mt_ux[offset_1],1,half_len*sizeof(double),0        ,false,0,STREAM0_OUTPUT_PRIORITY,18);
		Out1_no2 = raw_dma_p2p(part_uy+half_par,1,length  *sizeof(double),0,&mt_uy[offset_1],2,half_len*sizeof(double),row_step3,false,0,STREAM0_OUTPUT_PRIORITY,19);

		//free resource
		raw_dma_wait((int)Out0_no0);
		raw_dma_wait((int)Out0_no1);
		raw_dma_wait((int)Out0_no2);
		raw_dma_wait((int)Out1_no0);
		raw_dma_wait((int)Out1_no1);
		raw_dma_wait((int)Out1_no2);
	

}
