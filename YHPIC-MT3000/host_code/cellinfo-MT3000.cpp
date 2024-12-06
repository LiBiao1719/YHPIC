/*
 * This document only provides the code related to the MT-3000, primarily covering the configuration of the two-level storage structure, particle advection and electromagnetic field interpolation computation, and particle rearrangement
 */

/* Configure the cache size required for particle rearrangement based on the number of particles per process. */
unsigned long cell_particle = 4000;
unsigned long region_index  = 28000;
unsigned long region_buffer = 30000;

// two-level storage structure
void CellInfo::mt_malloc(int cluster_id)
{
	int i,j,k;
	struct particle* part;
	unsigned long index=0;
	unsigned long cell_len;
	unsigned long I_range,J_range,K_range;
	unsigned long size1,size2,size3,size4;
	unsigned long i_cell,j_cell,k_cell;
	unsigned long region_len;

	I_range = IMAX - IMIN;
	J_range = JMAX - JMIN;
	K_range = KMAX - KMIN;
	//cell_len = I_range*J_range*K_range;

	i_cell = I_range/4;
	j_cell = J_range/4;
	k_cell = K_range/4;
	region_len = i_cell*j_cell*k_cell;

	region_max_index = (unsigned long *)malloc(region_len * sizeof(long));

	// Particle data of a 4x4x4 grid is allocated to an accelerator array, treating this grid region as a unit.
	// size1 represents the particle size.
	// size2 represents the particle corresponding index size,
	// size3 represents the buffer data size for the region unit.

	size1 = region_len*(64*cell_particle)*sizeof(Scalar);
	size2 = region_len*(64*cell_particle)*4;
	size3 = region_len*region_buffer*sizeof(Scalar);

	mt_x      = (Scalar*)hthread_malloc(cluster_id,size1,HT_MEM_RW);
	mt_y      = (Scalar*)hthread_malloc(cluster_id,size1,HT_MEM_RW);
	mt_z      = (Scalar*)hthread_malloc(cluster_id,size1,HT_MEM_RW);
	mt_igamma = (Scalar*)hthread_malloc(cluster_id,size1,HT_MEM_RW);
	mt_ux     = (Scalar*)hthread_malloc(cluster_id,size1,HT_MEM_RW);
	mt_uy     = (Scalar*)hthread_malloc(cluster_id,size1,HT_MEM_RW);
	mt_uz     = (Scalar*)hthread_malloc(cluster_id,size1,HT_MEM_RW);
	mt_zm     = (Scalar*)hthread_malloc(cluster_id,size1,HT_MEM_RW);
	mt_zn     = (Scalar*)hthread_malloc(cluster_id,size1,HT_MEM_RW);
	mt_n      = (Scalar*)hthread_malloc(cluster_id,size1,HT_MEM_RW);
	mt_m      = (Scalar*)hthread_malloc(cluster_id,size1,HT_MEM_RW);

	memset(mt_x,  0, size1);
	memset(mt_y,  0, size1);
	memset(mt_z,  0, size1);
	memset(mt_igamma, 0, size1);
	memset(mt_ux, 0, size1);
	memset(mt_uy, 0, size1);
	memset(mt_uz, 0, size1);
	memset(mt_zm, 0, size1);
	memset(mt_zn, 0, size1);
	memset(mt_n,  0, size1);
	memset(mt_m,  0, size1);

	particle_id   = (int*)hthread_malloc(cluster_id,size2,HT_MEM_RW);
	memset(particle_id, 0, size2);

	buffer_x      = (Scalar*)hthread_malloc(cluster_id,size3,HT_MEM_RW);
	buffer_y      = (Scalar*)hthread_malloc(cluster_id,size3,HT_MEM_RW);
	buffer_z      = (Scalar*)hthread_malloc(cluster_id,size3,HT_MEM_RW);
	buffer_ux     = (Scalar*)hthread_malloc(cluster_id,size3,HT_MEM_RW);
	buffer_uy     = (Scalar*)hthread_malloc(cluster_id,size3,HT_MEM_RW);
	buffer_uz     = (Scalar*)hthread_malloc(cluster_id,size3,HT_MEM_RW);
	buffer_zm     = (Scalar*)hthread_malloc(cluster_id,size3,HT_MEM_RW);
	buffer_zn     = (Scalar*)hthread_malloc(cluster_id,size3,HT_MEM_RW);
	buffer_n      = (Scalar*)hthread_malloc(cluster_id,size3,HT_MEM_RW);
	buffer_m      = (Scalar*)hthread_malloc(cluster_id,size3,HT_MEM_RW);

	memset(buffer_x,  0, size3);
	memset(buffer_y,  0, size3);
	memset(buffer_z,  0, size3);
	memset(buffer_ux, 0, size3);
	memset(buffer_uy, 0, size3);
	memset(buffer_uz, 0, size3);
	memset(buffer_zm, 0, size3);
	memset(buffer_zn, 0, size3);
	memset(buffer_n,  0, size3);
	memset(buffer_m,  0, size3);

}

void CellInfo::mt_free()
{
	hthread_free(mt_x       );
	hthread_free(mt_y      );
	hthread_free(mt_z      );
	hthread_free(mt_igamma );
	hthread_free(mt_ux     );
	hthread_free(mt_uy     );
	hthread_free(mt_uz     );
	hthread_free(mt_zm     );
	hthread_free(mt_zn     );
	hthread_free(mt_n     );
	hthread_free(mt_m     );
	hthread_free(particle_id);
	hthread_free(region_max_index);

	hthread_free(buffer_x       );
	hthread_free(buffer_y      );
	hthread_free(buffer_z      );
	hthread_free(buffer_ux     );
	hthread_free(buffer_uy     );
	hthread_free(buffer_uz     );
	hthread_free(buffer_zm     );
	hthread_free(buffer_zn     );
	hthread_free(buffer_n     );
	hthread_free(buffer_m     );
}

// Convert the particle data represented by linked lists into contiguous storage data, implementing a two-level storage scheme. This operation is performed before entering the PIC loop.
void CellInfo::mt_read_particle()
{
	int i,j,k;
	struct particle* part;
	unsigned long index=0;
	int ii,jj,kk;

	unsigned long I_range,J_range,K_range;
	unsigned long size1,cell_len,size2;
	unsigned long i_cell,j_cell,k_cell;
	unsigned long index_region,region;
	unsigned long tmp;
	unsigned long region_len;
	unsigned long part_len;

  
	I_range = IMAX - IMIN;
	J_range = JMAX - JMIN;
	K_range = KMAX - KMIN;

	i_cell = I_range/4;
	j_cell = J_range/4;
	k_cell = K_range/4;
	region_len = i_cell*j_cell*k_cell;

	for(int i=0;i<region_len;i++)
		region_max_index[i] = 64*cell_particle;

	size1 = region_len*(64*cell_particle)*sizeof(Scalar);
	size2 = region_len*(64*cell_particle)*4;

	region = 0;
	for(ii=0;ii<i_cell;ii++)
	for(jj=0;jj<j_cell;jj++)
	for(kk=0;kk<k_cell;kk++)
	{
		index = region*(64*cell_particle);

		for(i=ii*4;i<(ii+1)*4;i++)	
		for(j=jj*4;j<(jj+1)*4;j++)	
		for(k=kk*4;k<(kk+1)*4;k++)	
		{
			
  			if(CellsM[i][j][k].first!=NULL){
    			for(part=CellsM[i][j][k].first; part!=NULL;part=part->next){
       	   		mt_x [index]     = part->x;
       	   		mt_y [index]     = part->y;
       	   		mt_z [index]     = part->z;
       	   		mt_ux[index]     = part->ux;
       	   		mt_uy[index]     = part->uy;
       	   		mt_uz[index]     = part->uz;
       	   		mt_igamma[index] = part->igamma;
       	   		mt_zm[index]     = para->parameter->specy[part->species].zm;
       	   		mt_zn[index]     = para->parameter->specy[part->species].zn;
       	   		mt_n [index]     = para->parameter->specy[part->species].n;
       	   		mt_m [index]     = para->parameter->specy[part->species].m;
            
				particle_id[index]		 = 1;
				index++;
            
				}
			}
			
		}// i,j,k
		region++;
	}//ii,jj,kk
}

/* PIC loop */
//mt_advance()：this function's core computation includes the following two functions.
//1. Particle advection and higher-order interpolation computation
//	mt_cpu_accelerate(mt_x,mt_y,mt_z,mt_ux,mt_uy,mt_uz,mt_zm,mt_igamma);
//
//2.  Particle sorting algorithm based on a 4x4x4 region
//	mt_change_region();


void CellInfo::mt_cpu_accelerate()
{
	int i,j,k,ii,jj,kk,im,jm,km;
	unsigned long index=0;
	int I_dist, J_dist, K_dist;
	unsigned long count;
	unsigned long i_cell,j_cell,k_cell;
	unsigned long region_len,region_dist,tmp1,tmp2,region;
	unsigned long t_size[6];

	unsigned long vector_each_iter = 970; // 742.5KB
	unsigned long particle_each_iter = 16*vector_each_iter; // each_iter_vec * 16;

	I_dist = IMAX-IMIN;
	J_dist = JMAX-JMIN;
	K_dist = KMAX-KMIN;

	t_size[0] = IMIN;
	t_size[1] = IMAX;
	t_size[2] = JMIN;
	t_size[3] = JMAX;
	t_size[4] = KMIN;
	t_size[5] = KMAX;

	i_cell = I_dist/4;
	j_cell = J_dist/4;
	k_cell = K_dist/4;

	// Store the electromagnetic field data in contiguous blocks, with the block sequence corresponding to the order of subregions. 
	// Each block contains the 7x7x7 grid data of a subregion.
	// Within each block, the arrays mt_EB1(/2/3) alternately store the electric field and magnetic field data. 
	// This is done to facilitate double-word data loading in assembly, improving transfer efficiency.
	region = 0;
	long tmp;
	for(ii=0;ii<i_cell;ii++)
	for(jj=0;jj<j_cell;jj++)
	for(kk=0;kk<k_cell;kk++)
	{
		// region 7*7*7 = 343
		// The storage of 6 arrays of data is approximately 16.07KB
		region_dist = 343*2; 
		tmp = 0;

    	for(i=4*ii-1;i<=4*(ii+1)+1;i++)
    	for(j=4*jj-1;j<=4*(jj+1)+1;j++)
    	for(k=4*kk-1;k<=4*(kk+1)+1;k++)
    	{
			mt_EB1[region+2*tmp+0] = field->ENode[i][j][k].e1();
			mt_EB1[region+2*tmp+1] = field->BNode[i][j][k].e1();
			mt_EB2[region+2*tmp+0] = field->ENode[i][j][k].e2();
			mt_EB2[region+2*tmp+1] = field->BNode[i][j][k].e2();
			mt_EB3[region+2*tmp+0] = field->ENode[i][j][k].e3();
			mt_EB3[region+2*tmp+1] = field->BNode[i][j][k].e3();
			tmp++;
		}
		region += region_dist;
	}

	// Two-level iteration computation, corresponding to the operation shown in Figure 4 of the paper.
	unsigned long extern_num = coreNum * cell_particle * 64;
	unsigned long inner_num  = cell_particle*64;

	region_len   = i_cell * j_cell * k_cell;
	int extern_iter = region_len / coreNum; 
	int dsp_res  = region_len % coreNum;

	unsigned long args[17];
	args[0]  = (unsigned long)coreNum;
	args[1]  = (unsigned long)extern_iter;    
	args[2]  = (unsigned long)extern_num;     
	args[3]  = (unsigned long)inner_num;   
	args[4]  = (unsigned long)particle_each_iter  
	args[5]  = (unsigned long)t_size;
	args[6]  = (unsigned long)region_max_index; 
	args[7]  = (unsigned long)mt_EB1;
	args[8]  = (unsigned long)mt_EB2;
	args[9]  = (unsigned long)mt_EB3;
	args[10] = (unsigned long)mt_x;
	args[11] = (unsigned long)mt_y;
	args[12] = (unsigned long)mt_z;
	args[13] = (unsigned long)mt_ux;
	args[14] = (unsigned long)mt_uy;
	args[15] = (unsigned long)mt_uz;
	args[16] = (unsigned long)mt_igamma;

	//Initiate computation operations on the accelerator domain, with the specific implementation in kernel code.
	hthread_group_exec(thread_id, "dsp_mesh_accelerate",5,12,args);
	hthread_group_wait(thread_id);
}



// Particle sort
void CellInfo::mt_change_region()
{
	int i,j,k,ii,jj,kk;
	unsigned long index=0;
	int I_dist, J_dist, K_dist;
	unsigned long count;
	unsigned long i_cell,j_cell,k_cell;
	unsigned long region_len,tmp1,tmp2,region;
	int XMAX,XMIN,YMAX,YMIN,ZMAX,ZMIN;

	I_dist = IMAX-IMIN;
	J_dist = JMAX-JMIN;
	K_dist = KMAX-KMIN;
  
	i_cell = I_dist/4;
	j_cell = J_dist/4;
	k_cell = K_dist/4;
	region_len = i_cell*j_cell*k_cell; 

	unsigned long *region_out_index = (unsigned long*)malloc(region_len*region_index*sizeof(long));
	XMIN = para->parameter->XMIN; 
	XMAX = para->parameter->XMAX; 
	YMIN = para->parameter->YMIN; 
	YMAX = para->parameter->YMAX; 
	ZMIN = para->parameter->ZMIN; 
	ZMAX = para->parameter->ZMAX; 

	region = 0;
	int in_reg; // region index
	int region_index_bound[region_len];
	int region_buffer_bound[region_len];
	
	memset(region_index_bound, 0,region_len*4);
	memset(region_buffer_bound,0,region_len*4);

	// Count the indices of the overflowed particles and store them in the region_out_index array.
	for(ii=0;ii<i_cell;ii++)
	for(jj=0;jj<j_cell;jj++)
	for(kk=0;kk<k_cell;kk++)
	{
		in_reg = ii*j_cell*k_cell + jj*k_cell+kk;

		tmp1 = in_reg*(64*cell_particle);
		tmp2 = tmp1 + region_max_index[in_reg];

    	for(index=tmp1;index<tmp2;index++)
    	{
			i = (int)FLOOR(mt_x[index]);
			j = (int)FLOOR(mt_y[index]);
			k = (int)FLOOR(mt_z[index]);

			if(i < ii*4 || i > ii*4+3 || j< jj*4 || j>jj*4+3 || k<kk*4 || k> kk*4+3)
			{
				/* tmp3 represents the indices of the particles that overflow in the in_reg-th subregion, with tmp3 starting from the 0th data of this subregion and incrementing sequentially. The initial value of region_index_bound[in_reg] is 0, and for each overflowed particle, region_index_bound[in_reg] is incremented by 1. */

				int tmp3 = in_reg*region_index + region_index_bound[in_reg];
				region_out_index[tmp3] = index;
				region_index_bound[in_reg]++;
				particle_id[index] = 0;
			} // if (i < ii*4)
		} // for index
	}// for ii,jj,kk

	//Process the boundary data and store the overflowed data in the buffer.
	for(ii=0;ii<i_cell;ii++)
	for(jj=0;jj<j_cell;jj++)
	for(kk=0;kk<k_cell;kk++)
	{
		in_reg = ii*j_cell*k_cell + jj*k_cell + kk;

		for(int s=0;s<region_index_bound[in_reg];s++)
		{
			index = region_out_index[in_reg*region_index + s];

			i = (int)FLOOR(mt_x[index]);
			j = (int)FLOOR(mt_y[index]);
			k = (int)FLOOR(mt_z[index]);

			//The x and y directions have periodic boundaries, while particles crossing the boundary in the z direction will disappear.
			if( k>= ZMIN && k < ZMAX)
			{
				if(i < XMIN ){
					i = XMAX + i;
					mt_x[index] = mt_x[index] + XMAX;
				}
				if(i == XMAX){
					i = i - XMAX;
					mt_x[index] = mt_x[index] - XMAX;
				}

				if(j < YMIN ){
					j = YMAX + j;
					mt_y[index] = mt_y[index] + YMAX;
				}
				if(j == YMAX){
					j = j - YMAX;
					mt_y[index] = mt_y[index] - YMAX;
				}
				
				int i_1 = i/4;
				int j_1 = j/4;
				int k_1 = k/4;
				int in_reg_ch = i_1*j_cell*k_cell + j_1*k_cell + k_1;

				int tmp4 = in_reg_ch*region_buffer + region_buffer_bound[in_reg_ch];
				buffer_x [tmp4] = mt_x [index];
				buffer_y [tmp4] = mt_y [index];
				buffer_z [tmp4] = mt_z [index];
				buffer_ux[tmp4] = mt_ux[index];
				buffer_uy[tmp4] = mt_uy[index];
				buffer_uz[tmp4] = mt_uz[index];
				buffer_zm[tmp4] = mt_zm[index];
				buffer_zn[tmp4] = mt_zn[index];
				buffer_n [tmp4] = mt_n [index];
				buffer_m [tmp4] = mt_m [index];

				region_buffer_bound[in_reg_ch]++;
			} //if k
		} // for s
	} // for ii

	//Insert the data from the buffer array into the mt_ array
	for(ii=0;ii<i_cell;ii++)
	for(jj=0;jj<j_cell;jj++)
	for(kk=0;kk<k_cell;kk++)
	{
		in_reg = ii*j_cell*k_cell + jj*k_cell+ kk;
		int tmp1 = in_reg*region_buffer; 
		int tmp2 = in_reg*region_index;
		for(int s=0;s<region_buffer_bound[in_reg];s++)
		{
			mt_x [region_out_index[tmp2+s]] = buffer_x [tmp1 + s];
			mt_y [region_out_index[tmp2+s]] = buffer_y [tmp1 + s];
			mt_z [region_out_index[tmp2+s]] = buffer_z [tmp1 + s];
			mt_ux[region_out_index[tmp2+s]] = buffer_ux[tmp1 + s];
			mt_uy[region_out_index[tmp2+s]] = buffer_uy[tmp1 + s];
			mt_uz[region_out_index[tmp2+s]] = buffer_uz[tmp1 + s];
			mt_zm[region_out_index[tmp2+s]] = buffer_zm[tmp1 + s];
			mt_zn[region_out_index[tmp2+s]] = buffer_zn[tmp1 + s];
			mt_n [region_out_index[tmp2+s]] = buffer_n [tmp1 + s];
			mt_m [region_out_index[tmp2+s]] = buffer_m [tmp1 + s];
			particle_id[region_out_index[tmp2+s]] = 1;
		}
	}

	// Place the later particles into the positions of the previously overflowed particles.
	for(ii=0;ii<i_cell;ii++)
	for(jj=0;jj<j_cell;jj++)
	for(kk=0;kk<k_cell;kk++)
	{
	  in_reg = ii*j_cell*k_cell + jj*k_cell+ kk;
	  unsigned long tmp1 = in_reg*region_index;
	  unsigned long tmp2 = in_reg*64*cell_particle + region_max_index[in_reg] -1;
	  for(int s=region_buffer_bound[in_reg];s<region_index_bound[in_reg];s++)
	  {
	    while(particle_id[tmp2] == 0)
	    {
	    	tmp2--;
			region_index_bound[in_reg]--;
	    }
	    mt_x [region_out_index[tmp1+s]] = mt_x [tmp2];
	    mt_y [region_out_index[tmp1+s]] = mt_y [tmp2];
	    mt_z [region_out_index[tmp1+s]] = mt_z [tmp2];
	    mt_ux[region_out_index[tmp1+s]] = mt_ux[tmp2];
	    mt_uy[region_out_index[tmp1+s]] = mt_uy[tmp2];
	    mt_uz[region_out_index[tmp1+s]] = mt_uz[tmp2];
	    mt_zm[region_out_index[tmp1+s]] = mt_zm[tmp2];
	    mt_zn[region_out_index[tmp1+s]] = mt_zn[tmp2];
	    mt_n [region_out_index[tmp1+s]] = mt_n [tmp2];
	    mt_m [region_out_index[tmp1+s]] = mt_m [tmp2];
		particle_id[region_out_index[tmp1+s]] = 1;
		tmp2--;
	  }
	  region_max_index[in_reg] = tmp2  - in_reg*64*cell_particle +1;
	}
}


