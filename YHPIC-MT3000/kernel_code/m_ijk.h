
// 50 pai  vm_rint()
__shared__ void m_ijk(lvector double *x, lvector double *y,lvector double *z, lvector double * mijk)
{
	__asm__ __volatile__(
	"   SMOVI24   0x1,R29           \n\t"
	"   SMULIU.M1 R29,%[px],R30   	\n\t" 
	"|  SMULIU.M2 R29,%[py],R31     \n\t" 

	"  	SMULIU.M1 R29,%[pz],R32    	\n\t"
	"   VMOVI.M1	0x3ff0000000000000,VR63\n\t" // 1.0 -> vr63

	"   SMVAGA.M1 R30,AR0           \n\t" //smvaga 2
	"|  SMVAGA.M2 R31,AR1           \n\t" //smvaga 2

	"   SMVAGA.M1 R32,AR2           \n\t" //smvaga 2

	"		VLDW	  *+AR0[0], VR0				\n\t" // vldw 9
	"|	VLDW	  *+AR1[0], VR1				\n\t" //
	"|	VMOVI.M2	0x3fe0000000000000,VR62\n\t"// 0.5 -> vr62	
	"|	VMOVI.M1	1,VR61\n\t"// 1 -> vr61	

	"		VLDW	  *+AR2[0], VR2				\n\t" 
	" 	SNOP    7       \n\t"

	" 	VFDTRU.M1 VR0,VR10    \n\t" //vfdtru 3 //qu zheng
	" 	VFDTRU.M1 VR1,VR11    \n\t" // 1.4 -> 1 
	" 	VFDTRU.M1 VR2,VR12    \n\t" 

	" 	VFINTD.M1 VR10, VR3   \n\t" // 1-> 1.0
	" 	VFINTD.M1 VR11, VR4   \n\t"
	" 	VFINTD.M1 VR12, VR5   \n\t"

	"		SNOP	2						\n\t"

	" 	VFMULBD.M1  VR0,VR63,VR3,VR6 \n\t" // 1.4 - 1.0 = 0.4
	"|  VFMULBD.M2  VR1,VR63,VR4,VR7\n\t"
	"|  VFMULBD.M3  VR2,VR63,VR5,VR8\n\t" //vr63 ->1.0

	"		SNOP	5						\n\t"

	" 	VFCMPLD.M1 VR6,VR62,VR9	   \n\t"  // < 0.5 取值为1
	"|  VFCMPLD.M2 VR7,VR62,VR10   \n\t"  //VFCMPED 1
	"|  VFCMPLD.M3 VR8,VR62,VR11   \n\t"  // vr62 -> 0.5

	" 	VXOR  		VR9 ,VR61,VR12    \n\t"   //vr61 -> 1
	" 	VXOR  		VR10,VR61,VR13    \n\t"   
	" 	VXOR  		VR11,VR61,VR14    \n\t"   

	"		VFINTD.M1 VR12,VR15			\n\t"
	"		VFINTD.M1 VR13,VR16			\n\t"
	"		VFINTD.M1 VR14,VR17			\n\t"

	" 	SNOP    1       \n\t"
	"		SMOV %[ijk],R30\n\t"
	
	"   VFMULAD.M1  VR15, VR63, VR3, VR18\n\t" // m_ii
	"|  VFMULAD.M2  VR16, VR63, VR4, VR19\n\t" // m_jj
	"|  VFMULAD.M3  VR17, VR63, VR5, VR20\n\t" //m_kk

	"   SMVAGA.M1 R30,AR0           \n\t" //smvaga 2
	" 	SNOP    4       \n\t"
	"		VSTW	VR18,*+AR0[0] \n\t"
	"|	VSTW	VR19,*+AR0[16] \n\t"
	"		VSTW	VR20,*+AR0[32] \n\t"
	" 	SNOP   3        \n\t"
	:
	: [px]   "r" (x),
	  [py]   "r" (y),
	  [pz]   "r" (z),
		[ijk]  "r" (mijk)
	: "R29","R30","R31","R32",\
	  "VR0","VR1","VR2","VR3","VR4","VR5","VR6","VR7","VR8",\
	  "VR9","VR10","VR11","VR12","VR13","VR14","VR15","VR16",\
	  "VR17","VR18","VR19","VR20","VR61","VR62","VR63"
	);
}
