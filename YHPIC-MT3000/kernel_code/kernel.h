// Implement the operation shown in Figure 6 of the paper using assembly instructions of the MT-3000 architecture
// The symbol '|' indicates that the instruction and the adjacent previous instruction are issued in the same clock cycle.
__shared__ void matrix_trans(lvector double *EB1, lvector double *EB2, lvector double *EB3,lvector double *vec, long *node_id)
{
			__asm__ __volatile__(
			"   SMOVI24   0x1,R29           \n\t"
			"   SMULIU.M1 R29,%[node],R32   \n\t" //smuliu 3
			"   SMULIU.M1 R29,%[eb1],R30     \n\t" //smuliu 3
			"|  SMULIU.M2 R29,%[Vec],R31    \n\t" //smuliu 3
			"   SNOP      1                 \n\t"
			"   SMVAGA.M1 R32,AR15          \n\t" //smvaga 2
			"   SMVAGA.M1 R30,AR0           \n\t" //smvaga 2
			"|  SMVAGA.M2 R31,AR3           \n\t" //smvaga 2

			"	SLDW	  *+AR15[0], R10	\n\t" //sldw 7
			"   SMVCGC    R29,VLR           \n\t"  //snop 4  
			"	SNOP	  5					\n\t"
			"	SSHFLL	  3, R10, R10       \n\t"  //SSHFLL 1 
			"   SADDA     R10,AR0,AR1		\n\t" //sadda AVAF 3 SVAF 2
			"	SNOP	  2					\n\t"
//EB1 index[0] 
			"   VLDDWM2 *-AR1[57],  VR1:VR0    \n\t" //vr0-0,vr1-1
			"|  VLDDWM2 *-AR1[56],  VR3:VR2    \n\t" //vr2-3,vr3-4

			"   VLDDWM2 *-AR1[55],  VR5:VR4    \n\t" //vr4-6,vr5-7
			"|  VLDDWM2 *-AR1[50],  VR7:VR6    \n\t" //vr6-9,vr7-10

			"   VLDDWM2 *-AR1[49],  VR9:VR8    \n\t" //vr4-12,vr5-13
			"|  VLDDWM2 *-AR1[48],  VR11:VR10	\n\t" //vr10-15,vr11-16

			"   VLDDWM2 *-AR1[43],  VR13:VR12  \n\t" //vr12-18,vr13-19
			"|  VLDDWM2 *-AR1[42],  VR15:VR14  \n\t" //vr14-21,vr15-22

			"   VLDDWM2 *-AR1[41],  VR17:VR16  \n\t" //vr16-24,vr17-25
			"|	VLDDWM2	*-AR1[8],	 VR19:VR18	\n\t" //v18-2   

			"   VLDDWM2 *-AR1[7],  VR21:VR20  \n\t" //vr16-24,vr17-25
			"|  VLDDWM2 *-AR1[6],  VR23:VR22  \n\t" //vr16-24,vr17-25

			"   VLDDWM2 *-AR1[1],  VR25:VR24  \n\t" //vr16-24,vr17-25
			"|  VLDDWM2 *+AR1[0],  VR27:VR26  \n\t" //vr16-24,vr17-25

			"   VLDDWM2 *+AR1[1],   VR29:VR28  \n\t" //vr16-24,vr17-25
			"|  VLDDWM2 *+AR1[6],  VR31:VR30  \n\t" //vr16-24,vr17-25

			"		SNOP	1\n\t"

			"   VSTDW   VR1:VR0, *+AR3[0]         \n\t"
			"|  VLDDWM2 *+AR1[7],  VR33:VR32  \n\t" //vr16-24,vr17-25
			"|	SLDW	  *+AR15[2], R10	\n\t" //sldw 7

			"   VSTDW   VR3:VR2, *+AR3[16]      \n\t"
			"|  VLDDWM2 *+AR1[8],  VR35:VR34  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR5:VR4, *+AR3[32]		\n\t"
			"|  VLDDWM2 *+AR1[41],  VR37:VR36  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR7:VR6, *+AR3[48]      \n\t"
			"|  VLDDWM2 *+AR1[42],  VR39:VR38  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR9:VR8, *+AR3[64]      \n\t"
			"|  VLDDWM2 *+AR1[43],  VR41:VR40  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR11:VR10, *+AR3[80]    \n\t"
			"|  VLDDWM2 *+AR1[48],  VR43:VR42  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR13:VR12, *+AR3[96]   \n\t"
			"|  VLDDWM2 *+AR1[49],  VR45:VR44  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR15:VR14, *+AR3[112]   \n\t"
			"|  VLDDWM2 *+AR1[50],  VR47:VR46  \n\t" //vr16-24,vr17-25
			"|	SSHFLL	  3, R10, R10       \n\t"  //SSHFLL 1  

			"   VSTDW   VR17:VR16, *+AR3[128]   \n\t"
			"|  VLDDWM2 *+AR1[55],  VR49:VR48  \n\t" //vr16-24,vr17-25
			"|  SADDA     R10,AR0,AR2		\n\t" //sadda AVAF 3 SVAF 2

			"   VSTDW   VR19:VR18, *+AR3[144]   \n\t"
			"|  VLDDWM2 *+AR1[56],  VR51:VR50  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR21:VR20, *+AR3[160]   \n\t"
			"|  VLDDWM2 *+AR1[57],  VR53:VR52  \n\t" //vr16-24,vr17-25
////////////////////////////////////////////////////////////////////////////
// index[2] read
			"	VSTDW   VR23:VR22, *+AR3[176]   \n\t"
			"|  VLDDWM2 *-AR2[57],  VR1:VR0    \n\t" //vr0-0,vr1-1

			"	VSTDW   VR25:VR24, *+AR3[192]   \n\t"
			"|  VLDDWM2 *-AR2[56],  VR3:VR2    \n\t" //vr2-3,vr3-4
			
			"	VSTDW   VR27:VR26, *+AR3[208]   \n\t"
			"|  VLDDWM2 *-AR2[55],  VR5:VR4    \n\t" //vr4-6,vr5-7

			"	VSTDW   VR29:VR28, *+AR3[224]   \n\t"
			"|  VLDDWM2 *-AR2[50],  VR7:VR6    \n\t" //vr6-9,vr7-10

			"	VSTDW   VR31:VR30, *+AR3[240]   \n\t"
			"|  VLDDWM2 *-AR2[49],  VR9:VR8    \n\t" //vr4-12,vr5-13

			"   VSTDW   VR33:VR32, *+AR3[256]         \n\t"
			"|  VLDDWM2 *-AR2[48],  VR11:VR10	\n\t" //vr10-15,vr11-16

			"   VSTDW   VR35:VR34, *+AR3[272]         \n\t"
			"|  VLDDWM2 *-AR2[43],  VR13:VR12  \n\t" //vr12-18,vr13-19

			"   VSTDW   VR37:VR36, *+AR3[288]         \n\t"
			"|  VLDDWM2 *-AR2[42],  VR15:VR14  \n\t" //vr14-21,vr15-22

			"   VSTDW   VR39:VR38, *+AR3[304]         \n\t"
			"|  VLDDWM2 *-AR2[41],  VR17:VR16  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR41:VR40, *+AR3[320]         \n\t"
			"|	VLDDWM2	*-AR2[8],	 VR19:VR18	\n\t" //v18-2   

			"   VSTDW   VR43:VR42, *+AR3[336]         \n\t"
			"|  VLDDWM2 *-AR2[7],  VR21:VR20  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR45:VR44, *+AR3[352]         \n\t"
			"|  VLDDWM2 *-AR2[6],  VR23:VR22  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR47:VR46, *+AR3[368]         \n\t"
			"|  VLDDWM2 *-AR2[1],  VR25:VR24  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR49:VR48, *+AR3[384]         \n\t"
			"|  VLDDWM2 *+AR2[0],  VR27:VR26  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR51:VR50, *+AR3[400]         \n\t"
			"|  VLDDWM2 *+AR2[1],   VR29:VR28  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR53:VR52, *+AR3[416]         \n\t"
			"|  VLDDWM2 *+AR2[6],  VR31:VR30  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR1:VR0, *+AR3[1]         \n\t"
			"|  VLDDWM2 *+AR2[7],  VR33:VR32  \n\t" //vr16-24,vr17-25
			"|	SLDW	  *+AR15[4], R10	\n\t" //sldw 7

			"   VSTDW   VR3:VR2, *+AR3[17]      \n\t"
			"|	VLDDWM2 *+AR2[8],  VR35:VR34  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR5:VR4, *+AR3[33]		\n\t"
			"|  VLDDWM2 *+AR2[41],  VR37:VR36  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR7:VR6, *+AR3[49]      \n\t"
			"|  VLDDWM2 *+AR2[42],  VR39:VR38  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR9:VR8, *+AR3[65]      \n\t"
			"|  VLDDWM2 *+AR2[43],  VR41:VR40  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR11:VR10, *+AR3[81]    \n\t"
			"|  VLDDWM2 *+AR2[48],  VR43:VR42  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR13:VR12, *+AR3[97]   \n\t"
			"|  VLDDWM2 *+AR2[49],  VR45:VR44  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR15:VR14, *+AR3[113]   \n\t"
			"|  VLDDWM2 *+AR2[50],  VR47:VR46  \n\t" //vr16-24,vr17-25
			"|	SSHFLL	  3, R10, R10       \n\t"  //SSHFLL 1  

			"   VSTDW   VR17:VR16, *+AR3[129]   \n\t"
			"|  VLDDWM2 *+AR2[55],  VR49:VR48  \n\t" //vr16-24,vr17-25
			"|  SADDA     R10,AR0,AR1		\n\t" //sadda AVAF 3 SVAF 2

			"   VSTDW   VR19:VR18, *+AR3[145]   \n\t"
			"|  VLDDWM2 *+AR2[56],  VR51:VR50  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR21:VR20, *+AR3[161]   \n\t"
			"|  VLDDWM2 *+AR2[57],  VR53:VR52  \n\t" //vr16-24,vr17-25
// index[4] read
			"	VSTDW   VR23:VR22, *+AR3[177]   \n\t"
			"|  VLDDWM2 *-AR1[57],  VR1:VR0    \n\t" //vr0-0,vr1-1

			"	VSTDW   VR25:VR24, *+AR3[193]   \n\t"
			"|  VLDDWM2 *-AR1[56],  VR3:VR2    \n\t" //vr2-3,vr3-4
			
			"	VSTDW   VR27:VR26, *+AR3[209]   \n\t"
			"|  VLDDWM2 *-AR1[55],  VR5:VR4    \n\t" //vr4-6,vr5-7

			"	VSTDW   VR29:VR28, *+AR3[225]   \n\t"
			"|  VLDDWM2 *-AR1[50],  VR7:VR6    \n\t" //vr6-9,vr7-10

			"	VSTDW   VR31:VR30, *+AR3[241]   \n\t"
			"|  VLDDWM2 *-AR1[49],  VR9:VR8    \n\t" //vr4-12,vr5-13

			"   VSTDW   VR33:VR32, *+AR3[257]         \n\t"
			"|  VLDDWM2 *-AR1[48],  VR11:VR10	\n\t" //vr10-15,vr11-16

			"   VSTDW   VR35:VR34, *+AR3[273]         \n\t"
			"|  VLDDWM2 *-AR1[43],  VR13:VR12  \n\t" //vr12-18,vr13-19

			"   VSTDW   VR37:VR36, *+AR3[289]         \n\t"
			"|  VLDDWM2 *-AR1[42],  VR15:VR14  \n\t" //vr14-21,vr15-22

			"   VSTDW   VR39:VR38, *+AR3[305]         \n\t"
			"|  VLDDWM2 *-AR1[41],  VR17:VR16  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR41:VR40, *+AR3[321]         \n\t"
			"|	VLDDWM2	*-AR1[8],	 VR19:VR18	\n\t" //v18-2   

			"   VSTDW   VR43:VR42, *+AR3[337]         \n\t"
			"|  VLDDWM2 *-AR1[7],  VR21:VR20  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR45:VR44, *+AR3[353]         \n\t"
			"|  VLDDWM2 *-AR1[6],  VR23:VR22  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR47:VR46, *+AR3[369]         \n\t"
			"|  VLDDWM2 *-AR1[1],  VR25:VR24  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR49:VR48, *+AR3[385]         \n\t"
			"|  VLDDWM2 *+AR1[0],  VR27:VR26  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR51:VR50, *+AR3[401]         \n\t"
			"|  VLDDWM2 *+AR1[1],   VR29:VR28  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR53:VR52, *+AR3[417]         \n\t"
			"|  VLDDWM2 *+AR1[6],  VR31:VR30  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR1:VR0, *+AR3[2]         \n\t"
			"|  VLDDWM2 *+AR1[7],  VR33:VR32  \n\t" //vr16-24,vr17-25
			"|	SLDW	  *+AR15[6], R10	\n\t" //sldw 7

			"   VSTDW   VR3:VR2, *+AR3[18]      \n\t"
			"|	VLDDWM2 *+AR1[8],  VR35:VR34  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR5:VR4, *+AR3[34]		\n\t"
			"|  VLDDWM2 *+AR1[41],  VR37:VR36  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR7:VR6, *+AR3[50]      \n\t"
			"|  VLDDWM2 *+AR1[42],  VR39:VR38  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR9:VR8, *+AR3[66]      \n\t"
			"|  VLDDWM2 *+AR1[43],  VR41:VR40  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR11:VR10, *+AR3[82]    \n\t"
			"|  VLDDWM2 *+AR1[48],  VR43:VR42  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR13:VR12, *+AR3[98]   \n\t"
			"|  VLDDWM2 *+AR1[49],  VR45:VR44  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR15:VR14, *+AR3[114]   \n\t"
			"|  VLDDWM2 *+AR1[50],  VR47:VR46  \n\t" //vr16-24,vr17-25
			"|	SSHFLL	  3, R10, R10       \n\t"  //SSHFLL 1  

			"   VSTDW   VR17:VR16, *+AR3[130]   \n\t"
			"|  VLDDWM2 *+AR1[55],  VR49:VR48  \n\t" //vr16-24,vr17-25
			"|  SADDA     R10,AR0,AR2		\n\t" //sadda AVAF 3 SVAF 2

			"   VSTDW   VR19:VR18, *+AR3[146]   \n\t"
			"|  VLDDWM2 *+AR1[56],  VR51:VR50  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR21:VR20, *+AR3[162]   \n\t"
			"|  VLDDWM2 *+AR1[57],  VR53:VR52  \n\t" //vr16-24,vr17-25
// index[6] read
			"	VSTDW   VR23:VR22, *+AR3[178]   \n\t"
			"|  VLDDWM2 *-AR2[57],  VR1:VR0    \n\t" //vr0-0,vr1-1

			"	VSTDW   VR25:VR24, *+AR3[194]   \n\t"
			"|  VLDDWM2 *-AR2[56],  VR3:VR2    \n\t" //vr2-3,vr3-4
			
			"	VSTDW   VR27:VR26, *+AR3[210]   \n\t"
			"|  VLDDWM2 *-AR2[55],  VR5:VR4    \n\t" //vr4-6,vr5-7

			"	VSTDW   VR29:VR28, *+AR3[226]   \n\t"
			"|  VLDDWM2 *-AR2[50],  VR7:VR6    \n\t" //vr6-9,vr7-10

			"	VSTDW   VR31:VR30, *+AR3[242]   \n\t"
			"|  VLDDWM2 *-AR2[49],  VR9:VR8    \n\t" //vr4-12,vr5-13

			"   VSTDW   VR33:VR32, *+AR3[258]         \n\t"
			"|  VLDDWM2 *-AR2[48],  VR11:VR10	\n\t" //vr10-15,vr11-16

			"   VSTDW   VR35:VR34, *+AR3[274]         \n\t"
			"|  VLDDWM2 *-AR2[43],  VR13:VR12  \n\t" //vr12-18,vr13-19

			"   VSTDW   VR37:VR36, *+AR3[290]         \n\t"
			"|  VLDDWM2 *-AR2[42],  VR15:VR14  \n\t" //vr14-21,vr15-22

			"   VSTDW   VR39:VR38, *+AR3[306]         \n\t"
			"|  VLDDWM2 *-AR2[41],  VR17:VR16  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR41:VR40, *+AR3[322]         \n\t"
			"|	VLDDWM2	*-AR2[8],	 VR19:VR18	\n\t" //v18-2   

			"   VSTDW   VR43:VR42, *+AR3[338]         \n\t"
			"|  VLDDWM2 *-AR2[7],  VR21:VR20  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR45:VR44, *+AR3[354]         \n\t"
			"|  VLDDWM2 *-AR2[6],  VR23:VR22  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR47:VR46, *+AR3[370]         \n\t"
			"|  VLDDWM2 *-AR2[1],  VR25:VR24  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR49:VR48, *+AR3[386]         \n\t"
			"|  VLDDWM2 *+AR2[0],  VR27:VR26  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR51:VR50, *+AR3[402]         \n\t"
			"|  VLDDWM2 *+AR2[1],   VR29:VR28  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR53:VR52, *+AR3[418]         \n\t"
			"|  VLDDWM2 *+AR2[6],  VR31:VR30  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR1:VR0, *+AR3[3]         \n\t"
			"|  VLDDWM2 *+AR2[7],  VR33:VR32  \n\t" //vr16-24,vr17-25
			"|	SLDW	  *+AR15[8], R10	\n\t" //sldw 7

			"   VSTDW   VR3:VR2, *+AR3[19]      \n\t"
			"|	VLDDWM2 *+AR2[8],  VR35:VR34  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR5:VR4, *+AR3[35]		\n\t"
			"|  VLDDWM2 *+AR2[41],  VR37:VR36  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR7:VR6, *+AR3[51]      \n\t"
			"|  VLDDWM2 *+AR2[42],  VR39:VR38  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR9:VR8, *+AR3[67]      \n\t"
			"|  VLDDWM2 *+AR2[43],  VR41:VR40  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR11:VR10, *+AR3[83]    \n\t"
			"|  VLDDWM2 *+AR2[48],  VR43:VR42  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR13:VR12, *+AR3[99]   \n\t"
			"|  VLDDWM2 *+AR2[49],  VR45:VR44  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR15:VR14, *+AR3[115]   \n\t"
			"|  VLDDWM2 *+AR2[50],  VR47:VR46  \n\t" //vr16-24,vr17-25
			"|	SSHFLL	  3, R10, R10       \n\t"  //SSHFLL 1  

			"   VSTDW   VR17:VR16, *+AR3[131]   \n\t"
			"|  VLDDWM2 *+AR2[55],  VR49:VR48  \n\t" //vr16-24,vr17-25
			"|  SADDA     R10,AR0,AR1		\n\t" //sadda AVAF 3 SVAF 2

			"   VSTDW   VR19:VR18, *+AR3[147]   \n\t"
			"|  VLDDWM2 *+AR2[56],  VR51:VR50  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR21:VR20, *+AR3[163]   \n\t"
			"|  VLDDWM2 *+AR2[57],  VR53:VR52  \n\t" //vr16-24,vr17-25
// index[8] read
			"	VSTDW   VR23:VR22, *+AR3[179]   \n\t"
			"|  VLDDWM2 *-AR1[57],  VR1:VR0    \n\t" //vr0-0,vr1-1

			"	VSTDW   VR25:VR24, *+AR3[195]   \n\t"
			"|  VLDDWM2 *-AR1[56],  VR3:VR2    \n\t" //vr2-3,vr3-4
			
			"	VSTDW   VR27:VR26, *+AR3[211]   \n\t"
			"|  VLDDWM2 *-AR1[55],  VR5:VR4    \n\t" //vr4-6,vr5-7

			"	VSTDW   VR29:VR28, *+AR3[227]   \n\t"
			"|  VLDDWM2 *-AR1[50],  VR7:VR6    \n\t" //vr6-9,vr7-10

			"	VSTDW   VR31:VR30, *+AR3[243]   \n\t"
			"|  VLDDWM2 *-AR1[49],  VR9:VR8    \n\t" //vr4-12,vr5-13

			"   VSTDW   VR33:VR32, *+AR3[259]         \n\t"
			"|  VLDDWM2 *-AR1[48],  VR11:VR10	\n\t" //vr10-15,vr11-16

			"   VSTDW   VR35:VR34, *+AR3[275]         \n\t"
			"|  VLDDWM2 *-AR1[43],  VR13:VR12  \n\t" //vr12-18,vr13-19

			"   VSTDW   VR37:VR36, *+AR3[291]         \n\t"
			"|  VLDDWM2 *-AR1[42],  VR15:VR14  \n\t" //vr14-21,vr15-22

			"   VSTDW   VR39:VR38, *+AR3[307]         \n\t"
			"|  VLDDWM2 *-AR1[41],  VR17:VR16  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR41:VR40, *+AR3[323]         \n\t"
			"|	VLDDWM2	*-AR1[8],	 VR19:VR18	\n\t" //v18-2   

			"   VSTDW   VR43:VR42, *+AR3[339]         \n\t"
			"|  VLDDWM2 *-AR1[7],  VR21:VR20  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR45:VR44, *+AR3[355]         \n\t"
			"|  VLDDWM2 *-AR1[6],  VR23:VR22  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR47:VR46, *+AR3[371]         \n\t"
			"|  VLDDWM2 *-AR1[1],  VR25:VR24  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR49:VR48, *+AR3[387]         \n\t"
			"|  VLDDWM2 *+AR1[0],  VR27:VR26  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR51:VR50, *+AR3[403]         \n\t"
			"|  VLDDWM2 *+AR1[1],   VR29:VR28  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR53:VR52, *+AR3[419]         \n\t"
			"|  VLDDWM2 *+AR1[6],  VR31:VR30  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR1:VR0, *+AR3[4]         \n\t"
			"|  VLDDWM2 *+AR1[7],  VR33:VR32  \n\t" //vr16-24,vr17-25
			"|	SLDW	  *+AR15[10], R10	\n\t" //sldw 7

			"   VSTDW   VR3:VR2, *+AR3[20]      \n\t"
			"|	VLDDWM2 *+AR1[8],  VR35:VR34  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR5:VR4, *+AR3[36]		\n\t"
			"|  VLDDWM2 *+AR1[41],  VR37:VR36  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR7:VR6, *+AR3[52]      \n\t"
			"|  VLDDWM2 *+AR1[42],  VR39:VR38  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR9:VR8, *+AR3[68]      \n\t"
			"|  VLDDWM2 *+AR1[43],  VR41:VR40  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR11:VR10, *+AR3[84]    \n\t"
			"|  VLDDWM2 *+AR1[48],  VR43:VR42  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR13:VR12, *+AR3[100]   \n\t"
			"|  VLDDWM2 *+AR1[49],  VR45:VR44  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR15:VR14, *+AR3[116]   \n\t"
			"|  VLDDWM2 *+AR1[50],  VR47:VR46  \n\t" //vr16-24,vr17-25
			"|	SSHFLL	  3, R10, R10       \n\t"  //SSHFLL 1  

			"   VSTDW   VR17:VR16, *+AR3[132]   \n\t"
			"|  VLDDWM2 *+AR1[55],  VR49:VR48  \n\t" //vr16-24,vr17-25
			"|  SADDA     R10,AR0,AR2		\n\t" //sadda AVAF 3 SVAF 2

			"   VSTDW   VR19:VR18, *+AR3[148]   \n\t"
			"|  VLDDWM2 *+AR1[56],  VR51:VR50  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR21:VR20, *+AR3[164]   \n\t"
			"|  VLDDWM2 *+AR1[57],  VR53:VR52  \n\t" //vr16-24,vr17-25
// index[10] read
			"	VSTDW   VR23:VR22, *+AR3[180]   \n\t"
			"|  VLDDWM2 *-AR2[57],  VR1:VR0    \n\t" //vr0-0,vr1-1

			"	VSTDW   VR25:VR24, *+AR3[196]   \n\t"
			"|  VLDDWM2 *-AR2[56],  VR3:VR2    \n\t" //vr2-3,vr3-4
			
			"	VSTDW   VR27:VR26, *+AR3[212]   \n\t"
			"|  VLDDWM2 *-AR2[55],  VR5:VR4    \n\t" //vr4-6,vr5-7

			"	VSTDW   VR29:VR28, *+AR3[228]   \n\t"
			"|  VLDDWM2 *-AR2[50],  VR7:VR6    \n\t" //vr6-9,vr7-10

			"	VSTDW   VR31:VR30, *+AR3[244]   \n\t"
			"|  VLDDWM2 *-AR2[49],  VR9:VR8    \n\t" //vr4-12,vr5-13

			"   VSTDW   VR33:VR32, *+AR3[260]         \n\t"
			"|  VLDDWM2 *-AR2[48],  VR11:VR10	\n\t" //vr10-15,vr11-16

			"   VSTDW   VR35:VR34, *+AR3[276]         \n\t"
			"|  VLDDWM2 *-AR2[43],  VR13:VR12  \n\t" //vr12-18,vr13-19

			"   VSTDW   VR37:VR36, *+AR3[292]         \n\t"
			"|  VLDDWM2 *-AR2[42],  VR15:VR14  \n\t" //vr14-21,vr15-22

			"   VSTDW   VR39:VR38, *+AR3[308]         \n\t"
			"|  VLDDWM2 *-AR2[41],  VR17:VR16  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR41:VR40, *+AR3[324]         \n\t"
			"|	VLDDWM2	*-AR2[8],	 VR19:VR18	\n\t" //v18-2   

			"   VSTDW   VR43:VR42, *+AR3[340]         \n\t"
			"|  VLDDWM2 *-AR2[7],  VR21:VR20  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR45:VR44, *+AR3[356]         \n\t"
			"|  VLDDWM2 *-AR2[6],  VR23:VR22  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR47:VR46, *+AR3[372]         \n\t"
			"|  VLDDWM2 *-AR2[1],  VR25:VR24  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR49:VR48, *+AR3[388]         \n\t"
			"|  VLDDWM2 *+AR2[0],  VR27:VR26  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR51:VR50, *+AR3[404]         \n\t"
			"|  VLDDWM2 *+AR2[1],   VR29:VR28  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR53:VR52, *+AR3[420]         \n\t"
			"|  VLDDWM2 *+AR2[6],  VR31:VR30  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR1:VR0, *+AR3[5]         \n\t"
			"|  VLDDWM2 *+AR2[7],  VR33:VR32  \n\t" //vr16-24,vr17-25
			"|	SLDW	  *+AR15[12], R10	\n\t" //sldw 7

			"   VSTDW   VR3:VR2, *+AR3[21]      \n\t"
			"|	VLDDWM2 *+AR2[8],  VR35:VR34  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR5:VR4, *+AR3[37]		\n\t"
			"|  VLDDWM2 *+AR2[41],  VR37:VR36  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR7:VR6, *+AR3[53]      \n\t"
			"|  VLDDWM2 *+AR2[42],  VR39:VR38  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR9:VR8, *+AR3[69]      \n\t"
			"|  VLDDWM2 *+AR2[43],  VR41:VR40  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR11:VR10, *+AR3[85]    \n\t"
			"|  VLDDWM2 *+AR2[48],  VR43:VR42  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR13:VR12, *+AR3[101]   \n\t"
			"|  VLDDWM2 *+AR2[49],  VR45:VR44  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR15:VR14, *+AR3[117]   \n\t"
			"|  VLDDWM2 *+AR2[50],  VR47:VR46  \n\t" //vr16-24,vr17-25
			"|	SSHFLL	  3, R10, R10       \n\t"  //SSHFLL 1  

			"   VSTDW   VR17:VR16, *+AR3[133]   \n\t"
			"|  VLDDWM2 *+AR2[55],  VR49:VR48  \n\t" //vr16-24,vr17-25
			"|  SADDA     R10,AR0,AR1		\n\t" //sadda AVAF 3 SVAF 2

			"   VSTDW   VR19:VR18, *+AR3[149]   \n\t"
			"|  VLDDWM2 *+AR2[56],  VR51:VR50  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR21:VR20, *+AR3[165]   \n\t"
			"|  VLDDWM2 *+AR2[57],  VR53:VR52  \n\t" //vr16-24,vr17-25
// index[12] read
			"	VSTDW   VR23:VR22, *+AR3[181]   \n\t"
			"|  VLDDWM2 *-AR1[57],  VR1:VR0    \n\t" //vr0-0,vr1-1

			"	VSTDW   VR25:VR24, *+AR3[197]   \n\t"
			"|  VLDDWM2 *-AR1[56],  VR3:VR2    \n\t" //vr2-3,vr3-4
			
			"	VSTDW   VR27:VR26, *+AR3[213]   \n\t"
			"|  VLDDWM2 *-AR1[55],  VR5:VR4    \n\t" //vr4-6,vr5-7

			"	VSTDW   VR29:VR28, *+AR3[229]   \n\t"
			"|  VLDDWM2 *-AR1[50],  VR7:VR6    \n\t" //vr6-9,vr7-10

			"	VSTDW   VR31:VR30, *+AR3[245]   \n\t"
			"|  VLDDWM2 *-AR1[49],  VR9:VR8    \n\t" //vr4-12,vr5-13

			"   VSTDW   VR33:VR32, *+AR3[261]         \n\t"
			"|  VLDDWM2 *-AR1[48],  VR11:VR10	\n\t" //vr10-15,vr11-16

			"   VSTDW   VR35:VR34, *+AR3[277]         \n\t"
			"|  VLDDWM2 *-AR1[43],  VR13:VR12  \n\t" //vr12-18,vr13-19

			"   VSTDW   VR37:VR36, *+AR3[293]         \n\t"
			"|  VLDDWM2 *-AR1[42],  VR15:VR14  \n\t" //vr14-21,vr15-22

			"   VSTDW   VR39:VR38, *+AR3[309]         \n\t"
			"|  VLDDWM2 *-AR1[41],  VR17:VR16  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR41:VR40, *+AR3[325]         \n\t"
			"|	VLDDWM2	*-AR1[8],	 VR19:VR18	\n\t" //v18-2   

			"   VSTDW   VR43:VR42, *+AR3[341]         \n\t"
			"|  VLDDWM2 *-AR1[7],  VR21:VR20  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR45:VR44, *+AR3[357]         \n\t"
			"|  VLDDWM2 *-AR1[6],  VR23:VR22  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR47:VR46, *+AR3[373]         \n\t"
			"|  VLDDWM2 *-AR1[1],  VR25:VR24  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR49:VR48, *+AR3[389]         \n\t"
			"|  VLDDWM2 *+AR1[0],  VR27:VR26  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR51:VR50, *+AR3[405]         \n\t"
			"|  VLDDWM2 *+AR1[1],   VR29:VR28  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR53:VR52, *+AR3[421]         \n\t"
			"|  VLDDWM2 *+AR1[6],  VR31:VR30  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR1:VR0, *+AR3[6]         \n\t"
			"|  VLDDWM2 *+AR1[7],  VR33:VR32  \n\t" //vr16-24,vr17-25
			"|	SLDW	  *+AR15[14], R10	\n\t" //sldw 7

			"   VSTDW   VR3:VR2, *+AR3[22]      \n\t"
			"|	VLDDWM2 *+AR1[8],  VR35:VR34  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR5:VR4, *+AR3[38]		\n\t"
			"|  VLDDWM2 *+AR1[41],  VR37:VR36  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR7:VR6, *+AR3[54]      \n\t"
			"|  VLDDWM2 *+AR1[42],  VR39:VR38  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR9:VR8, *+AR3[70]      \n\t"
			"|  VLDDWM2 *+AR1[43],  VR41:VR40  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR11:VR10, *+AR3[86]    \n\t"
			"|  VLDDWM2 *+AR1[48],  VR43:VR42  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR13:VR12, *+AR3[102]   \n\t"
			"|  VLDDWM2 *+AR1[49],  VR45:VR44  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR15:VR14, *+AR3[118]   \n\t"
			"|  VLDDWM2 *+AR1[50],  VR47:VR46  \n\t" //vr16-24,vr17-25
			"|	SSHFLL	  3, R10, R10       \n\t"  //SSHFLL 1  

			"   VSTDW   VR17:VR16, *+AR3[134]   \n\t"
			"|  VLDDWM2 *+AR1[55],  VR49:VR48  \n\t" //vr16-24,vr17-25
			"|  SADDA     R10,AR0,AR2		\n\t" //sadda AVAF 3 SVAF 2

			"   VSTDW   VR19:VR18, *+AR3[150]   \n\t"
			"|  VLDDWM2 *+AR1[56],  VR51:VR50  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR21:VR20, *+AR3[166]   \n\t"
			"|  VLDDWM2 *+AR1[57],  VR53:VR52  \n\t" //vr16-24,vr17-25
// index[14] read
			"	VSTDW   VR23:VR22, *+AR3[182]   \n\t"
			"|  VLDDWM2 *-AR2[57],  VR1:VR0    \n\t" //vr0-0,vr1-1

			"	VSTDW   VR25:VR24, *+AR3[198]   \n\t"
			"|  VLDDWM2 *-AR2[56],  VR3:VR2    \n\t" //vr2-3,vr3-4
			
			"	VSTDW   VR27:VR26, *+AR3[214]   \n\t"
			"|  VLDDWM2 *-AR2[55],  VR5:VR4    \n\t" //vr4-6,vr5-7

			"	VSTDW   VR29:VR28, *+AR3[230]   \n\t"
			"|  VLDDWM2 *-AR2[50],  VR7:VR6    \n\t" //vr6-9,vr7-10

			"	VSTDW   VR31:VR30, *+AR3[246]   \n\t"
			"|  VLDDWM2 *-AR2[49],  VR9:VR8    \n\t" //vr4-12,vr5-13

			"   VSTDW   VR33:VR32, *+AR3[262]         \n\t"
			"|  VLDDWM2 *-AR2[48],  VR11:VR10	\n\t" //vr10-15,vr11-16

			"   VSTDW   VR35:VR34, *+AR3[278]         \n\t"
			"|  VLDDWM2 *-AR2[43],  VR13:VR12  \n\t" //vr12-18,vr13-19

			"   VSTDW   VR37:VR36, *+AR3[294]         \n\t"
			"|  VLDDWM2 *-AR2[42],  VR15:VR14  \n\t" //vr14-21,vr15-22

			"   VSTDW   VR39:VR38, *+AR3[310]         \n\t"
			"|  VLDDWM2 *-AR2[41],  VR17:VR16  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR41:VR40, *+AR3[326]         \n\t"
			"|	VLDDWM2	*-AR2[8],	 VR19:VR18	\n\t" //v18-2   

			"   VSTDW   VR43:VR42, *+AR3[342]         \n\t"
			"|  VLDDWM2 *-AR2[7],  VR21:VR20  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR45:VR44, *+AR3[358]         \n\t"
			"|  VLDDWM2 *-AR2[6],  VR23:VR22  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR47:VR46, *+AR3[374]         \n\t"
			"|  VLDDWM2 *-AR2[1],  VR25:VR24  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR49:VR48, *+AR3[390]         \n\t"
			"|  VLDDWM2 *+AR2[0],  VR27:VR26  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR51:VR50, *+AR3[406]         \n\t"
			"|  VLDDWM2 *+AR2[1],   VR29:VR28  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR53:VR52, *+AR3[422]         \n\t"
			"|  VLDDWM2 *+AR2[6],  VR31:VR30  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR1:VR0, *+AR3[7]         \n\t"
			"|  VLDDWM2 *+AR2[7],  VR33:VR32  \n\t" //vr16-24,vr17-25
			"|	SLDW	  *+AR15[1], R10	\n\t" //sldw 7

			"   VSTDW   VR3:VR2, *+AR3[23]      \n\t"
			"|	VLDDWM2 *+AR2[8],  VR35:VR34  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR5:VR4, *+AR3[39]		\n\t"
			"|  VLDDWM2 *+AR2[41],  VR37:VR36  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR7:VR6, *+AR3[55]      \n\t"
			"|  VLDDWM2 *+AR2[42],  VR39:VR38  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR9:VR8, *+AR3[71]      \n\t"
			"|  VLDDWM2 *+AR2[43],  VR41:VR40  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR11:VR10, *+AR3[87]    \n\t"
			"|  VLDDWM2 *+AR2[48],  VR43:VR42  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR13:VR12, *+AR3[103]   \n\t"
			"|  VLDDWM2 *+AR2[49],  VR45:VR44  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR15:VR14, *+AR3[119]   \n\t"
			"|  VLDDWM2 *+AR2[50],  VR47:VR46  \n\t" //vr16-24,vr17-25
			"|	SSHFLL	  3, R10, R10       \n\t"  //SSHFLL 1  

			"   VSTDW   VR17:VR16, *+AR3[135]   \n\t"
			"|  VLDDWM2 *+AR2[55],  VR49:VR48  \n\t" //vr16-24,vr17-25
			"|  SADDA     R10,AR0,AR1		\n\t" //sadda AVAF 3 SVAF 2

			"   VSTDW   VR19:VR18, *+AR3[151]   \n\t"
			"|  VLDDWM2 *+AR2[56],  VR51:VR50  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR21:VR20, *+AR3[167]   \n\t"
			"|  VLDDWM2 *+AR2[57],  VR53:VR52  \n\t" //vr16-24,vr17-25

			"	SNOP	1\n\t"

			"	VSTDW   VR23:VR22, *+AR3[183]   \n\t"
			"|	VSTDW   VR25:VR24, *+AR3[199]   \n\t"

			"	VSTDW   VR27:VR26, *+AR3[215]   \n\t"
			"|	VSTDW   VR29:VR28, *+AR3[231]   \n\t"

			"	VSTDW   VR31:VR30, *+AR3[247]   \n\t"
			"|  VSTDW   VR33:VR32, *+AR3[263]         \n\t"

			"   VSTDW   VR35:VR34, *+AR3[279]         \n\t"
			"|  VSTDW   VR37:VR36, *+AR3[295]         \n\t"

			"   VSTDW   VR39:VR38, *+AR3[311]         \n\t"
			"|  VSTDW   VR41:VR40, *+AR3[327]         \n\t"

			"   VSTDW   VR43:VR42, *+AR3[343]         \n\t"
			"|  VSTDW   VR45:VR44, *+AR3[359]         \n\t"

			"   VSTDW   VR47:VR46, *+AR3[375]         \n\t"
			"|  VSTDW   VR49:VR48, *+AR3[391]         \n\t"

			"   VSTDW   VR51:VR50, *+AR3[407]         \n\t"
			"|  VSTDW   VR53:VR52, *+AR3[423]         \n\t"

//			"|  SMVAGA.M2 R31,AR3           \n\t" //R31 //index[1] read

			"   VLDDWM2 *-AR1[57],  VR1:VR0    \n\t" //vr0-0,vr1-1
			"|  VLDDWM2 *-AR1[56],  VR3:VR2    \n\t" //vr2-3,vr3-4

			"   VLDDWM2 *-AR1[55],  VR5:VR4    \n\t" //vr4-6,vr5-7
			"|  VLDDWM2 *-AR1[50],  VR7:VR6    \n\t" //vr6-9,vr7-10

			"   VLDDWM2 *-AR1[49],  VR9:VR8    \n\t" //vr4-12,vr5-13
			"|  VLDDWM2 *-AR1[48],  VR11:VR10	\n\t" //vr10-15,vr11-16

			"   VLDDWM2 *-AR1[43],  VR13:VR12  \n\t" //vr12-18,vr13-19
			"|  VLDDWM2 *-AR1[42],  VR15:VR14  \n\t" //vr14-21,vr15-22

			"   VLDDWM2 *-AR1[41],  VR17:VR16  \n\t" //vr16-24,vr17-25
			"|	VLDDWM2	*-AR1[8],	 VR19:VR18	\n\t" //v18-2   

			"   VLDDWM2 *-AR1[7],  VR21:VR20  \n\t" //vr16-24,vr17-25
			"|  VLDDWM2 *-AR1[6],  VR23:VR22  \n\t" //vr16-24,vr17-25

			"   VLDDWM2 *-AR1[1],  VR25:VR24  \n\t" //vr16-24,vr17-25
			"|  VLDDWM2 *+AR1[0],  VR27:VR26  \n\t" //vr16-24,vr17-25

			"   VLDDWM2 *+AR1[1],   VR29:VR28  \n\t" //vr16-24,vr17-25
			"|  VLDDWM2 *+AR1[6],  VR31:VR30  \n\t" //vr16-24,vr17-25

			"	SNOP	1\n\t"

			"   VSTW   VR0, *++AR3[1]         \n\t"  
			"|  VLDDWM2 *+AR1[7],  VR33:VR32  \n\t" //vr16-24,vr17-25
			"|	SLDW	  *+AR15[3], R10	\n\t" //sldw 7

			"   VSTW   VR1, *++AR3[16]      \n\t"
			"|  VLDDWM2 *+AR1[8],  VR35:VR34  \n\t" //vr16-24,vr17-25

			"   VSTW   VR2, *++AR3[16]		\n\t"
			"|  VLDDWM2 *+AR1[41],  VR37:VR36  \n\t" //vr16-24,vr17-25

			"   VSTW   VR3, *++AR3[16]      \n\t"
			"|  VLDDWM2 *+AR1[42],  VR39:VR38  \n\t" //vr16-24,vr17-25

			"   VSTW   VR4, *++AR3[16]      \n\t"
			"|  VLDDWM2 *+AR1[43],  VR41:VR40  \n\t" //vr16-24,vr17-25

			"   VSTW   VR5, *++AR3[16]    \n\t"
			"|  VLDDWM2 *+AR1[48],  VR43:VR42  \n\t" //vr16-24,vr17-25

			"   VSTW   VR6, *++AR3[16]   \n\t"
			"|  VLDDWM2 *+AR1[49],  VR45:VR44  \n\t" //vr16-24,vr17-25

			"   VSTW   VR7, *++AR3[16]   \n\t"
			"|  VLDDWM2 *+AR1[50],  VR47:VR46  \n\t" //vr16-24,vr17-25
			"|	SSHFLL	  3, R10, R10       \n\t"  //SSHFLL 1  

			"   VSTW   VR8, *++AR3[16]   \n\t"
			"|  VLDDWM2 *+AR1[55],  VR49:VR48  \n\t" //vr16-24,vr17-25
			"|  SADDA     R10,AR0,AR2		\n\t" //sadda AVAF 3 SVAF 2

			"   VSTW   VR9, *++AR3[16]   \n\t"
			"|  VLDDWM2 *+AR1[56],  VR51:VR50  \n\t" //vr16-24,vr17-25

			"   VSTW   VR10, *++AR3[16]   \n\t"
			"|  VLDDWM2 *+AR1[57],  VR53:VR52  \n\t" //vr16-24,vr17-25

			"   VSTW   VR11, *+AR3[16]   \n\t"
			"|  VSTW   VR12, *+AR3[32]   \n\t"

			"   VSTW   VR13, *+AR3[48]   \n\t"
			"|  VSTW   VR14, *+AR3[64]   \n\t"

			"   VSTW   VR15, *+AR3[80]   \n\t"
			"|  VSTW   VR16, *+AR3[96]   \n\t"

			"   VSTW   VR17, *+AR3[112]   \n\t"
			"|  VSTW   VR18, *+AR3[128]   \n\t"

			"   VSTW   VR19, *+AR3[144]   \n\t"
			"|  VSTW   VR20, *+AR3[160]   \n\t"

			"   VSTW   VR21, *+AR3[176]   \n\t"
			"|  VSTW   VR22, *+AR3[192]   \n\t"

			"   VSTW   VR23, *+AR3[208]   \n\t"
			"|  VSTW   VR24, *+AR3[224]   \n\t"

			"   VSTW   VR25, *+AR3[240]   \n\t"
			"|  VSTW   VR26, *+AR3[256]   \n\t"

			"   VSTW   VR27, *+AR3[272]   \n\t"
			"|  VSTW   VR28, *+AR3[288]   \n\t"
			"|  SADDA     R10,AR0,AR1		\n\t" //sadda AVAF 3 SVAF 2

			"   VSTW   VR29, *+AR3[304]   \n\t"
			"|  VSTW   VR30, *+AR3[320]   \n\t"

			"   VSTW   VR31, *+AR3[336]   \n\t"
			"|  VSTW   VR32, *+AR3[352]   \n\t"

			"   VSTW   VR33, *+AR3[368]   \n\t"
			"|  VSTW   VR34, *+AR3[384]   \n\t"

			"   VSTW   VR35, *+AR3[400]   \n\t"
			"|  VSTW   VR36, *+AR3[416]   \n\t"

			"   VSTW   VR37, *++AR3[432]   \n\t"
			"|  SMVAGA.M2 R31,AR4           \n\t" //R31 
// index[3] read
			"	VSTW   VR38, *++AR3[16]   \n\t"
			"|  VLDDWM2 *-AR1[57],  VR1:VR0    \n\t" //vr0-0,vr1-1

			"	VSTW   VR39, *++AR3[16]   \n\t"
			"|  VLDDWM2 *-AR1[56],  VR3:VR2    \n\t" //vr2-3,vr3-4
			
			"	VSTW   VR40, *++AR3[16]   \n\t"
			"|  VLDDWM2 *-AR1[55],  VR5:VR4    \n\t" //vr4-6,vr5-7

			"	VSTW   VR41, *++AR3[16]   \n\t"
			"|  VLDDWM2 *-AR1[50],  VR7:VR6    \n\t" //vr6-9,vr7-10

			"	VSTW   VR42, *++AR3[16]   \n\t"
			"|  VLDDWM2 *-AR1[49],  VR9:VR8    \n\t" //vr4-12,vr5-13

			"   VSTW   VR43, *++AR3[16]         \n\t"
			"|  VLDDWM2 *-AR1[48],  VR11:VR10	\n\t" //vr10-15,vr11-16

			"   VSTW   VR44, *++AR3[16]         \n\t"
			"|  VLDDWM2 *-AR1[43],  VR13:VR12  \n\t" //vr12-18,vr13-19

			"   VSTW   VR45, *++AR3[16]         \n\t"
			"|  VLDDWM2 *-AR1[42],  VR15:VR14  \n\t" //vr14-21,vr15-22

			"   VSTW   VR46, *++AR3[16]         \n\t"
			"|  VLDDWM2 *-AR1[41],  VR17:VR16  \n\t" //vr16-24,vr17-25

			"   VSTW   VR47, *++AR3[16]         \n\t"
			"|	VLDDWM2	*-AR1[8],	 VR19:VR18	\n\t" //v18-2   

			"   VSTW   VR48, *++AR3[16]         \n\t"
			"|  VLDDWM2 *-AR1[7],  VR21:VR20  \n\t" //vr16-24,vr17-25

			"   VSTW   VR49, *++AR3[16]         \n\t"
			"|  VLDDWM2 *-AR1[6],  VR23:VR22  \n\t" //vr16-24,vr17-25

			"   VSTW   VR50, *++AR3[16]         \n\t"
			"|  VLDDWM2 *-AR1[1],  VR25:VR24  \n\t" //vr16-24,vr17-25

			"   VSTW   VR51, *++AR3[16]         \n\t"
			"|  VLDDWM2 *+AR1[0],  VR27:VR26  \n\t" //vr16-24,vr17-25

			"   VSTW   VR52, *++AR3[16]         \n\t"
			"|  VLDDWM2 *+AR1[1],   VR29:VR28  \n\t" //vr16-24,vr17-25

			"   VSTW   VR53, *++AR3[16]         \n\t"
			"|  VLDDWM2 *+AR1[6],  VR31:VR30  \n\t" //vr16-24,vr17-25

			"   VSTW   VR0, *++AR4[3]         \n\t"
			"|  VLDDWM2 *+AR1[7],  VR33:VR32  \n\t" //vr16-24,vr17-25
			"|	SLDW	  *+AR15[5], R10	\n\t" //sldw 7

			"   VSTW   VR1, *++AR4[16]      \n\t"
			"|	VLDDWM2 *+AR1[8],  VR35:VR34  \n\t" //vr16-24,vr17-25

			"   VSTW   VR2, *++AR4[16]		\n\t"
			"|  VLDDWM2 *+AR1[41],  VR37:VR36  \n\t" //vr16-24,vr17-25

			"   VSTW   VR3, *++AR4[16]      \n\t"
			"|  VLDDWM2 *+AR1[42],  VR39:VR38  \n\t" //vr16-24,vr17-25

			"   VSTW   VR4, *++AR4[16]      \n\t"
			"|  VLDDWM2 *+AR1[43],  VR41:VR40  \n\t" //vr16-24,vr17-25

			"   VSTW   VR5, *++AR4[16]    \n\t"
			"|  VLDDWM2 *+AR1[48],  VR43:VR42  \n\t" //vr16-24,vr17-25

			"   VSTW   VR6, *++AR4[16]   \n\t"
			"|  VLDDWM2 *+AR1[49],  VR45:VR44  \n\t" //vr16-24,vr17-25

			"   VSTW   VR7, *++AR4[16]   \n\t"
			"|  VLDDWM2 *+AR1[50],  VR47:VR46  \n\t" //vr16-24,vr17-25
			"|	SSHFLL	  3, R10, R10       \n\t"  //SSHFLL 1  

			"   VSTW   VR8, *++AR4[16]   \n\t"
			"|  VLDDWM2 *+AR1[55],  VR49:VR48  \n\t" //vr16-24,vr17-25

			"   VSTW   VR9, *++AR4[16]   \n\t"
			"|  VLDDWM2 *+AR1[56],  VR51:VR50  \n\t" //vr16-24,vr17-25

			"   VSTW   VR10, *++AR4[16]   \n\t"
			"|  VLDDWM2 *+AR1[57],  VR53:VR52  \n\t" //vr16-24,vr17-25

			"   VSTW   VR11, *+AR4[16]   \n\t"
			"|  VSTW   VR12, *+AR4[32]   \n\t"

			"   VSTW   VR13, *+AR4[48]   \n\t"
			"|  VSTW   VR14, *+AR4[64]   \n\t"

			"   VSTW   VR15, *+AR4[80]   \n\t"
			"|  VSTW   VR16, *+AR4[96]   \n\t"

			"   VSTW   VR17, *+AR4[112]   \n\t"
			"|  VSTW   VR18, *+AR4[128]   \n\t"

			"   VSTW   VR19, *+AR4[144]   \n\t"
			"|  VSTW   VR20, *+AR4[160]   \n\t"

			"   VSTW   VR21, *+AR4[176]   \n\t"
			"|  VSTW   VR22, *+AR4[192]   \n\t"

			"   VSTW   VR23, *+AR4[208]   \n\t"
			"|  VSTW   VR24, *+AR4[224]   \n\t"

			"   VSTW   VR25, *+AR4[240]   \n\t"
			"|  VSTW   VR26, *+AR4[256]   \n\t"

			"   VSTW   VR27, *+AR4[272]   \n\t"
			"|  VSTW   VR28, *+AR4[288]   \n\t"
			"|  SADDA     R10,AR0,AR1		\n\t" //sadda AVAF 3 SVAF 2

			"   VSTW   VR29, *+AR4[304]   \n\t"
			"|  VSTW   VR30, *+AR4[320]   \n\t"

			"   VSTW   VR31, *+AR4[336]   \n\t"
			"|  VSTW   VR32, *+AR4[352]   \n\t"

			"   VSTW   VR33, *+AR4[368]   \n\t"
			"|  VSTW   VR34, *+AR4[384]   \n\t"

			"   VSTW   VR35, *+AR4[400]   \n\t"
			"|  VSTW   VR36, *+AR4[416]   \n\t"

			"   VSTW   VR37, *++AR4[432]   \n\t"
			"|  SMVAGA.M2 R31,AR3           \n\t" //R31 // index[5] read
			"	VSTW   VR38, *++AR4[16]   \n\t"
			"|  VLDDWM2 *-AR1[57],  VR1:VR0    \n\t" //vr0-0,vr1-1

			"	VSTW   VR39, *++AR4[16]   \n\t"
			"|  VLDDWM2 *-AR1[56],  VR3:VR2    \n\t" //vr2-3,vr3-4
			
			"	VSTW   VR40, *++AR4[16]   \n\t"
			"|  VLDDWM2 *-AR1[55],  VR5:VR4    \n\t" //vr4-6,vr5-7

			"	VSTW   VR41, *++AR4[16]   \n\t"
			"|  VLDDWM2 *-AR1[50],  VR7:VR6    \n\t" //vr6-9,vr7-10

			"	VSTW   VR42, *++AR4[16]   \n\t"
			"|  VLDDWM2 *-AR1[49],  VR9:VR8    \n\t" //vr4-12,vr5-13

			"   VSTW   VR43, *++AR4[16]         \n\t"
			"|  VLDDWM2 *-AR1[48],  VR11:VR10	\n\t" //vr10-15,vr11-16

			"   VSTW   VR44, *++AR4[16]         \n\t"
			"|  VLDDWM2 *-AR1[43],  VR13:VR12  \n\t" //vr12-18,vr13-19

			"   VSTW   VR45, *++AR4[16]         \n\t"
			"|  VLDDWM2 *-AR1[42],  VR15:VR14  \n\t" //vr14-21,vr15-22

			"   VSTW   VR46, *++AR4[16]         \n\t"
			"|  VLDDWM2 *-AR1[41],  VR17:VR16  \n\t" //vr16-24,vr17-25

			"   VSTW   VR47, *++AR4[16]         \n\t"
			"|	VLDDWM2	*-AR1[8],	 VR19:VR18	\n\t" //v18-2   

			"   VSTW   VR48, *++AR4[16]         \n\t"
			"|  VLDDWM2 *-AR1[7],  VR21:VR20  \n\t" //vr16-24,vr17-25

			"   VSTW   VR49, *++AR4[16]         \n\t"
			"|  VLDDWM2 *-AR1[6],  VR23:VR22  \n\t" //vr16-24,vr17-25

			"   VSTW   VR50, *++AR4[16]         \n\t"
			"|  VLDDWM2 *-AR1[1],  VR25:VR24  \n\t" //vr16-24,vr17-25

			"   VSTW   VR51, *++AR4[16]         \n\t"
			"|  VLDDWM2 *+AR1[0],  VR27:VR26  \n\t" //vr16-24,vr17-25

			"   VSTW   VR52, *++AR4[16]         \n\t"
			"|  VLDDWM2 *+AR1[1],   VR29:VR28  \n\t" //vr16-24,vr17-25

			"   VSTW   VR53, *++AR4[16]         \n\t"
			"|  VLDDWM2 *+AR1[6],  VR31:VR30  \n\t" //vr16-24,vr17-25

			"   VSTW   VR0, *++AR3[5]         \n\t"
			"|  VLDDWM2 *+AR1[7],  VR33:VR32  \n\t" //vr16-24,vr17-25
			"|	SLDW	  *+AR15[7], R10	\n\t" //sldw 7

			"   VSTW   VR1, *++AR3[16]      \n\t"
			"|	VLDDWM2 *+AR1[8],  VR35:VR34  \n\t" //vr16-24,vr17-25

			"   VSTW   VR2, *++AR3[16]		\n\t"
			"|  VLDDWM2 *+AR1[41],  VR37:VR36  \n\t" //vr16-24,vr17-25

			"   VSTW   VR3, *++AR3[16]      \n\t"
			"|  VLDDWM2 *+AR1[42],  VR39:VR38  \n\t" //vr16-24,vr17-25

			"   VSTW   VR4, *++AR3[16]      \n\t"
			"|  VLDDWM2 *+AR1[43],  VR41:VR40  \n\t" //vr16-24,vr17-25

			"   VSTW   VR5, *++AR3[16]    \n\t"
			"|  VLDDWM2 *+AR1[48],  VR43:VR42  \n\t" //vr16-24,vr17-25

			"   VSTW   VR6, *++AR3[16]   \n\t"
			"|  VLDDWM2 *+AR1[49],  VR45:VR44  \n\t" //vr16-24,vr17-25

			"   VSTW   VR7, *++AR3[16]   \n\t"
			"|  VLDDWM2 *+AR1[50],  VR47:VR46  \n\t" //vr16-24,vr17-25
			"|	SSHFLL	  3, R10, R10       \n\t"  //SSHFLL 1  

			"   VSTW   VR8, *++AR3[16]   \n\t"
			"|  VLDDWM2 *+AR1[55],  VR49:VR48  \n\t" //vr16-24,vr17-25

			"   VSTW   VR9, *++AR3[16]   \n\t"
			"|  VLDDWM2 *+AR1[56],  VR51:VR50  \n\t" //vr16-24,vr17-25

			"   VSTW   VR10, *++AR3[16]   \n\t"
			"|  VLDDWM2 *+AR1[57],  VR53:VR52  \n\t" //vr16-24,vr17-25

			"   VSTW   VR11, *+AR3[16]   \n\t"
			"|  VSTW   VR12, *+AR3[32]   \n\t"

			"   VSTW   VR13, *+AR3[48]   \n\t"
			"|  VSTW   VR14, *+AR3[64]   \n\t"

			"   VSTW   VR15, *+AR3[80]   \n\t"
			"|  VSTW   VR16, *+AR3[96]   \n\t"

			"   VSTW   VR17, *+AR3[112]   \n\t"
			"|  VSTW   VR18, *+AR3[128]   \n\t"

			"   VSTW   VR19, *+AR3[144]   \n\t"
			"|  VSTW   VR20, *+AR3[160]   \n\t"

			"   VSTW   VR21, *+AR3[176]   \n\t"
			"|  VSTW   VR22, *+AR3[192]   \n\t"

			"   VSTW   VR23, *+AR3[208]   \n\t"
			"|  VSTW   VR24, *+AR3[224]   \n\t"

			"   VSTW   VR25, *+AR3[240]   \n\t"
			"|  VSTW   VR26, *+AR3[256]   \n\t"

			"   VSTW   VR27, *+AR3[272]   \n\t"
			"|  VSTW   VR28, *+AR3[288]   \n\t"
			"|  SADDA     R10,AR0,AR1		\n\t" //sadda AVAF 3 SVAF 2

			"   VSTW   VR29, *+AR3[304]   \n\t"
			"|  VSTW   VR30, *+AR3[320]   \n\t"

			"   VSTW   VR31, *+AR3[336]   \n\t"
			"|  VSTW   VR32, *+AR3[352]   \n\t"

			"   VSTW   VR33, *+AR3[368]   \n\t"
			"|  VSTW   VR34, *+AR3[384]   \n\t"

			"   VSTW   VR35, *+AR3[400]   \n\t"
			"|  VSTW   VR36, *+AR3[416]   \n\t"

			"   VSTW   VR37, *++AR3[432]   \n\t"
			"|  SMVAGA.M2 R31,AR4           \n\t" //R31 
// index[7] read
			"	VSTW   VR38, *++AR3[16]   \n\t"
			"|  VLDDWM2 *-AR1[57],  VR1:VR0    \n\t" //vr0-0,vr1-1

			"	VSTW   VR39, *++AR3[16]   \n\t"
			"|  VLDDWM2 *-AR1[56],  VR3:VR2    \n\t" //vr2-3,vr3-4
			
			"	VSTW   VR40, *++AR3[16]   \n\t"
			"|  VLDDWM2 *-AR1[55],  VR5:VR4    \n\t" //vr4-6,vr5-7

			"	VSTW   VR41, *++AR3[16]   \n\t"
			"|  VLDDWM2 *-AR1[50],  VR7:VR6    \n\t" //vr6-9,vr7-10

			"	VSTW   VR42, *++AR3[16]   \n\t"
			"|  VLDDWM2 *-AR1[49],  VR9:VR8    \n\t" //vr4-12,vr5-13

			"   VSTW   VR43, *++AR3[16]         \n\t"
			"|  VLDDWM2 *-AR1[48],  VR11:VR10	\n\t" //vr10-15,vr11-16

			"   VSTW   VR44, *++AR3[16]         \n\t"
			"|  VLDDWM2 *-AR1[43],  VR13:VR12  \n\t" //vr12-18,vr13-19

			"   VSTW   VR45, *++AR3[16]         \n\t"
			"|  VLDDWM2 *-AR1[42],  VR15:VR14  \n\t" //vr14-21,vr15-22

			"   VSTW   VR46, *++AR3[16]         \n\t"
			"|  VLDDWM2 *-AR1[41],  VR17:VR16  \n\t" //vr16-24,vr17-25

			"   VSTW   VR47, *++AR3[16]         \n\t"
			"|	VLDDWM2	*-AR1[8],	 VR19:VR18	\n\t" //v18-2   

			"   VSTW   VR48, *++AR3[16]         \n\t"
			"|  VLDDWM2 *-AR1[7],  VR21:VR20  \n\t" //vr16-24,vr17-25

			"   VSTW   VR49, *++AR3[16]         \n\t"
			"|  VLDDWM2 *-AR1[6],  VR23:VR22  \n\t" //vr16-24,vr17-25

			"   VSTW   VR50, *++AR3[16]         \n\t"
			"|  VLDDWM2 *-AR1[1],  VR25:VR24  \n\t" //vr16-24,vr17-25

			"   VSTW   VR51, *++AR3[16]         \n\t"
			"|  VLDDWM2 *+AR1[0],  VR27:VR26  \n\t" //vr16-24,vr17-25

			"   VSTW   VR52, *++AR3[16]         \n\t"
			"|  VLDDWM2 *+AR1[1],   VR29:VR28  \n\t" //vr16-24,vr17-25

			"   VSTW   VR53, *++AR3[16]         \n\t"
			"|  VLDDWM2 *+AR1[6],  VR31:VR30  \n\t" //vr16-24,vr17-25

			"   VSTW   VR0, *++AR4[7]         \n\t"
			"|  VLDDWM2 *+AR1[7],  VR33:VR32  \n\t" //vr16-24,vr17-25
			"|	SLDW	  *+AR15[9], R10	\n\t" //sldw 7

			"   VSTW   VR1, *++AR4[16]      \n\t"
			"|	VLDDWM2 *+AR1[8],  VR35:VR34  \n\t" //vr16-24,vr17-25

			"   VSTW   VR2, *++AR4[16]		\n\t"
			"|  VLDDWM2 *+AR1[41],  VR37:VR36  \n\t" //vr16-24,vr17-25

			"   VSTW   VR3, *++AR4[16]      \n\t"
			"|  VLDDWM2 *+AR1[42],  VR39:VR38  \n\t" //vr16-24,vr17-25

			"   VSTW   VR4, *++AR4[16]      \n\t"
			"|  VLDDWM2 *+AR1[43],  VR41:VR40  \n\t" //vr16-24,vr17-25

			"   VSTW   VR5, *++AR4[16]    \n\t"
			"|  VLDDWM2 *+AR1[48],  VR43:VR42  \n\t" //vr16-24,vr17-25

			"   VSTW   VR6, *++AR4[16]   \n\t"
			"|  VLDDWM2 *+AR1[49],  VR45:VR44  \n\t" //vr16-24,vr17-25

			"   VSTW   VR7, *++AR4[16]   \n\t"
			"|  VLDDWM2 *+AR1[50],  VR47:VR46  \n\t" //vr16-24,vr17-25
			"|	SSHFLL	  3, R10, R10       \n\t"  //SSHFLL 1  

			"   VSTW   VR8, *++AR4[16]   \n\t"
			"|  VLDDWM2 *+AR1[55],  VR49:VR48  \n\t" //vr16-24,vr17-25

			"   VSTW   VR9, *++AR4[16]   \n\t"
			"|  VLDDWM2 *+AR1[56],  VR51:VR50  \n\t" //vr16-24,vr17-25

			"   VSTW   VR10, *++AR4[16]   \n\t"
			"|  VLDDWM2 *+AR1[57],  VR53:VR52  \n\t" //vr16-24,vr17-25

			"   VSTW   VR11, *+AR4[16]   \n\t"
			"|  VSTW   VR12, *+AR4[32]   \n\t"

			"   VSTW   VR13, *+AR4[48]   \n\t"
			"|  VSTW   VR14, *+AR4[64]   \n\t"

			"   VSTW   VR15, *+AR4[80]   \n\t"
			"|  VSTW   VR16, *+AR4[96]   \n\t"

			"   VSTW   VR17, *+AR4[112]   \n\t"
			"|  VSTW   VR18, *+AR4[128]   \n\t"

			"   VSTW   VR19, *+AR4[144]   \n\t"
			"|  VSTW   VR20, *+AR4[160]   \n\t"

			"   VSTW   VR21, *+AR4[176]   \n\t"
			"|  VSTW   VR22, *+AR4[192]   \n\t"

			"   VSTW   VR23, *+AR4[208]   \n\t"
			"|  VSTW   VR24, *+AR4[224]   \n\t"

			"   VSTW   VR25, *+AR4[240]   \n\t"
			"|  VSTW   VR26, *+AR4[256]   \n\t"

			"   VSTW   VR27, *+AR4[272]   \n\t"
			"|  VSTW   VR28, *+AR4[288]   \n\t"
			"|  SADDA     R10,AR0,AR1		\n\t" //sadda AVAF 3 SVAF 2

			"   VSTW   VR29, *+AR4[304]   \n\t"
			"|  VSTW   VR30, *+AR4[320]   \n\t"

			"   VSTW   VR31, *+AR4[336]   \n\t"
			"|  VSTW   VR32, *+AR4[352]   \n\t"

			"   VSTW   VR33, *+AR4[368]   \n\t"
			"|  VSTW   VR34, *+AR4[384]   \n\t"

			"   VSTW   VR35, *+AR4[400]   \n\t"
			"|  VSTW   VR36, *+AR4[416]   \n\t"

			"   VSTW   VR37, *++AR4[432]   \n\t"
			"|  SMVAGA.M2 R31,AR3           \n\t" //R31 
// index[9] read
			"	VSTW   VR38, *++AR4[16]   \n\t"
			"|  VLDDWM2 *-AR1[57],  VR1:VR0    \n\t" //vr0-0,vr1-1

			"	VSTW   VR39, *++AR4[16]   \n\t"
			"|  VLDDWM2 *-AR1[56],  VR3:VR2    \n\t" //vr2-3,vr3-4
			
			"	VSTW   VR40, *++AR4[16]   \n\t"
			"|  VLDDWM2 *-AR1[55],  VR5:VR4    \n\t" //vr4-6,vr5-7

			"	VSTW   VR41, *++AR4[16]   \n\t"
			"|  VLDDWM2 *-AR1[50],  VR7:VR6    \n\t" //vr6-9,vr7-10

			"	VSTW   VR42, *++AR4[16]   \n\t"
			"|  VLDDWM2 *-AR1[49],  VR9:VR8    \n\t" //vr4-12,vr5-13

			"   VSTW   VR43, *++AR4[16]         \n\t"
			"|  VLDDWM2 *-AR1[48],  VR11:VR10	\n\t" //vr10-15,vr11-16

			"   VSTW   VR44, *++AR4[16]         \n\t"
			"|  VLDDWM2 *-AR1[43],  VR13:VR12  \n\t" //vr12-18,vr13-19

			"   VSTW   VR45, *++AR4[16]         \n\t"
			"|  VLDDWM2 *-AR1[42],  VR15:VR14  \n\t" //vr14-21,vr15-22

			"   VSTW   VR46, *++AR4[16]         \n\t"
			"|  VLDDWM2 *-AR1[41],  VR17:VR16  \n\t" //vr16-24,vr17-25

			"   VSTW   VR47, *++AR4[16]         \n\t"
			"|	VLDDWM2	*-AR1[8],	 VR19:VR18	\n\t" //v18-2   

			"   VSTW   VR48, *++AR4[16]         \n\t"
			"|  VLDDWM2 *-AR1[7],  VR21:VR20  \n\t" //vr16-24,vr17-25

			"   VSTW   VR49, *++AR4[16]         \n\t"
			"|  VLDDWM2 *-AR1[6],  VR23:VR22  \n\t" //vr16-24,vr17-25

			"   VSTW   VR50, *++AR4[16]         \n\t"
			"|  VLDDWM2 *-AR1[1],  VR25:VR24  \n\t" //vr16-24,vr17-25

			"   VSTW   VR51, *++AR4[16]         \n\t"
			"|  VLDDWM2 *+AR1[0],  VR27:VR26  \n\t" //vr16-24,vr17-25

			"   VSTW   VR52, *++AR4[16]         \n\t"
			"|  VLDDWM2 *+AR1[1],   VR29:VR28  \n\t" //vr16-24,vr17-25

			"   VSTW   VR53, *++AR4[16]         \n\t"
			"|  VLDDWM2 *+AR1[6],  VR31:VR30  \n\t" //vr16-24,vr17-25

			"   VSTW   VR0, *++AR3[9]         \n\t"
			"|  VLDDWM2 *+AR1[7],  VR33:VR32  \n\t" //vr16-24,vr17-25
			"|	SLDW	  *+AR15[11], R10	\n\t" //sldw 7

			"   VSTW   VR1, *++AR3[16]      \n\t"
			"|	VLDDWM2 *+AR1[8],  VR35:VR34  \n\t" //vr16-24,vr17-25

			"   VSTW   VR2, *++AR3[16]		\n\t"
			"|  VLDDWM2 *+AR1[41],  VR37:VR36  \n\t" //vr16-24,vr17-25

			"   VSTW   VR3, *++AR3[16]      \n\t"
			"|  VLDDWM2 *+AR1[42],  VR39:VR38  \n\t" //vr16-24,vr17-25

			"   VSTW   VR4, *++AR3[16]      \n\t"
			"|  VLDDWM2 *+AR1[43],  VR41:VR40  \n\t" //vr16-24,vr17-25

			"   VSTW   VR5, *++AR3[16]    \n\t"
			"|  VLDDWM2 *+AR1[48],  VR43:VR42  \n\t" //vr16-24,vr17-25

			"   VSTW   VR6, *++AR3[16]   \n\t"
			"|  VLDDWM2 *+AR1[49],  VR45:VR44  \n\t" //vr16-24,vr17-25

			"   VSTW   VR7, *++AR3[16]   \n\t"
			"|  VLDDWM2 *+AR1[50],  VR47:VR46  \n\t" //vr16-24,vr17-25
			"|	SSHFLL	  3, R10, R10       \n\t"  //SSHFLL 1  

			"   VSTW   VR8, *++AR3[16]   \n\t"
			"|  VLDDWM2 *+AR1[55],  VR49:VR48  \n\t" //vr16-24,vr17-25

			"   VSTW   VR9, *++AR3[16]   \n\t"
			"|  VLDDWM2 *+AR1[56],  VR51:VR50  \n\t" //vr16-24,vr17-25

			"   VSTW   VR10, *++AR3[16]   \n\t"
			"|  VLDDWM2 *+AR1[57],  VR53:VR52  \n\t" //vr16-24,vr17-25

			"   VSTW   VR11, *+AR3[16]   \n\t"
			"|  VSTW   VR12, *+AR3[32]   \n\t"

			"   VSTW   VR13, *+AR3[48]   \n\t"
			"|  VSTW   VR14, *+AR3[64]   \n\t"

			"   VSTW   VR15, *+AR3[80]   \n\t"
			"|  VSTW   VR16, *+AR3[96]   \n\t"

			"   VSTW   VR17, *+AR3[112]   \n\t"
			"|  VSTW   VR18, *+AR3[128]   \n\t"

			"   VSTW   VR19, *+AR3[144]   \n\t"
			"|  VSTW   VR20, *+AR3[160]   \n\t"

			"   VSTW   VR21, *+AR3[176]   \n\t"
			"|  VSTW   VR22, *+AR3[192]   \n\t"

			"   VSTW   VR23, *+AR3[208]   \n\t"
			"|  VSTW   VR24, *+AR3[224]   \n\t"

			"   VSTW   VR25, *+AR3[240]   \n\t"
			"|  VSTW   VR26, *+AR3[256]   \n\t"

			"   VSTW   VR27, *+AR3[272]   \n\t"
			"|  VSTW   VR28, *+AR3[288]   \n\t"
			"|  SADDA     R10,AR0,AR1		\n\t" //sadda AVAF 3 SVAF 2

			"   VSTW   VR29, *+AR3[304]   \n\t"
			"|  VSTW   VR30, *+AR3[320]   \n\t"

			"   VSTW   VR31, *+AR3[336]   \n\t"
			"|  VSTW   VR32, *+AR3[352]   \n\t"

			"   VSTW   VR33, *+AR3[368]   \n\t"
			"|  VSTW   VR34, *+AR3[384]   \n\t"

			"   VSTW   VR35, *+AR3[400]   \n\t"
			"|  VSTW   VR36, *+AR3[416]   \n\t"

			"   VSTW   VR37, *++AR3[432]   \n\t"
			"|  SMVAGA.M2 R31,AR4           \n\t" //R31 
// index[11] read
			"	VSTW   VR38, *++AR3[16]   \n\t"
			"|  VLDDWM2 *-AR1[57],  VR1:VR0    \n\t" //vr0-0,vr1-1

			"	VSTW   VR39, *++AR3[16]   \n\t"
			"|  VLDDWM2 *-AR1[56],  VR3:VR2    \n\t" //vr2-3,vr3-4
			
			"	VSTW   VR40, *++AR3[16]   \n\t"
			"|  VLDDWM2 *-AR1[55],  VR5:VR4    \n\t" //vr4-6,vr5-7

			"	VSTW   VR41, *++AR3[16]   \n\t"
			"|  VLDDWM2 *-AR1[50],  VR7:VR6    \n\t" //vr6-9,vr7-10

			"	VSTW   VR42, *++AR3[16]   \n\t"
			"|  VLDDWM2 *-AR1[49],  VR9:VR8    \n\t" //vr4-12,vr5-13

			"   VSTW   VR43, *++AR3[16]         \n\t"
			"|  VLDDWM2 *-AR1[48],  VR11:VR10	\n\t" //vr10-15,vr11-16

			"   VSTW   VR44, *++AR3[16]         \n\t"
			"|  VLDDWM2 *-AR1[43],  VR13:VR12  \n\t" //vr12-18,vr13-19

			"   VSTW   VR45, *++AR3[16]         \n\t"
			"|  VLDDWM2 *-AR1[42],  VR15:VR14  \n\t" //vr14-21,vr15-22

			"   VSTW   VR46, *++AR3[16]         \n\t"
			"|  VLDDWM2 *-AR1[41],  VR17:VR16  \n\t" //vr16-24,vr17-25

			"   VSTW   VR47, *++AR3[16]         \n\t"
			"|	VLDDWM2	*-AR1[8],	 VR19:VR18	\n\t" //v18-2   

			"   VSTW   VR48, *++AR3[16]         \n\t"
			"|  VLDDWM2 *-AR1[7],  VR21:VR20  \n\t" //vr16-24,vr17-25

			"   VSTW   VR49, *++AR3[16]         \n\t"
			"|  VLDDWM2 *-AR1[6],  VR23:VR22  \n\t" //vr16-24,vr17-25

			"   VSTW   VR50, *++AR3[16]         \n\t"
			"|  VLDDWM2 *-AR1[1],  VR25:VR24  \n\t" //vr16-24,vr17-25

			"   VSTW   VR51, *++AR3[16]         \n\t"
			"|  VLDDWM2 *+AR1[0],  VR27:VR26  \n\t" //vr16-24,vr17-25

			"   VSTW   VR52, *++AR3[16]         \n\t"
			"|  VLDDWM2 *+AR1[1],   VR29:VR28  \n\t" //vr16-24,vr17-25

			"   VSTW   VR53, *++AR3[16]         \n\t"
			"|  VLDDWM2 *+AR1[6],  VR31:VR30  \n\t" //vr16-24,vr17-25

			"   VSTW   VR0, *++AR4[11]         \n\t"
			"|  VLDDWM2 *+AR1[7],  VR33:VR32  \n\t" //vr16-24,vr17-25
			"|	SLDW	  *+AR15[13], R10	\n\t" //sldw 7

			"   VSTW   VR1, *++AR4[16]      \n\t"
			"|	VLDDWM2 *+AR1[8],  VR35:VR34  \n\t" //vr16-24,vr17-25

			"   VSTW   VR2, *++AR4[16]		\n\t"
			"|  VLDDWM2 *+AR1[41],  VR37:VR36  \n\t" //vr16-24,vr17-25

			"   VSTW   VR3, *++AR4[16]      \n\t"
			"|  VLDDWM2 *+AR1[42],  VR39:VR38  \n\t" //vr16-24,vr17-25

			"   VSTW   VR4, *++AR4[16]      \n\t"
			"|  VLDDWM2 *+AR1[43],  VR41:VR40  \n\t" //vr16-24,vr17-25

			"   VSTW   VR5, *++AR4[16]    \n\t"
			"|  VLDDWM2 *+AR1[48],  VR43:VR42  \n\t" //vr16-24,vr17-25

			"   VSTW   VR6, *++AR4[16]   \n\t"
			"|  VLDDWM2 *+AR1[49],  VR45:VR44  \n\t" //vr16-24,vr17-25

			"   VSTW   VR7, *++AR4[16]   \n\t"
			"|  VLDDWM2 *+AR1[50],  VR47:VR46  \n\t" //vr16-24,vr17-25
			"|	SSHFLL	  3, R10, R10       \n\t"  //SSHFLL 1  

			"   VSTW   VR8, *++AR4[16]   \n\t"
			"|  VLDDWM2 *+AR1[55],  VR49:VR48  \n\t" //vr16-24,vr17-25

			"   VSTW   VR9, *++AR4[16]   \n\t"
			"|  VLDDWM2 *+AR1[56],  VR51:VR50  \n\t" //vr16-24,vr17-25

			"   VSTW   VR10, *++AR4[16]   \n\t"
			"|  VLDDWM2 *+AR1[57],  VR53:VR52  \n\t" //vr16-24,vr17-25

			"   VSTW   VR11, *+AR4[16]   \n\t"
			"|  VSTW   VR12, *+AR4[32]   \n\t"

			"   VSTW   VR13, *+AR4[48]   \n\t"
			"|  VSTW   VR14, *+AR4[64]   \n\t"

			"   VSTW   VR15, *+AR4[80]   \n\t"
			"|  VSTW   VR16, *+AR4[96]   \n\t"

			"   VSTW   VR17, *+AR4[112]   \n\t"
			"|  VSTW   VR18, *+AR4[128]   \n\t"

			"   VSTW   VR19, *+AR4[144]   \n\t"
			"|  VSTW   VR20, *+AR4[160]   \n\t"

			"   VSTW   VR21, *+AR4[176]   \n\t"
			"|  VSTW   VR22, *+AR4[192]   \n\t"

			"   VSTW   VR23, *+AR4[208]   \n\t"
			"|  VSTW   VR24, *+AR4[224]   \n\t"

			"   VSTW   VR25, *+AR4[240]   \n\t"
			"|  VSTW   VR26, *+AR4[256]   \n\t"

			"   VSTW   VR27, *+AR4[272]   \n\t"
			"|  VSTW   VR28, *+AR4[288]   \n\t"
			"|  SADDA     R10,AR0,AR1		\n\t" //sadda AVAF 3 SVAF 2

			"   VSTW   VR29, *+AR4[304]   \n\t"
			"|  VSTW   VR30, *+AR4[320]   \n\t"

			"   VSTW   VR31, *+AR4[336]   \n\t"
			"|  VSTW   VR32, *+AR4[352]   \n\t"

			"   VSTW   VR33, *+AR4[368]   \n\t"
			"|  VSTW   VR34, *+AR4[384]   \n\t"

			"   VSTW   VR35, *+AR4[400]   \n\t"
			"|  VSTW   VR36, *+AR4[416]   \n\t"

			"   VSTW   VR37, *++AR4[432]   \n\t"
			"|  SMVAGA.M2 R31,AR3           \n\t" //R31 
// index[13] read
			"	VSTW   VR38, *++AR4[16]   \n\t"
			"|  VLDDWM2 *-AR1[57],  VR1:VR0    \n\t" //vr0-0,vr1-1

			"	VSTW   VR39, *++AR4[16]   \n\t"
			"|  VLDDWM2 *-AR1[56],  VR3:VR2    \n\t" //vr2-3,vr3-4
			
			"	VSTW   VR40, *++AR4[16]   \n\t"
			"|  VLDDWM2 *-AR1[55],  VR5:VR4    \n\t" //vr4-6,vr5-7

			"	VSTW   VR41, *++AR4[16]   \n\t"
			"|  VLDDWM2 *-AR1[50],  VR7:VR6    \n\t" //vr6-9,vr7-10

			"	VSTW   VR42, *++AR4[16]   \n\t"
			"|  VLDDWM2 *-AR1[49],  VR9:VR8    \n\t" //vr4-12,vr5-13

			"   VSTW   VR43, *++AR4[16]         \n\t"
			"|  VLDDWM2 *-AR1[48],  VR11:VR10	\n\t" //vr10-15,vr11-16

			"   VSTW   VR44, *++AR4[16]         \n\t"
			"|  VLDDWM2 *-AR1[43],  VR13:VR12  \n\t" //vr12-18,vr13-19

			"   VSTW   VR45, *++AR4[16]         \n\t"
			"|  VLDDWM2 *-AR1[42],  VR15:VR14  \n\t" //vr14-21,vr15-22

			"   VSTW   VR46, *++AR4[16]         \n\t"
			"|  VLDDWM2 *-AR1[41],  VR17:VR16  \n\t" //vr16-24,vr17-25

			"   VSTW   VR47, *++AR4[16]         \n\t"
			"|	VLDDWM2	*-AR1[8],	 VR19:VR18	\n\t" //v18-2   

			"   VSTW   VR48, *++AR4[16]         \n\t"
			"|  VLDDWM2 *-AR1[7],  VR21:VR20  \n\t" //vr16-24,vr17-25

			"   VSTW   VR49, *++AR4[16]         \n\t"
			"|  VLDDWM2 *-AR1[6],  VR23:VR22  \n\t" //vr16-24,vr17-25

			"   VSTW   VR50, *++AR4[16]         \n\t"
			"|  VLDDWM2 *-AR1[1],  VR25:VR24  \n\t" //vr16-24,vr17-25

			"   VSTW   VR51, *++AR4[16]         \n\t"
			"|  VLDDWM2 *+AR1[0],  VR27:VR26  \n\t" //vr16-24,vr17-25

			"   VSTW   VR52, *++AR4[16]         \n\t"
			"|  VLDDWM2 *+AR1[1],   VR29:VR28  \n\t" //vr16-24,vr17-25

			"   VSTW   VR53, *++AR4[16]         \n\t"
			"|  VLDDWM2 *+AR1[6],  VR31:VR30  \n\t" //vr16-24,vr17-25

			"   VSTW   VR0, *++AR3[13]         \n\t"
			"|  VLDDWM2 *+AR1[7],  VR33:VR32  \n\t" //vr16-24,vr17-25
			"|	SLDW	  *+AR15[15], R10	\n\t" //sldw 7

			"   VSTW   VR1, *++AR3[16]      \n\t"
			"|	VLDDWM2 *+AR1[8],  VR35:VR34  \n\t" //vr16-24,vr17-25

			"   VSTW   VR2, *++AR3[16]		\n\t"
			"|  VLDDWM2 *+AR1[41],  VR37:VR36  \n\t" //vr16-24,vr17-25

			"   VSTW   VR3, *++AR3[16]      \n\t"
			"|  VLDDWM2 *+AR1[42],  VR39:VR38  \n\t" //vr16-24,vr17-25

			"   VSTW   VR4, *++AR3[16]      \n\t"
			"|  VLDDWM2 *+AR1[43],  VR41:VR40  \n\t" //vr16-24,vr17-25

			"   VSTW   VR5, *++AR3[16]    \n\t"
			"|  VLDDWM2 *+AR1[48],  VR43:VR42  \n\t" //vr16-24,vr17-25

			"   VSTW   VR6, *++AR3[16]   \n\t"
			"|  VLDDWM2 *+AR1[49],  VR45:VR44  \n\t" //vr16-24,vr17-25

			"   VSTW   VR7, *++AR3[16]   \n\t"
			"|  VLDDWM2 *+AR1[50],  VR47:VR46  \n\t" //vr16-24,vr17-25
			"|	SSHFLL	  3, R10, R10       \n\t"  //SSHFLL 1  

			"   VSTW   VR8, *++AR3[16]   \n\t"
			"|  VLDDWM2 *+AR1[55],  VR49:VR48  \n\t" //vr16-24,vr17-25

			"   VSTW   VR9, *++AR3[16]   \n\t"
			"|  VLDDWM2 *+AR1[56],  VR51:VR50  \n\t" //vr16-24,vr17-25

			"   VSTW   VR10, *++AR3[16]   \n\t"
			"|  VLDDWM2 *+AR1[57],  VR53:VR52  \n\t" //vr16-24,vr17-25

			"   VSTW   VR11, *+AR3[16]   \n\t"
			"|  VSTW   VR12, *+AR3[32]   \n\t"

			"   VSTW   VR13, *+AR3[48]   \n\t"
			"|  VSTW   VR14, *+AR3[64]   \n\t"

			"   VSTW   VR15, *+AR3[80]   \n\t"
			"|  VSTW   VR16, *+AR3[96]   \n\t"

			"   VSTW   VR17, *+AR3[112]   \n\t"
			"|  VSTW   VR18, *+AR3[128]   \n\t"

			"   VSTW   VR19, *+AR3[144]   \n\t"
			"|  VSTW   VR20, *+AR3[160]   \n\t"

			"   VSTW   VR21, *+AR3[176]   \n\t"
			"|  VSTW   VR22, *+AR3[192]   \n\t"

			"   VSTW   VR23, *+AR3[208]   \n\t"
			"|  VSTW   VR24, *+AR3[224]   \n\t"

			"   VSTW   VR25, *+AR3[240]   \n\t"
			"|  VSTW   VR26, *+AR3[256]   \n\t"

			"   VSTW   VR27, *+AR3[272]   \n\t"
			"|  VSTW   VR28, *+AR3[288]   \n\t"
			"|  SADDA     R10,AR0,AR1		\n\t" //sadda AVAF 3 SVAF 2

			"   VSTW   VR29, *+AR3[304]   \n\t"
			"|  VSTW   VR30, *+AR3[320]   \n\t"

			"   VSTW   VR31, *+AR3[336]   \n\t"
			"|  VSTW   VR32, *+AR3[352]   \n\t"

			"   VSTW   VR33, *+AR3[368]   \n\t"
			"|  VSTW   VR34, *+AR3[384]   \n\t"

			"   VSTW   VR35, *+AR3[400]   \n\t"
			"|  VSTW   VR36, *+AR3[416]   \n\t"

			"   VSTW   VR37, *++AR3[432]   \n\t"
			"|  SMVAGA.M2 R31,AR4           \n\t" //R31 // index[15] read
			"	VSTW   VR38, *++AR3[16]   \n\t"
			"|  VLDDWM2 *-AR1[57],  VR1:VR0    \n\t" //vr0-0,vr1-1

			"	VSTW   VR39, *++AR3[16]   \n\t"
			"|  VLDDWM2 *-AR1[56],  VR3:VR2    \n\t" //vr2-3,vr3-4
			
			"	VSTW   VR40, *++AR3[16]   \n\t"
			"|  VLDDWM2 *-AR1[55],  VR5:VR4    \n\t" //vr4-6,vr5-7

			"	VSTW   VR41, *++AR3[16]   \n\t"
			"|  VLDDWM2 *-AR1[50],  VR7:VR6    \n\t" //vr6-9,vr7-10

			"	VSTW   VR42, *++AR3[16]   \n\t"
			"|  VLDDWM2 *-AR1[49],  VR9:VR8    \n\t" //vr4-12,vr5-13

			"   VSTW   VR43, *++AR3[16]         \n\t"
			"|  VLDDWM2 *-AR1[48],  VR11:VR10	\n\t" //vr10-15,vr11-16

			"   VSTW   VR44, *++AR3[16]         \n\t"
			"|  VLDDWM2 *-AR1[43],  VR13:VR12  \n\t" //vr12-18,vr13-19

			"   VSTW   VR45, *++AR3[16]         \n\t"
			"|  VLDDWM2 *-AR1[42],  VR15:VR14  \n\t" //vr14-21,vr15-22

			"   VSTW   VR46, *++AR3[16]         \n\t"
			"|  VLDDWM2 *-AR1[41],  VR17:VR16  \n\t" //vr16-24,vr17-25

			"   VSTW   VR47, *++AR3[16]         \n\t"
			"|	VLDDWM2	*-AR1[8],	 VR19:VR18	\n\t" //v18-2   

			"   VSTW   VR48, *++AR3[16]         \n\t"
			"|  VLDDWM2 *-AR1[7],  VR21:VR20  \n\t" //vr16-24,vr17-25

			"   VSTW   VR49, *++AR3[16]         \n\t"
			"|  VLDDWM2 *-AR1[6],  VR23:VR22  \n\t" //vr16-24,vr17-25

			"   VSTW   VR50, *++AR3[16]         \n\t"
			"|  VLDDWM2 *-AR1[1],  VR25:VR24  \n\t" //vr16-24,vr17-25

			"   VSTW   VR51, *++AR3[16]         \n\t"
			"|  VLDDWM2 *+AR1[0],  VR27:VR26  \n\t" //vr16-24,vr17-25

			"   VSTW   VR52, *++AR3[16]         \n\t"
			"|  VLDDWM2 *+AR1[1],   VR29:VR28  \n\t" //vr16-24,vr17-25

			"   VSTW   VR53, *++AR3[16]         \n\t"
			"|  VLDDWM2 *+AR1[6],  VR31:VR30  \n\t" //vr16-24,vr17-25

			"   VSTW   VR0, *++AR4[15]         \n\t"
			"|  VLDDWM2 *+AR1[7],  VR33:VR32  \n\t" //vr16-24,vr17-25
			"|	SLDW	  *+AR15[0], R10	\n\t" //sldw 7

			"   VSTW   VR1, *++AR4[16]      \n\t"
			"|	VLDDWM2 *+AR1[8],  VR35:VR34  \n\t" //vr16-24,vr17-25

			"   VSTW   VR2, *++AR4[16]		\n\t"
			"|  VLDDWM2 *+AR1[41],  VR37:VR36  \n\t" //vr16-24,vr17-25

			"   VSTW   VR3, *++AR4[16]      \n\t"
			"|  VLDDWM2 *+AR1[42],  VR39:VR38  \n\t" //vr16-24,vr17-25

			"   VSTW   VR4, *++AR4[16]      \n\t"
			"|  VLDDWM2 *+AR1[43],  VR41:VR40  \n\t" //vr16-24,vr17-25

			"   VSTW   VR5, *++AR4[16]    \n\t"
			"|  VLDDWM2 *+AR1[48],  VR43:VR42  \n\t" //vr16-24,vr17-25

			"   VSTW   VR6, *++AR4[16]   \n\t"
			"|  VLDDWM2 *+AR1[49],  VR45:VR44  \n\t" //vr16-24,vr17-25

			"   VSTW   VR7, *++AR4[16]   \n\t"
			"|  VLDDWM2 *+AR1[50],  VR47:VR46  \n\t" //vr16-24,vr17-25
			"|	SSHFLL	  3, R10, R10       \n\t"  //SSHFLL 1  

			"   VSTW   VR8, *++AR4[16]   \n\t"
			"|  VLDDWM2 *+AR1[55],  VR49:VR48  \n\t" //vr16-24,vr17-25

			"   VSTW   VR9, *++AR4[16]   \n\t"
			"|  VLDDWM2 *+AR1[56],  VR51:VR50  \n\t" //vr16-24,vr17-25

			"   VSTW   VR10, *++AR4[16]   \n\t"
			"|  VLDDWM2 *+AR1[57],  VR53:VR52  \n\t" //vr16-24,vr17-25

			"   VSTW   VR11, *+AR4[16]   \n\t"
			"|  VSTW   VR12, *+AR4[32]   \n\t"
			"|  SMULIU.M1 R29,%[eb2],R30     \n\t" //smuliu 3

			"   VSTW   VR13, *+AR4[48]   \n\t"
			"|  VSTW   VR14, *+AR4[64]   \n\t"

			"   VSTW   VR15, *+AR4[80]   \n\t"
			"|  VSTW   VR16, *+AR4[96]   \n\t"

			"   VSTW   VR17, *+AR4[112]   \n\t"
			"|  VSTW   VR18, *+AR4[128]   \n\t"
			"|  SMVAGA.M1 R30,AR0           \n\t" //smvaga 2

			"   VSTW   VR19, *+AR4[144]   \n\t"
			"|  VSTW   VR20, *+AR4[160]   \n\t"

			"   VSTW   VR21, *+AR4[176]   \n\t"
			"|  VSTW   VR22, *+AR4[192]   \n\t"

			"   VSTW   VR23, *+AR4[208]   \n\t"
			"|  VSTW   VR24, *+AR4[224]   \n\t"

			"   VSTW   VR25, *+AR4[240]   \n\t"
			"|  VSTW   VR26, *+AR4[256]   \n\t"

			"   VSTW   VR27, *+AR4[272]   \n\t"
			"|  VSTW   VR28, *+AR4[288]   \n\t"
			"|  SADDA     R10,AR0,AR1		\n\t" //sadda AVAF 3 SVAF 2

			"   VSTW   VR29, *+AR4[304]   \n\t"
			"|  VSTW   VR30, *+AR4[320]   \n\t"
			"|	SMOVI24	0x1B00,R11		  \n\t" // offset 54*16*8 bytes

			"   VSTW   VR31, *+AR4[336]   \n\t"
			"|  VSTW   VR32, *+AR4[352]   \n\t"
			"|	SADD	R11,R31,R31		  \n\t" //R31 EB2 addr for write

			"   VSTW   VR33, *+AR4[368]   \n\t"
			"|  VSTW   VR34, *+AR4[384]   \n\t"

			"   VSTW   VR35, *+AR4[400]   \n\t"
			"|  VSTW   VR36, *+AR4[416]   \n\t"

			"   VSTW   VR37, *+AR4[432]   \n\t"
			"|  SMVAGA.M2 R31,AR3           \n\t" //R31 
			"	VSTW   VR38, *+AR4[448]   \n\t"
			"|	VSTW   VR39, *+AR4[464]   \n\t"

			"	VSTW   VR40, *+AR4[480]   \n\t"
			"|	VSTW   VR41, *+AR4[496]   \n\t"

			"	VSTW   VR42, *+AR4[512]   \n\t"
			"|  VSTW   VR43, *+AR4[528]         \n\t"

			"   VSTW   VR44, *+AR4[544]         \n\t"
			"|  VSTW   VR45, *+AR4[560]         \n\t"

			"   VSTW   VR46, *+AR4[576]         \n\t"
			"|  VSTW   VR47, *+AR4[592]         \n\t"

			"   VSTW   VR48, *+AR4[608]         \n\t"
			"|  VSTW   VR49, *+AR4[624]         \n\t"

			"   VSTW   VR50, *+AR4[640]         \n\t"
			"|  VSTW   VR51, *+AR4[656]         \n\t"

			"   VSTW   VR52, *+AR4[672]         \n\t"
			"|  VSTW   VR53, *+AR4[688]         \n\t"

// EB2 index[0] // index[15] write
			"   VLDDWM2 *-AR1[57],  VR1:VR0    \n\t" //vr0-0,vr1-1
			"|  VLDDWM2 *-AR1[56],  VR3:VR2    \n\t" //vr2-3,vr3-4

			"   VLDDWM2 *-AR1[55],  VR5:VR4    \n\t" //vr4-6,vr5-7
			"|  VLDDWM2 *-AR1[50],  VR7:VR6    \n\t" //vr6-9,vr7-10

			"   VLDDWM2 *-AR1[49],  VR9:VR8    \n\t" //vr4-12,vr5-13
			"|  VLDDWM2 *-AR1[48],  VR11:VR10	\n\t" //vr10-15,vr11-16

			"   VLDDWM2 *-AR1[43],  VR13:VR12  \n\t" //vr12-18,vr13-19
			"|  VLDDWM2 *-AR1[42],  VR15:VR14  \n\t" //vr14-21,vr15-22

			"   VLDDWM2 *-AR1[41],  VR17:VR16  \n\t" //vr16-24,vr17-25
			"|	VLDDWM2	*-AR1[8],	 VR19:VR18	\n\t" //v18-2   

			"   VLDDWM2 *-AR1[7],  VR21:VR20  \n\t" //vr16-24,vr17-25
			"|  VLDDWM2 *-AR1[6],  VR23:VR22  \n\t" //vr16-24,vr17-25

			"   VLDDWM2 *-AR1[1],  VR25:VR24  \n\t" //vr16-24,vr17-25
			"|  VLDDWM2 *+AR1[0],  VR27:VR26  \n\t" //vr16-24,vr17-25

			"   VLDDWM2 *+AR1[1],   VR29:VR28  \n\t" //vr16-24,vr17-25
			"|  VLDDWM2 *+AR1[6],  VR31:VR30  \n\t" //vr16-24,vr17-25

			"	SNOP	1\n\t"

			"   VSTDW   VR1:VR0, *+AR3[0]         \n\t"
			"|  VLDDWM2 *+AR1[7],  VR33:VR32  \n\t" //vr16-24,vr17-25
			"|	SLDW	  *+AR15[2], R10	\n\t" //sldw 7

			"   VSTDW   VR3:VR2, *+AR3[16]      \n\t"
			"|  VLDDWM2 *+AR1[8],  VR35:VR34  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR5:VR4, *+AR3[32]		\n\t"
			"|  VLDDWM2 *+AR1[41],  VR37:VR36  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR7:VR6, *+AR3[48]      \n\t"
			"|  VLDDWM2 *+AR1[42],  VR39:VR38  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR9:VR8, *+AR3[64]      \n\t"
			"|  VLDDWM2 *+AR1[43],  VR41:VR40  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR11:VR10, *+AR3[80]    \n\t"
			"|  VLDDWM2 *+AR1[48],  VR43:VR42  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR13:VR12, *+AR3[96]   \n\t"
			"|  VLDDWM2 *+AR1[49],  VR45:VR44  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR15:VR14, *+AR3[112]   \n\t"
			"|  VLDDWM2 *+AR1[50],  VR47:VR46  \n\t" //vr16-24,vr17-25
			"|	SSHFLL	  3, R10, R10       \n\t"  //SSHFLL 1  

			"   VSTDW   VR17:VR16, *+AR3[128]   \n\t"
			"|  VLDDWM2 *+AR1[55],  VR49:VR48  \n\t" //vr16-24,vr17-25
			"|  SADDA     R10,AR0,AR2		\n\t" //sadda AVAF 3 SVAF 2

			"   VSTDW   VR19:VR18, *+AR3[144]   \n\t"
			"|  VLDDWM2 *+AR1[56],  VR51:VR50  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR21:VR20, *+AR3[160]   \n\t"
			"|  VLDDWM2 *+AR1[57],  VR53:VR52  \n\t" //vr16-24,vr17-25
////////////////////////////////////////////////////////////////////////////
// index[2] read
			"	VSTDW   VR23:VR22, *+AR3[176]   \n\t"
			"|  VLDDWM2 *-AR2[57],  VR1:VR0    \n\t" //vr0-0,vr1-1

			"	VSTDW   VR25:VR24, *+AR3[192]   \n\t"
			"|  VLDDWM2 *-AR2[56],  VR3:VR2    \n\t" //vr2-3,vr3-4
			
			"	VSTDW   VR27:VR26, *+AR3[208]   \n\t"
			"|  VLDDWM2 *-AR2[55],  VR5:VR4    \n\t" //vr4-6,vr5-7

			"	VSTDW   VR29:VR28, *+AR3[224]   \n\t"
			"|  VLDDWM2 *-AR2[50],  VR7:VR6    \n\t" //vr6-9,vr7-10

			"	VSTDW   VR31:VR30, *+AR3[240]   \n\t"
			"|  VLDDWM2 *-AR2[49],  VR9:VR8    \n\t" //vr4-12,vr5-13

			"   VSTDW   VR33:VR32, *+AR3[256]         \n\t"
			"|  VLDDWM2 *-AR2[48],  VR11:VR10	\n\t" //vr10-15,vr11-16

			"   VSTDW   VR35:VR34, *+AR3[272]         \n\t"
			"|  VLDDWM2 *-AR2[43],  VR13:VR12  \n\t" //vr12-18,vr13-19

			"   VSTDW   VR37:VR36, *+AR3[288]         \n\t"
			"|  VLDDWM2 *-AR2[42],  VR15:VR14  \n\t" //vr14-21,vr15-22

			"   VSTDW   VR39:VR38, *+AR3[304]         \n\t"
			"|  VLDDWM2 *-AR2[41],  VR17:VR16  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR41:VR40, *+AR3[320]         \n\t"
			"|	VLDDWM2	*-AR2[8],	 VR19:VR18	\n\t" //v18-2   

			"   VSTDW   VR43:VR42, *+AR3[336]         \n\t"
			"|  VLDDWM2 *-AR2[7],  VR21:VR20  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR45:VR44, *+AR3[352]         \n\t"
			"|  VLDDWM2 *-AR2[6],  VR23:VR22  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR47:VR46, *+AR3[368]         \n\t"
			"|  VLDDWM2 *-AR2[1],  VR25:VR24  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR49:VR48, *+AR3[384]         \n\t"
			"|  VLDDWM2 *+AR2[0],  VR27:VR26  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR51:VR50, *+AR3[400]         \n\t"
			"|  VLDDWM2 *+AR2[1],   VR29:VR28  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR53:VR52, *+AR3[416]         \n\t"
			"|  VLDDWM2 *+AR2[6],  VR31:VR30  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR1:VR0, *+AR3[1]         \n\t"
			"|  VLDDWM2 *+AR2[7],  VR33:VR32  \n\t" //vr16-24,vr17-25
			"|	SLDW	  *+AR15[4], R10	\n\t" //sldw 7

			"   VSTDW   VR3:VR2, *+AR3[17]      \n\t"
			"|	VLDDWM2 *+AR2[8],  VR35:VR34  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR5:VR4, *+AR3[33]		\n\t"
			"|  VLDDWM2 *+AR2[41],  VR37:VR36  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR7:VR6, *+AR3[49]      \n\t"
			"|  VLDDWM2 *+AR2[42],  VR39:VR38  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR9:VR8, *+AR3[65]      \n\t"
			"|  VLDDWM2 *+AR2[43],  VR41:VR40  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR11:VR10, *+AR3[81]    \n\t"
			"|  VLDDWM2 *+AR2[48],  VR43:VR42  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR13:VR12, *+AR3[97]   \n\t"
			"|  VLDDWM2 *+AR2[49],  VR45:VR44  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR15:VR14, *+AR3[113]   \n\t"
			"|  VLDDWM2 *+AR2[50],  VR47:VR46  \n\t" //vr16-24,vr17-25
			"|	SSHFLL	  3, R10, R10       \n\t"  //SSHFLL 1  

			"   VSTDW   VR17:VR16, *+AR3[129]   \n\t"
			"|  VLDDWM2 *+AR2[55],  VR49:VR48  \n\t" //vr16-24,vr17-25
			"|  SADDA     R10,AR0,AR1		\n\t" //sadda AVAF 3 SVAF 2

			"   VSTDW   VR19:VR18, *+AR3[145]   \n\t"
			"|  VLDDWM2 *+AR2[56],  VR51:VR50  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR21:VR20, *+AR3[161]   \n\t"
			"|  VLDDWM2 *+AR2[57],  VR53:VR52  \n\t" //vr16-24,vr17-25
// index[4] read
			"	VSTDW   VR23:VR22, *+AR3[177]   \n\t"
			"|  VLDDWM2 *-AR1[57],  VR1:VR0    \n\t" //vr0-0,vr1-1

			"	VSTDW   VR25:VR24, *+AR3[193]   \n\t"
			"|  VLDDWM2 *-AR1[56],  VR3:VR2    \n\t" //vr2-3,vr3-4
			
			"	VSTDW   VR27:VR26, *+AR3[209]   \n\t"
			"|  VLDDWM2 *-AR1[55],  VR5:VR4    \n\t" //vr4-6,vr5-7

			"	VSTDW   VR29:VR28, *+AR3[225]   \n\t"
			"|  VLDDWM2 *-AR1[50],  VR7:VR6    \n\t" //vr6-9,vr7-10

			"	VSTDW   VR31:VR30, *+AR3[241]   \n\t"
			"|  VLDDWM2 *-AR1[49],  VR9:VR8    \n\t" //vr4-12,vr5-13

			"   VSTDW   VR33:VR32, *+AR3[257]         \n\t"
			"|  VLDDWM2 *-AR1[48],  VR11:VR10	\n\t" //vr10-15,vr11-16

			"   VSTDW   VR35:VR34, *+AR3[273]         \n\t"
			"|  VLDDWM2 *-AR1[43],  VR13:VR12  \n\t" //vr12-18,vr13-19

			"   VSTDW   VR37:VR36, *+AR3[289]         \n\t"
			"|  VLDDWM2 *-AR1[42],  VR15:VR14  \n\t" //vr14-21,vr15-22

			"   VSTDW   VR39:VR38, *+AR3[305]         \n\t"
			"|  VLDDWM2 *-AR1[41],  VR17:VR16  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR41:VR40, *+AR3[321]         \n\t"
			"|	VLDDWM2	*-AR1[8],	 VR19:VR18	\n\t" //v18-2   

			"   VSTDW   VR43:VR42, *+AR3[337]         \n\t"
			"|  VLDDWM2 *-AR1[7],  VR21:VR20  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR45:VR44, *+AR3[353]         \n\t"
			"|  VLDDWM2 *-AR1[6],  VR23:VR22  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR47:VR46, *+AR3[369]         \n\t"
			"|  VLDDWM2 *-AR1[1],  VR25:VR24  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR49:VR48, *+AR3[385]         \n\t"
			"|  VLDDWM2 *+AR1[0],  VR27:VR26  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR51:VR50, *+AR3[401]         \n\t"
			"|  VLDDWM2 *+AR1[1],   VR29:VR28  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR53:VR52, *+AR3[417]         \n\t"
			"|  VLDDWM2 *+AR1[6],  VR31:VR30  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR1:VR0, *+AR3[2]         \n\t"
			"|  VLDDWM2 *+AR1[7],  VR33:VR32  \n\t" //vr16-24,vr17-25
			"|	SLDW	  *+AR15[6], R10	\n\t" //sldw 7

			"   VSTDW   VR3:VR2, *+AR3[18]      \n\t"
			"|	VLDDWM2 *+AR1[8],  VR35:VR34  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR5:VR4, *+AR3[34]		\n\t"
			"|  VLDDWM2 *+AR1[41],  VR37:VR36  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR7:VR6, *+AR3[50]      \n\t"
			"|  VLDDWM2 *+AR1[42],  VR39:VR38  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR9:VR8, *+AR3[66]      \n\t"
			"|  VLDDWM2 *+AR1[43],  VR41:VR40  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR11:VR10, *+AR3[82]    \n\t"
			"|  VLDDWM2 *+AR1[48],  VR43:VR42  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR13:VR12, *+AR3[98]   \n\t"
			"|  VLDDWM2 *+AR1[49],  VR45:VR44  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR15:VR14, *+AR3[114]   \n\t"
			"|  VLDDWM2 *+AR1[50],  VR47:VR46  \n\t" //vr16-24,vr17-25
			"|	SSHFLL	  3, R10, R10       \n\t"  //SSHFLL 1  

			"   VSTDW   VR17:VR16, *+AR3[130]   \n\t"
			"|  VLDDWM2 *+AR1[55],  VR49:VR48  \n\t" //vr16-24,vr17-25
			"|  SADDA     R10,AR0,AR2		\n\t" //sadda AVAF 3 SVAF 2

			"   VSTDW   VR19:VR18, *+AR3[146]   \n\t"
			"|  VLDDWM2 *+AR1[56],  VR51:VR50  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR21:VR20, *+AR3[162]   \n\t"
			"|  VLDDWM2 *+AR1[57],  VR53:VR52  \n\t" //vr16-24,vr17-25
// index[6] read
			"	VSTDW   VR23:VR22, *+AR3[178]   \n\t"
			"|  VLDDWM2 *-AR2[57],  VR1:VR0    \n\t" //vr0-0,vr1-1

			"	VSTDW   VR25:VR24, *+AR3[194]   \n\t"
			"|  VLDDWM2 *-AR2[56],  VR3:VR2    \n\t" //vr2-3,vr3-4
			
			"	VSTDW   VR27:VR26, *+AR3[210]   \n\t"
			"|  VLDDWM2 *-AR2[55],  VR5:VR4    \n\t" //vr4-6,vr5-7

			"	VSTDW   VR29:VR28, *+AR3[226]   \n\t"
			"|  VLDDWM2 *-AR2[50],  VR7:VR6    \n\t" //vr6-9,vr7-10

			"	VSTDW   VR31:VR30, *+AR3[242]   \n\t"
			"|  VLDDWM2 *-AR2[49],  VR9:VR8    \n\t" //vr4-12,vr5-13

			"   VSTDW   VR33:VR32, *+AR3[258]         \n\t"
			"|  VLDDWM2 *-AR2[48],  VR11:VR10	\n\t" //vr10-15,vr11-16

			"   VSTDW   VR35:VR34, *+AR3[274]         \n\t"
			"|  VLDDWM2 *-AR2[43],  VR13:VR12  \n\t" //vr12-18,vr13-19

			"   VSTDW   VR37:VR36, *+AR3[290]         \n\t"
			"|  VLDDWM2 *-AR2[42],  VR15:VR14  \n\t" //vr14-21,vr15-22

			"   VSTDW   VR39:VR38, *+AR3[306]         \n\t"
			"|  VLDDWM2 *-AR2[41],  VR17:VR16  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR41:VR40, *+AR3[322]         \n\t"
			"|	VLDDWM2	*-AR2[8],	 VR19:VR18	\n\t" //v18-2   

			"   VSTDW   VR43:VR42, *+AR3[338]         \n\t"
			"|  VLDDWM2 *-AR2[7],  VR21:VR20  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR45:VR44, *+AR3[354]         \n\t"
			"|  VLDDWM2 *-AR2[6],  VR23:VR22  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR47:VR46, *+AR3[370]         \n\t"
			"|  VLDDWM2 *-AR2[1],  VR25:VR24  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR49:VR48, *+AR3[386]         \n\t"
			"|  VLDDWM2 *+AR2[0],  VR27:VR26  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR51:VR50, *+AR3[402]         \n\t"
			"|  VLDDWM2 *+AR2[1],   VR29:VR28  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR53:VR52, *+AR3[418]         \n\t"
			"|  VLDDWM2 *+AR2[6],  VR31:VR30  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR1:VR0, *+AR3[3]         \n\t"
			"|  VLDDWM2 *+AR2[7],  VR33:VR32  \n\t" //vr16-24,vr17-25
			"|	SLDW	  *+AR15[8], R10	\n\t" //sldw 7

			"   VSTDW   VR3:VR2, *+AR3[19]      \n\t"
			"|	VLDDWM2 *+AR2[8],  VR35:VR34  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR5:VR4, *+AR3[35]		\n\t"
			"|  VLDDWM2 *+AR2[41],  VR37:VR36  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR7:VR6, *+AR3[51]      \n\t"
			"|  VLDDWM2 *+AR2[42],  VR39:VR38  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR9:VR8, *+AR3[67]      \n\t"
			"|  VLDDWM2 *+AR2[43],  VR41:VR40  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR11:VR10, *+AR3[83]    \n\t"
			"|  VLDDWM2 *+AR2[48],  VR43:VR42  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR13:VR12, *+AR3[99]   \n\t"
			"|  VLDDWM2 *+AR2[49],  VR45:VR44  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR15:VR14, *+AR3[115]   \n\t"
			"|  VLDDWM2 *+AR2[50],  VR47:VR46  \n\t" //vr16-24,vr17-25
			"|	SSHFLL	  3, R10, R10       \n\t"  //SSHFLL 1  

			"   VSTDW   VR17:VR16, *+AR3[131]   \n\t"
			"|  VLDDWM2 *+AR2[55],  VR49:VR48  \n\t" //vr16-24,vr17-25
			"|  SADDA     R10,AR0,AR1		\n\t" //sadda AVAF 3 SVAF 2

			"   VSTDW   VR19:VR18, *+AR3[147]   \n\t"
			"|  VLDDWM2 *+AR2[56],  VR51:VR50  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR21:VR20, *+AR3[163]   \n\t"
			"|  VLDDWM2 *+AR2[57],  VR53:VR52  \n\t" //vr16-24,vr17-25
// index[8] read
			"	VSTDW   VR23:VR22, *+AR3[179]   \n\t"
			"|  VLDDWM2 *-AR1[57],  VR1:VR0    \n\t" //vr0-0,vr1-1

			"	VSTDW   VR25:VR24, *+AR3[195]   \n\t"
			"|  VLDDWM2 *-AR1[56],  VR3:VR2    \n\t" //vr2-3,vr3-4
			
			"	VSTDW   VR27:VR26, *+AR3[211]   \n\t"
			"|  VLDDWM2 *-AR1[55],  VR5:VR4    \n\t" //vr4-6,vr5-7

			"	VSTDW   VR29:VR28, *+AR3[227]   \n\t"
			"|  VLDDWM2 *-AR1[50],  VR7:VR6    \n\t" //vr6-9,vr7-10

			"	VSTDW   VR31:VR30, *+AR3[243]   \n\t"
			"|  VLDDWM2 *-AR1[49],  VR9:VR8    \n\t" //vr4-12,vr5-13

			"   VSTDW   VR33:VR32, *+AR3[259]         \n\t"
			"|  VLDDWM2 *-AR1[48],  VR11:VR10	\n\t" //vr10-15,vr11-16

			"   VSTDW   VR35:VR34, *+AR3[275]         \n\t"
			"|  VLDDWM2 *-AR1[43],  VR13:VR12  \n\t" //vr12-18,vr13-19

			"   VSTDW   VR37:VR36, *+AR3[291]         \n\t"
			"|  VLDDWM2 *-AR1[42],  VR15:VR14  \n\t" //vr14-21,vr15-22

			"   VSTDW   VR39:VR38, *+AR3[307]         \n\t"
			"|  VLDDWM2 *-AR1[41],  VR17:VR16  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR41:VR40, *+AR3[323]         \n\t"
			"|	VLDDWM2	*-AR1[8],	 VR19:VR18	\n\t" //v18-2   

			"   VSTDW   VR43:VR42, *+AR3[339]         \n\t"
			"|  VLDDWM2 *-AR1[7],  VR21:VR20  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR45:VR44, *+AR3[355]         \n\t"
			"|  VLDDWM2 *-AR1[6],  VR23:VR22  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR47:VR46, *+AR3[371]         \n\t"
			"|  VLDDWM2 *-AR1[1],  VR25:VR24  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR49:VR48, *+AR3[387]         \n\t"
			"|  VLDDWM2 *+AR1[0],  VR27:VR26  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR51:VR50, *+AR3[403]         \n\t"
			"|  VLDDWM2 *+AR1[1],   VR29:VR28  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR53:VR52, *+AR3[419]         \n\t"
			"|  VLDDWM2 *+AR1[6],  VR31:VR30  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR1:VR0, *+AR3[4]         \n\t"
			"|  VLDDWM2 *+AR1[7],  VR33:VR32  \n\t" //vr16-24,vr17-25
			"|	SLDW	  *+AR15[10], R10	\n\t" //sldw 7

			"   VSTDW   VR3:VR2, *+AR3[20]      \n\t"
			"|	VLDDWM2 *+AR1[8],  VR35:VR34  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR5:VR4, *+AR3[36]		\n\t"
			"|  VLDDWM2 *+AR1[41],  VR37:VR36  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR7:VR6, *+AR3[52]      \n\t"
			"|  VLDDWM2 *+AR1[42],  VR39:VR38  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR9:VR8, *+AR3[68]      \n\t"
			"|  VLDDWM2 *+AR1[43],  VR41:VR40  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR11:VR10, *+AR3[84]    \n\t"
			"|  VLDDWM2 *+AR1[48],  VR43:VR42  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR13:VR12, *+AR3[100]   \n\t"
			"|  VLDDWM2 *+AR1[49],  VR45:VR44  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR15:VR14, *+AR3[116]   \n\t"
			"|  VLDDWM2 *+AR1[50],  VR47:VR46  \n\t" //vr16-24,vr17-25
			"|	SSHFLL	  3, R10, R10       \n\t"  //SSHFLL 1  

			"   VSTDW   VR17:VR16, *+AR3[132]   \n\t"
			"|  VLDDWM2 *+AR1[55],  VR49:VR48  \n\t" //vr16-24,vr17-25
			"|  SADDA     R10,AR0,AR2		\n\t" //sadda AVAF 3 SVAF 2

			"   VSTDW   VR19:VR18, *+AR3[148]   \n\t"
			"|  VLDDWM2 *+AR1[56],  VR51:VR50  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR21:VR20, *+AR3[164]   \n\t"
			"|  VLDDWM2 *+AR1[57],  VR53:VR52  \n\t" //vr16-24,vr17-25
// index[10] read
			"	VSTDW   VR23:VR22, *+AR3[180]   \n\t"
			"|  VLDDWM2 *-AR2[57],  VR1:VR0    \n\t" //vr0-0,vr1-1

			"	VSTDW   VR25:VR24, *+AR3[196]   \n\t"
			"|  VLDDWM2 *-AR2[56],  VR3:VR2    \n\t" //vr2-3,vr3-4
			
			"	VSTDW   VR27:VR26, *+AR3[212]   \n\t"
			"|  VLDDWM2 *-AR2[55],  VR5:VR4    \n\t" //vr4-6,vr5-7

			"	VSTDW   VR29:VR28, *+AR3[228]   \n\t"
			"|  VLDDWM2 *-AR2[50],  VR7:VR6    \n\t" //vr6-9,vr7-10

			"	VSTDW   VR31:VR30, *+AR3[244]   \n\t"
			"|  VLDDWM2 *-AR2[49],  VR9:VR8    \n\t" //vr4-12,vr5-13

			"   VSTDW   VR33:VR32, *+AR3[260]         \n\t"
			"|  VLDDWM2 *-AR2[48],  VR11:VR10	\n\t" //vr10-15,vr11-16

			"   VSTDW   VR35:VR34, *+AR3[276]         \n\t"
			"|  VLDDWM2 *-AR2[43],  VR13:VR12  \n\t" //vr12-18,vr13-19

			"   VSTDW   VR37:VR36, *+AR3[292]         \n\t"
			"|  VLDDWM2 *-AR2[42],  VR15:VR14  \n\t" //vr14-21,vr15-22

			"   VSTDW   VR39:VR38, *+AR3[308]         \n\t"
			"|  VLDDWM2 *-AR2[41],  VR17:VR16  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR41:VR40, *+AR3[324]         \n\t"
			"|	VLDDWM2	*-AR2[8],	 VR19:VR18	\n\t" //v18-2   

			"   VSTDW   VR43:VR42, *+AR3[340]         \n\t"
			"|  VLDDWM2 *-AR2[7],  VR21:VR20  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR45:VR44, *+AR3[356]         \n\t"
			"|  VLDDWM2 *-AR2[6],  VR23:VR22  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR47:VR46, *+AR3[372]         \n\t"
			"|  VLDDWM2 *-AR2[1],  VR25:VR24  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR49:VR48, *+AR3[388]         \n\t"
			"|  VLDDWM2 *+AR2[0],  VR27:VR26  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR51:VR50, *+AR3[404]         \n\t"
			"|  VLDDWM2 *+AR2[1],   VR29:VR28  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR53:VR52, *+AR3[420]         \n\t"
			"|  VLDDWM2 *+AR2[6],  VR31:VR30  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR1:VR0, *+AR3[5]         \n\t"
			"|  VLDDWM2 *+AR2[7],  VR33:VR32  \n\t" //vr16-24,vr17-25
			"|	SLDW	  *+AR15[12], R10	\n\t" //sldw 7

			"   VSTDW   VR3:VR2, *+AR3[21]      \n\t"
			"|	VLDDWM2 *+AR2[8],  VR35:VR34  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR5:VR4, *+AR3[37]		\n\t"
			"|  VLDDWM2 *+AR2[41],  VR37:VR36  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR7:VR6, *+AR3[53]      \n\t"
			"|  VLDDWM2 *+AR2[42],  VR39:VR38  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR9:VR8, *+AR3[69]      \n\t"
			"|  VLDDWM2 *+AR2[43],  VR41:VR40  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR11:VR10, *+AR3[85]    \n\t"
			"|  VLDDWM2 *+AR2[48],  VR43:VR42  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR13:VR12, *+AR3[101]   \n\t"
			"|  VLDDWM2 *+AR2[49],  VR45:VR44  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR15:VR14, *+AR3[117]   \n\t"
			"|  VLDDWM2 *+AR2[50],  VR47:VR46  \n\t" //vr16-24,vr17-25
			"|	SSHFLL	  3, R10, R10       \n\t"  //SSHFLL 1  

			"   VSTDW   VR17:VR16, *+AR3[133]   \n\t"
			"|  VLDDWM2 *+AR2[55],  VR49:VR48  \n\t" //vr16-24,vr17-25
			"|  SADDA     R10,AR0,AR1		\n\t" //sadda AVAF 3 SVAF 2

			"   VSTDW   VR19:VR18, *+AR3[149]   \n\t"
			"|  VLDDWM2 *+AR2[56],  VR51:VR50  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR21:VR20, *+AR3[165]   \n\t"
			"|  VLDDWM2 *+AR2[57],  VR53:VR52  \n\t" //vr16-24,vr17-25
// index[12] read
			"	VSTDW   VR23:VR22, *+AR3[181]   \n\t"
			"|  VLDDWM2 *-AR1[57],  VR1:VR0    \n\t" //vr0-0,vr1-1

			"	VSTDW   VR25:VR24, *+AR3[197]   \n\t"
			"|  VLDDWM2 *-AR1[56],  VR3:VR2    \n\t" //vr2-3,vr3-4
			
			"	VSTDW   VR27:VR26, *+AR3[213]   \n\t"
			"|  VLDDWM2 *-AR1[55],  VR5:VR4    \n\t" //vr4-6,vr5-7

			"	VSTDW   VR29:VR28, *+AR3[229]   \n\t"
			"|  VLDDWM2 *-AR1[50],  VR7:VR6    \n\t" //vr6-9,vr7-10

			"	VSTDW   VR31:VR30, *+AR3[245]   \n\t"
			"|  VLDDWM2 *-AR1[49],  VR9:VR8    \n\t" //vr4-12,vr5-13

			"   VSTDW   VR33:VR32, *+AR3[261]         \n\t"
			"|  VLDDWM2 *-AR1[48],  VR11:VR10	\n\t" //vr10-15,vr11-16

			"   VSTDW   VR35:VR34, *+AR3[277]         \n\t"
			"|  VLDDWM2 *-AR1[43],  VR13:VR12  \n\t" //vr12-18,vr13-19

			"   VSTDW   VR37:VR36, *+AR3[293]         \n\t"
			"|  VLDDWM2 *-AR1[42],  VR15:VR14  \n\t" //vr14-21,vr15-22

			"   VSTDW   VR39:VR38, *+AR3[309]         \n\t"
			"|  VLDDWM2 *-AR1[41],  VR17:VR16  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR41:VR40, *+AR3[325]         \n\t"
			"|	VLDDWM2	*-AR1[8],	 VR19:VR18	\n\t" //v18-2   

			"   VSTDW   VR43:VR42, *+AR3[341]         \n\t"
			"|  VLDDWM2 *-AR1[7],  VR21:VR20  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR45:VR44, *+AR3[357]         \n\t"
			"|  VLDDWM2 *-AR1[6],  VR23:VR22  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR47:VR46, *+AR3[373]         \n\t"
			"|  VLDDWM2 *-AR1[1],  VR25:VR24  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR49:VR48, *+AR3[389]         \n\t"
			"|  VLDDWM2 *+AR1[0],  VR27:VR26  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR51:VR50, *+AR3[405]         \n\t"
			"|  VLDDWM2 *+AR1[1],   VR29:VR28  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR53:VR52, *+AR3[421]         \n\t"
			"|  VLDDWM2 *+AR1[6],  VR31:VR30  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR1:VR0, *+AR3[6]         \n\t"
			"|  VLDDWM2 *+AR1[7],  VR33:VR32  \n\t" //vr16-24,vr17-25
			"|	SLDW	  *+AR15[14], R10	\n\t" //sldw 7

			"   VSTDW   VR3:VR2, *+AR3[22]      \n\t"
			"|	VLDDWM2 *+AR1[8],  VR35:VR34  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR5:VR4, *+AR3[38]		\n\t"
			"|  VLDDWM2 *+AR1[41],  VR37:VR36  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR7:VR6, *+AR3[54]      \n\t"
			"|  VLDDWM2 *+AR1[42],  VR39:VR38  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR9:VR8, *+AR3[70]      \n\t"
			"|  VLDDWM2 *+AR1[43],  VR41:VR40  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR11:VR10, *+AR3[86]    \n\t"
			"|  VLDDWM2 *+AR1[48],  VR43:VR42  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR13:VR12, *+AR3[102]   \n\t"
			"|  VLDDWM2 *+AR1[49],  VR45:VR44  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR15:VR14, *+AR3[118]   \n\t"
			"|  VLDDWM2 *+AR1[50],  VR47:VR46  \n\t" //vr16-24,vr17-25
			"|	SSHFLL	  3, R10, R10       \n\t"  //SSHFLL 1  

			"   VSTDW   VR17:VR16, *+AR3[134]   \n\t"
			"|  VLDDWM2 *+AR1[55],  VR49:VR48  \n\t" //vr16-24,vr17-25
			"|  SADDA     R10,AR0,AR2		\n\t" //sadda AVAF 3 SVAF 2

			"   VSTDW   VR19:VR18, *+AR3[150]   \n\t"
			"|  VLDDWM2 *+AR1[56],  VR51:VR50  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR21:VR20, *+AR3[166]   \n\t"
			"|  VLDDWM2 *+AR1[57],  VR53:VR52  \n\t" //vr16-24,vr17-25
// index[14] read
			"	VSTDW   VR23:VR22, *+AR3[182]   \n\t"
			"|  VLDDWM2 *-AR2[57],  VR1:VR0    \n\t" //vr0-0,vr1-1

			"	VSTDW   VR25:VR24, *+AR3[198]   \n\t"
			"|  VLDDWM2 *-AR2[56],  VR3:VR2    \n\t" //vr2-3,vr3-4
			
			"	VSTDW   VR27:VR26, *+AR3[214]   \n\t"
			"|  VLDDWM2 *-AR2[55],  VR5:VR4    \n\t" //vr4-6,vr5-7

			"	VSTDW   VR29:VR28, *+AR3[230]   \n\t"
			"|  VLDDWM2 *-AR2[50],  VR7:VR6    \n\t" //vr6-9,vr7-10

			"	VSTDW   VR31:VR30, *+AR3[246]   \n\t"
			"|  VLDDWM2 *-AR2[49],  VR9:VR8    \n\t" //vr4-12,vr5-13

			"   VSTDW   VR33:VR32, *+AR3[262]         \n\t"
			"|  VLDDWM2 *-AR2[48],  VR11:VR10	\n\t" //vr10-15,vr11-16

			"   VSTDW   VR35:VR34, *+AR3[278]         \n\t"
			"|  VLDDWM2 *-AR2[43],  VR13:VR12  \n\t" //vr12-18,vr13-19

			"   VSTDW   VR37:VR36, *+AR3[294]         \n\t"
			"|  VLDDWM2 *-AR2[42],  VR15:VR14  \n\t" //vr14-21,vr15-22

			"   VSTDW   VR39:VR38, *+AR3[310]         \n\t"
			"|  VLDDWM2 *-AR2[41],  VR17:VR16  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR41:VR40, *+AR3[326]         \n\t"
			"|	VLDDWM2	*-AR2[8],	 VR19:VR18	\n\t" //v18-2   

			"   VSTDW   VR43:VR42, *+AR3[342]         \n\t"
			"|  VLDDWM2 *-AR2[7],  VR21:VR20  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR45:VR44, *+AR3[358]         \n\t"
			"|  VLDDWM2 *-AR2[6],  VR23:VR22  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR47:VR46, *+AR3[374]         \n\t"
			"|  VLDDWM2 *-AR2[1],  VR25:VR24  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR49:VR48, *+AR3[390]         \n\t"
			"|  VLDDWM2 *+AR2[0],  VR27:VR26  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR51:VR50, *+AR3[406]         \n\t"
			"|  VLDDWM2 *+AR2[1],   VR29:VR28  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR53:VR52, *+AR3[422]         \n\t"
			"|  VLDDWM2 *+AR2[6],  VR31:VR30  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR1:VR0, *+AR3[7]         \n\t"
			"|  VLDDWM2 *+AR2[7],  VR33:VR32  \n\t" //vr16-24,vr17-25
			"|	SLDW	  *+AR15[1], R10	\n\t" //sldw 7

			"   VSTDW   VR3:VR2, *+AR3[23]      \n\t"
			"|	VLDDWM2 *+AR2[8],  VR35:VR34  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR5:VR4, *+AR3[39]		\n\t"
			"|  VLDDWM2 *+AR2[41],  VR37:VR36  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR7:VR6, *+AR3[55]      \n\t"
			"|  VLDDWM2 *+AR2[42],  VR39:VR38  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR9:VR8, *+AR3[71]      \n\t"
			"|  VLDDWM2 *+AR2[43],  VR41:VR40  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR11:VR10, *+AR3[87]    \n\t"
			"|  VLDDWM2 *+AR2[48],  VR43:VR42  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR13:VR12, *+AR3[103]   \n\t"
			"|  VLDDWM2 *+AR2[49],  VR45:VR44  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR15:VR14, *+AR3[119]   \n\t"
			"|  VLDDWM2 *+AR2[50],  VR47:VR46  \n\t" //vr16-24,vr17-25
			"|	SSHFLL	  3, R10, R10       \n\t"  //SSHFLL 1  

			"   VSTDW   VR17:VR16, *+AR3[135]   \n\t"
			"|  VLDDWM2 *+AR2[55],  VR49:VR48  \n\t" //vr16-24,vr17-25
			"|  SADDA     R10,AR0,AR1		\n\t" //sadda AVAF 3 SVAF 2

			"   VSTDW   VR19:VR18, *+AR3[151]   \n\t"
			"|  VLDDWM2 *+AR2[56],  VR51:VR50  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR21:VR20, *+AR3[167]   \n\t"
			"|  VLDDWM2 *+AR2[57],  VR53:VR52  \n\t" //vr16-24,vr17-25

			"	SNOP	1\n\t"

			"	VSTDW   VR23:VR22, *+AR3[183]   \n\t"
			"|	VSTDW   VR25:VR24, *+AR3[199]   \n\t"

			"	VSTDW   VR27:VR26, *+AR3[215]   \n\t"
			"|	VSTDW   VR29:VR28, *+AR3[231]   \n\t"

			"	VSTDW   VR31:VR30, *+AR3[247]   \n\t"
			"|  VSTDW   VR33:VR32, *+AR3[263]         \n\t"

			"   VSTDW   VR35:VR34, *+AR3[279]         \n\t"
			"|  VSTDW   VR37:VR36, *+AR3[295]         \n\t"

			"   VSTDW   VR39:VR38, *+AR3[311]         \n\t"
			"|  VSTDW   VR41:VR40, *+AR3[327]         \n\t"

			"   VSTDW   VR43:VR42, *+AR3[343]         \n\t"
			"|  VSTDW   VR45:VR44, *+AR3[359]         \n\t"

			"   VSTDW   VR47:VR46, *+AR3[375]         \n\t"
			"|  VSTDW   VR49:VR48, *+AR3[391]         \n\t"

			"   VSTDW   VR51:VR50, *+AR3[407]         \n\t"
			"|  VSTDW   VR53:VR52, *+AR3[423]         \n\t"

//			"|  SMVAGA.M2 R31,AR3           \n\t" //R31 
//index[1] read
			"   VLDDWM2 *-AR1[57],  VR1:VR0    \n\t" //vr0-0,vr1-1
			"|  VLDDWM2 *-AR1[56],  VR3:VR2    \n\t" //vr2-3,vr3-4

			"   VLDDWM2 *-AR1[55],  VR5:VR4    \n\t" //vr4-6,vr5-7
			"|  VLDDWM2 *-AR1[50],  VR7:VR6    \n\t" //vr6-9,vr7-10

			"   VLDDWM2 *-AR1[49],  VR9:VR8    \n\t" //vr4-12,vr5-13
			"|  VLDDWM2 *-AR1[48],  VR11:VR10	\n\t" //vr10-15,vr11-16

			"   VLDDWM2 *-AR1[43],  VR13:VR12  \n\t" //vr12-18,vr13-19
			"|  VLDDWM2 *-AR1[42],  VR15:VR14  \n\t" //vr14-21,vr15-22

			"   VLDDWM2 *-AR1[41],  VR17:VR16  \n\t" //vr16-24,vr17-25
			"|	VLDDWM2	*-AR1[8],	 VR19:VR18	\n\t" //v18-2   

			"   VLDDWM2 *-AR1[7],  VR21:VR20  \n\t" //vr16-24,vr17-25
			"|  VLDDWM2 *-AR1[6],  VR23:VR22  \n\t" //vr16-24,vr17-25

			"   VLDDWM2 *-AR1[1],  VR25:VR24  \n\t" //vr16-24,vr17-25
			"|  VLDDWM2 *+AR1[0],  VR27:VR26  \n\t" //vr16-24,vr17-25

			"   VLDDWM2 *+AR1[1],   VR29:VR28  \n\t" //vr16-24,vr17-25
			"|  VLDDWM2 *+AR1[6],  VR31:VR30  \n\t" //vr16-24,vr17-25

			"	SNOP	1\n\t"

			"   VSTW   VR0, *++AR3[1]         \n\t"  
			"|  VLDDWM2 *+AR1[7],  VR33:VR32  \n\t" //vr16-24,vr17-25
			"|	SLDW	  *+AR15[3], R10	\n\t" //sldw 7

			"   VSTW   VR1, *++AR3[16]      \n\t"
			"|  VLDDWM2 *+AR1[8],  VR35:VR34  \n\t" //vr16-24,vr17-25

			"   VSTW   VR2, *++AR3[16]		\n\t"
			"|  VLDDWM2 *+AR1[41],  VR37:VR36  \n\t" //vr16-24,vr17-25

			"   VSTW   VR3, *++AR3[16]      \n\t"
			"|  VLDDWM2 *+AR1[42],  VR39:VR38  \n\t" //vr16-24,vr17-25

			"   VSTW   VR4, *++AR3[16]      \n\t"
			"|  VLDDWM2 *+AR1[43],  VR41:VR40  \n\t" //vr16-24,vr17-25

			"   VSTW   VR5, *++AR3[16]    \n\t"
			"|  VLDDWM2 *+AR1[48],  VR43:VR42  \n\t" //vr16-24,vr17-25

			"   VSTW   VR6, *++AR3[16]   \n\t"
			"|  VLDDWM2 *+AR1[49],  VR45:VR44  \n\t" //vr16-24,vr17-25

			"   VSTW   VR7, *++AR3[16]   \n\t"
			"|  VLDDWM2 *+AR1[50],  VR47:VR46  \n\t" //vr16-24,vr17-25
			"|	SSHFLL	  3, R10, R10       \n\t"  //SSHFLL 1  

			"   VSTW   VR8, *++AR3[16]   \n\t"
			"|  VLDDWM2 *+AR1[55],  VR49:VR48  \n\t" //vr16-24,vr17-25
			"|  SADDA     R10,AR0,AR2		\n\t" //sadda AVAF 3 SVAF 2

			"   VSTW   VR9, *++AR3[16]   \n\t"
			"|  VLDDWM2 *+AR1[56],  VR51:VR50  \n\t" //vr16-24,vr17-25

			"   VSTW   VR10, *++AR3[16]   \n\t"
			"|  VLDDWM2 *+AR1[57],  VR53:VR52  \n\t" //vr16-24,vr17-25

			"   VSTW   VR11, *+AR3[16]   \n\t"
			"|  VSTW   VR12, *+AR3[32]   \n\t"

			"   VSTW   VR13, *+AR3[48]   \n\t"
			"|  VSTW   VR14, *+AR3[64]   \n\t"

			"   VSTW   VR15, *+AR3[80]   \n\t"
			"|  VSTW   VR16, *+AR3[96]   \n\t"

			"   VSTW   VR17, *+AR3[112]   \n\t"
			"|  VSTW   VR18, *+AR3[128]   \n\t"

			"   VSTW   VR19, *+AR3[144]   \n\t"
			"|  VSTW   VR20, *+AR3[160]   \n\t"

			"   VSTW   VR21, *+AR3[176]   \n\t"
			"|  VSTW   VR22, *+AR3[192]   \n\t"

			"   VSTW   VR23, *+AR3[208]   \n\t"
			"|  VSTW   VR24, *+AR3[224]   \n\t"

			"   VSTW   VR25, *+AR3[240]   \n\t"
			"|  VSTW   VR26, *+AR3[256]   \n\t"

			"   VSTW   VR27, *+AR3[272]   \n\t"
			"|  VSTW   VR28, *+AR3[288]   \n\t"
			"|  SADDA     R10,AR0,AR1		\n\t" //sadda AVAF 3 SVAF 2

			"   VSTW   VR29, *+AR3[304]   \n\t"
			"|  VSTW   VR30, *+AR3[320]   \n\t"

			"   VSTW   VR31, *+AR3[336]   \n\t"
			"|  VSTW   VR32, *+AR3[352]   \n\t"

			"   VSTW   VR33, *+AR3[368]   \n\t"
			"|  VSTW   VR34, *+AR3[384]   \n\t"

			"   VSTW   VR35, *+AR3[400]   \n\t"
			"|  VSTW   VR36, *+AR3[416]   \n\t"

			"   VSTW   VR37, *++AR3[432]   \n\t"
			"|  SMVAGA.M2 R31,AR4           \n\t" //R31 ////////////////////////////////////////////////////////////////////////////
// index[3] read
			"	VSTW   VR38, *++AR3[16]   \n\t"
			"|  VLDDWM2 *-AR1[57],  VR1:VR0    \n\t" //vr0-0,vr1-1

			"	VSTW   VR39, *++AR3[16]   \n\t"
			"|  VLDDWM2 *-AR1[56],  VR3:VR2    \n\t" //vr2-3,vr3-4
			
			"	VSTW   VR40, *++AR3[16]   \n\t"
			"|  VLDDWM2 *-AR1[55],  VR5:VR4    \n\t" //vr4-6,vr5-7

			"	VSTW   VR41, *++AR3[16]   \n\t"
			"|  VLDDWM2 *-AR1[50],  VR7:VR6    \n\t" //vr6-9,vr7-10

			"	VSTW   VR42, *++AR3[16]   \n\t"
			"|  VLDDWM2 *-AR1[49],  VR9:VR8    \n\t" //vr4-12,vr5-13

			"   VSTW   VR43, *++AR3[16]         \n\t"
			"|  VLDDWM2 *-AR1[48],  VR11:VR10	\n\t" //vr10-15,vr11-16

			"   VSTW   VR44, *++AR3[16]         \n\t"
			"|  VLDDWM2 *-AR1[43],  VR13:VR12  \n\t" //vr12-18,vr13-19

			"   VSTW   VR45, *++AR3[16]         \n\t"
			"|  VLDDWM2 *-AR1[42],  VR15:VR14  \n\t" //vr14-21,vr15-22

			"   VSTW   VR46, *++AR3[16]         \n\t"
			"|  VLDDWM2 *-AR1[41],  VR17:VR16  \n\t" //vr16-24,vr17-25

			"   VSTW   VR47, *++AR3[16]         \n\t"
			"|	VLDDWM2	*-AR1[8],	 VR19:VR18	\n\t" //v18-2   

			"   VSTW   VR48, *++AR3[16]         \n\t"
			"|  VLDDWM2 *-AR1[7],  VR21:VR20  \n\t" //vr16-24,vr17-25

			"   VSTW   VR49, *++AR3[16]         \n\t"
			"|  VLDDWM2 *-AR1[6],  VR23:VR22  \n\t" //vr16-24,vr17-25

			"   VSTW   VR50, *++AR3[16]         \n\t"
			"|  VLDDWM2 *-AR1[1],  VR25:VR24  \n\t" //vr16-24,vr17-25

			"   VSTW   VR51, *++AR3[16]         \n\t"
			"|  VLDDWM2 *+AR1[0],  VR27:VR26  \n\t" //vr16-24,vr17-25

			"   VSTW   VR52, *++AR3[16]         \n\t"
			"|  VLDDWM2 *+AR1[1],   VR29:VR28  \n\t" //vr16-24,vr17-25

			"   VSTW   VR53, *++AR3[16]         \n\t"
			"|  VLDDWM2 *+AR1[6],  VR31:VR30  \n\t" //vr16-24,vr17-25

			"   VSTW   VR0, *++AR4[3]         \n\t"
			"|  VLDDWM2 *+AR1[7],  VR33:VR32  \n\t" //vr16-24,vr17-25
			"|	SLDW	  *+AR15[5], R10	\n\t" //sldw 7

			"   VSTW   VR1, *++AR4[16]      \n\t"
			"|	VLDDWM2 *+AR1[8],  VR35:VR34  \n\t" //vr16-24,vr17-25

			"   VSTW   VR2, *++AR4[16]		\n\t"
			"|  VLDDWM2 *+AR1[41],  VR37:VR36  \n\t" //vr16-24,vr17-25

			"   VSTW   VR3, *++AR4[16]      \n\t"
			"|  VLDDWM2 *+AR1[42],  VR39:VR38  \n\t" //vr16-24,vr17-25

			"   VSTW   VR4, *++AR4[16]      \n\t"
			"|  VLDDWM2 *+AR1[43],  VR41:VR40  \n\t" //vr16-24,vr17-25

			"   VSTW   VR5, *++AR4[16]    \n\t"
			"|  VLDDWM2 *+AR1[48],  VR43:VR42  \n\t" //vr16-24,vr17-25

			"   VSTW   VR6, *++AR4[16]   \n\t"
			"|  VLDDWM2 *+AR1[49],  VR45:VR44  \n\t" //vr16-24,vr17-25

			"   VSTW   VR7, *++AR4[16]   \n\t"
			"|  VLDDWM2 *+AR1[50],  VR47:VR46  \n\t" //vr16-24,vr17-25
			"|	SSHFLL	  3, R10, R10       \n\t"  //SSHFLL 1  

			"   VSTW   VR8, *++AR4[16]   \n\t"
			"|  VLDDWM2 *+AR1[55],  VR49:VR48  \n\t" //vr16-24,vr17-25

			"   VSTW   VR9, *++AR4[16]   \n\t"
			"|  VLDDWM2 *+AR1[56],  VR51:VR50  \n\t" //vr16-24,vr17-25

			"   VSTW   VR10, *++AR4[16]   \n\t"
			"|  VLDDWM2 *+AR1[57],  VR53:VR52  \n\t" //vr16-24,vr17-25

			"   VSTW   VR11, *+AR4[16]   \n\t"
			"|  VSTW   VR12, *+AR4[32]   \n\t"

			"   VSTW   VR13, *+AR4[48]   \n\t"
			"|  VSTW   VR14, *+AR4[64]   \n\t"

			"   VSTW   VR15, *+AR4[80]   \n\t"
			"|  VSTW   VR16, *+AR4[96]   \n\t"

			"   VSTW   VR17, *+AR4[112]   \n\t"
			"|  VSTW   VR18, *+AR4[128]   \n\t"

			"   VSTW   VR19, *+AR4[144]   \n\t"
			"|  VSTW   VR20, *+AR4[160]   \n\t"

			"   VSTW   VR21, *+AR4[176]   \n\t"
			"|  VSTW   VR22, *+AR4[192]   \n\t"

			"   VSTW   VR23, *+AR4[208]   \n\t"
			"|  VSTW   VR24, *+AR4[224]   \n\t"

			"   VSTW   VR25, *+AR4[240]   \n\t"
			"|  VSTW   VR26, *+AR4[256]   \n\t"

			"   VSTW   VR27, *+AR4[272]   \n\t"
			"|  VSTW   VR28, *+AR4[288]   \n\t"
			"|  SADDA     R10,AR0,AR1		\n\t" //sadda AVAF 3 SVAF 2

			"   VSTW   VR29, *+AR4[304]   \n\t"
			"|  VSTW   VR30, *+AR4[320]   \n\t"

			"   VSTW   VR31, *+AR4[336]   \n\t"
			"|  VSTW   VR32, *+AR4[352]   \n\t"

			"   VSTW   VR33, *+AR4[368]   \n\t"
			"|  VSTW   VR34, *+AR4[384]   \n\t"

			"   VSTW   VR35, *+AR4[400]   \n\t"
			"|  VSTW   VR36, *+AR4[416]   \n\t"

			"   VSTW   VR37, *++AR4[432]   \n\t"
			"|  SMVAGA.M2 R31,AR3           \n\t" //R31 // index[5] read
			"	VSTW   VR38, *++AR4[16]   \n\t"
			"|  VLDDWM2 *-AR1[57],  VR1:VR0    \n\t" //vr0-0,vr1-1

			"	VSTW   VR39, *++AR4[16]   \n\t"
			"|  VLDDWM2 *-AR1[56],  VR3:VR2    \n\t" //vr2-3,vr3-4
			
			"	VSTW   VR40, *++AR4[16]   \n\t"
			"|  VLDDWM2 *-AR1[55],  VR5:VR4    \n\t" //vr4-6,vr5-7

			"	VSTW   VR41, *++AR4[16]   \n\t"
			"|  VLDDWM2 *-AR1[50],  VR7:VR6    \n\t" //vr6-9,vr7-10

			"	VSTW   VR42, *++AR4[16]   \n\t"
			"|  VLDDWM2 *-AR1[49],  VR9:VR8    \n\t" //vr4-12,vr5-13

			"   VSTW   VR43, *++AR4[16]         \n\t"
			"|  VLDDWM2 *-AR1[48],  VR11:VR10	\n\t" //vr10-15,vr11-16

			"   VSTW   VR44, *++AR4[16]         \n\t"
			"|  VLDDWM2 *-AR1[43],  VR13:VR12  \n\t" //vr12-18,vr13-19

			"   VSTW   VR45, *++AR4[16]         \n\t"
			"|  VLDDWM2 *-AR1[42],  VR15:VR14  \n\t" //vr14-21,vr15-22

			"   VSTW   VR46, *++AR4[16]         \n\t"
			"|  VLDDWM2 *-AR1[41],  VR17:VR16  \n\t" //vr16-24,vr17-25

			"   VSTW   VR47, *++AR4[16]         \n\t"
			"|	VLDDWM2	*-AR1[8],	 VR19:VR18	\n\t" //v18-2   

			"   VSTW   VR48, *++AR4[16]         \n\t"
			"|  VLDDWM2 *-AR1[7],  VR21:VR20  \n\t" //vr16-24,vr17-25

			"   VSTW   VR49, *++AR4[16]         \n\t"
			"|  VLDDWM2 *-AR1[6],  VR23:VR22  \n\t" //vr16-24,vr17-25

			"   VSTW   VR50, *++AR4[16]         \n\t"
			"|  VLDDWM2 *-AR1[1],  VR25:VR24  \n\t" //vr16-24,vr17-25

			"   VSTW   VR51, *++AR4[16]         \n\t"
			"|  VLDDWM2 *+AR1[0],  VR27:VR26  \n\t" //vr16-24,vr17-25

			"   VSTW   VR52, *++AR4[16]         \n\t"
			"|  VLDDWM2 *+AR1[1],   VR29:VR28  \n\t" //vr16-24,vr17-25

			"   VSTW   VR53, *++AR4[16]         \n\t"
			"|  VLDDWM2 *+AR1[6],  VR31:VR30  \n\t" //vr16-24,vr17-25

			"   VSTW   VR0, *++AR3[5]         \n\t"
			"|  VLDDWM2 *+AR1[7],  VR33:VR32  \n\t" //vr16-24,vr17-25
			"|	SLDW	  *+AR15[7], R10	\n\t" //sldw 7

			"   VSTW   VR1, *++AR3[16]      \n\t"
			"|	VLDDWM2 *+AR1[8],  VR35:VR34  \n\t" //vr16-24,vr17-25

			"   VSTW   VR2, *++AR3[16]		\n\t"
			"|  VLDDWM2 *+AR1[41],  VR37:VR36  \n\t" //vr16-24,vr17-25

			"   VSTW   VR3, *++AR3[16]      \n\t"
			"|  VLDDWM2 *+AR1[42],  VR39:VR38  \n\t" //vr16-24,vr17-25

			"   VSTW   VR4, *++AR3[16]      \n\t"
			"|  VLDDWM2 *+AR1[43],  VR41:VR40  \n\t" //vr16-24,vr17-25

			"   VSTW   VR5, *++AR3[16]    \n\t"
			"|  VLDDWM2 *+AR1[48],  VR43:VR42  \n\t" //vr16-24,vr17-25

			"   VSTW   VR6, *++AR3[16]   \n\t"
			"|  VLDDWM2 *+AR1[49],  VR45:VR44  \n\t" //vr16-24,vr17-25

			"   VSTW   VR7, *++AR3[16]   \n\t"
			"|  VLDDWM2 *+AR1[50],  VR47:VR46  \n\t" //vr16-24,vr17-25
			"|	SSHFLL	  3, R10, R10       \n\t"  //SSHFLL 1  

			"   VSTW   VR8, *++AR3[16]   \n\t"
			"|  VLDDWM2 *+AR1[55],  VR49:VR48  \n\t" //vr16-24,vr17-25

			"   VSTW   VR9, *++AR3[16]   \n\t"
			"|  VLDDWM2 *+AR1[56],  VR51:VR50  \n\t" //vr16-24,vr17-25

			"   VSTW   VR10, *++AR3[16]   \n\t"
			"|  VLDDWM2 *+AR1[57],  VR53:VR52  \n\t" //vr16-24,vr17-25

			"   VSTW   VR11, *+AR3[16]   \n\t"
			"|  VSTW   VR12, *+AR3[32]   \n\t"

			"   VSTW   VR13, *+AR3[48]   \n\t"
			"|  VSTW   VR14, *+AR3[64]   \n\t"

			"   VSTW   VR15, *+AR3[80]   \n\t"
			"|  VSTW   VR16, *+AR3[96]   \n\t"

			"   VSTW   VR17, *+AR3[112]   \n\t"
			"|  VSTW   VR18, *+AR3[128]   \n\t"

			"   VSTW   VR19, *+AR3[144]   \n\t"
			"|  VSTW   VR20, *+AR3[160]   \n\t"

			"   VSTW   VR21, *+AR3[176]   \n\t"
			"|  VSTW   VR22, *+AR3[192]   \n\t"

			"   VSTW   VR23, *+AR3[208]   \n\t"
			"|  VSTW   VR24, *+AR3[224]   \n\t"

			"   VSTW   VR25, *+AR3[240]   \n\t"
			"|  VSTW   VR26, *+AR3[256]   \n\t"

			"   VSTW   VR27, *+AR3[272]   \n\t"
			"|  VSTW   VR28, *+AR3[288]   \n\t"
			"|  SADDA     R10,AR0,AR1		\n\t" //sadda AVAF 3 SVAF 2

			"   VSTW   VR29, *+AR3[304]   \n\t"
			"|  VSTW   VR30, *+AR3[320]   \n\t"

			"   VSTW   VR31, *+AR3[336]   \n\t"
			"|  VSTW   VR32, *+AR3[352]   \n\t"

			"   VSTW   VR33, *+AR3[368]   \n\t"
			"|  VSTW   VR34, *+AR3[384]   \n\t"

			"   VSTW   VR35, *+AR3[400]   \n\t"
			"|  VSTW   VR36, *+AR3[416]   \n\t"

			"   VSTW   VR37, *++AR3[432]   \n\t"
			"|  SMVAGA.M2 R31,AR4           \n\t" //R31 
// index[7] read
			"	VSTW   VR38, *++AR3[16]   \n\t"
			"|  VLDDWM2 *-AR1[57],  VR1:VR0    \n\t" //vr0-0,vr1-1

			"	VSTW   VR39, *++AR3[16]   \n\t"
			"|  VLDDWM2 *-AR1[56],  VR3:VR2    \n\t" //vr2-3,vr3-4
			
			"	VSTW   VR40, *++AR3[16]   \n\t"
			"|  VLDDWM2 *-AR1[55],  VR5:VR4    \n\t" //vr4-6,vr5-7

			"	VSTW   VR41, *++AR3[16]   \n\t"
			"|  VLDDWM2 *-AR1[50],  VR7:VR6    \n\t" //vr6-9,vr7-10

			"	VSTW   VR42, *++AR3[16]   \n\t"
			"|  VLDDWM2 *-AR1[49],  VR9:VR8    \n\t" //vr4-12,vr5-13

			"   VSTW   VR43, *++AR3[16]         \n\t"
			"|  VLDDWM2 *-AR1[48],  VR11:VR10	\n\t" //vr10-15,vr11-16

			"   VSTW   VR44, *++AR3[16]         \n\t"
			"|  VLDDWM2 *-AR1[43],  VR13:VR12  \n\t" //vr12-18,vr13-19

			"   VSTW   VR45, *++AR3[16]         \n\t"
			"|  VLDDWM2 *-AR1[42],  VR15:VR14  \n\t" //vr14-21,vr15-22

			"   VSTW   VR46, *++AR3[16]         \n\t"
			"|  VLDDWM2 *-AR1[41],  VR17:VR16  \n\t" //vr16-24,vr17-25

			"   VSTW   VR47, *++AR3[16]         \n\t"
			"|	VLDDWM2	*-AR1[8],	 VR19:VR18	\n\t" //v18-2   

			"   VSTW   VR48, *++AR3[16]         \n\t"
			"|  VLDDWM2 *-AR1[7],  VR21:VR20  \n\t" //vr16-24,vr17-25

			"   VSTW   VR49, *++AR3[16]         \n\t"
			"|  VLDDWM2 *-AR1[6],  VR23:VR22  \n\t" //vr16-24,vr17-25

			"   VSTW   VR50, *++AR3[16]         \n\t"
			"|  VLDDWM2 *-AR1[1],  VR25:VR24  \n\t" //vr16-24,vr17-25

			"   VSTW   VR51, *++AR3[16]         \n\t"
			"|  VLDDWM2 *+AR1[0],  VR27:VR26  \n\t" //vr16-24,vr17-25

			"   VSTW   VR52, *++AR3[16]         \n\t"
			"|  VLDDWM2 *+AR1[1],   VR29:VR28  \n\t" //vr16-24,vr17-25

			"   VSTW   VR53, *++AR3[16]         \n\t"
			"|  VLDDWM2 *+AR1[6],  VR31:VR30  \n\t" //vr16-24,vr17-25

			"   VSTW   VR0, *++AR4[7]         \n\t"
			"|  VLDDWM2 *+AR1[7],  VR33:VR32  \n\t" //vr16-24,vr17-25
			"|	SLDW	  *+AR15[9], R10	\n\t" //sldw 7

			"   VSTW   VR1, *++AR4[16]      \n\t"
			"|	VLDDWM2 *+AR1[8],  VR35:VR34  \n\t" //vr16-24,vr17-25

			"   VSTW   VR2, *++AR4[16]		\n\t"
			"|  VLDDWM2 *+AR1[41],  VR37:VR36  \n\t" //vr16-24,vr17-25

			"   VSTW   VR3, *++AR4[16]      \n\t"
			"|  VLDDWM2 *+AR1[42],  VR39:VR38  \n\t" //vr16-24,vr17-25

			"   VSTW   VR4, *++AR4[16]      \n\t"
			"|  VLDDWM2 *+AR1[43],  VR41:VR40  \n\t" //vr16-24,vr17-25

			"   VSTW   VR5, *++AR4[16]    \n\t"
			"|  VLDDWM2 *+AR1[48],  VR43:VR42  \n\t" //vr16-24,vr17-25

			"   VSTW   VR6, *++AR4[16]   \n\t"
			"|  VLDDWM2 *+AR1[49],  VR45:VR44  \n\t" //vr16-24,vr17-25

			"   VSTW   VR7, *++AR4[16]   \n\t"
			"|  VLDDWM2 *+AR1[50],  VR47:VR46  \n\t" //vr16-24,vr17-25
			"|	SSHFLL	  3, R10, R10       \n\t"  //SSHFLL 1  

			"   VSTW   VR8, *++AR4[16]   \n\t"
			"|  VLDDWM2 *+AR1[55],  VR49:VR48  \n\t" //vr16-24,vr17-25

			"   VSTW   VR9, *++AR4[16]   \n\t"
			"|  VLDDWM2 *+AR1[56],  VR51:VR50  \n\t" //vr16-24,vr17-25

			"   VSTW   VR10, *++AR4[16]   \n\t"
			"|  VLDDWM2 *+AR1[57],  VR53:VR52  \n\t" //vr16-24,vr17-25

			"   VSTW   VR11, *+AR4[16]   \n\t"
			"|  VSTW   VR12, *+AR4[32]   \n\t"

			"   VSTW   VR13, *+AR4[48]   \n\t"
			"|  VSTW   VR14, *+AR4[64]   \n\t"

			"   VSTW   VR15, *+AR4[80]   \n\t"
			"|  VSTW   VR16, *+AR4[96]   \n\t"

			"   VSTW   VR17, *+AR4[112]   \n\t"
			"|  VSTW   VR18, *+AR4[128]   \n\t"

			"   VSTW   VR19, *+AR4[144]   \n\t"
			"|  VSTW   VR20, *+AR4[160]   \n\t"

			"   VSTW   VR21, *+AR4[176]   \n\t"
			"|  VSTW   VR22, *+AR4[192]   \n\t"

			"   VSTW   VR23, *+AR4[208]   \n\t"
			"|  VSTW   VR24, *+AR4[224]   \n\t"

			"   VSTW   VR25, *+AR4[240]   \n\t"
			"|  VSTW   VR26, *+AR4[256]   \n\t"

			"   VSTW   VR27, *+AR4[272]   \n\t"
			"|  VSTW   VR28, *+AR4[288]   \n\t"
			"|  SADDA     R10,AR0,AR1		\n\t" //sadda AVAF 3 SVAF 2

			"   VSTW   VR29, *+AR4[304]   \n\t"
			"|  VSTW   VR30, *+AR4[320]   \n\t"

			"   VSTW   VR31, *+AR4[336]   \n\t"
			"|  VSTW   VR32, *+AR4[352]   \n\t"

			"   VSTW   VR33, *+AR4[368]   \n\t"
			"|  VSTW   VR34, *+AR4[384]   \n\t"

			"   VSTW   VR35, *+AR4[400]   \n\t"
			"|  VSTW   VR36, *+AR4[416]   \n\t"

			"   VSTW   VR37, *++AR4[432]   \n\t"
			"|  SMVAGA.M2 R31,AR3           \n\t" //R31 
// index[9] read
			"	VSTW   VR38, *++AR4[16]   \n\t"
			"|  VLDDWM2 *-AR1[57],  VR1:VR0    \n\t" //vr0-0,vr1-1

			"	VSTW   VR39, *++AR4[16]   \n\t"
			"|  VLDDWM2 *-AR1[56],  VR3:VR2    \n\t" //vr2-3,vr3-4
			
			"	VSTW   VR40, *++AR4[16]   \n\t"
			"|  VLDDWM2 *-AR1[55],  VR5:VR4    \n\t" //vr4-6,vr5-7

			"	VSTW   VR41, *++AR4[16]   \n\t"
			"|  VLDDWM2 *-AR1[50],  VR7:VR6    \n\t" //vr6-9,vr7-10

			"	VSTW   VR42, *++AR4[16]   \n\t"
			"|  VLDDWM2 *-AR1[49],  VR9:VR8    \n\t" //vr4-12,vr5-13

			"   VSTW   VR43, *++AR4[16]         \n\t"
			"|  VLDDWM2 *-AR1[48],  VR11:VR10	\n\t" //vr10-15,vr11-16

			"   VSTW   VR44, *++AR4[16]         \n\t"
			"|  VLDDWM2 *-AR1[43],  VR13:VR12  \n\t" //vr12-18,vr13-19

			"   VSTW   VR45, *++AR4[16]         \n\t"
			"|  VLDDWM2 *-AR1[42],  VR15:VR14  \n\t" //vr14-21,vr15-22

			"   VSTW   VR46, *++AR4[16]         \n\t"
			"|  VLDDWM2 *-AR1[41],  VR17:VR16  \n\t" //vr16-24,vr17-25

			"   VSTW   VR47, *++AR4[16]         \n\t"
			"|	VLDDWM2	*-AR1[8],	 VR19:VR18	\n\t" //v18-2   

			"   VSTW   VR48, *++AR4[16]         \n\t"
			"|  VLDDWM2 *-AR1[7],  VR21:VR20  \n\t" //vr16-24,vr17-25

			"   VSTW   VR49, *++AR4[16]         \n\t"
			"|  VLDDWM2 *-AR1[6],  VR23:VR22  \n\t" //vr16-24,vr17-25

			"   VSTW   VR50, *++AR4[16]         \n\t"
			"|  VLDDWM2 *-AR1[1],  VR25:VR24  \n\t" //vr16-24,vr17-25

			"   VSTW   VR51, *++AR4[16]         \n\t"
			"|  VLDDWM2 *+AR1[0],  VR27:VR26  \n\t" //vr16-24,vr17-25

			"   VSTW   VR52, *++AR4[16]         \n\t"
			"|  VLDDWM2 *+AR1[1],   VR29:VR28  \n\t" //vr16-24,vr17-25

			"   VSTW   VR53, *++AR4[16]         \n\t"
			"|  VLDDWM2 *+AR1[6],  VR31:VR30  \n\t" //vr16-24,vr17-25

			"   VSTW   VR0, *++AR3[9]         \n\t"
			"|  VLDDWM2 *+AR1[7],  VR33:VR32  \n\t" //vr16-24,vr17-25
			"|	SLDW	  *+AR15[11], R10	\n\t" //sldw 7

			"   VSTW   VR1, *++AR3[16]      \n\t"
			"|	VLDDWM2 *+AR1[8],  VR35:VR34  \n\t" //vr16-24,vr17-25

			"   VSTW   VR2, *++AR3[16]		\n\t"
			"|  VLDDWM2 *+AR1[41],  VR37:VR36  \n\t" //vr16-24,vr17-25

			"   VSTW   VR3, *++AR3[16]      \n\t"
			"|  VLDDWM2 *+AR1[42],  VR39:VR38  \n\t" //vr16-24,vr17-25

			"   VSTW   VR4, *++AR3[16]      \n\t"
			"|  VLDDWM2 *+AR1[43],  VR41:VR40  \n\t" //vr16-24,vr17-25

			"   VSTW   VR5, *++AR3[16]    \n\t"
			"|  VLDDWM2 *+AR1[48],  VR43:VR42  \n\t" //vr16-24,vr17-25

			"   VSTW   VR6, *++AR3[16]   \n\t"
			"|  VLDDWM2 *+AR1[49],  VR45:VR44  \n\t" //vr16-24,vr17-25

			"   VSTW   VR7, *++AR3[16]   \n\t"
			"|  VLDDWM2 *+AR1[50],  VR47:VR46  \n\t" //vr16-24,vr17-25
			"|	SSHFLL	  3, R10, R10       \n\t"  //SSHFLL 1  

			"   VSTW   VR8, *++AR3[16]   \n\t"
			"|  VLDDWM2 *+AR1[55],  VR49:VR48  \n\t" //vr16-24,vr17-25

			"   VSTW   VR9, *++AR3[16]   \n\t"
			"|  VLDDWM2 *+AR1[56],  VR51:VR50  \n\t" //vr16-24,vr17-25

			"   VSTW   VR10, *++AR3[16]   \n\t"
			"|  VLDDWM2 *+AR1[57],  VR53:VR52  \n\t" //vr16-24,vr17-25

			"   VSTW   VR11, *+AR3[16]   \n\t"
			"|  VSTW   VR12, *+AR3[32]   \n\t"

			"   VSTW   VR13, *+AR3[48]   \n\t"
			"|  VSTW   VR14, *+AR3[64]   \n\t"

			"   VSTW   VR15, *+AR3[80]   \n\t"
			"|  VSTW   VR16, *+AR3[96]   \n\t"

			"   VSTW   VR17, *+AR3[112]   \n\t"
			"|  VSTW   VR18, *+AR3[128]   \n\t"

			"   VSTW   VR19, *+AR3[144]   \n\t"
			"|  VSTW   VR20, *+AR3[160]   \n\t"

			"   VSTW   VR21, *+AR3[176]   \n\t"
			"|  VSTW   VR22, *+AR3[192]   \n\t"

			"   VSTW   VR23, *+AR3[208]   \n\t"
			"|  VSTW   VR24, *+AR3[224]   \n\t"

			"   VSTW   VR25, *+AR3[240]   \n\t"
			"|  VSTW   VR26, *+AR3[256]   \n\t"

			"   VSTW   VR27, *+AR3[272]   \n\t"
			"|  VSTW   VR28, *+AR3[288]   \n\t"
			"|  SADDA     R10,AR0,AR1		\n\t" //sadda AVAF 3 SVAF 2

			"   VSTW   VR29, *+AR3[304]   \n\t"
			"|  VSTW   VR30, *+AR3[320]   \n\t"

			"   VSTW   VR31, *+AR3[336]   \n\t"
			"|  VSTW   VR32, *+AR3[352]   \n\t"

			"   VSTW   VR33, *+AR3[368]   \n\t"
			"|  VSTW   VR34, *+AR3[384]   \n\t"

			"   VSTW   VR35, *+AR3[400]   \n\t"
			"|  VSTW   VR36, *+AR3[416]   \n\t"

			"   VSTW   VR37, *++AR3[432]   \n\t"
			"|  SMVAGA.M2 R31,AR4           \n\t" //R31 
// index[11] read
			"	VSTW   VR38, *++AR3[16]   \n\t"
			"|  VLDDWM2 *-AR1[57],  VR1:VR0    \n\t" //vr0-0,vr1-1

			"	VSTW   VR39, *++AR3[16]   \n\t"
			"|  VLDDWM2 *-AR1[56],  VR3:VR2    \n\t" //vr2-3,vr3-4
			
			"	VSTW   VR40, *++AR3[16]   \n\t"
			"|  VLDDWM2 *-AR1[55],  VR5:VR4    \n\t" //vr4-6,vr5-7

			"	VSTW   VR41, *++AR3[16]   \n\t"
			"|  VLDDWM2 *-AR1[50],  VR7:VR6    \n\t" //vr6-9,vr7-10

			"	VSTW   VR42, *++AR3[16]   \n\t"
			"|  VLDDWM2 *-AR1[49],  VR9:VR8    \n\t" //vr4-12,vr5-13

			"   VSTW   VR43, *++AR3[16]         \n\t"
			"|  VLDDWM2 *-AR1[48],  VR11:VR10	\n\t" //vr10-15,vr11-16

			"   VSTW   VR44, *++AR3[16]         \n\t"
			"|  VLDDWM2 *-AR1[43],  VR13:VR12  \n\t" //vr12-18,vr13-19

			"   VSTW   VR45, *++AR3[16]         \n\t"
			"|  VLDDWM2 *-AR1[42],  VR15:VR14  \n\t" //vr14-21,vr15-22

			"   VSTW   VR46, *++AR3[16]         \n\t"
			"|  VLDDWM2 *-AR1[41],  VR17:VR16  \n\t" //vr16-24,vr17-25

			"   VSTW   VR47, *++AR3[16]         \n\t"
			"|	VLDDWM2	*-AR1[8],	 VR19:VR18	\n\t" //v18-2   

			"   VSTW   VR48, *++AR3[16]         \n\t"
			"|  VLDDWM2 *-AR1[7],  VR21:VR20  \n\t" //vr16-24,vr17-25

			"   VSTW   VR49, *++AR3[16]         \n\t"
			"|  VLDDWM2 *-AR1[6],  VR23:VR22  \n\t" //vr16-24,vr17-25

			"   VSTW   VR50, *++AR3[16]         \n\t"
			"|  VLDDWM2 *-AR1[1],  VR25:VR24  \n\t" //vr16-24,vr17-25

			"   VSTW   VR51, *++AR3[16]         \n\t"
			"|  VLDDWM2 *+AR1[0],  VR27:VR26  \n\t" //vr16-24,vr17-25

			"   VSTW   VR52, *++AR3[16]         \n\t"
			"|  VLDDWM2 *+AR1[1],   VR29:VR28  \n\t" //vr16-24,vr17-25

			"   VSTW   VR53, *++AR3[16]         \n\t"
			"|  VLDDWM2 *+AR1[6],  VR31:VR30  \n\t" //vr16-24,vr17-25

			"   VSTW   VR0, *++AR4[11]         \n\t"
			"|  VLDDWM2 *+AR1[7],  VR33:VR32  \n\t" //vr16-24,vr17-25
			"|	SLDW	  *+AR15[13], R10	\n\t" //sldw 7

			"   VSTW   VR1, *++AR4[16]      \n\t"
			"|	VLDDWM2 *+AR1[8],  VR35:VR34  \n\t" //vr16-24,vr17-25

			"   VSTW   VR2, *++AR4[16]		\n\t"
			"|  VLDDWM2 *+AR1[41],  VR37:VR36  \n\t" //vr16-24,vr17-25

			"   VSTW   VR3, *++AR4[16]      \n\t"
			"|  VLDDWM2 *+AR1[42],  VR39:VR38  \n\t" //vr16-24,vr17-25

			"   VSTW   VR4, *++AR4[16]      \n\t"
			"|  VLDDWM2 *+AR1[43],  VR41:VR40  \n\t" //vr16-24,vr17-25

			"   VSTW   VR5, *++AR4[16]    \n\t"
			"|  VLDDWM2 *+AR1[48],  VR43:VR42  \n\t" //vr16-24,vr17-25

			"   VSTW   VR6, *++AR4[16]   \n\t"
			"|  VLDDWM2 *+AR1[49],  VR45:VR44  \n\t" //vr16-24,vr17-25

			"   VSTW   VR7, *++AR4[16]   \n\t"
			"|  VLDDWM2 *+AR1[50],  VR47:VR46  \n\t" //vr16-24,vr17-25
			"|	SSHFLL	  3, R10, R10       \n\t"  //SSHFLL 1  

			"   VSTW   VR8, *++AR4[16]   \n\t"
			"|  VLDDWM2 *+AR1[55],  VR49:VR48  \n\t" //vr16-24,vr17-25

			"   VSTW   VR9, *++AR4[16]   \n\t"
			"|  VLDDWM2 *+AR1[56],  VR51:VR50  \n\t" //vr16-24,vr17-25

			"   VSTW   VR10, *++AR4[16]   \n\t"
			"|  VLDDWM2 *+AR1[57],  VR53:VR52  \n\t" //vr16-24,vr17-25

			"   VSTW   VR11, *+AR4[16]   \n\t"
			"|  VSTW   VR12, *+AR4[32]   \n\t"

			"   VSTW   VR13, *+AR4[48]   \n\t"
			"|  VSTW   VR14, *+AR4[64]   \n\t"

			"   VSTW   VR15, *+AR4[80]   \n\t"
			"|  VSTW   VR16, *+AR4[96]   \n\t"

			"   VSTW   VR17, *+AR4[112]   \n\t"
			"|  VSTW   VR18, *+AR4[128]   \n\t"

			"   VSTW   VR19, *+AR4[144]   \n\t"
			"|  VSTW   VR20, *+AR4[160]   \n\t"

			"   VSTW   VR21, *+AR4[176]   \n\t"
			"|  VSTW   VR22, *+AR4[192]   \n\t"

			"   VSTW   VR23, *+AR4[208]   \n\t"
			"|  VSTW   VR24, *+AR4[224]   \n\t"

			"   VSTW   VR25, *+AR4[240]   \n\t"
			"|  VSTW   VR26, *+AR4[256]   \n\t"

			"   VSTW   VR27, *+AR4[272]   \n\t"
			"|  VSTW   VR28, *+AR4[288]   \n\t"
			"|  SADDA     R10,AR0,AR1		\n\t" //sadda AVAF 3 SVAF 2

			"   VSTW   VR29, *+AR4[304]   \n\t"
			"|  VSTW   VR30, *+AR4[320]   \n\t"

			"   VSTW   VR31, *+AR4[336]   \n\t"
			"|  VSTW   VR32, *+AR4[352]   \n\t"

			"   VSTW   VR33, *+AR4[368]   \n\t"
			"|  VSTW   VR34, *+AR4[384]   \n\t"

			"   VSTW   VR35, *+AR4[400]   \n\t"
			"|  VSTW   VR36, *+AR4[416]   \n\t"

			"   VSTW   VR37, *++AR4[432]   \n\t"
			"|  SMVAGA.M2 R31,AR3           \n\t" //R31 
// index[13] read
			"	VSTW   VR38, *++AR4[16]   \n\t"
			"|  VLDDWM2 *-AR1[57],  VR1:VR0    \n\t" //vr0-0,vr1-1

			"	VSTW   VR39, *++AR4[16]   \n\t"
			"|  VLDDWM2 *-AR1[56],  VR3:VR2    \n\t" //vr2-3,vr3-4
			
			"	VSTW   VR40, *++AR4[16]   \n\t"
			"|  VLDDWM2 *-AR1[55],  VR5:VR4    \n\t" //vr4-6,vr5-7

			"	VSTW   VR41, *++AR4[16]   \n\t"
			"|  VLDDWM2 *-AR1[50],  VR7:VR6    \n\t" //vr6-9,vr7-10

			"	VSTW   VR42, *++AR4[16]   \n\t"
			"|  VLDDWM2 *-AR1[49],  VR9:VR8    \n\t" //vr4-12,vr5-13

			"   VSTW   VR43, *++AR4[16]         \n\t"
			"|  VLDDWM2 *-AR1[48],  VR11:VR10	\n\t" //vr10-15,vr11-16

			"   VSTW   VR44, *++AR4[16]         \n\t"
			"|  VLDDWM2 *-AR1[43],  VR13:VR12  \n\t" //vr12-18,vr13-19

			"   VSTW   VR45, *++AR4[16]         \n\t"
			"|  VLDDWM2 *-AR1[42],  VR15:VR14  \n\t" //vr14-21,vr15-22

			"   VSTW   VR46, *++AR4[16]         \n\t"
			"|  VLDDWM2 *-AR1[41],  VR17:VR16  \n\t" //vr16-24,vr17-25

			"   VSTW   VR47, *++AR4[16]         \n\t"
			"|	VLDDWM2	*-AR1[8],	 VR19:VR18	\n\t" //v18-2   

			"   VSTW   VR48, *++AR4[16]         \n\t"
			"|  VLDDWM2 *-AR1[7],  VR21:VR20  \n\t" //vr16-24,vr17-25

			"   VSTW   VR49, *++AR4[16]         \n\t"
			"|  VLDDWM2 *-AR1[6],  VR23:VR22  \n\t" //vr16-24,vr17-25

			"   VSTW   VR50, *++AR4[16]         \n\t"
			"|  VLDDWM2 *-AR1[1],  VR25:VR24  \n\t" //vr16-24,vr17-25

			"   VSTW   VR51, *++AR4[16]         \n\t"
			"|  VLDDWM2 *+AR1[0],  VR27:VR26  \n\t" //vr16-24,vr17-25

			"   VSTW   VR52, *++AR4[16]         \n\t"
			"|  VLDDWM2 *+AR1[1],   VR29:VR28  \n\t" //vr16-24,vr17-25

			"   VSTW   VR53, *++AR4[16]         \n\t"
			"|  VLDDWM2 *+AR1[6],  VR31:VR30  \n\t" //vr16-24,vr17-25

			"   VSTW   VR0, *++AR3[13]         \n\t"
			"|  VLDDWM2 *+AR1[7],  VR33:VR32  \n\t" //vr16-24,vr17-25
			"|	SLDW	  *+AR15[15], R10	\n\t" //sldw 7

			"   VSTW   VR1, *++AR3[16]      \n\t"
			"|	VLDDWM2 *+AR1[8],  VR35:VR34  \n\t" //vr16-24,vr17-25

			"   VSTW   VR2, *++AR3[16]		\n\t"
			"|  VLDDWM2 *+AR1[41],  VR37:VR36  \n\t" //vr16-24,vr17-25

			"   VSTW   VR3, *++AR3[16]      \n\t"
			"|  VLDDWM2 *+AR1[42],  VR39:VR38  \n\t" //vr16-24,vr17-25

			"   VSTW   VR4, *++AR3[16]      \n\t"
			"|  VLDDWM2 *+AR1[43],  VR41:VR40  \n\t" //vr16-24,vr17-25

			"   VSTW   VR5, *++AR3[16]    \n\t"
			"|  VLDDWM2 *+AR1[48],  VR43:VR42  \n\t" //vr16-24,vr17-25

			"   VSTW   VR6, *++AR3[16]   \n\t"
			"|  VLDDWM2 *+AR1[49],  VR45:VR44  \n\t" //vr16-24,vr17-25

			"   VSTW   VR7, *++AR3[16]   \n\t"
			"|  VLDDWM2 *+AR1[50],  VR47:VR46  \n\t" //vr16-24,vr17-25
			"|	SSHFLL	  3, R10, R10       \n\t"  //SSHFLL 1  

			"   VSTW   VR8, *++AR3[16]   \n\t"
			"|  VLDDWM2 *+AR1[55],  VR49:VR48  \n\t" //vr16-24,vr17-25

			"   VSTW   VR9, *++AR3[16]   \n\t"
			"|  VLDDWM2 *+AR1[56],  VR51:VR50  \n\t" //vr16-24,vr17-25

			"   VSTW   VR10, *++AR3[16]   \n\t"
			"|  VLDDWM2 *+AR1[57],  VR53:VR52  \n\t" //vr16-24,vr17-25

			"   VSTW   VR11, *+AR3[16]   \n\t"
			"|  VSTW   VR12, *+AR3[32]   \n\t"

			"   VSTW   VR13, *+AR3[48]   \n\t"
			"|  VSTW   VR14, *+AR3[64]   \n\t"

			"   VSTW   VR15, *+AR3[80]   \n\t"
			"|  VSTW   VR16, *+AR3[96]   \n\t"

			"   VSTW   VR17, *+AR3[112]   \n\t"
			"|  VSTW   VR18, *+AR3[128]   \n\t"

			"   VSTW   VR19, *+AR3[144]   \n\t"
			"|  VSTW   VR20, *+AR3[160]   \n\t"

			"   VSTW   VR21, *+AR3[176]   \n\t"
			"|  VSTW   VR22, *+AR3[192]   \n\t"

			"   VSTW   VR23, *+AR3[208]   \n\t"
			"|  VSTW   VR24, *+AR3[224]   \n\t"

			"   VSTW   VR25, *+AR3[240]   \n\t"
			"|  VSTW   VR26, *+AR3[256]   \n\t"

			"   VSTW   VR27, *+AR3[272]   \n\t"
			"|  VSTW   VR28, *+AR3[288]   \n\t"
			"|  SADDA     R10,AR0,AR1		\n\t" //sadda AVAF 3 SVAF 2

			"   VSTW   VR29, *+AR3[304]   \n\t"
			"|  VSTW   VR30, *+AR3[320]   \n\t"

			"   VSTW   VR31, *+AR3[336]   \n\t"
			"|  VSTW   VR32, *+AR3[352]   \n\t"

			"   VSTW   VR33, *+AR3[368]   \n\t"
			"|  VSTW   VR34, *+AR3[384]   \n\t"

			"   VSTW   VR35, *+AR3[400]   \n\t"
			"|  VSTW   VR36, *+AR3[416]   \n\t"

			"   VSTW   VR37, *++AR3[432]   \n\t"
			"|  SMVAGA.M2 R31,AR4           \n\t" //R31 // index[15] read
			"	VSTW   VR38, *++AR3[16]   \n\t"
			"|  VLDDWM2 *-AR1[57],  VR1:VR0    \n\t" //vr0-0,vr1-1

			"	VSTW   VR39, *++AR3[16]   \n\t"
			"|  VLDDWM2 *-AR1[56],  VR3:VR2    \n\t" //vr2-3,vr3-4
			
			"	VSTW   VR40, *++AR3[16]   \n\t"
			"|  VLDDWM2 *-AR1[55],  VR5:VR4    \n\t" //vr4-6,vr5-7

			"	VSTW   VR41, *++AR3[16]   \n\t"
			"|  VLDDWM2 *-AR1[50],  VR7:VR6    \n\t" //vr6-9,vr7-10

			"	VSTW   VR42, *++AR3[16]   \n\t"
			"|  VLDDWM2 *-AR1[49],  VR9:VR8    \n\t" //vr4-12,vr5-13

			"   VSTW   VR43, *++AR3[16]         \n\t"
			"|  VLDDWM2 *-AR1[48],  VR11:VR10	\n\t" //vr10-15,vr11-16

			"   VSTW   VR44, *++AR3[16]         \n\t"
			"|  VLDDWM2 *-AR1[43],  VR13:VR12  \n\t" //vr12-18,vr13-19

			"   VSTW   VR45, *++AR3[16]         \n\t"
			"|  VLDDWM2 *-AR1[42],  VR15:VR14  \n\t" //vr14-21,vr15-22

			"   VSTW   VR46, *++AR3[16]         \n\t"
			"|  VLDDWM2 *-AR1[41],  VR17:VR16  \n\t" //vr16-24,vr17-25

			"   VSTW   VR47, *++AR3[16]         \n\t"
			"|	VLDDWM2	*-AR1[8],	 VR19:VR18	\n\t" //v18-2   

			"   VSTW   VR48, *++AR3[16]         \n\t"
			"|  VLDDWM2 *-AR1[7],  VR21:VR20  \n\t" //vr16-24,vr17-25

			"   VSTW   VR49, *++AR3[16]         \n\t"
			"|  VLDDWM2 *-AR1[6],  VR23:VR22  \n\t" //vr16-24,vr17-25

			"   VSTW   VR50, *++AR3[16]         \n\t"
			"|  VLDDWM2 *-AR1[1],  VR25:VR24  \n\t" //vr16-24,vr17-25

			"   VSTW   VR51, *++AR3[16]         \n\t"
			"|  VLDDWM2 *+AR1[0],  VR27:VR26  \n\t" //vr16-24,vr17-25

			"   VSTW   VR52, *++AR3[16]         \n\t"
			"|  VLDDWM2 *+AR1[1],   VR29:VR28  \n\t" //vr16-24,vr17-25

			"   VSTW   VR53, *++AR3[16]         \n\t"
			"|  VLDDWM2 *+AR1[6],  VR31:VR30  \n\t" //vr16-24,vr17-25

			"   VSTW   VR0, *++AR4[15]         \n\t"
			"|  VLDDWM2 *+AR1[7],  VR33:VR32  \n\t" //vr16-24,vr17-25
			"|	SLDW	  *+AR15[0], R10	\n\t" //sldw 7

			"   VSTW   VR1, *++AR4[16]      \n\t"
			"|	VLDDWM2 *+AR1[8],  VR35:VR34  \n\t" //vr16-24,vr17-25

			"   VSTW   VR2, *++AR4[16]		\n\t"
			"|  VLDDWM2 *+AR1[41],  VR37:VR36  \n\t" //vr16-24,vr17-25

			"   VSTW   VR3, *++AR4[16]      \n\t"
			"|  VLDDWM2 *+AR1[42],  VR39:VR38  \n\t" //vr16-24,vr17-25

			"   VSTW   VR4, *++AR4[16]      \n\t"
			"|  VLDDWM2 *+AR1[43],  VR41:VR40  \n\t" //vr16-24,vr17-25

			"   VSTW   VR5, *++AR4[16]    \n\t"
			"|  VLDDWM2 *+AR1[48],  VR43:VR42  \n\t" //vr16-24,vr17-25

			"   VSTW   VR6, *++AR4[16]   \n\t"
			"|  VLDDWM2 *+AR1[49],  VR45:VR44  \n\t" //vr16-24,vr17-25

			"   VSTW   VR7, *++AR4[16]   \n\t"
			"|  VLDDWM2 *+AR1[50],  VR47:VR46  \n\t" //vr16-24,vr17-25
			"|	SSHFLL	  3, R10, R10       \n\t"  //SSHFLL 1  

			"   VSTW   VR8, *++AR4[16]   \n\t"
			"|  VLDDWM2 *+AR1[55],  VR49:VR48  \n\t" //vr16-24,vr17-25

			"   VSTW   VR9, *++AR4[16]   \n\t"
			"|  VLDDWM2 *+AR1[56],  VR51:VR50  \n\t" //vr16-24,vr17-25

			"   VSTW   VR10, *++AR4[16]   \n\t"
			"|  VLDDWM2 *+AR1[57],  VR53:VR52  \n\t" //vr16-24,vr17-25

			"   VSTW   VR11, *+AR4[16]   \n\t"
			"|  VSTW   VR12, *+AR4[32]   \n\t"
			"|  SMULIU.M1 R29,%[eb3],R30     \n\t" //smuliu 3

			"   VSTW   VR13, *+AR4[48]   \n\t"
			"|  VSTW   VR14, *+AR4[64]   \n\t"

			"   VSTW   VR15, *+AR4[80]   \n\t"
			"|  VSTW   VR16, *+AR4[96]   \n\t"

			"   VSTW   VR17, *+AR4[112]   \n\t"
			"|  VSTW   VR18, *+AR4[128]   \n\t"
			"|  SMVAGA.M1 R30,AR0           \n\t" //smvaga 2

			"   VSTW   VR19, *+AR4[144]   \n\t"
			"|  VSTW   VR20, *+AR4[160]   \n\t"

			"   VSTW   VR21, *+AR4[176]   \n\t"
			"|  VSTW   VR22, *+AR4[192]   \n\t"

			"   VSTW   VR23, *+AR4[208]   \n\t"
			"|  VSTW   VR24, *+AR4[224]   \n\t"

			"   VSTW   VR25, *+AR4[240]   \n\t"
			"|  VSTW   VR26, *+AR4[256]   \n\t"

			"   VSTW   VR27, *+AR4[272]   \n\t"
			"|  VSTW   VR28, *+AR4[288]   \n\t"
			"|  SADDA     R10,AR0,AR1		\n\t" //sadda AVAF 3 SVAF 2

			"   VSTW   VR29, *+AR4[304]   \n\t"
			"|  VSTW   VR30, *+AR4[320]   \n\t"
			"|	SMOVI24	0x1B00,R11		  \n\t" 

			"   VSTW   VR31, *+AR4[336]   \n\t"
			"|  VSTW   VR32, *+AR4[352]   \n\t"
			"|	SADD	R11,R31,R31		  \n\t" 

			"   VSTW   VR33, *+AR4[368]   \n\t"
			"|  VSTW   VR34, *+AR4[384]   \n\t"

			"   VSTW   VR35, *+AR4[400]   \n\t"
			"|  VSTW   VR36, *+AR4[416]   \n\t"

			"   VSTW   VR37, *+AR4[432]   \n\t"
			"|  SMVAGA.M2 R31,AR3           \n\t" //R31 
			"	VSTW   VR38, *+AR4[448]   \n\t"
			"|	VSTW   VR39, *+AR4[464]   \n\t"

			"	VSTW   VR40, *+AR4[480]   \n\t"
			"|	VSTW   VR41, *+AR4[496]   \n\t"

			"	VSTW   VR42, *+AR4[512]   \n\t"
			"|  VSTW   VR43, *+AR4[528]         \n\t"

			"   VSTW   VR44, *+AR4[544]         \n\t"
			"|  VSTW   VR45, *+AR4[560]         \n\t"

			"   VSTW   VR46, *+AR4[576]         \n\t"
			"|  VSTW   VR47, *+AR4[592]         \n\t"

			"   VSTW   VR48, *+AR4[608]         \n\t"
			"|  VSTW   VR49, *+AR4[624]         \n\t"

			"   VSTW   VR50, *+AR4[640]         \n\t"
			"|  VSTW   VR51, *+AR4[656]         \n\t"

			"   VSTW   VR52, *+AR4[672]         \n\t"
			"|  VSTW   VR53, *+AR4[688]         \n\t"

//EB3 //index[0] read
			"   VLDDWM2 *-AR1[57],  VR1:VR0    \n\t" //vr0-0,vr1-1
			"|  VLDDWM2 *-AR1[56],  VR3:VR2    \n\t" //vr2-3,vr3-4

			"   VLDDWM2 *-AR1[55],  VR5:VR4    \n\t" //vr4-6,vr5-7
			"|  VLDDWM2 *-AR1[50],  VR7:VR6    \n\t" //vr6-9,vr7-10

			"   VLDDWM2 *-AR1[49],  VR9:VR8    \n\t" //vr4-12,vr5-13
			"|  VLDDWM2 *-AR1[48],  VR11:VR10	\n\t" //vr10-15,vr11-16

			"   VLDDWM2 *-AR1[43],  VR13:VR12  \n\t" //vr12-18,vr13-19
			"|  VLDDWM2 *-AR1[42],  VR15:VR14  \n\t" //vr14-21,vr15-22

			"   VLDDWM2 *-AR1[41],  VR17:VR16  \n\t" //vr16-24,vr17-25
			"|	VLDDWM2	*-AR1[8],	 VR19:VR18	\n\t" //v18-2   

			"   VLDDWM2 *-AR1[7],  VR21:VR20  \n\t" //vr16-24,vr17-25
			"|  VLDDWM2 *-AR1[6],  VR23:VR22  \n\t" //vr16-24,vr17-25

			"   VLDDWM2 *-AR1[1],  VR25:VR24  \n\t" //vr16-24,vr17-25
			"|  VLDDWM2 *+AR1[0],  VR27:VR26  \n\t" //vr16-24,vr17-25

			"   VLDDWM2 *+AR1[1],   VR29:VR28  \n\t" //vr16-24,vr17-25
			"|  VLDDWM2 *+AR1[6],  VR31:VR30  \n\t" //vr16-24,vr17-25

			"	SNOP	1\n\t"

			"   VSTDW   VR1:VR0, *+AR3[0]         \n\t"
			"|  VLDDWM2 *+AR1[7],  VR33:VR32  \n\t" //vr16-24,vr17-25
			"|	SLDW	  *+AR15[2], R10	\n\t" //sldw 7

			"   VSTDW   VR3:VR2, *+AR3[16]      \n\t"
			"|  VLDDWM2 *+AR1[8],  VR35:VR34  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR5:VR4, *+AR3[32]		\n\t"
			"|  VLDDWM2 *+AR1[41],  VR37:VR36  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR7:VR6, *+AR3[48]      \n\t"
			"|  VLDDWM2 *+AR1[42],  VR39:VR38  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR9:VR8, *+AR3[64]      \n\t"
			"|  VLDDWM2 *+AR1[43],  VR41:VR40  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR11:VR10, *+AR3[80]    \n\t"
			"|  VLDDWM2 *+AR1[48],  VR43:VR42  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR13:VR12, *+AR3[96]   \n\t"
			"|  VLDDWM2 *+AR1[49],  VR45:VR44  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR15:VR14, *+AR3[112]   \n\t"
			"|  VLDDWM2 *+AR1[50],  VR47:VR46  \n\t" //vr16-24,vr17-25
			"|	SSHFLL	  3, R10, R10       \n\t"  //SSHFLL 1  

			"   VSTDW   VR17:VR16, *+AR3[128]   \n\t"
			"|  VLDDWM2 *+AR1[55],  VR49:VR48  \n\t" //vr16-24,vr17-25
			"|  SADDA     R10,AR0,AR2		\n\t" //sadda AVAF 3 SVAF 2

			"   VSTDW   VR19:VR18, *+AR3[144]   \n\t"
			"|  VLDDWM2 *+AR1[56],  VR51:VR50  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR21:VR20, *+AR3[160]   \n\t"
			"|  VLDDWM2 *+AR1[57],  VR53:VR52  \n\t" //vr16-24,vr17-25
////////////////////////////////////////////////////////////////////////////
// index[2] read
			"	VSTDW   VR23:VR22, *+AR3[176]   \n\t"
			"|  VLDDWM2 *-AR2[57],  VR1:VR0    \n\t" //vr0-0,vr1-1

			"	VSTDW   VR25:VR24, *+AR3[192]   \n\t"
			"|  VLDDWM2 *-AR2[56],  VR3:VR2    \n\t" //vr2-3,vr3-4
			
			"	VSTDW   VR27:VR26, *+AR3[208]   \n\t"
			"|  VLDDWM2 *-AR2[55],  VR5:VR4    \n\t" //vr4-6,vr5-7

			"	VSTDW   VR29:VR28, *+AR3[224]   \n\t"
			"|  VLDDWM2 *-AR2[50],  VR7:VR6    \n\t" //vr6-9,vr7-10

			"	VSTDW   VR31:VR30, *+AR3[240]   \n\t"
			"|  VLDDWM2 *-AR2[49],  VR9:VR8    \n\t" //vr4-12,vr5-13

			"   VSTDW   VR33:VR32, *+AR3[256]         \n\t"
			"|  VLDDWM2 *-AR2[48],  VR11:VR10	\n\t" //vr10-15,vr11-16

			"   VSTDW   VR35:VR34, *+AR3[272]         \n\t"
			"|  VLDDWM2 *-AR2[43],  VR13:VR12  \n\t" //vr12-18,vr13-19

			"   VSTDW   VR37:VR36, *+AR3[288]         \n\t"
			"|  VLDDWM2 *-AR2[42],  VR15:VR14  \n\t" //vr14-21,vr15-22

			"   VSTDW   VR39:VR38, *+AR3[304]         \n\t"
			"|  VLDDWM2 *-AR2[41],  VR17:VR16  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR41:VR40, *+AR3[320]         \n\t"
			"|	VLDDWM2	*-AR2[8],	 VR19:VR18	\n\t" //v18-2   

			"   VSTDW   VR43:VR42, *+AR3[336]         \n\t"
			"|  VLDDWM2 *-AR2[7],  VR21:VR20  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR45:VR44, *+AR3[352]         \n\t"
			"|  VLDDWM2 *-AR2[6],  VR23:VR22  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR47:VR46, *+AR3[368]         \n\t"
			"|  VLDDWM2 *-AR2[1],  VR25:VR24  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR49:VR48, *+AR3[384]         \n\t"
			"|  VLDDWM2 *+AR2[0],  VR27:VR26  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR51:VR50, *+AR3[400]         \n\t"
			"|  VLDDWM2 *+AR2[1],   VR29:VR28  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR53:VR52, *+AR3[416]         \n\t"
			"|  VLDDWM2 *+AR2[6],  VR31:VR30  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR1:VR0, *+AR3[1]         \n\t"
			"|  VLDDWM2 *+AR2[7],  VR33:VR32  \n\t" //vr16-24,vr17-25
			"|	SLDW	  *+AR15[4], R10	\n\t" //sldw 7

			"   VSTDW   VR3:VR2, *+AR3[17]      \n\t"
			"|	VLDDWM2 *+AR2[8],  VR35:VR34  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR5:VR4, *+AR3[33]		\n\t"
			"|  VLDDWM2 *+AR2[41],  VR37:VR36  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR7:VR6, *+AR3[49]      \n\t"
			"|  VLDDWM2 *+AR2[42],  VR39:VR38  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR9:VR8, *+AR3[65]      \n\t"
			"|  VLDDWM2 *+AR2[43],  VR41:VR40  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR11:VR10, *+AR3[81]    \n\t"
			"|  VLDDWM2 *+AR2[48],  VR43:VR42  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR13:VR12, *+AR3[97]   \n\t"
			"|  VLDDWM2 *+AR2[49],  VR45:VR44  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR15:VR14, *+AR3[113]   \n\t"
			"|  VLDDWM2 *+AR2[50],  VR47:VR46  \n\t" //vr16-24,vr17-25
			"|	SSHFLL	  3, R10, R10       \n\t"  //SSHFLL 1  

			"   VSTDW   VR17:VR16, *+AR3[129]   \n\t"
			"|  VLDDWM2 *+AR2[55],  VR49:VR48  \n\t" //vr16-24,vr17-25
			"|  SADDA     R10,AR0,AR1		\n\t" //sadda AVAF 3 SVAF 2

			"   VSTDW   VR19:VR18, *+AR3[145]   \n\t"
			"|  VLDDWM2 *+AR2[56],  VR51:VR50  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR21:VR20, *+AR3[161]   \n\t"
			"|  VLDDWM2 *+AR2[57],  VR53:VR52  \n\t" //vr16-24,vr17-25
// index[4] read
			"	VSTDW   VR23:VR22, *+AR3[177]   \n\t"
			"|  VLDDWM2 *-AR1[57],  VR1:VR0    \n\t" //vr0-0,vr1-1

			"	VSTDW   VR25:VR24, *+AR3[193]   \n\t"
			"|  VLDDWM2 *-AR1[56],  VR3:VR2    \n\t" //vr2-3,vr3-4
			
			"	VSTDW   VR27:VR26, *+AR3[209]   \n\t"
			"|  VLDDWM2 *-AR1[55],  VR5:VR4    \n\t" //vr4-6,vr5-7

			"	VSTDW   VR29:VR28, *+AR3[225]   \n\t"
			"|  VLDDWM2 *-AR1[50],  VR7:VR6    \n\t" //vr6-9,vr7-10

			"	VSTDW   VR31:VR30, *+AR3[241]   \n\t"
			"|  VLDDWM2 *-AR1[49],  VR9:VR8    \n\t" //vr4-12,vr5-13

			"   VSTDW   VR33:VR32, *+AR3[257]         \n\t"
			"|  VLDDWM2 *-AR1[48],  VR11:VR10	\n\t" //vr10-15,vr11-16

			"   VSTDW   VR35:VR34, *+AR3[273]         \n\t"
			"|  VLDDWM2 *-AR1[43],  VR13:VR12  \n\t" //vr12-18,vr13-19

			"   VSTDW   VR37:VR36, *+AR3[289]         \n\t"
			"|  VLDDWM2 *-AR1[42],  VR15:VR14  \n\t" //vr14-21,vr15-22

			"   VSTDW   VR39:VR38, *+AR3[305]         \n\t"
			"|  VLDDWM2 *-AR1[41],  VR17:VR16  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR41:VR40, *+AR3[321]         \n\t"
			"|	VLDDWM2	*-AR1[8],	 VR19:VR18	\n\t" //v18-2   

			"   VSTDW   VR43:VR42, *+AR3[337]         \n\t"
			"|  VLDDWM2 *-AR1[7],  VR21:VR20  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR45:VR44, *+AR3[353]         \n\t"
			"|  VLDDWM2 *-AR1[6],  VR23:VR22  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR47:VR46, *+AR3[369]         \n\t"
			"|  VLDDWM2 *-AR1[1],  VR25:VR24  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR49:VR48, *+AR3[385]         \n\t"
			"|  VLDDWM2 *+AR1[0],  VR27:VR26  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR51:VR50, *+AR3[401]         \n\t"
			"|  VLDDWM2 *+AR1[1],   VR29:VR28  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR53:VR52, *+AR3[417]         \n\t"
			"|  VLDDWM2 *+AR1[6],  VR31:VR30  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR1:VR0, *+AR3[2]         \n\t"
			"|  VLDDWM2 *+AR1[7],  VR33:VR32  \n\t" //vr16-24,vr17-25
			"|	SLDW	  *+AR15[6], R10	\n\t" //sldw 7

			"   VSTDW   VR3:VR2, *+AR3[18]      \n\t"
			"|	VLDDWM2 *+AR1[8],  VR35:VR34  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR5:VR4, *+AR3[34]		\n\t"
			"|  VLDDWM2 *+AR1[41],  VR37:VR36  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR7:VR6, *+AR3[50]      \n\t"
			"|  VLDDWM2 *+AR1[42],  VR39:VR38  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR9:VR8, *+AR3[66]      \n\t"
			"|  VLDDWM2 *+AR1[43],  VR41:VR40  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR11:VR10, *+AR3[82]    \n\t"
			"|  VLDDWM2 *+AR1[48],  VR43:VR42  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR13:VR12, *+AR3[98]   \n\t"
			"|  VLDDWM2 *+AR1[49],  VR45:VR44  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR15:VR14, *+AR3[114]   \n\t"
			"|  VLDDWM2 *+AR1[50],  VR47:VR46  \n\t" //vr16-24,vr17-25
			"|	SSHFLL	  3, R10, R10       \n\t"  //SSHFLL 1  

			"   VSTDW   VR17:VR16, *+AR3[130]   \n\t"
			"|  VLDDWM2 *+AR1[55],  VR49:VR48  \n\t" //vr16-24,vr17-25
			"|  SADDA     R10,AR0,AR2		\n\t" //sadda AVAF 3 SVAF 2

			"   VSTDW   VR19:VR18, *+AR3[146]   \n\t"
			"|  VLDDWM2 *+AR1[56],  VR51:VR50  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR21:VR20, *+AR3[162]   \n\t"
			"|  VLDDWM2 *+AR1[57],  VR53:VR52  \n\t" //vr16-24,vr17-25
// index[6] read
			"	VSTDW   VR23:VR22, *+AR3[178]   \n\t"
			"|  VLDDWM2 *-AR2[57],  VR1:VR0    \n\t" //vr0-0,vr1-1

			"	VSTDW   VR25:VR24, *+AR3[194]   \n\t"
			"|  VLDDWM2 *-AR2[56],  VR3:VR2    \n\t" //vr2-3,vr3-4
			
			"	VSTDW   VR27:VR26, *+AR3[210]   \n\t"
			"|  VLDDWM2 *-AR2[55],  VR5:VR4    \n\t" //vr4-6,vr5-7

			"	VSTDW   VR29:VR28, *+AR3[226]   \n\t"
			"|  VLDDWM2 *-AR2[50],  VR7:VR6    \n\t" //vr6-9,vr7-10

			"	VSTDW   VR31:VR30, *+AR3[242]   \n\t"
			"|  VLDDWM2 *-AR2[49],  VR9:VR8    \n\t" //vr4-12,vr5-13

			"   VSTDW   VR33:VR32, *+AR3[258]         \n\t"
			"|  VLDDWM2 *-AR2[48],  VR11:VR10	\n\t" //vr10-15,vr11-16

			"   VSTDW   VR35:VR34, *+AR3[274]         \n\t"
			"|  VLDDWM2 *-AR2[43],  VR13:VR12  \n\t" //vr12-18,vr13-19

			"   VSTDW   VR37:VR36, *+AR3[290]         \n\t"
			"|  VLDDWM2 *-AR2[42],  VR15:VR14  \n\t" //vr14-21,vr15-22

			"   VSTDW   VR39:VR38, *+AR3[306]         \n\t"
			"|  VLDDWM2 *-AR2[41],  VR17:VR16  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR41:VR40, *+AR3[322]         \n\t"
			"|	VLDDWM2	*-AR2[8],	 VR19:VR18	\n\t" //v18-2   

			"   VSTDW   VR43:VR42, *+AR3[338]         \n\t"
			"|  VLDDWM2 *-AR2[7],  VR21:VR20  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR45:VR44, *+AR3[354]         \n\t"
			"|  VLDDWM2 *-AR2[6],  VR23:VR22  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR47:VR46, *+AR3[370]         \n\t"
			"|  VLDDWM2 *-AR2[1],  VR25:VR24  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR49:VR48, *+AR3[386]         \n\t"
			"|  VLDDWM2 *+AR2[0],  VR27:VR26  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR51:VR50, *+AR3[402]         \n\t"
			"|  VLDDWM2 *+AR2[1],   VR29:VR28  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR53:VR52, *+AR3[418]         \n\t"
			"|  VLDDWM2 *+AR2[6],  VR31:VR30  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR1:VR0, *+AR3[3]         \n\t"
			"|  VLDDWM2 *+AR2[7],  VR33:VR32  \n\t" //vr16-24,vr17-25
			"|	SLDW	  *+AR15[8], R10	\n\t" //sldw 7

			"   VSTDW   VR3:VR2, *+AR3[19]      \n\t"
			"|	VLDDWM2 *+AR2[8],  VR35:VR34  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR5:VR4, *+AR3[35]		\n\t"
			"|  VLDDWM2 *+AR2[41],  VR37:VR36  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR7:VR6, *+AR3[51]      \n\t"
			"|  VLDDWM2 *+AR2[42],  VR39:VR38  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR9:VR8, *+AR3[67]      \n\t"
			"|  VLDDWM2 *+AR2[43],  VR41:VR40  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR11:VR10, *+AR3[83]    \n\t"
			"|  VLDDWM2 *+AR2[48],  VR43:VR42  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR13:VR12, *+AR3[99]   \n\t"
			"|  VLDDWM2 *+AR2[49],  VR45:VR44  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR15:VR14, *+AR3[115]   \n\t"
			"|  VLDDWM2 *+AR2[50],  VR47:VR46  \n\t" //vr16-24,vr17-25
			"|	SSHFLL	  3, R10, R10       \n\t"  //SSHFLL 1  

			"   VSTDW   VR17:VR16, *+AR3[131]   \n\t"
			"|  VLDDWM2 *+AR2[55],  VR49:VR48  \n\t" //vr16-24,vr17-25
			"|  SADDA     R10,AR0,AR1		\n\t" //sadda AVAF 3 SVAF 2

			"   VSTDW   VR19:VR18, *+AR3[147]   \n\t"
			"|  VLDDWM2 *+AR2[56],  VR51:VR50  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR21:VR20, *+AR3[163]   \n\t"
			"|  VLDDWM2 *+AR2[57],  VR53:VR52  \n\t" //vr16-24,vr17-25
// index[8] read
			"	VSTDW   VR23:VR22, *+AR3[179]   \n\t"
			"|  VLDDWM2 *-AR1[57],  VR1:VR0    \n\t" //vr0-0,vr1-1

			"	VSTDW   VR25:VR24, *+AR3[195]   \n\t"
			"|  VLDDWM2 *-AR1[56],  VR3:VR2    \n\t" //vr2-3,vr3-4
			
			"	VSTDW   VR27:VR26, *+AR3[211]   \n\t"
			"|  VLDDWM2 *-AR1[55],  VR5:VR4    \n\t" //vr4-6,vr5-7

			"	VSTDW   VR29:VR28, *+AR3[227]   \n\t"
			"|  VLDDWM2 *-AR1[50],  VR7:VR6    \n\t" //vr6-9,vr7-10

			"	VSTDW   VR31:VR30, *+AR3[243]   \n\t"
			"|  VLDDWM2 *-AR1[49],  VR9:VR8    \n\t" //vr4-12,vr5-13

			"   VSTDW   VR33:VR32, *+AR3[259]         \n\t"
			"|  VLDDWM2 *-AR1[48],  VR11:VR10	\n\t" //vr10-15,vr11-16

			"   VSTDW   VR35:VR34, *+AR3[275]         \n\t"
			"|  VLDDWM2 *-AR1[43],  VR13:VR12  \n\t" //vr12-18,vr13-19

			"   VSTDW   VR37:VR36, *+AR3[291]         \n\t"
			"|  VLDDWM2 *-AR1[42],  VR15:VR14  \n\t" //vr14-21,vr15-22

			"   VSTDW   VR39:VR38, *+AR3[307]         \n\t"
			"|  VLDDWM2 *-AR1[41],  VR17:VR16  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR41:VR40, *+AR3[323]         \n\t"
			"|	VLDDWM2	*-AR1[8],	 VR19:VR18	\n\t" //v18-2   

			"   VSTDW   VR43:VR42, *+AR3[339]         \n\t"
			"|  VLDDWM2 *-AR1[7],  VR21:VR20  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR45:VR44, *+AR3[355]         \n\t"
			"|  VLDDWM2 *-AR1[6],  VR23:VR22  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR47:VR46, *+AR3[371]         \n\t"
			"|  VLDDWM2 *-AR1[1],  VR25:VR24  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR49:VR48, *+AR3[387]         \n\t"
			"|  VLDDWM2 *+AR1[0],  VR27:VR26  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR51:VR50, *+AR3[403]         \n\t"
			"|  VLDDWM2 *+AR1[1],   VR29:VR28  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR53:VR52, *+AR3[419]         \n\t"
			"|  VLDDWM2 *+AR1[6],  VR31:VR30  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR1:VR0, *+AR3[4]         \n\t"
			"|  VLDDWM2 *+AR1[7],  VR33:VR32  \n\t" //vr16-24,vr17-25
			"|	SLDW	  *+AR15[10], R10	\n\t" //sldw 7

			"   VSTDW   VR3:VR2, *+AR3[20]      \n\t"
			"|	VLDDWM2 *+AR1[8],  VR35:VR34  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR5:VR4, *+AR3[36]		\n\t"
			"|  VLDDWM2 *+AR1[41],  VR37:VR36  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR7:VR6, *+AR3[52]      \n\t"
			"|  VLDDWM2 *+AR1[42],  VR39:VR38  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR9:VR8, *+AR3[68]      \n\t"
			"|  VLDDWM2 *+AR1[43],  VR41:VR40  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR11:VR10, *+AR3[84]    \n\t"
			"|  VLDDWM2 *+AR1[48],  VR43:VR42  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR13:VR12, *+AR3[100]   \n\t"
			"|  VLDDWM2 *+AR1[49],  VR45:VR44  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR15:VR14, *+AR3[116]   \n\t"
			"|  VLDDWM2 *+AR1[50],  VR47:VR46  \n\t" //vr16-24,vr17-25
			"|	SSHFLL	  3, R10, R10       \n\t"  //SSHFLL 1  

			"   VSTDW   VR17:VR16, *+AR3[132]   \n\t"
			"|  VLDDWM2 *+AR1[55],  VR49:VR48  \n\t" //vr16-24,vr17-25
			"|  SADDA     R10,AR0,AR2		\n\t" //sadda AVAF 3 SVAF 2

			"   VSTDW   VR19:VR18, *+AR3[148]   \n\t"
			"|  VLDDWM2 *+AR1[56],  VR51:VR50  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR21:VR20, *+AR3[164]   \n\t"
			"|  VLDDWM2 *+AR1[57],  VR53:VR52  \n\t" //vr16-24,vr17-25
// index[10] read
			"	VSTDW   VR23:VR22, *+AR3[180]   \n\t"
			"|  VLDDWM2 *-AR2[57],  VR1:VR0    \n\t" //vr0-0,vr1-1

			"	VSTDW   VR25:VR24, *+AR3[196]   \n\t"
			"|  VLDDWM2 *-AR2[56],  VR3:VR2    \n\t" //vr2-3,vr3-4
			
			"	VSTDW   VR27:VR26, *+AR3[212]   \n\t"
			"|  VLDDWM2 *-AR2[55],  VR5:VR4    \n\t" //vr4-6,vr5-7

			"	VSTDW   VR29:VR28, *+AR3[228]   \n\t"
			"|  VLDDWM2 *-AR2[50],  VR7:VR6    \n\t" //vr6-9,vr7-10

			"	VSTDW   VR31:VR30, *+AR3[244]   \n\t"
			"|  VLDDWM2 *-AR2[49],  VR9:VR8    \n\t" //vr4-12,vr5-13

			"   VSTDW   VR33:VR32, *+AR3[260]         \n\t"
			"|  VLDDWM2 *-AR2[48],  VR11:VR10	\n\t" //vr10-15,vr11-16

			"   VSTDW   VR35:VR34, *+AR3[276]         \n\t"
			"|  VLDDWM2 *-AR2[43],  VR13:VR12  \n\t" //vr12-18,vr13-19

			"   VSTDW   VR37:VR36, *+AR3[292]         \n\t"
			"|  VLDDWM2 *-AR2[42],  VR15:VR14  \n\t" //vr14-21,vr15-22

			"   VSTDW   VR39:VR38, *+AR3[308]         \n\t"
			"|  VLDDWM2 *-AR2[41],  VR17:VR16  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR41:VR40, *+AR3[324]         \n\t"
			"|	VLDDWM2	*-AR2[8],	 VR19:VR18	\n\t" //v18-2   

			"   VSTDW   VR43:VR42, *+AR3[340]         \n\t"
			"|  VLDDWM2 *-AR2[7],  VR21:VR20  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR45:VR44, *+AR3[356]         \n\t"
			"|  VLDDWM2 *-AR2[6],  VR23:VR22  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR47:VR46, *+AR3[372]         \n\t"
			"|  VLDDWM2 *-AR2[1],  VR25:VR24  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR49:VR48, *+AR3[388]         \n\t"
			"|  VLDDWM2 *+AR2[0],  VR27:VR26  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR51:VR50, *+AR3[404]         \n\t"
			"|  VLDDWM2 *+AR2[1],   VR29:VR28  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR53:VR52, *+AR3[420]         \n\t"
			"|  VLDDWM2 *+AR2[6],  VR31:VR30  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR1:VR0, *+AR3[5]         \n\t"
			"|  VLDDWM2 *+AR2[7],  VR33:VR32  \n\t" //vr16-24,vr17-25
			"|	SLDW	  *+AR15[12], R10	\n\t" //sldw 7

			"   VSTDW   VR3:VR2, *+AR3[21]      \n\t"
			"|	VLDDWM2 *+AR2[8],  VR35:VR34  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR5:VR4, *+AR3[37]		\n\t"
			"|  VLDDWM2 *+AR2[41],  VR37:VR36  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR7:VR6, *+AR3[53]      \n\t"
			"|  VLDDWM2 *+AR2[42],  VR39:VR38  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR9:VR8, *+AR3[69]      \n\t"
			"|  VLDDWM2 *+AR2[43],  VR41:VR40  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR11:VR10, *+AR3[85]    \n\t"
			"|  VLDDWM2 *+AR2[48],  VR43:VR42  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR13:VR12, *+AR3[101]   \n\t"
			"|  VLDDWM2 *+AR2[49],  VR45:VR44  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR15:VR14, *+AR3[117]   \n\t"
			"|  VLDDWM2 *+AR2[50],  VR47:VR46  \n\t" //vr16-24,vr17-25
			"|	SSHFLL	  3, R10, R10       \n\t"  //SSHFLL 1  

			"   VSTDW   VR17:VR16, *+AR3[133]   \n\t"
			"|  VLDDWM2 *+AR2[55],  VR49:VR48  \n\t" //vr16-24,vr17-25
			"|  SADDA     R10,AR0,AR1		\n\t" //sadda AVAF 3 SVAF 2

			"   VSTDW   VR19:VR18, *+AR3[149]   \n\t"
			"|  VLDDWM2 *+AR2[56],  VR51:VR50  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR21:VR20, *+AR3[165]   \n\t"
			"|  VLDDWM2 *+AR2[57],  VR53:VR52  \n\t" //vr16-24,vr17-25
// index[12] read
			"	VSTDW   VR23:VR22, *+AR3[181]   \n\t"
			"|  VLDDWM2 *-AR1[57],  VR1:VR0    \n\t" //vr0-0,vr1-1

			"	VSTDW   VR25:VR24, *+AR3[197]   \n\t"
			"|  VLDDWM2 *-AR1[56],  VR3:VR2    \n\t" //vr2-3,vr3-4
			
			"	VSTDW   VR27:VR26, *+AR3[213]   \n\t"
			"|  VLDDWM2 *-AR1[55],  VR5:VR4    \n\t" //vr4-6,vr5-7

			"	VSTDW   VR29:VR28, *+AR3[229]   \n\t"
			"|  VLDDWM2 *-AR1[50],  VR7:VR6    \n\t" //vr6-9,vr7-10

			"	VSTDW   VR31:VR30, *+AR3[245]   \n\t"
			"|  VLDDWM2 *-AR1[49],  VR9:VR8    \n\t" //vr4-12,vr5-13

			"   VSTDW   VR33:VR32, *+AR3[261]         \n\t"
			"|  VLDDWM2 *-AR1[48],  VR11:VR10	\n\t" //vr10-15,vr11-16

			"   VSTDW   VR35:VR34, *+AR3[277]         \n\t"
			"|  VLDDWM2 *-AR1[43],  VR13:VR12  \n\t" //vr12-18,vr13-19

			"   VSTDW   VR37:VR36, *+AR3[293]         \n\t"
			"|  VLDDWM2 *-AR1[42],  VR15:VR14  \n\t" //vr14-21,vr15-22

			"   VSTDW   VR39:VR38, *+AR3[309]         \n\t"
			"|  VLDDWM2 *-AR1[41],  VR17:VR16  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR41:VR40, *+AR3[325]         \n\t"
			"|	VLDDWM2	*-AR1[8],	 VR19:VR18	\n\t" //v18-2   

			"   VSTDW   VR43:VR42, *+AR3[341]         \n\t"
			"|  VLDDWM2 *-AR1[7],  VR21:VR20  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR45:VR44, *+AR3[357]         \n\t"
			"|  VLDDWM2 *-AR1[6],  VR23:VR22  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR47:VR46, *+AR3[373]         \n\t"
			"|  VLDDWM2 *-AR1[1],  VR25:VR24  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR49:VR48, *+AR3[389]         \n\t"
			"|  VLDDWM2 *+AR1[0],  VR27:VR26  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR51:VR50, *+AR3[405]         \n\t"
			"|  VLDDWM2 *+AR1[1],   VR29:VR28  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR53:VR52, *+AR3[421]         \n\t"
			"|  VLDDWM2 *+AR1[6],  VR31:VR30  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR1:VR0, *+AR3[6]         \n\t"
			"|  VLDDWM2 *+AR1[7],  VR33:VR32  \n\t" //vr16-24,vr17-25
			"|	SLDW	  *+AR15[14], R10	\n\t" //sldw 7

			"   VSTDW   VR3:VR2, *+AR3[22]      \n\t"
			"|	VLDDWM2 *+AR1[8],  VR35:VR34  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR5:VR4, *+AR3[38]		\n\t"
			"|  VLDDWM2 *+AR1[41],  VR37:VR36  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR7:VR6, *+AR3[54]      \n\t"
			"|  VLDDWM2 *+AR1[42],  VR39:VR38  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR9:VR8, *+AR3[70]      \n\t"
			"|  VLDDWM2 *+AR1[43],  VR41:VR40  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR11:VR10, *+AR3[86]    \n\t"
			"|  VLDDWM2 *+AR1[48],  VR43:VR42  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR13:VR12, *+AR3[102]   \n\t"
			"|  VLDDWM2 *+AR1[49],  VR45:VR44  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR15:VR14, *+AR3[118]   \n\t"
			"|  VLDDWM2 *+AR1[50],  VR47:VR46  \n\t" //vr16-24,vr17-25
			"|	SSHFLL	  3, R10, R10       \n\t"  //SSHFLL 1  

			"   VSTDW   VR17:VR16, *+AR3[134]   \n\t"
			"|  VLDDWM2 *+AR1[55],  VR49:VR48  \n\t" //vr16-24,vr17-25
			"|  SADDA     R10,AR0,AR2		\n\t" //sadda AVAF 3 SVAF 2

			"   VSTDW   VR19:VR18, *+AR3[150]   \n\t"
			"|  VLDDWM2 *+AR1[56],  VR51:VR50  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR21:VR20, *+AR3[166]   \n\t"
			"|  VLDDWM2 *+AR1[57],  VR53:VR52  \n\t" //vr16-24,vr17-25
// index[14] read
			"	VSTDW   VR23:VR22, *+AR3[182]   \n\t"
			"|  VLDDWM2 *-AR2[57],  VR1:VR0    \n\t" //vr0-0,vr1-1

			"	VSTDW   VR25:VR24, *+AR3[198]   \n\t"
			"|  VLDDWM2 *-AR2[56],  VR3:VR2    \n\t" //vr2-3,vr3-4
			
			"	VSTDW   VR27:VR26, *+AR3[214]   \n\t"
			"|  VLDDWM2 *-AR2[55],  VR5:VR4    \n\t" //vr4-6,vr5-7

			"	VSTDW   VR29:VR28, *+AR3[230]   \n\t"
			"|  VLDDWM2 *-AR2[50],  VR7:VR6    \n\t" //vr6-9,vr7-10

			"	VSTDW   VR31:VR30, *+AR3[246]   \n\t"
			"|  VLDDWM2 *-AR2[49],  VR9:VR8    \n\t" //vr4-12,vr5-13

			"   VSTDW   VR33:VR32, *+AR3[262]         \n\t"
			"|  VLDDWM2 *-AR2[48],  VR11:VR10	\n\t" //vr10-15,vr11-16

			"   VSTDW   VR35:VR34, *+AR3[278]         \n\t"
			"|  VLDDWM2 *-AR2[43],  VR13:VR12  \n\t" //vr12-18,vr13-19

			"   VSTDW   VR37:VR36, *+AR3[294]         \n\t"
			"|  VLDDWM2 *-AR2[42],  VR15:VR14  \n\t" //vr14-21,vr15-22

			"   VSTDW   VR39:VR38, *+AR3[310]         \n\t"
			"|  VLDDWM2 *-AR2[41],  VR17:VR16  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR41:VR40, *+AR3[326]         \n\t"
			"|	VLDDWM2	*-AR2[8],	 VR19:VR18	\n\t" //v18-2   

			"   VSTDW   VR43:VR42, *+AR3[342]         \n\t"
			"|  VLDDWM2 *-AR2[7],  VR21:VR20  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR45:VR44, *+AR3[358]         \n\t"
			"|  VLDDWM2 *-AR2[6],  VR23:VR22  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR47:VR46, *+AR3[374]         \n\t"
			"|  VLDDWM2 *-AR2[1],  VR25:VR24  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR49:VR48, *+AR3[390]         \n\t"
			"|  VLDDWM2 *+AR2[0],  VR27:VR26  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR51:VR50, *+AR3[406]         \n\t"
			"|  VLDDWM2 *+AR2[1],   VR29:VR28  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR53:VR52, *+AR3[422]         \n\t"
			"|  VLDDWM2 *+AR2[6],  VR31:VR30  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR1:VR0, *+AR3[7]         \n\t"
			"|  VLDDWM2 *+AR2[7],  VR33:VR32  \n\t" //vr16-24,vr17-25
			"|	SLDW	  *+AR15[1], R10	\n\t" //sldw 7

			"   VSTDW   VR3:VR2, *+AR3[23]      \n\t"
			"|	VLDDWM2 *+AR2[8],  VR35:VR34  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR5:VR4, *+AR3[39]		\n\t"
			"|  VLDDWM2 *+AR2[41],  VR37:VR36  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR7:VR6, *+AR3[55]      \n\t"
			"|  VLDDWM2 *+AR2[42],  VR39:VR38  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR9:VR8, *+AR3[71]      \n\t"
			"|  VLDDWM2 *+AR2[43],  VR41:VR40  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR11:VR10, *+AR3[87]    \n\t"
			"|  VLDDWM2 *+AR2[48],  VR43:VR42  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR13:VR12, *+AR3[103]   \n\t"
			"|  VLDDWM2 *+AR2[49],  VR45:VR44  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR15:VR14, *+AR3[119]   \n\t"
			"|  VLDDWM2 *+AR2[50],  VR47:VR46  \n\t" //vr16-24,vr17-25
			"|	SSHFLL	  3, R10, R10       \n\t"  //SSHFLL 1  

			"   VSTDW   VR17:VR16, *+AR3[135]   \n\t"
			"|  VLDDWM2 *+AR2[55],  VR49:VR48  \n\t" //vr16-24,vr17-25
			"|  SADDA     R10,AR0,AR1		\n\t" //sadda AVAF 3 SVAF 2

			"   VSTDW   VR19:VR18, *+AR3[151]   \n\t"
			"|  VLDDWM2 *+AR2[56],  VR51:VR50  \n\t" //vr16-24,vr17-25

			"   VSTDW   VR21:VR20, *+AR3[167]   \n\t"
			"|  VLDDWM2 *+AR2[57],  VR53:VR52  \n\t" //vr16-24,vr17-25

			"	SNOP	1\n\t"

			"	VSTDW   VR23:VR22, *+AR3[183]   \n\t"
			"|	VSTDW   VR25:VR24, *+AR3[199]   \n\t"

			"	VSTDW   VR27:VR26, *+AR3[215]   \n\t"
			"|	VSTDW   VR29:VR28, *+AR3[231]   \n\t"

			"	VSTDW   VR31:VR30, *+AR3[247]   \n\t"
			"|  VSTDW   VR33:VR32, *+AR3[263]         \n\t"

			"   VSTDW   VR35:VR34, *+AR3[279]         \n\t"
			"|  VSTDW   VR37:VR36, *+AR3[295]         \n\t"

			"   VSTDW   VR39:VR38, *+AR3[311]         \n\t"
			"|  VSTDW   VR41:VR40, *+AR3[327]         \n\t"

			"   VSTDW   VR43:VR42, *+AR3[343]         \n\t"
			"|  VSTDW   VR45:VR44, *+AR3[359]         \n\t"

			"   VSTDW   VR47:VR46, *+AR3[375]         \n\t"
			"|  VSTDW   VR49:VR48, *+AR3[391]         \n\t"

			"   VSTDW   VR51:VR50, *+AR3[407]         \n\t"
			"|  VSTDW   VR53:VR52, *+AR3[423]         \n\t"

//			"|  SMVAGA.M2 R31,AR3           \n\t" //R31 //index[1] read

			"   VLDDWM2 *-AR1[57],  VR1:VR0    \n\t" //vr0-0,vr1-1
			"|  VLDDWM2 *-AR1[56],  VR3:VR2    \n\t" //vr2-3,vr3-4

			"   VLDDWM2 *-AR1[55],  VR5:VR4    \n\t" //vr4-6,vr5-7
			"|  VLDDWM2 *-AR1[50],  VR7:VR6    \n\t" //vr6-9,vr7-10

			"   VLDDWM2 *-AR1[49],  VR9:VR8    \n\t" //vr4-12,vr5-13
			"|  VLDDWM2 *-AR1[48],  VR11:VR10	\n\t" //vr10-15,vr11-16

			"   VLDDWM2 *-AR1[43],  VR13:VR12  \n\t" //vr12-18,vr13-19
			"|  VLDDWM2 *-AR1[42],  VR15:VR14  \n\t" //vr14-21,vr15-22

			"   VLDDWM2 *-AR1[41],  VR17:VR16  \n\t" //vr16-24,vr17-25
			"|	VLDDWM2	*-AR1[8],	 VR19:VR18	\n\t" //v18-2   

			"   VLDDWM2 *-AR1[7],  VR21:VR20  \n\t" //vr16-24,vr17-25
			"|  VLDDWM2 *-AR1[6],  VR23:VR22  \n\t" //vr16-24,vr17-25

			"   VLDDWM2 *-AR1[1],  VR25:VR24  \n\t" //vr16-24,vr17-25
			"|  VLDDWM2 *+AR1[0],  VR27:VR26  \n\t" //vr16-24,vr17-25

			"   VLDDWM2 *+AR1[1],   VR29:VR28  \n\t" //vr16-24,vr17-25
			"|  VLDDWM2 *+AR1[6],  VR31:VR30  \n\t" //vr16-24,vr17-25

			"	SNOP	1\n\t"

			"   VSTW   VR0, *++AR3[1]         \n\t"  
			"|  VLDDWM2 *+AR1[7],  VR33:VR32  \n\t" //vr16-24,vr17-25
			"|	SLDW	  *+AR15[3], R10	\n\t" //sldw 7

			"   VSTW   VR1, *++AR3[16]      \n\t"
			"|  VLDDWM2 *+AR1[8],  VR35:VR34  \n\t" //vr16-24,vr17-25

			"   VSTW   VR2, *++AR3[16]		\n\t"
			"|  VLDDWM2 *+AR1[41],  VR37:VR36  \n\t" //vr16-24,vr17-25

			"   VSTW   VR3, *++AR3[16]      \n\t"
			"|  VLDDWM2 *+AR1[42],  VR39:VR38  \n\t" //vr16-24,vr17-25

			"   VSTW   VR4, *++AR3[16]      \n\t"
			"|  VLDDWM2 *+AR1[43],  VR41:VR40  \n\t" //vr16-24,vr17-25

			"   VSTW   VR5, *++AR3[16]    \n\t"
			"|  VLDDWM2 *+AR1[48],  VR43:VR42  \n\t" //vr16-24,vr17-25

			"   VSTW   VR6, *++AR3[16]   \n\t"
			"|  VLDDWM2 *+AR1[49],  VR45:VR44  \n\t" //vr16-24,vr17-25

			"   VSTW   VR7, *++AR3[16]   \n\t"
			"|  VLDDWM2 *+AR1[50],  VR47:VR46  \n\t" //vr16-24,vr17-25
			"|	SSHFLL	  3, R10, R10       \n\t"  //SSHFLL 1  

			"   VSTW   VR8, *++AR3[16]   \n\t"
			"|  VLDDWM2 *+AR1[55],  VR49:VR48  \n\t" //vr16-24,vr17-25
			"|  SADDA     R10,AR0,AR2		\n\t" //sadda AVAF 3 SVAF 2

			"   VSTW   VR9, *++AR3[16]   \n\t"
			"|  VLDDWM2 *+AR1[56],  VR51:VR50  \n\t" //vr16-24,vr17-25

			"   VSTW   VR10, *++AR3[16]   \n\t"
			"|  VLDDWM2 *+AR1[57],  VR53:VR52  \n\t" //vr16-24,vr17-25

			"   VSTW   VR11, *+AR3[16]   \n\t"
			"|  VSTW   VR12, *+AR3[32]   \n\t"

			"   VSTW   VR13, *+AR3[48]   \n\t"
			"|  VSTW   VR14, *+AR3[64]   \n\t"

			"   VSTW   VR15, *+AR3[80]   \n\t"
			"|  VSTW   VR16, *+AR3[96]   \n\t"

			"   VSTW   VR17, *+AR3[112]   \n\t"
			"|  VSTW   VR18, *+AR3[128]   \n\t"

			"   VSTW   VR19, *+AR3[144]   \n\t"
			"|  VSTW   VR20, *+AR3[160]   \n\t"

			"   VSTW   VR21, *+AR3[176]   \n\t"
			"|  VSTW   VR22, *+AR3[192]   \n\t"

			"   VSTW   VR23, *+AR3[208]   \n\t"
			"|  VSTW   VR24, *+AR3[224]   \n\t"

			"   VSTW   VR25, *+AR3[240]   \n\t"
			"|  VSTW   VR26, *+AR3[256]   \n\t"

			"   VSTW   VR27, *+AR3[272]   \n\t"
			"|  VSTW   VR28, *+AR3[288]   \n\t"
			"|  SADDA     R10,AR0,AR1		\n\t" //sadda AVAF 3 SVAF 2

			"   VSTW   VR29, *+AR3[304]   \n\t"
			"|  VSTW   VR30, *+AR3[320]   \n\t"

			"   VSTW   VR31, *+AR3[336]   \n\t"
			"|  VSTW   VR32, *+AR3[352]   \n\t"

			"   VSTW   VR33, *+AR3[368]   \n\t"
			"|  VSTW   VR34, *+AR3[384]   \n\t"

			"   VSTW   VR35, *+AR3[400]   \n\t"
			"|  VSTW   VR36, *+AR3[416]   \n\t"

			"   VSTW   VR37, *++AR3[432]   \n\t"
			"|  SMVAGA.M2 R31,AR4           \n\t" //R31 ////////////////////////////////////////////////////////////////////////////
// index[3] read
			"	VSTW   VR38, *++AR3[16]   \n\t"
			"|  VLDDWM2 *-AR1[57],  VR1:VR0    \n\t" //vr0-0,vr1-1

			"	VSTW   VR39, *++AR3[16]   \n\t"
			"|  VLDDWM2 *-AR1[56],  VR3:VR2    \n\t" //vr2-3,vr3-4
			
			"	VSTW   VR40, *++AR3[16]   \n\t"
			"|  VLDDWM2 *-AR1[55],  VR5:VR4    \n\t" //vr4-6,vr5-7

			"	VSTW   VR41, *++AR3[16]   \n\t"
			"|  VLDDWM2 *-AR1[50],  VR7:VR6    \n\t" //vr6-9,vr7-10

			"	VSTW   VR42, *++AR3[16]   \n\t"
			"|  VLDDWM2 *-AR1[49],  VR9:VR8    \n\t" //vr4-12,vr5-13

			"   VSTW   VR43, *++AR3[16]         \n\t"
			"|  VLDDWM2 *-AR1[48],  VR11:VR10	\n\t" //vr10-15,vr11-16

			"   VSTW   VR44, *++AR3[16]         \n\t"
			"|  VLDDWM2 *-AR1[43],  VR13:VR12  \n\t" //vr12-18,vr13-19

			"   VSTW   VR45, *++AR3[16]         \n\t"
			"|  VLDDWM2 *-AR1[42],  VR15:VR14  \n\t" //vr14-21,vr15-22

			"   VSTW   VR46, *++AR3[16]         \n\t"
			"|  VLDDWM2 *-AR1[41],  VR17:VR16  \n\t" //vr16-24,vr17-25

			"   VSTW   VR47, *++AR3[16]         \n\t"
			"|	VLDDWM2	*-AR1[8],	 VR19:VR18	\n\t" //v18-2   

			"   VSTW   VR48, *++AR3[16]         \n\t"
			"|  VLDDWM2 *-AR1[7],  VR21:VR20  \n\t" //vr16-24,vr17-25

			"   VSTW   VR49, *++AR3[16]         \n\t"
			"|  VLDDWM2 *-AR1[6],  VR23:VR22  \n\t" //vr16-24,vr17-25

			"   VSTW   VR50, *++AR3[16]         \n\t"
			"|  VLDDWM2 *-AR1[1],  VR25:VR24  \n\t" //vr16-24,vr17-25

			"   VSTW   VR51, *++AR3[16]         \n\t"
			"|  VLDDWM2 *+AR1[0],  VR27:VR26  \n\t" //vr16-24,vr17-25

			"   VSTW   VR52, *++AR3[16]         \n\t"
			"|  VLDDWM2 *+AR1[1],   VR29:VR28  \n\t" //vr16-24,vr17-25

			"   VSTW   VR53, *++AR3[16]         \n\t"
			"|  VLDDWM2 *+AR1[6],  VR31:VR30  \n\t" //vr16-24,vr17-25

			"   VSTW   VR0, *++AR4[3]         \n\t"
			"|  VLDDWM2 *+AR1[7],  VR33:VR32  \n\t" //vr16-24,vr17-25
			"|	SLDW	  *+AR15[5], R10	\n\t" //sldw 7

			"   VSTW   VR1, *++AR4[16]      \n\t"
			"|	VLDDWM2 *+AR1[8],  VR35:VR34  \n\t" //vr16-24,vr17-25

			"   VSTW   VR2, *++AR4[16]		\n\t"
			"|  VLDDWM2 *+AR1[41],  VR37:VR36  \n\t" //vr16-24,vr17-25

			"   VSTW   VR3, *++AR4[16]      \n\t"
			"|  VLDDWM2 *+AR1[42],  VR39:VR38  \n\t" //vr16-24,vr17-25

			"   VSTW   VR4, *++AR4[16]      \n\t"
			"|  VLDDWM2 *+AR1[43],  VR41:VR40  \n\t" //vr16-24,vr17-25

			"   VSTW   VR5, *++AR4[16]    \n\t"
			"|  VLDDWM2 *+AR1[48],  VR43:VR42  \n\t" //vr16-24,vr17-25

			"   VSTW   VR6, *++AR4[16]   \n\t"
			"|  VLDDWM2 *+AR1[49],  VR45:VR44  \n\t" //vr16-24,vr17-25

			"   VSTW   VR7, *++AR4[16]   \n\t"
			"|  VLDDWM2 *+AR1[50],  VR47:VR46  \n\t" //vr16-24,vr17-25
			"|	SSHFLL	  3, R10, R10       \n\t"  //SSHFLL 1  

			"   VSTW   VR8, *++AR4[16]   \n\t"
			"|  VLDDWM2 *+AR1[55],  VR49:VR48  \n\t" //vr16-24,vr17-25

			"   VSTW   VR9, *++AR4[16]   \n\t"
			"|  VLDDWM2 *+AR1[56],  VR51:VR50  \n\t" //vr16-24,vr17-25

			"   VSTW   VR10, *++AR4[16]   \n\t"
			"|  VLDDWM2 *+AR1[57],  VR53:VR52  \n\t" //vr16-24,vr17-25

			"   VSTW   VR11, *+AR4[16]   \n\t"
			"|  VSTW   VR12, *+AR4[32]   \n\t"

			"   VSTW   VR13, *+AR4[48]   \n\t"
			"|  VSTW   VR14, *+AR4[64]   \n\t"

			"   VSTW   VR15, *+AR4[80]   \n\t"
			"|  VSTW   VR16, *+AR4[96]   \n\t"

			"   VSTW   VR17, *+AR4[112]   \n\t"
			"|  VSTW   VR18, *+AR4[128]   \n\t"

			"   VSTW   VR19, *+AR4[144]   \n\t"
			"|  VSTW   VR20, *+AR4[160]   \n\t"

			"   VSTW   VR21, *+AR4[176]   \n\t"
			"|  VSTW   VR22, *+AR4[192]   \n\t"

			"   VSTW   VR23, *+AR4[208]   \n\t"
			"|  VSTW   VR24, *+AR4[224]   \n\t"

			"   VSTW   VR25, *+AR4[240]   \n\t"
			"|  VSTW   VR26, *+AR4[256]   \n\t"

			"   VSTW   VR27, *+AR4[272]   \n\t"
			"|  VSTW   VR28, *+AR4[288]   \n\t"
			"|  SADDA     R10,AR0,AR1		\n\t" //sadda AVAF 3 SVAF 2

			"   VSTW   VR29, *+AR4[304]   \n\t"
			"|  VSTW   VR30, *+AR4[320]   \n\t"

			"   VSTW   VR31, *+AR4[336]   \n\t"
			"|  VSTW   VR32, *+AR4[352]   \n\t"

			"   VSTW   VR33, *+AR4[368]   \n\t"
			"|  VSTW   VR34, *+AR4[384]   \n\t"

			"   VSTW   VR35, *+AR4[400]   \n\t"
			"|  VSTW   VR36, *+AR4[416]   \n\t"

			"   VSTW   VR37, *++AR4[432]   \n\t"
			"|  SMVAGA.M2 R31,AR3           \n\t" //R31 // index[5] read
			"	VSTW   VR38, *++AR4[16]   \n\t"
			"|  VLDDWM2 *-AR1[57],  VR1:VR0    \n\t" //vr0-0,vr1-1

			"	VSTW   VR39, *++AR4[16]   \n\t"
			"|  VLDDWM2 *-AR1[56],  VR3:VR2    \n\t" //vr2-3,vr3-4
			
			"	VSTW   VR40, *++AR4[16]   \n\t"
			"|  VLDDWM2 *-AR1[55],  VR5:VR4    \n\t" //vr4-6,vr5-7

			"	VSTW   VR41, *++AR4[16]   \n\t"
			"|  VLDDWM2 *-AR1[50],  VR7:VR6    \n\t" //vr6-9,vr7-10

			"	VSTW   VR42, *++AR4[16]   \n\t"
			"|  VLDDWM2 *-AR1[49],  VR9:VR8    \n\t" //vr4-12,vr5-13

			"   VSTW   VR43, *++AR4[16]         \n\t"
			"|  VLDDWM2 *-AR1[48],  VR11:VR10	\n\t" //vr10-15,vr11-16

			"   VSTW   VR44, *++AR4[16]         \n\t"
			"|  VLDDWM2 *-AR1[43],  VR13:VR12  \n\t" //vr12-18,vr13-19

			"   VSTW   VR45, *++AR4[16]         \n\t"
			"|  VLDDWM2 *-AR1[42],  VR15:VR14  \n\t" //vr14-21,vr15-22

			"   VSTW   VR46, *++AR4[16]         \n\t"
			"|  VLDDWM2 *-AR1[41],  VR17:VR16  \n\t" //vr16-24,vr17-25

			"   VSTW   VR47, *++AR4[16]         \n\t"
			"|	VLDDWM2	*-AR1[8],	 VR19:VR18	\n\t" //v18-2   

			"   VSTW   VR48, *++AR4[16]         \n\t"
			"|  VLDDWM2 *-AR1[7],  VR21:VR20  \n\t" //vr16-24,vr17-25

			"   VSTW   VR49, *++AR4[16]         \n\t"
			"|  VLDDWM2 *-AR1[6],  VR23:VR22  \n\t" //vr16-24,vr17-25

			"   VSTW   VR50, *++AR4[16]         \n\t"
			"|  VLDDWM2 *-AR1[1],  VR25:VR24  \n\t" //vr16-24,vr17-25

			"   VSTW   VR51, *++AR4[16]         \n\t"
			"|  VLDDWM2 *+AR1[0],  VR27:VR26  \n\t" //vr16-24,vr17-25

			"   VSTW   VR52, *++AR4[16]         \n\t"
			"|  VLDDWM2 *+AR1[1],   VR29:VR28  \n\t" //vr16-24,vr17-25

			"   VSTW   VR53, *++AR4[16]         \n\t"
			"|  VLDDWM2 *+AR1[6],  VR31:VR30  \n\t" //vr16-24,vr17-25

			"   VSTW   VR0, *++AR3[5]         \n\t"
			"|  VLDDWM2 *+AR1[7],  VR33:VR32  \n\t" //vr16-24,vr17-25
			"|	SLDW	  *+AR15[7], R10	\n\t" //sldw 7

			"   VSTW   VR1, *++AR3[16]      \n\t"
			"|	VLDDWM2 *+AR1[8],  VR35:VR34  \n\t" //vr16-24,vr17-25

			"   VSTW   VR2, *++AR3[16]		\n\t"
			"|  VLDDWM2 *+AR1[41],  VR37:VR36  \n\t" //vr16-24,vr17-25

			"   VSTW   VR3, *++AR3[16]      \n\t"
			"|  VLDDWM2 *+AR1[42],  VR39:VR38  \n\t" //vr16-24,vr17-25

			"   VSTW   VR4, *++AR3[16]      \n\t"
			"|  VLDDWM2 *+AR1[43],  VR41:VR40  \n\t" //vr16-24,vr17-25

			"   VSTW   VR5, *++AR3[16]    \n\t"
			"|  VLDDWM2 *+AR1[48],  VR43:VR42  \n\t" //vr16-24,vr17-25

			"   VSTW   VR6, *++AR3[16]   \n\t"
			"|  VLDDWM2 *+AR1[49],  VR45:VR44  \n\t" //vr16-24,vr17-25

			"   VSTW   VR7, *++AR3[16]   \n\t"
			"|  VLDDWM2 *+AR1[50],  VR47:VR46  \n\t" //vr16-24,vr17-25
			"|	SSHFLL	  3, R10, R10       \n\t"  //SSHFLL 1  

			"   VSTW   VR8, *++AR3[16]   \n\t"
			"|  VLDDWM2 *+AR1[55],  VR49:VR48  \n\t" //vr16-24,vr17-25

			"   VSTW   VR9, *++AR3[16]   \n\t"
			"|  VLDDWM2 *+AR1[56],  VR51:VR50  \n\t" //vr16-24,vr17-25

			"   VSTW   VR10, *++AR3[16]   \n\t"
			"|  VLDDWM2 *+AR1[57],  VR53:VR52  \n\t" //vr16-24,vr17-25

			"   VSTW   VR11, *+AR3[16]   \n\t"
			"|  VSTW   VR12, *+AR3[32]   \n\t"

			"   VSTW   VR13, *+AR3[48]   \n\t"
			"|  VSTW   VR14, *+AR3[64]   \n\t"

			"   VSTW   VR15, *+AR3[80]   \n\t"
			"|  VSTW   VR16, *+AR3[96]   \n\t"

			"   VSTW   VR17, *+AR3[112]   \n\t"
			"|  VSTW   VR18, *+AR3[128]   \n\t"

			"   VSTW   VR19, *+AR3[144]   \n\t"
			"|  VSTW   VR20, *+AR3[160]   \n\t"

			"   VSTW   VR21, *+AR3[176]   \n\t"
			"|  VSTW   VR22, *+AR3[192]   \n\t"

			"   VSTW   VR23, *+AR3[208]   \n\t"
			"|  VSTW   VR24, *+AR3[224]   \n\t"

			"   VSTW   VR25, *+AR3[240]   \n\t"
			"|  VSTW   VR26, *+AR3[256]   \n\t"

			"   VSTW   VR27, *+AR3[272]   \n\t"
			"|  VSTW   VR28, *+AR3[288]   \n\t"
			"|  SADDA     R10,AR0,AR1		\n\t" //sadda AVAF 3 SVAF 2

			"   VSTW   VR29, *+AR3[304]   \n\t"
			"|  VSTW   VR30, *+AR3[320]   \n\t"

			"   VSTW   VR31, *+AR3[336]   \n\t"
			"|  VSTW   VR32, *+AR3[352]   \n\t"

			"   VSTW   VR33, *+AR3[368]   \n\t"
			"|  VSTW   VR34, *+AR3[384]   \n\t"

			"   VSTW   VR35, *+AR3[400]   \n\t"
			"|  VSTW   VR36, *+AR3[416]   \n\t"

			"   VSTW   VR37, *++AR3[432]   \n\t"
			"|  SMVAGA.M2 R31,AR4           \n\t" //R31 
// index[7] read
			"	VSTW   VR38, *++AR3[16]   \n\t"
			"|  VLDDWM2 *-AR1[57],  VR1:VR0    \n\t" //vr0-0,vr1-1

			"	VSTW   VR39, *++AR3[16]   \n\t"
			"|  VLDDWM2 *-AR1[56],  VR3:VR2    \n\t" //vr2-3,vr3-4
			
			"	VSTW   VR40, *++AR3[16]   \n\t"
			"|  VLDDWM2 *-AR1[55],  VR5:VR4    \n\t" //vr4-6,vr5-7

			"	VSTW   VR41, *++AR3[16]   \n\t"
			"|  VLDDWM2 *-AR1[50],  VR7:VR6    \n\t" //vr6-9,vr7-10

			"	VSTW   VR42, *++AR3[16]   \n\t"
			"|  VLDDWM2 *-AR1[49],  VR9:VR8    \n\t" //vr4-12,vr5-13

			"   VSTW   VR43, *++AR3[16]         \n\t"
			"|  VLDDWM2 *-AR1[48],  VR11:VR10	\n\t" //vr10-15,vr11-16

			"   VSTW   VR44, *++AR3[16]         \n\t"
			"|  VLDDWM2 *-AR1[43],  VR13:VR12  \n\t" //vr12-18,vr13-19

			"   VSTW   VR45, *++AR3[16]         \n\t"
			"|  VLDDWM2 *-AR1[42],  VR15:VR14  \n\t" //vr14-21,vr15-22

			"   VSTW   VR46, *++AR3[16]         \n\t"
			"|  VLDDWM2 *-AR1[41],  VR17:VR16  \n\t" //vr16-24,vr17-25

			"   VSTW   VR47, *++AR3[16]         \n\t"
			"|	VLDDWM2	*-AR1[8],	 VR19:VR18	\n\t" //v18-2   

			"   VSTW   VR48, *++AR3[16]         \n\t"
			"|  VLDDWM2 *-AR1[7],  VR21:VR20  \n\t" //vr16-24,vr17-25

			"   VSTW   VR49, *++AR3[16]         \n\t"
			"|  VLDDWM2 *-AR1[6],  VR23:VR22  \n\t" //vr16-24,vr17-25

			"   VSTW   VR50, *++AR3[16]         \n\t"
			"|  VLDDWM2 *-AR1[1],  VR25:VR24  \n\t" //vr16-24,vr17-25

			"   VSTW   VR51, *++AR3[16]         \n\t"
			"|  VLDDWM2 *+AR1[0],  VR27:VR26  \n\t" //vr16-24,vr17-25

			"   VSTW   VR52, *++AR3[16]         \n\t"
			"|  VLDDWM2 *+AR1[1],   VR29:VR28  \n\t" //vr16-24,vr17-25

			"   VSTW   VR53, *++AR3[16]         \n\t"
			"|  VLDDWM2 *+AR1[6],  VR31:VR30  \n\t" //vr16-24,vr17-25

			"   VSTW   VR0, *++AR4[7]         \n\t"
			"|  VLDDWM2 *+AR1[7],  VR33:VR32  \n\t" //vr16-24,vr17-25
			"|	SLDW	  *+AR15[9], R10	\n\t" //sldw 7

			"   VSTW   VR1, *++AR4[16]      \n\t"
			"|	VLDDWM2 *+AR1[8],  VR35:VR34  \n\t" //vr16-24,vr17-25

			"   VSTW   VR2, *++AR4[16]		\n\t"
			"|  VLDDWM2 *+AR1[41],  VR37:VR36  \n\t" //vr16-24,vr17-25

			"   VSTW   VR3, *++AR4[16]      \n\t"
			"|  VLDDWM2 *+AR1[42],  VR39:VR38  \n\t" //vr16-24,vr17-25

			"   VSTW   VR4, *++AR4[16]      \n\t"
			"|  VLDDWM2 *+AR1[43],  VR41:VR40  \n\t" //vr16-24,vr17-25

			"   VSTW   VR5, *++AR4[16]    \n\t"
			"|  VLDDWM2 *+AR1[48],  VR43:VR42  \n\t" //vr16-24,vr17-25

			"   VSTW   VR6, *++AR4[16]   \n\t"
			"|  VLDDWM2 *+AR1[49],  VR45:VR44  \n\t" //vr16-24,vr17-25

			"   VSTW   VR7, *++AR4[16]   \n\t"
			"|  VLDDWM2 *+AR1[50],  VR47:VR46  \n\t" //vr16-24,vr17-25
			"|	SSHFLL	  3, R10, R10       \n\t"  //SSHFLL 1  

			"   VSTW   VR8, *++AR4[16]   \n\t"
			"|  VLDDWM2 *+AR1[55],  VR49:VR48  \n\t" //vr16-24,vr17-25

			"   VSTW   VR9, *++AR4[16]   \n\t"
			"|  VLDDWM2 *+AR1[56],  VR51:VR50  \n\t" //vr16-24,vr17-25

			"   VSTW   VR10, *++AR4[16]   \n\t"
			"|  VLDDWM2 *+AR1[57],  VR53:VR52  \n\t" //vr16-24,vr17-25

			"   VSTW   VR11, *+AR4[16]   \n\t"
			"|  VSTW   VR12, *+AR4[32]   \n\t"

			"   VSTW   VR13, *+AR4[48]   \n\t"
			"|  VSTW   VR14, *+AR4[64]   \n\t"

			"   VSTW   VR15, *+AR4[80]   \n\t"
			"|  VSTW   VR16, *+AR4[96]   \n\t"

			"   VSTW   VR17, *+AR4[112]   \n\t"
			"|  VSTW   VR18, *+AR4[128]   \n\t"

			"   VSTW   VR19, *+AR4[144]   \n\t"
			"|  VSTW   VR20, *+AR4[160]   \n\t"

			"   VSTW   VR21, *+AR4[176]   \n\t"
			"|  VSTW   VR22, *+AR4[192]   \n\t"

			"   VSTW   VR23, *+AR4[208]   \n\t"
			"|  VSTW   VR24, *+AR4[224]   \n\t"

			"   VSTW   VR25, *+AR4[240]   \n\t"
			"|  VSTW   VR26, *+AR4[256]   \n\t"

			"   VSTW   VR27, *+AR4[272]   \n\t"
			"|  VSTW   VR28, *+AR4[288]   \n\t"
			"|  SADDA     R10,AR0,AR1		\n\t" //sadda AVAF 3 SVAF 2

			"   VSTW   VR29, *+AR4[304]   \n\t"
			"|  VSTW   VR30, *+AR4[320]   \n\t"

			"   VSTW   VR31, *+AR4[336]   \n\t"
			"|  VSTW   VR32, *+AR4[352]   \n\t"

			"   VSTW   VR33, *+AR4[368]   \n\t"
			"|  VSTW   VR34, *+AR4[384]   \n\t"

			"   VSTW   VR35, *+AR4[400]   \n\t"
			"|  VSTW   VR36, *+AR4[416]   \n\t"

			"   VSTW   VR37, *++AR4[432]   \n\t"
			"|  SMVAGA.M2 R31,AR3           \n\t" //R31 
// index[9] read
			"	VSTW   VR38, *++AR4[16]   \n\t"
			"|  VLDDWM2 *-AR1[57],  VR1:VR0    \n\t" //vr0-0,vr1-1

			"	VSTW   VR39, *++AR4[16]   \n\t"
			"|  VLDDWM2 *-AR1[56],  VR3:VR2    \n\t" //vr2-3,vr3-4
			
			"	VSTW   VR40, *++AR4[16]   \n\t"
			"|  VLDDWM2 *-AR1[55],  VR5:VR4    \n\t" //vr4-6,vr5-7

			"	VSTW   VR41, *++AR4[16]   \n\t"
			"|  VLDDWM2 *-AR1[50],  VR7:VR6    \n\t" //vr6-9,vr7-10

			"	VSTW   VR42, *++AR4[16]   \n\t"
			"|  VLDDWM2 *-AR1[49],  VR9:VR8    \n\t" //vr4-12,vr5-13

			"   VSTW   VR43, *++AR4[16]         \n\t"
			"|  VLDDWM2 *-AR1[48],  VR11:VR10	\n\t" //vr10-15,vr11-16

			"   VSTW   VR44, *++AR4[16]         \n\t"
			"|  VLDDWM2 *-AR1[43],  VR13:VR12  \n\t" //vr12-18,vr13-19

			"   VSTW   VR45, *++AR4[16]         \n\t"
			"|  VLDDWM2 *-AR1[42],  VR15:VR14  \n\t" //vr14-21,vr15-22

			"   VSTW   VR46, *++AR4[16]         \n\t"
			"|  VLDDWM2 *-AR1[41],  VR17:VR16  \n\t" //vr16-24,vr17-25

			"   VSTW   VR47, *++AR4[16]         \n\t"
			"|	VLDDWM2	*-AR1[8],	 VR19:VR18	\n\t" //v18-2   

			"   VSTW   VR48, *++AR4[16]         \n\t"
			"|  VLDDWM2 *-AR1[7],  VR21:VR20  \n\t" //vr16-24,vr17-25

			"   VSTW   VR49, *++AR4[16]         \n\t"
			"|  VLDDWM2 *-AR1[6],  VR23:VR22  \n\t" //vr16-24,vr17-25

			"   VSTW   VR50, *++AR4[16]         \n\t"
			"|  VLDDWM2 *-AR1[1],  VR25:VR24  \n\t" //vr16-24,vr17-25

			"   VSTW   VR51, *++AR4[16]         \n\t"
			"|  VLDDWM2 *+AR1[0],  VR27:VR26  \n\t" //vr16-24,vr17-25

			"   VSTW   VR52, *++AR4[16]         \n\t"
			"|  VLDDWM2 *+AR1[1],   VR29:VR28  \n\t" //vr16-24,vr17-25

			"   VSTW   VR53, *++AR4[16]         \n\t"
			"|  VLDDWM2 *+AR1[6],  VR31:VR30  \n\t" //vr16-24,vr17-25

			"   VSTW   VR0, *++AR3[9]         \n\t"
			"|  VLDDWM2 *+AR1[7],  VR33:VR32  \n\t" //vr16-24,vr17-25
			"|	SLDW	  *+AR15[11], R10	\n\t" //sldw 7

			"   VSTW   VR1, *++AR3[16]      \n\t"
			"|	VLDDWM2 *+AR1[8],  VR35:VR34  \n\t" //vr16-24,vr17-25

			"   VSTW   VR2, *++AR3[16]		\n\t"
			"|  VLDDWM2 *+AR1[41],  VR37:VR36  \n\t" //vr16-24,vr17-25

			"   VSTW   VR3, *++AR3[16]      \n\t"
			"|  VLDDWM2 *+AR1[42],  VR39:VR38  \n\t" //vr16-24,vr17-25

			"   VSTW   VR4, *++AR3[16]      \n\t"
			"|  VLDDWM2 *+AR1[43],  VR41:VR40  \n\t" //vr16-24,vr17-25

			"   VSTW   VR5, *++AR3[16]    \n\t"
			"|  VLDDWM2 *+AR1[48],  VR43:VR42  \n\t" //vr16-24,vr17-25

			"   VSTW   VR6, *++AR3[16]   \n\t"
			"|  VLDDWM2 *+AR1[49],  VR45:VR44  \n\t" //vr16-24,vr17-25

			"   VSTW   VR7, *++AR3[16]   \n\t"
			"|  VLDDWM2 *+AR1[50],  VR47:VR46  \n\t" //vr16-24,vr17-25
			"|	SSHFLL	  3, R10, R10       \n\t"  //SSHFLL 1  

			"   VSTW   VR8, *++AR3[16]   \n\t"
			"|  VLDDWM2 *+AR1[55],  VR49:VR48  \n\t" //vr16-24,vr17-25

			"   VSTW   VR9, *++AR3[16]   \n\t"
			"|  VLDDWM2 *+AR1[56],  VR51:VR50  \n\t" //vr16-24,vr17-25

			"   VSTW   VR10, *++AR3[16]   \n\t"
			"|  VLDDWM2 *+AR1[57],  VR53:VR52  \n\t" //vr16-24,vr17-25

			"   VSTW   VR11, *+AR3[16]   \n\t"
			"|  VSTW   VR12, *+AR3[32]   \n\t"

			"   VSTW   VR13, *+AR3[48]   \n\t"
			"|  VSTW   VR14, *+AR3[64]   \n\t"

			"   VSTW   VR15, *+AR3[80]   \n\t"
			"|  VSTW   VR16, *+AR3[96]   \n\t"

			"   VSTW   VR17, *+AR3[112]   \n\t"
			"|  VSTW   VR18, *+AR3[128]   \n\t"

			"   VSTW   VR19, *+AR3[144]   \n\t"
			"|  VSTW   VR20, *+AR3[160]   \n\t"

			"   VSTW   VR21, *+AR3[176]   \n\t"
			"|  VSTW   VR22, *+AR3[192]   \n\t"

			"   VSTW   VR23, *+AR3[208]   \n\t"
			"|  VSTW   VR24, *+AR3[224]   \n\t"

			"   VSTW   VR25, *+AR3[240]   \n\t"
			"|  VSTW   VR26, *+AR3[256]   \n\t"

			"   VSTW   VR27, *+AR3[272]   \n\t"
			"|  VSTW   VR28, *+AR3[288]   \n\t"
			"|  SADDA     R10,AR0,AR1		\n\t" //sadda AVAF 3 SVAF 2

			"   VSTW   VR29, *+AR3[304]   \n\t"
			"|  VSTW   VR30, *+AR3[320]   \n\t"

			"   VSTW   VR31, *+AR3[336]   \n\t"
			"|  VSTW   VR32, *+AR3[352]   \n\t"

			"   VSTW   VR33, *+AR3[368]   \n\t"
			"|  VSTW   VR34, *+AR3[384]   \n\t"

			"   VSTW   VR35, *+AR3[400]   \n\t"
			"|  VSTW   VR36, *+AR3[416]   \n\t"

			"   VSTW   VR37, *++AR3[432]   \n\t"
			"|  SMVAGA.M2 R31,AR4           \n\t" //R31 
// index[11] read
			"	VSTW   VR38, *++AR3[16]   \n\t"
			"|  VLDDWM2 *-AR1[57],  VR1:VR0    \n\t" //vr0-0,vr1-1

			"	VSTW   VR39, *++AR3[16]   \n\t"
			"|  VLDDWM2 *-AR1[56],  VR3:VR2    \n\t" //vr2-3,vr3-4
			
			"	VSTW   VR40, *++AR3[16]   \n\t"
			"|  VLDDWM2 *-AR1[55],  VR5:VR4    \n\t" //vr4-6,vr5-7

			"	VSTW   VR41, *++AR3[16]   \n\t"
			"|  VLDDWM2 *-AR1[50],  VR7:VR6    \n\t" //vr6-9,vr7-10

			"	VSTW   VR42, *++AR3[16]   \n\t"
			"|  VLDDWM2 *-AR1[49],  VR9:VR8    \n\t" //vr4-12,vr5-13

			"   VSTW   VR43, *++AR3[16]         \n\t"
			"|  VLDDWM2 *-AR1[48],  VR11:VR10	\n\t" //vr10-15,vr11-16

			"   VSTW   VR44, *++AR3[16]         \n\t"
			"|  VLDDWM2 *-AR1[43],  VR13:VR12  \n\t" //vr12-18,vr13-19

			"   VSTW   VR45, *++AR3[16]         \n\t"
			"|  VLDDWM2 *-AR1[42],  VR15:VR14  \n\t" //vr14-21,vr15-22

			"   VSTW   VR46, *++AR3[16]         \n\t"
			"|  VLDDWM2 *-AR1[41],  VR17:VR16  \n\t" //vr16-24,vr17-25

			"   VSTW   VR47, *++AR3[16]         \n\t"
			"|	VLDDWM2	*-AR1[8],	 VR19:VR18	\n\t" //v18-2   

			"   VSTW   VR48, *++AR3[16]         \n\t"
			"|  VLDDWM2 *-AR1[7],  VR21:VR20  \n\t" //vr16-24,vr17-25

			"   VSTW   VR49, *++AR3[16]         \n\t"
			"|  VLDDWM2 *-AR1[6],  VR23:VR22  \n\t" //vr16-24,vr17-25

			"   VSTW   VR50, *++AR3[16]         \n\t"
			"|  VLDDWM2 *-AR1[1],  VR25:VR24  \n\t" //vr16-24,vr17-25

			"   VSTW   VR51, *++AR3[16]         \n\t"
			"|  VLDDWM2 *+AR1[0],  VR27:VR26  \n\t" //vr16-24,vr17-25

			"   VSTW   VR52, *++AR3[16]         \n\t"
			"|  VLDDWM2 *+AR1[1],   VR29:VR28  \n\t" //vr16-24,vr17-25

			"   VSTW   VR53, *++AR3[16]         \n\t"
			"|  VLDDWM2 *+AR1[6],  VR31:VR30  \n\t" //vr16-24,vr17-25

			"   VSTW   VR0, *++AR4[11]         \n\t"
			"|  VLDDWM2 *+AR1[7],  VR33:VR32  \n\t" //vr16-24,vr17-25
			"|	SLDW	  *+AR15[13], R10	\n\t" //sldw 7

			"   VSTW   VR1, *++AR4[16]      \n\t"
			"|	VLDDWM2 *+AR1[8],  VR35:VR34  \n\t" //vr16-24,vr17-25

			"   VSTW   VR2, *++AR4[16]		\n\t"
			"|  VLDDWM2 *+AR1[41],  VR37:VR36  \n\t" //vr16-24,vr17-25

			"   VSTW   VR3, *++AR4[16]      \n\t"
			"|  VLDDWM2 *+AR1[42],  VR39:VR38  \n\t" //vr16-24,vr17-25

			"   VSTW   VR4, *++AR4[16]      \n\t"
			"|  VLDDWM2 *+AR1[43],  VR41:VR40  \n\t" //vr16-24,vr17-25

			"   VSTW   VR5, *++AR4[16]    \n\t"
			"|  VLDDWM2 *+AR1[48],  VR43:VR42  \n\t" //vr16-24,vr17-25

			"   VSTW   VR6, *++AR4[16]   \n\t"
			"|  VLDDWM2 *+AR1[49],  VR45:VR44  \n\t" //vr16-24,vr17-25

			"   VSTW   VR7, *++AR4[16]   \n\t"
			"|  VLDDWM2 *+AR1[50],  VR47:VR46  \n\t" //vr16-24,vr17-25
			"|	SSHFLL	  3, R10, R10       \n\t"  //SSHFLL 1  

			"   VSTW   VR8, *++AR4[16]   \n\t"
			"|  VLDDWM2 *+AR1[55],  VR49:VR48  \n\t" //vr16-24,vr17-25

			"   VSTW   VR9, *++AR4[16]   \n\t"
			"|  VLDDWM2 *+AR1[56],  VR51:VR50  \n\t" //vr16-24,vr17-25

			"   VSTW   VR10, *++AR4[16]   \n\t"
			"|  VLDDWM2 *+AR1[57],  VR53:VR52  \n\t" //vr16-24,vr17-25

			"   VSTW   VR11, *+AR4[16]   \n\t"
			"|  VSTW   VR12, *+AR4[32]   \n\t"

			"   VSTW   VR13, *+AR4[48]   \n\t"
			"|  VSTW   VR14, *+AR4[64]   \n\t"

			"   VSTW   VR15, *+AR4[80]   \n\t"
			"|  VSTW   VR16, *+AR4[96]   \n\t"

			"   VSTW   VR17, *+AR4[112]   \n\t"
			"|  VSTW   VR18, *+AR4[128]   \n\t"

			"   VSTW   VR19, *+AR4[144]   \n\t"
			"|  VSTW   VR20, *+AR4[160]   \n\t"

			"   VSTW   VR21, *+AR4[176]   \n\t"
			"|  VSTW   VR22, *+AR4[192]   \n\t"

			"   VSTW   VR23, *+AR4[208]   \n\t"
			"|  VSTW   VR24, *+AR4[224]   \n\t"

			"   VSTW   VR25, *+AR4[240]   \n\t"
			"|  VSTW   VR26, *+AR4[256]   \n\t"

			"   VSTW   VR27, *+AR4[272]   \n\t"
			"|  VSTW   VR28, *+AR4[288]   \n\t"
			"|  SADDA     R10,AR0,AR1		\n\t" //sadda AVAF 3 SVAF 2

			"   VSTW   VR29, *+AR4[304]   \n\t"
			"|  VSTW   VR30, *+AR4[320]   \n\t"

			"   VSTW   VR31, *+AR4[336]   \n\t"
			"|  VSTW   VR32, *+AR4[352]   \n\t"

			"   VSTW   VR33, *+AR4[368]   \n\t"
			"|  VSTW   VR34, *+AR4[384]   \n\t"

			"   VSTW   VR35, *+AR4[400]   \n\t"
			"|  VSTW   VR36, *+AR4[416]   \n\t"

			"   VSTW   VR37, *++AR4[432]   \n\t"
			"|  SMVAGA.M2 R31,AR3           \n\t" //R31 
// index[13] read
			"	VSTW   VR38, *++AR4[16]   \n\t"
			"|  VLDDWM2 *-AR1[57],  VR1:VR0    \n\t" //vr0-0,vr1-1

			"	VSTW   VR39, *++AR4[16]   \n\t"
			"|  VLDDWM2 *-AR1[56],  VR3:VR2    \n\t" //vr2-3,vr3-4
			
			"	VSTW   VR40, *++AR4[16]   \n\t"
			"|  VLDDWM2 *-AR1[55],  VR5:VR4    \n\t" //vr4-6,vr5-7

			"	VSTW   VR41, *++AR4[16]   \n\t"
			"|  VLDDWM2 *-AR1[50],  VR7:VR6    \n\t" //vr6-9,vr7-10

			"	VSTW   VR42, *++AR4[16]   \n\t"
			"|  VLDDWM2 *-AR1[49],  VR9:VR8    \n\t" //vr4-12,vr5-13

			"   VSTW   VR43, *++AR4[16]         \n\t"
			"|  VLDDWM2 *-AR1[48],  VR11:VR10	\n\t" //vr10-15,vr11-16

			"   VSTW   VR44, *++AR4[16]         \n\t"
			"|  VLDDWM2 *-AR1[43],  VR13:VR12  \n\t" //vr12-18,vr13-19

			"   VSTW   VR45, *++AR4[16]         \n\t"
			"|  VLDDWM2 *-AR1[42],  VR15:VR14  \n\t" //vr14-21,vr15-22

			"   VSTW   VR46, *++AR4[16]         \n\t"
			"|  VLDDWM2 *-AR1[41],  VR17:VR16  \n\t" //vr16-24,vr17-25

			"   VSTW   VR47, *++AR4[16]         \n\t"
			"|	VLDDWM2	*-AR1[8],	 VR19:VR18	\n\t" //v18-2   

			"   VSTW   VR48, *++AR4[16]         \n\t"
			"|  VLDDWM2 *-AR1[7],  VR21:VR20  \n\t" //vr16-24,vr17-25

			"   VSTW   VR49, *++AR4[16]         \n\t"
			"|  VLDDWM2 *-AR1[6],  VR23:VR22  \n\t" //vr16-24,vr17-25

			"   VSTW   VR50, *++AR4[16]         \n\t"
			"|  VLDDWM2 *-AR1[1],  VR25:VR24  \n\t" //vr16-24,vr17-25

			"   VSTW   VR51, *++AR4[16]         \n\t"
			"|  VLDDWM2 *+AR1[0],  VR27:VR26  \n\t" //vr16-24,vr17-25

			"   VSTW   VR52, *++AR4[16]         \n\t"
			"|  VLDDWM2 *+AR1[1],   VR29:VR28  \n\t" //vr16-24,vr17-25

			"   VSTW   VR53, *++AR4[16]         \n\t"
			"|  VLDDWM2 *+AR1[6],  VR31:VR30  \n\t" //vr16-24,vr17-25

			"   VSTW   VR0, *++AR3[13]         \n\t"
			"|  VLDDWM2 *+AR1[7],  VR33:VR32  \n\t" //vr16-24,vr17-25
			"|	SLDW	  *+AR15[15], R10	\n\t" //sldw 7

			"   VSTW   VR1, *++AR3[16]      \n\t"
			"|	VLDDWM2 *+AR1[8],  VR35:VR34  \n\t" //vr16-24,vr17-25

			"   VSTW   VR2, *++AR3[16]		\n\t"
			"|  VLDDWM2 *+AR1[41],  VR37:VR36  \n\t" //vr16-24,vr17-25

			"   VSTW   VR3, *++AR3[16]      \n\t"
			"|  VLDDWM2 *+AR1[42],  VR39:VR38  \n\t" //vr16-24,vr17-25

			"   VSTW   VR4, *++AR3[16]      \n\t"
			"|  VLDDWM2 *+AR1[43],  VR41:VR40  \n\t" //vr16-24,vr17-25

			"   VSTW   VR5, *++AR3[16]    \n\t"
			"|  VLDDWM2 *+AR1[48],  VR43:VR42  \n\t" //vr16-24,vr17-25

			"   VSTW   VR6, *++AR3[16]   \n\t"
			"|  VLDDWM2 *+AR1[49],  VR45:VR44  \n\t" //vr16-24,vr17-25

			"   VSTW   VR7, *++AR3[16]   \n\t"
			"|  VLDDWM2 *+AR1[50],  VR47:VR46  \n\t" //vr16-24,vr17-25
			"|	SSHFLL	  3, R10, R10       \n\t"  //SSHFLL 1  

			"   VSTW   VR8, *++AR3[16]   \n\t"
			"|  VLDDWM2 *+AR1[55],  VR49:VR48  \n\t" //vr16-24,vr17-25

			"   VSTW   VR9, *++AR3[16]   \n\t"
			"|  VLDDWM2 *+AR1[56],  VR51:VR50  \n\t" //vr16-24,vr17-25

			"   VSTW   VR10, *++AR3[16]   \n\t"
			"|  VLDDWM2 *+AR1[57],  VR53:VR52  \n\t" //vr16-24,vr17-25

			"   VSTW   VR11, *+AR3[16]   \n\t"
			"|  VSTW   VR12, *+AR3[32]   \n\t"

			"   VSTW   VR13, *+AR3[48]   \n\t"
			"|  VSTW   VR14, *+AR3[64]   \n\t"

			"   VSTW   VR15, *+AR3[80]   \n\t"
			"|  VSTW   VR16, *+AR3[96]   \n\t"

			"   VSTW   VR17, *+AR3[112]   \n\t"
			"|  VSTW   VR18, *+AR3[128]   \n\t"

			"   VSTW   VR19, *+AR3[144]   \n\t"
			"|  VSTW   VR20, *+AR3[160]   \n\t"

			"   VSTW   VR21, *+AR3[176]   \n\t"
			"|  VSTW   VR22, *+AR3[192]   \n\t"

			"   VSTW   VR23, *+AR3[208]   \n\t"
			"|  VSTW   VR24, *+AR3[224]   \n\t"

			"   VSTW   VR25, *+AR3[240]   \n\t"
			"|  VSTW   VR26, *+AR3[256]   \n\t"

			"   VSTW   VR27, *+AR3[272]   \n\t"
			"|  VSTW   VR28, *+AR3[288]   \n\t"
			"|  SADDA     R10,AR0,AR1		\n\t" //sadda AVAF 3 SVAF 2

			"   VSTW   VR29, *+AR3[304]   \n\t"
			"|  VSTW   VR30, *+AR3[320]   \n\t"

			"   VSTW   VR31, *+AR3[336]   \n\t"
			"|  VSTW   VR32, *+AR3[352]   \n\t"

			"   VSTW   VR33, *+AR3[368]   \n\t"
			"|  VSTW   VR34, *+AR3[384]   \n\t"

			"   VSTW   VR35, *+AR3[400]   \n\t"
			"|  VSTW   VR36, *+AR3[416]   \n\t"

			"   VSTW   VR37, *++AR3[432]   \n\t"
			"|  SMVAGA.M2 R31,AR4           \n\t" //R31 // index[15] read
			"	VSTW   VR38, *++AR3[16]   \n\t"
			"|  VLDDWM2 *-AR1[57],  VR1:VR0    \n\t" //vr0-0,vr1-1

			"	VSTW   VR39, *++AR3[16]   \n\t"
			"|  VLDDWM2 *-AR1[56],  VR3:VR2    \n\t" //vr2-3,vr3-4
			
			"	VSTW   VR40, *++AR3[16]   \n\t"
			"|  VLDDWM2 *-AR1[55],  VR5:VR4    \n\t" //vr4-6,vr5-7

			"	VSTW   VR41, *++AR3[16]   \n\t"
			"|  VLDDWM2 *-AR1[50],  VR7:VR6    \n\t" //vr6-9,vr7-10

			"	VSTW   VR42, *++AR3[16]   \n\t"
			"|  VLDDWM2 *-AR1[49],  VR9:VR8    \n\t" //vr4-12,vr5-13

			"   VSTW   VR43, *++AR3[16]         \n\t"
			"|  VLDDWM2 *-AR1[48],  VR11:VR10	\n\t" //vr10-15,vr11-16

			"   VSTW   VR44, *++AR3[16]         \n\t"
			"|  VLDDWM2 *-AR1[43],  VR13:VR12  \n\t" //vr12-18,vr13-19

			"   VSTW   VR45, *++AR3[16]         \n\t"
			"|  VLDDWM2 *-AR1[42],  VR15:VR14  \n\t" //vr14-21,vr15-22

			"   VSTW   VR46, *++AR3[16]         \n\t"
			"|  VLDDWM2 *-AR1[41],  VR17:VR16  \n\t" //vr16-24,vr17-25

			"   VSTW   VR47, *++AR3[16]         \n\t"
			"|	VLDDWM2	*-AR1[8],	 VR19:VR18	\n\t" //v18-2   

			"   VSTW   VR48, *++AR3[16]         \n\t"
			"|  VLDDWM2 *-AR1[7],  VR21:VR20  \n\t" //vr16-24,vr17-25

			"   VSTW   VR49, *++AR3[16]         \n\t"
			"|  VLDDWM2 *-AR1[6],  VR23:VR22  \n\t" //vr16-24,vr17-25

			"   VSTW   VR50, *++AR3[16]         \n\t"
			"|  VLDDWM2 *-AR1[1],  VR25:VR24  \n\t" //vr16-24,vr17-25

			"   VSTW   VR51, *++AR3[16]         \n\t"
			"|  VLDDWM2 *+AR1[0],  VR27:VR26  \n\t" //vr16-24,vr17-25

			"   VSTW   VR52, *++AR3[16]         \n\t"
			"|  VLDDWM2 *+AR1[1],   VR29:VR28  \n\t" //vr16-24,vr17-25

			"   VSTW   VR53, *++AR3[16]         \n\t"
			"|  VLDDWM2 *+AR1[6],  VR31:VR30  \n\t" //vr16-24,vr17-25

			"   VSTW   VR0, *++AR4[15]         \n\t"
			"|  VLDDWM2 *+AR1[7],  VR33:VR32  \n\t" //vr16-24,vr17-25
			"|	SLDW	  *+AR15[11], R10	\n\t" //sldw 7

			"   VSTW   VR1, *++AR4[16]      \n\t"
			"|	VLDDWM2 *+AR1[8],  VR35:VR34  \n\t" //vr16-24,vr17-25

			"   VSTW   VR2, *++AR4[16]		\n\t"
			"|  VLDDWM2 *+AR1[41],  VR37:VR36  \n\t" //vr16-24,vr17-25

			"   VSTW   VR3, *++AR4[16]      \n\t"
			"|  VLDDWM2 *+AR1[42],  VR39:VR38  \n\t" //vr16-24,vr17-25

			"   VSTW   VR4, *++AR4[16]      \n\t"
			"|  VLDDWM2 *+AR1[43],  VR41:VR40  \n\t" //vr16-24,vr17-25

			"   VSTW   VR5, *++AR4[16]    \n\t"
			"|  VLDDWM2 *+AR1[48],  VR43:VR42  \n\t" //vr16-24,vr17-25

			"   VSTW   VR6, *++AR4[16]   \n\t"
			"|  VLDDWM2 *+AR1[49],  VR45:VR44  \n\t" //vr16-24,vr17-25

			"   VSTW   VR7, *++AR4[16]   \n\t"
			"|  VLDDWM2 *+AR1[50],  VR47:VR46  \n\t" //vr16-24,vr17-25
			"|	SSHFLL	  3, R10, R10       \n\t"  //SSHFLL 1  

			"   VSTW   VR8, *++AR4[16]   \n\t"
			"|  VLDDWM2 *+AR1[55],  VR49:VR48  \n\t" //vr16-24,vr17-25

			"   VSTW   VR9, *++AR4[16]   \n\t"
			"|  VLDDWM2 *+AR1[56],  VR51:VR50  \n\t" //vr16-24,vr17-25

			"   VSTW   VR10, *++AR4[16]   \n\t"
			"|  VLDDWM2 *+AR1[57],  VR53:VR52  \n\t" //vr16-24,vr17-25

			"   VSTW   VR11, *+AR4[16]   \n\t"
			"|  VSTW   VR12, *+AR4[32]   \n\t"

			"   VSTW   VR13, *+AR4[48]   \n\t"
			"|  VSTW   VR14, *+AR4[64]   \n\t"

			"   VSTW   VR15, *+AR4[80]   \n\t"
			"|  VSTW   VR16, *+AR4[96]   \n\t"

			"   VSTW   VR17, *+AR4[112]   \n\t"
			"|  VSTW   VR18, *+AR4[128]   \n\t"

			"   VSTW   VR19, *+AR4[144]   \n\t"
			"|  VSTW   VR20, *+AR4[160]   \n\t"

			"   VSTW   VR21, *+AR4[176]   \n\t"
			"|  VSTW   VR22, *+AR4[192]   \n\t"

			"   VSTW   VR23, *+AR4[208]   \n\t"
			"|  VSTW   VR24, *+AR4[224]   \n\t"

			"   VSTW   VR25, *+AR4[240]   \n\t"
			"|  VSTW   VR26, *+AR4[256]   \n\t"

			"   VSTW   VR27, *+AR4[272]   \n\t"
			"|  VSTW   VR28, *+AR4[288]   \n\t"

			"   VSTW   VR29, *+AR4[304]   \n\t"
			"|  VSTW   VR30, *+AR4[320]   \n\t"

			"   VSTW   VR31, *+AR4[336]   \n\t"
			"|  VSTW   VR32, *+AR4[352]   \n\t"

			"   VSTW   VR33, *+AR4[368]   \n\t"
			"|  VSTW   VR34, *+AR4[384]   \n\t"

			"   VSTW   VR35, *+AR4[400]   \n\t"
			"|  VSTW   VR36, *+AR4[416]   \n\t"

			"   VSTW   VR37, *+AR4[432]   \n\t"

			"	VSTW   VR38, *+AR4[448]   \n\t"
			"|	VSTW   VR39, *+AR4[464]   \n\t"

			"	VSTW   VR40, *+AR4[480]   \n\t"
			"|	VSTW   VR41, *+AR4[496]   \n\t"

			"	VSTW   VR42, *+AR4[512]   \n\t"
			"|  VSTW   VR43, *+AR4[528]         \n\t"

			"   VSTW   VR44, *+AR4[544]         \n\t"
			"|  VSTW   VR45, *+AR4[560]         \n\t"

			"   VSTW   VR46, *+AR4[576]         \n\t"
			"|  VSTW   VR47, *+AR4[592]         \n\t"

			"   VSTW   VR48, *+AR4[608]         \n\t"
			"|  VSTW   VR49, *+AR4[624]         \n\t"

			"   VSTW   VR50, *+AR4[640]         \n\t"
			"|  VSTW   VR51, *+AR4[656]         \n\t"

			"   VSTW   VR52, *+AR4[672]         \n\t"
			"|  VSTW   VR53, *+AR4[688]         \n\t"
			"|  SMOVI24   0xFFFF,R29           \n\t"

			"	SNOP	3						\n\t"
			"   SMVCGC  R29,VLR           \n\t" 
			"	SNOP	4						\n\t"
			:
			: [eb1]   "r" (EB1),
			  [eb2]   "r" (EB2),
			  [eb3]   "r" (EB3),
			  [Vec]  "r" (vec),
			  [node] "r" (node_id)
			: "R10","R11","R29","R30","R31","R32",\
			  "VR0","VR1","VR2","VR3","VR4","VR5","VR6","VR7","VR8",\
			  "VR9","VR10","VR11","VR12","VR13","VR14","VR15","VR16",\
			  "VR17","VR18","VR19","VR20","VR21","VR22","VR23","VR24",\
			  "VR25","VR26","VR27","VR28","VR29","VR30","VR31","VR32",\
			  "VR33","VR34","VR35","VR36","VR37","VR38","VR39","VR40",\
			  "VR41","VR42","VR43","VR44","VR45","VR46","VR47","VR48",\
			  "VR49","VR50","VR51","VR52","VR53"
			);
			
}


// Higher-order interpolation vector computation based on assembly instructions of the MT-3000 architecture.
__shared__ void S0_multiply_eb(lvector double *vec, lvector double *BE_node,lvector double *S0)
{
	
			__asm__ __volatile__(
			"   SMOVI24   0x1,R29           \n\t"
			"   SMULIU.M1 R29,%[Vec],R31    \n\t" //smuliu 3
			"|  SMULIU.M2 R29,%[s0], R30    \n\t" //smuliu 3
			"   SMULIU.M1 R29,%[BE],R32    \n\t" //smuliu 3
			"   SNOP      1                 \n\t"
			"   SMVAGA.M1 R31,AR0           \n\t" //smvaga 2
			"|  SMVAGA.M2 R30,AR1           \n\t" //smvaga 2
			"   SMVAGA.M1 R32,AR2           \n\t" //smvaga 2

			"   VLDDW     *+AR1[0], VR1:VR0 \n\t"
			"|  VLDDW     *+AR1[16],VR3:VR2 \n\t"
			"|  VMOVI       1,VR29          \n\t"

			"   VLDDW     *+AR1[32],VR5:VR4 \n\t"
			"|  VLDDW     *+AR1[48],VR7:VR6 \n\t"
			"|  VFINTD.M1 VR29, VR63        \n\t" //vfintd 3

			"   VLDW     *+AR1[128],VR8     \n\t"
			"|  VLDW     *+AR0[640],VR20    \n\t" //vr16-24,vr17-25

			"   VLDW	*+AR0[672], VR21   \n\t" 
			"|  VLDW	*+AR0[704], VR22   \n\t" 

			"   VLDW	*+AR0[736], VR23   \n\t" 
			"|  VLDW	*+AR0[768], VR24   \n\t" 

			"   VLDW	*+AR0[800], VR25  \n\t" 
			"|  VLDW	*+AR0[832], VR26  \n\t"

			"   VLDW *+AR0[448], VR14   \n\t"
			"|  VLDW *+AR0[480], VR15   \n\t" 

			"   VLDW *+AR0[512], VR16   \n\t"
			"|  VLDW *+AR0[544], VR17   \n\t"

			"   VLDW *+AR0[576], VR18   \n\t"
			"|  VLDW *+AR0[608], VR19   \n\t"
//			"   SNOP    3                   \n\t"

			"   VFMULD.M1 VR0,VR3,VR9      \n\t"   //vfmuld 4
			"|  VFMULD.M2 VR1,VR3,VR10      \n\t"
			"|  VFMULD.M3 VR2,VR3,VR11      \n\t"

			"   VFMULD.M1 VR0,VR4,VR57      \n\t"   //vfmuld 4 VR 57-VR 14
			"|  VFMULD.M2 VR1,VR4,VR58      \n\t"  //vr58-vr15
			"|  VFMULD.M3 VR2,VR4,VR59      \n\t" //vr59-vr16

			"   VFMULD.M1 VR0,VR5,VR60      \n\t"   //vfmuld 4
			"|  VFMULD.M2 VR1,VR5,VR61      \n\t" // vr17-19 - vr60-62
			"|  VFMULD.M3 VR2,VR5,VR62      \n\t"

			"   SNOP    1                   \n\t"

			"   VFMULD.M1 VR9,VR6,VR30     \n\t"   //vfmuld 4
			"|  VFMULD.M2 VR9,VR7,VR31     \n\t"
			"|  VFMULD.M3 VR9,VR8,VR32     \n\t"
			"|  VLDW *+AR0[96],  VR3		\n\t" 
			"|  VLDW *+AR0[384], VR12   \n\t" //vr16-24,vr17-25

			"   VFMULD.M1 VR10,VR6,VR39     \n\t"   //vfmuld 4
			"|  VFMULD.M2 VR10,VR7,VR40     \n\t"
			"|  VFMULD.M3 VR10,VR8,VR41     \n\t"
			"|  VLDW *+AR0[128], VR4		\n\t" 
			"|  VLDW *+AR0[416], VR13   \n\t" //vr16-24,vr17-25

			"   VFMULD.M1 VR11,VR6,VR48     \n\t"   //vfmuld 4
			"|  VFMULD.M2 VR11,VR7,VR49     \n\t"
			"|  VFMULD.M3 VR11,VR8,VR50     \n\t"
			"|  VLDW *+AR0[160], VR5    \n\t" 
			"|  VLDW *+AR0[0],   VR0    \n\t" 

			"   VFMULD.M1 VR57,VR6,VR33     \n\t"   //vfmuld 4
			"|  VFMULD.M2 VR57,VR7,VR34     \n\t"
			"|  VFMULD.M3 VR57,VR8,VR35     \n\t"
			"|  VLDW *+AR0[32],  VR1    \n\t" 
			"|  VLDW *+AR0[64],  VR2    \n\t" 

			"   VFMULD.M1 VR58,VR6,VR42     \n\t"   //vfmuld 4
			"|  VFMULD.M2 VR58,VR7,VR43     \n\t"
			"|  VFMULD.M3 VR58,VR8,VR44     \n\t"
			"|  VLDW *+AR0[288], VR9    \n\t" 

			"   VFMULD.M1 VR59,VR6,VR51     \n\t"   //vfmuld 4
			"|  VFMULD.M2 VR59,VR7,VR52     \n\t"
			"|  VFMULD.M3 VR59,VR8,VR53     \n\t"
			"|  VLDW *+AR0[320], VR10   \n\t" 

			"   VFMULD.M1 VR60,VR6,VR36     \n\t"   //vfmuld 4
			"|  VFMULD.M2 VR60,VR7,VR37     \n\t"
			"|  VFMULD.M3 VR60,VR8,VR38     \n\t"
			"|  VLDW *+AR0[352], VR11   \n\t" //vr16-24,vr17-25

			"   VFMULD.M1 VR61,VR6,VR45     \n\t"   //vfmuld 4
			"|  VFMULD.M2 VR61,VR7,VR46     \n\t"
			"|  VFMULD.M3 VR61,VR8,VR47     \n\t"

			"   VFMULD.M1 VR62,VR6,VR54     \n\t"   //vfmuld 4
			"|  VFMULD.M2 VR62,VR7,VR55     \n\t"
			"|  VFMULD.M3 VR62,VR8,VR56     \n\t"
//E1
			"   VFMULD.M1   VR21,VR51,VR21  \n\t"
			"|  VFMULD.M2   VR22,VR52,VR22  \n\t"
			"|  VFMULD.M3   VR23,VR53,VR23  \n\t"

			"   VFMULD.M1   VR12,VR42,VR12  \n\t"
			"|  VFMULD.M2   VR13,VR43,VR13  \n\t"
			"|  VFMULD.M3   VR14,VR44,VR14  \n\t"

			"   VFMULD.M1   VR15,VR45,VR15  \n\t"
			"|  VFMULD.M2   VR16,VR46,VR16  \n\t"
			"|  VFMULD.M3   VR17,VR47,VR17  \n\t"

			"   VFMULD.M1   VR18,VR48,VR18  \n\t" 
			"|  VFMULD.M2   VR19,VR49,VR19  \n\t"
			"|  VFMULD.M3   VR20,VR50,VR20  \n\t"
			"|  VLDW *+AR0[192], VR6    \n\t" 
			"|  VLDW *+AR0[224], VR7    \n\t" 

			"   VFMULAD.M1   VR0,VR30,VR21,VR21	  \n\t" //vfmuld 4
			"|  VFMULAD.M2   VR1,VR31,VR22,VR22   \n\t"
			"|  VFMULAD.M3   VR2,VR32,VR23,VR23   \n\t"
			"|  VLDW *+AR0[256], VR8    \n\t"

			"   VFMULAD.M1   VR3,VR33,VR12,VR12    \n\t"
			"|  VFMULAD.M2   VR4,VR34,VR13,VR13    \n\t"
			"|  VFMULAD.M3   VR5,VR35,VR14,VR14    \n\t"

			"   VFMULAD.M1   VR9, VR39,VR15,VR15  \n\t" 
			"|  VFMULAD.M2   VR10,VR40,VR16,VR16  \n\t"
			"|  VFMULAD.M3   VR11,VR41,VR17,VR17  \n\t"

			"   VFMULAD.M1   VR24,VR54,VR18,VR24  \n\t"
			"|  VFMULAD.M2   VR25,VR55,VR19,VR25  \n\t"
			"|  VFMULAD.M3   VR26,VR56,VR20,VR26  \n\t"

			"	SNOP 2					\n\t"

			"   VLDW *+AR0[16],  VR0    \n\t"
			"|  VLDW *+AR0[48],  VR1    \n\t" 
//			"	SNOP 3					\n\t"

			"   VFMULAD.M1   VR12,VR63,VR21,VR21	  \n\t" //vfmuld 4
			"|  VFMULAD.M2   VR13,VR63,VR22,VR22   \n\t"
			"|  VFMULAD.M3   VR14,VR63,VR23,VR23   \n\t"
			"|  VLDW *+AR0[80],  VR2    \n\t" 
			"|  VLDW *+AR0[112], VR3	\n\t"

			"   VLDW *+AR0[144], VR4    \n\t"
			"|  VLDW *+AR0[176], VR5    \n\t"
//			"	SNOP 1					\n\t"

			"   VFMULAD.M1   VR6,VR36,VR15,VR27	  \n\t" //vfmuld 4
			"|  VFMULAD.M2   VR7,VR37,VR16,VR28   \n\t"
			"|  VFMULAD.M3   VR8,VR38,VR17,VR29   \n\t"
			"|  VLDW *+AR0[304], VR9    \n\t"
			"|  VLDW *+AR0[336], VR10   \n\t"
			
			"   VLDW *+AR0[368], VR11   \n\t" 			
			"	SNOP 2					\n\t"
//			"	SNOP 3					\n\t"

			"   VFMULAD.M1   VR24,VR63,VR21,VR24  \n\t"
			"|  VFMULAD.M2   VR25,VR63,VR22,VR25  \n\t"
			"|  VFMULAD.M3   VR26,VR63,VR23,VR26  \n\t"
			"|  VLDW *+AR0[400], VR12   \n\t"
			"|  VLDW *+AR0[432], VR13   \n\t"

			"   VLDW *+AR0[464], VR14   \n\t"

			"   VLDW *+AR0[208], VR6    \n\t"
			"|  VLDW *+AR0[240], VR7    \n\t"
// B1
			"   VLDW *+AR0[272], VR8    \n\t"
			"|  VLDW *+AR0[496], VR15   \n\t"
			"|  VFMULD.M1   VR0,VR30,VR0    \n\t"
			"|  VFMULD.M2   VR1,VR31,VR1    \n\t"
			"|  VFMULD.M3   VR2,VR32,VR2    \n\t"

			"   VLDW *+AR0[528], VR16   \n\t"
			"|  VLDW *+AR0[560], VR17   \n\t"
			"|  VFMULD.M1   VR3,VR33,VR3    \n\t"
			"|  VFMULD.M2   VR4,VR34,VR4    \n\t"
			"|  VFMULD.M3   VR5,VR35,VR5    \n\t"

			"   VLDW *+AR0[592], VR18   \n\t"
			"|  VLDW *+AR0[624], VR19   \n\t"


//			"	SNOP 5					\n\t"

			"   VFMULAD.M1   VR24,VR63,VR27,VR27  \n\t"
			"|  VFMULAD.M2   VR25,VR63,VR28,VR28  \n\t"
			"|  VFMULAD.M3   VR26,VR63,VR29,VR29  \n\t"
			"|  VLDW *+AR0[656], VR20   \n\t"
			"|  VLDW *+AR0[688], VR21   \n\t"

			"   VLDW *+AR0[720], VR22   \n\t"
			"|  VLDW *+AR0[752], VR23   \n\t"

			"   SNOP    4                   \n\t"

			"   VFMULAD.M1   VR27,VR63,VR28,VR28    \n\t"
			"|  VLDW *+AR0[784], VR24   \n\t"
			"|  VLDW *+AR0[816], VR25   \n\t"

			"   VLDW *+AR0[848], VR26   \n\t"

			"   SNOP    1                   \n\t"

			"   VFMULD.M1   VR6,VR36,VR6    \n\t"
			"|  VFMULD.M2   VR7,VR37,VR7    \n\t"
			"|  VFMULD.M3   VR8,VR38,VR8    \n\t"
			
			"   VFMULD.M1   VR9, VR39,VR9    \n\t"
			"|  VFMULD.M2   VR10,VR40,VR10    \n\t"
			"|  VFMULD.M3   VR11,VR41,VR11    \n\t"

			"   VFMULAD.M1   VR12,VR42,VR0,VR12    \n\t"
			"|  VFMULAD.M2   VR13,VR43,VR1,VR13    \n\t"
			"|  VFMULAD.M3   VR14,VR44,VR2,VR14    \n\t"

//			"   SNOP    5                   \n\t"
			"   VFMULAD.M1   VR28,VR63,VR29,VR57    \n\t"
			
			"   VFMULAD.M1   VR15,VR45,VR3,VR15    \n\t"
			"|  VFMULAD.M2   VR16,VR46,VR4,VR16    \n\t"
			"|  VFMULAD.M3   VR17,VR47,VR5,VR17    \n\t"

			"   VFMULAD.M1   VR18,VR48,VR6,VR18    \n\t"
			"|  VFMULAD.M2   VR19,VR49,VR7,VR19    \n\t"
			"|  VFMULAD.M3   VR20,VR50,VR8,VR20    \n\t"

			"   VFMULAD.M1   VR21,VR51,VR9, VR21    \n\t"
			"|  VFMULAD.M2   VR22,VR52,VR10,VR22    \n\t"
			"|  VFMULAD.M3   VR23,VR53,VR11,VR23    \n\t"

			"   SMOVI24     0x1B00, R10             \n\t"

			"   VFMULAD.M1   VR24,VR54,VR12,VR24    \n\t"
			"|  VFMULAD.M2   VR25,VR55,VR13,VR25    \n\t"
			"|  VFMULAD.M3   VR26,VR56,VR14,VR26    \n\t"
			"|  SADDA     R10,AR0,AR0       \n\t" //sadda AVAF 3 SVAF 2

//			"   SNOP    5                   \n\t"
			"   VSTW        VR57,*+AR2[0]       \n\t"

			"   SNOP      1                 \n\t"

			"   VFMULAD.M1   VR15,VR63,VR18,VR27    \n\t"
			"|  VFMULAD.M2   VR16,VR63,VR19,VR28    \n\t"
			"|  VFMULAD.M3   VR17,VR63,VR20,VR29    \n\t"

			"   VLDW *+AR0[0],   VR0    \n\t"
			"|  VLDW *+AR0[32],  VR1    \n\t"

			"   VLDW *+AR0[64],  VR2    \n\t"
			"|  VLDW *+AR0[96],  VR3    \n\t"

			"   VFMULAD.M1   VR24,VR63,VR21,VR24    \n\t"
			"|  VFMULAD.M2   VR25,VR63,VR22,VR25    \n\t"
			"|  VFMULAD.M3   VR26,VR63,VR23,VR26    \n\t"
			"|  VLDW *+AR0[128], VR4    \n\t"
			"|  VLDW *+AR0[160], VR5    \n\t"

			"   VLDW *+AR0[192],  VR6   \n\t"
			"|  VLDW *+AR0[224],  VR7   \n\t"

			"   VLDW *+AR0[256],  VR8   \n\t"
			"|  VLDW *+AR0[288],  VR9   \n\t"

			"   VLDW *+AR0[320],  VR10  \n\t"
			"|  VLDW *+AR0[352],  VR11  \n\t"

			"   VLDW *+AR0[384],  VR12  \n\t"
			"|  VLDW *+AR0[416],  VR13  \n\t"

			"   VLDW *+AR0[448],  VR14  \n\t"
			"|  VLDW *+AR0[480],  VR15  \n\t"

			"   VFMULAD.M1   VR24,VR63,VR27,VR27    \n\t"
			"|  VFMULAD.M2   VR25,VR63,VR28,VR28    \n\t"
			"|  VFMULAD.M3   VR26,VR63,VR29,VR29    \n\t"
			"|  VLDW *+AR0[512],  VR16  \n\t"
			"|  VLDW *+AR0[544],  VR17  \n\t"
				
			"   VLDW *+AR0[576], VR18   \n\t"
			"|  VLDW *+AR0[608], VR19   \n\t"

			"   VLDW *+AR0[640], VR20   \n\t"
			"|  VLDW *+AR0[672], VR21   \n\t" 
//E2
			"   VLDW *+AR0[704],  VR22  \n\t"
			"   VLDW *+AR0[736],  VR23  \n\t"
			"|  VFMULD.M1   VR0,VR30,VR0    \n\t"
			"|  VFMULD.M2   VR1,VR31,VR1    \n\t"
			"|  VFMULD.M3   VR2,VR32,VR2    \n\t"

			"   VFMULD.M1   VR3,VR33,VR3    \n\t"
			"|  VFMULD.M2   VR4,VR34,VR4    \n\t"
			"|  VFMULD.M3   VR5,VR35,VR5    \n\t"

			"   VFMULD.M1   VR6,VR36,VR6    \n\t"
			"|  VFMULD.M2   VR7,VR37,VR7    \n\t"
			"|  VFMULD.M3   VR8,VR38,VR8    \n\t"

			"   VFMULAD.M1   VR27,VR63,VR28,VR28    \n\t"
			"|  VLDW *+AR0[768], VR24   \n\t"
			"|  VLDW *+AR0[800], VR25   \n\t"

			"   VLDW *+AR0[832], VR26   \n\t"
			"	SNOP 1\n\t"

			"   VFMULD.M1   VR9, VR39,VR9    \n\t"
			"|  VFMULD.M2   VR10,VR40,VR10    \n\t"
			"|  VFMULD.M3   VR11,VR41,VR11    \n\t"

			"   VFMULAD.M1   VR12,VR42,VR0,VR12    \n\t"
			"|  VFMULAD.M2   VR13,VR43,VR1,VR13    \n\t"
			"|  VFMULAD.M3   VR14,VR44,VR2,VR14    \n\t"

			"   VFMULAD.M1   VR15,VR45,VR3,VR15    \n\t"
			"|  VFMULAD.M2   VR16,VR46,VR4,VR16    \n\t"
			"|  VFMULAD.M3   VR17,VR47,VR5,VR17    \n\t"

			"   VFMULAD.M1   VR28,VR63,VR29,VR57    \n\t"

			"   VFMULAD.M1   VR18,VR48,VR6,VR18    \n\t"
			"|  VFMULAD.M2   VR19,VR49,VR7,VR19    \n\t"
			"|  VFMULAD.M3   VR20,VR50,VR8,VR20    \n\t"

			"   VFMULAD.M1   VR21,VR51,VR9, VR21    \n\t"
			"|  VFMULAD.M2   VR22,VR52,VR10,VR22    \n\t"
			"|  VFMULAD.M3   VR23,VR53,VR11,VR23    \n\t"

			"	SNOP 1		\n\t"  //vload vr24 9 

			"   VFMULAD.M1   VR24,VR54,VR12,VR24    \n\t"
			"|  VFMULAD.M2   VR25,VR55,VR13,VR25    \n\t"
			"|  VFMULAD.M3   VR26,VR56,VR14,VR26    \n\t"

			"	SNOP 1 \n\t"

			"   VSTW        VR57,*+AR2[16]       \n\t"

			"   VFMULAD.M1   VR15,VR63,VR18,VR27    \n\t"
			"|  VFMULAD.M2   VR16,VR63,VR19,VR28    \n\t"
			"|  VFMULAD.M3   VR17,VR63,VR20,VR29    \n\t"

			"   VLDW *+AR0[16],   VR0    \n\t"
			"|  VLDW *+AR0[48],  VR1    \n\t"

			"   VLDW *+AR0[80],  VR2    \n\t"
			"|  VLDW *+AR0[112],  VR3    \n\t"


			"   VFMULAD.M1   VR24,VR63,VR21,VR24    \n\t"
			"|  VFMULAD.M2   VR25,VR63,VR22,VR25    \n\t"
			"|  VFMULAD.M3   VR26,VR63,VR23,VR26    \n\t"
			"|  VLDW *+AR0[144], VR4    \n\t"
			"|  VLDW *+AR0[176], VR5    \n\t"

			"   VLDW *+AR0[208],  VR6   \n\t"
			"|  VLDW *+AR0[240],  VR7   \n\t"

			"   VLDW *+AR0[272],  VR8   \n\t"
			"|  VLDW *+AR0[304],  VR9   \n\t"

			"   VLDW *+AR0[336],  VR10  \n\t"
			"|  VLDW *+AR0[368],  VR11  \n\t"

			"   VLDW *+AR0[400],  VR12  \n\t"
			"|  VLDW *+AR0[432],  VR13  \n\t"

			"   VLDW *+AR0[464],  VR14  \n\t"
			"|  VLDW *+AR0[496],  VR15  \n\t"


			"   VFMULAD.M1   VR24,VR63,VR27,VR27    \n\t"
			"|  VFMULAD.M2   VR25,VR63,VR28,VR28    \n\t"
			"|  VFMULAD.M3   VR26,VR63,VR29,VR29    \n\t"
			"|  VLDW *+AR0[528],  VR16  \n\t"
			"|  VLDW *+AR0[560],  VR17  \n\t"
				
			"   VLDW *+AR0[592], VR18   \n\t"
			"|  VLDW *+AR0[624], VR19   \n\t"

			"   VLDW *+AR0[656], VR20   \n\t"
			"|  VLDW *+AR0[688], VR21   \n\t" 

//B2
			"   VLDW *+AR0[720],  VR22  \n\t"
			"   VLDW *+AR0[752],  VR23  \n\t"
			"|  VFMULD.M1   VR0,VR30,VR0    \n\t"
			"|  VFMULD.M2   VR1,VR31,VR1    \n\t"
			"|  VFMULD.M3   VR2,VR32,VR2    \n\t"

			"   VFMULD.M1   VR3,VR33,VR3    \n\t"
			"|  VFMULD.M2   VR4,VR34,VR4    \n\t"
			"|  VFMULD.M3   VR5,VR35,VR5    \n\t"

			"   VFMULD.M1   VR6,VR36,VR6    \n\t"
			"|  VFMULD.M2   VR7,VR37,VR7    \n\t"
			"|  VFMULD.M3   VR8,VR38,VR8    \n\t"

			"   VFMULAD.M1   VR27,VR63,VR28,VR28    \n\t"
			"|  VLDW *+AR0[784], VR24   \n\t"
			"|  VLDW *+AR0[816], VR25   \n\t"

			"   VLDW *+AR0[848], VR26   \n\t"
			"	SNOP 1\n\t"

			"   VFMULD.M1   VR9, VR39,VR9    \n\t"
			"|  VFMULD.M2   VR10,VR40,VR10    \n\t"
			"|  VFMULD.M3   VR11,VR41,VR11    \n\t"

			"   VFMULAD.M1   VR12,VR42,VR0,VR12    \n\t"
			"|  VFMULAD.M2   VR13,VR43,VR1,VR13    \n\t"
			"|  VFMULAD.M3   VR14,VR44,VR2,VR14    \n\t"

			"   VFMULAD.M1   VR15,VR45,VR3,VR15    \n\t"
			"|  VFMULAD.M2   VR16,VR46,VR4,VR16    \n\t"
			"|  VFMULAD.M3   VR17,VR47,VR5,VR17    \n\t"

			"   VFMULAD.M1   VR28,VR63,VR29,VR57    \n\t"

			"   VFMULAD.M1   VR18,VR48,VR6,VR18    \n\t"
			"|  VFMULAD.M2   VR19,VR49,VR7,VR19    \n\t"
			"|  VFMULAD.M3   VR20,VR50,VR8,VR20    \n\t"

			"   VFMULAD.M1   VR21,VR51,VR9, VR21    \n\t"
			"|  VFMULAD.M2   VR22,VR52,VR10,VR22    \n\t"
			"|  VFMULAD.M3   VR23,VR53,VR11,VR23    \n\t"

			"	SNOP 1		\n\t"  //vload vr24 9 

			"   VFMULAD.M1   VR24,VR54,VR12,VR24    \n\t"
			"|  VFMULAD.M2   VR25,VR55,VR13,VR25    \n\t"
			"|  VFMULAD.M3   VR26,VR56,VR14,VR26    \n\t"
			"|  SMOVI24     0x1B00, R10             \n\t"

			"   SADDA     R10,AR0,AR0       \n\t" //sadda AVAF 3 SVAF 2

			"   VSTW        VR57,*+AR2[32]       \n\t"

			"   VFMULAD.M1   VR15,VR63,VR18,VR27    \n\t"
			"|  VFMULAD.M2   VR16,VR63,VR19,VR28    \n\t"
			"|  VFMULAD.M3   VR17,VR63,VR20,VR29    \n\t"

			"   VLDW *+AR0[0],   VR0    \n\t"
			"|  VLDW *+AR0[32],  VR1    \n\t"

			"   VLDW *+AR0[64],  VR2    \n\t"
			"|  VLDW *+AR0[96],  VR3    \n\t"

			"   VFMULAD.M1   VR24,VR63,VR21,VR24    \n\t"
			"|  VFMULAD.M2   VR25,VR63,VR22,VR25    \n\t"
			"|  VFMULAD.M3   VR26,VR63,VR23,VR26    \n\t"
			"|  VLDW *+AR0[128], VR4    \n\t"
			"|  VLDW *+AR0[160], VR5    \n\t"

			"   VLDW *+AR0[192],  VR6   \n\t"
			"|  VLDW *+AR0[224],  VR7   \n\t"

			"   VLDW *+AR0[256],  VR8   \n\t"
			"|  VLDW *+AR0[288],  VR9   \n\t"

			"   VLDW *+AR0[320],  VR10  \n\t"
			"|  VLDW *+AR0[352],  VR11  \n\t"

			"   VLDW *+AR0[384],  VR12  \n\t"
			"|  VLDW *+AR0[416],  VR13  \n\t"

			"   VLDW *+AR0[448],  VR14  \n\t"
			"|  VLDW *+AR0[480],  VR15  \n\t"

			"   VFMULAD.M1   VR24,VR63,VR27,VR27    \n\t"
			"|  VFMULAD.M2   VR25,VR63,VR28,VR28    \n\t"
			"|  VFMULAD.M3   VR26,VR63,VR29,VR29    \n\t"
			"|  VLDW *+AR0[512],  VR16  \n\t"
			"|  VLDW *+AR0[544],  VR17  \n\t"
				
			"   VLDW *+AR0[576], VR18   \n\t"
			"|  VLDW *+AR0[608], VR19   \n\t"

			"   VLDW *+AR0[640], VR20   \n\t"
			"|  VLDW *+AR0[672], VR21   \n\t" 
//E3
			"   VLDW *+AR0[704],  VR22  \n\t"
			"   VLDW *+AR0[736],  VR23  \n\t"
			"|  VFMULD.M1   VR0,VR30,VR0    \n\t"
			"|  VFMULD.M2   VR1,VR31,VR1    \n\t"
			"|  VFMULD.M3   VR2,VR32,VR2    \n\t"

			"   VFMULD.M1   VR3,VR33,VR3    \n\t"
			"|  VFMULD.M2   VR4,VR34,VR4    \n\t"
			"|  VFMULD.M3   VR5,VR35,VR5    \n\t"

			"   VFMULD.M1   VR6,VR36,VR6    \n\t"
			"|  VFMULD.M2   VR7,VR37,VR7    \n\t"
			"|  VFMULD.M3   VR8,VR38,VR8    \n\t"

			"   VFMULAD.M1   VR27,VR63,VR28,VR28    \n\t"
			"|  VLDW *+AR0[768], VR24   \n\t"
			"|  VLDW *+AR0[800], VR25   \n\t"

			"   VLDW *+AR0[832], VR26   \n\t"
			"	SNOP 1\n\t"

			"   VFMULD.M1   VR9, VR39,VR9    \n\t"
			"|  VFMULD.M2   VR10,VR40,VR10    \n\t"
			"|  VFMULD.M3   VR11,VR41,VR11    \n\t"

			"   VFMULAD.M1   VR12,VR42,VR0,VR12    \n\t"
			"|  VFMULAD.M2   VR13,VR43,VR1,VR13    \n\t"
			"|  VFMULAD.M3   VR14,VR44,VR2,VR14    \n\t"

			"   VFMULAD.M1   VR15,VR45,VR3,VR15    \n\t"
			"|  VFMULAD.M2   VR16,VR46,VR4,VR16    \n\t"
			"|  VFMULAD.M3   VR17,VR47,VR5,VR17    \n\t"

			"   VFMULAD.M1   VR28,VR63,VR29,VR57    \n\t"

			"   VFMULAD.M1   VR18,VR48,VR6,VR18    \n\t"
			"|  VFMULAD.M2   VR19,VR49,VR7,VR19    \n\t"
			"|  VFMULAD.M3   VR20,VR50,VR8,VR20    \n\t"

			"   VFMULAD.M1   VR21,VR51,VR9, VR21    \n\t"
			"|  VFMULAD.M2   VR22,VR52,VR10,VR22    \n\t"
			"|  VFMULAD.M3   VR23,VR53,VR11,VR23    \n\t"

			"	SNOP 1		\n\t"  //vload vr24 9 

			"   VFMULAD.M1   VR24,VR54,VR12,VR24    \n\t"
			"|  VFMULAD.M2   VR25,VR55,VR13,VR25    \n\t"
			"|  VFMULAD.M3   VR26,VR56,VR14,VR26    \n\t"

			"	SNOP 1 \n\t"

			"   VSTW        VR57,*+AR2[48]       \n\t"

			"   VFMULAD.M1   VR15,VR63,VR18,VR27    \n\t"
			"|  VFMULAD.M2   VR16,VR63,VR19,VR28    \n\t"
			"|  VFMULAD.M3   VR17,VR63,VR20,VR29    \n\t"
			

			"   VLDW *+AR0[16],  VR0    \n\t"
			"|  VLDW *+AR0[48],  VR1    \n\t"

			"   VLDW *+AR0[80],  VR2    \n\t"
			"|  VLDW *+AR0[112], VR3    \n\t"


			"   VFMULAD.M1   VR24,VR63,VR21,VR24    \n\t"
			"|  VFMULAD.M2   VR25,VR63,VR22,VR25    \n\t"
			"|  VFMULAD.M3   VR26,VR63,VR23,VR26    \n\t"
			"|  VLDW *+AR0[144], VR4    \n\t"
			"|  VLDW *+AR0[176], VR5    \n\t"

			"   VLDW *+AR0[208],  VR6   \n\t"
			"|  VLDW *+AR0[240],  VR7   \n\t"

			"   VLDW *+AR0[272],  VR8   \n\t"
			"|  VLDW *+AR0[304],  VR9   \n\t"

			"   VLDW *+AR0[336],  VR10  \n\t"
			"|  VLDW *+AR0[368],  VR11  \n\t"

			"   VLDW *+AR0[400],  VR12  \n\t"
			"|  VLDW *+AR0[432],  VR13  \n\t"

			"   VLDW *+AR0[464],  VR14  \n\t"
			"|  VLDW *+AR0[496],  VR15  \n\t"


			"   VFMULAD.M1   VR24,VR63,VR27,VR27    \n\t"
			"|  VFMULAD.M2   VR25,VR63,VR28,VR28    \n\t"
			"|  VFMULAD.M3   VR26,VR63,VR29,VR29    \n\t"
			"|  VLDW *+AR0[528],  VR16  \n\t"
			"|  VLDW *+AR0[560],  VR17  \n\t"
				
			"   VLDW *+AR0[592], VR18   \n\t"
			"|  VLDW *+AR0[624], VR19   \n\t"

			"   VLDW *+AR0[656], VR20   \n\t"
			"|  VLDW *+AR0[688], VR21   \n\t" 

//B3
			"   VLDW *+AR0[720],  VR22  \n\t"
			"|  VLDW *+AR0[752],  VR23  \n\t"
			"|  VFMULD.M1   VR0,VR30,VR0    \n\t"
			"|  VFMULD.M2   VR1,VR31,VR1    \n\t"
			"|  VFMULD.M3   VR2,VR32,VR2    \n\t"

			"   VFMULD.M1   VR3,VR33,VR3    \n\t"
			"|  VFMULD.M2   VR4,VR34,VR4    \n\t"
			"|  VFMULD.M3   VR5,VR35,VR5    \n\t"

			"   VFMULD.M1   VR6,VR36,VR6    \n\t"
			"|  VFMULD.M2   VR7,VR37,VR7    \n\t"
			"|  VFMULD.M3   VR8,VR38,VR8    \n\t"

			"   VFMULAD.M1   VR27,VR63,VR28,VR28    \n\t"
			"|  VLDW *+AR0[784], VR24   \n\t"
			"|  VLDW *+AR0[816], VR25   \n\t"

			"   VLDW *+AR0[848], VR26   \n\t"
			"		SNOP 1\n\t"

			"   VFMULD.M1   VR9, VR39,VR9    \n\t"
			"|  VFMULD.M2   VR10,VR40,VR10    \n\t"
			"|  VFMULD.M3   VR11,VR41,VR11    \n\t"

			"   VFMULAD.M1   VR12,VR42,VR0,VR12    \n\t"
			"|  VFMULAD.M2   VR13,VR43,VR1,VR13    \n\t"
			"|  VFMULAD.M3   VR14,VR44,VR2,VR14    \n\t"

			"   VFMULAD.M1   VR15,VR45,VR3,VR15    \n\t"
			"|  VFMULAD.M2   VR16,VR46,VR4,VR16    \n\t"
			"|  VFMULAD.M3   VR17,VR47,VR5,VR17    \n\t"

			"   VFMULAD.M1   VR28,VR63,VR29,VR57    \n\t"

			"   VFMULAD.M1   VR18,VR48,VR6,VR18    \n\t"
			"|  VFMULAD.M2   VR19,VR49,VR7,VR19    \n\t"
			"|  VFMULAD.M3   VR20,VR50,VR8,VR20    \n\t"

			"   VFMULAD.M1   VR21,VR51,VR9, VR21    \n\t"
			"|  VFMULAD.M2   VR22,VR52,VR10,VR22    \n\t"
			"|  VFMULAD.M3   VR23,VR53,VR11,VR23    \n\t"

			"		SNOP 1		\n\t"  //vload vr24 9 

			"   VFMULAD.M1   VR24,VR54,VR12,VR24    \n\t"
			"|  VFMULAD.M2   VR25,VR55,VR13,VR25    \n\t"
			"|  VFMULAD.M3   VR26,VR56,VR14,VR26    \n\t"

			"		SNOP 1		\n\t"  //vload vr24 9 

			"   VSTW        VR57,*+AR2[64]       \n\t"

			"   VFMULAD.M1   VR15,VR63,VR18,VR27    \n\t"
			"|  VFMULAD.M2   VR16,VR63,VR19,VR28    \n\t"
			"|  VFMULAD.M3   VR17,VR63,VR20,VR29    \n\t"

			"		SNOP	2\n\t"

			"   VFMULAD.M1   VR24,VR63,VR21,VR24    \n\t"
			"|  VFMULAD.M2   VR25,VR63,VR22,VR25    \n\t"
			"|  VFMULAD.M3   VR26,VR63,VR23,VR26    \n\t"

			"		SNOP	5\n\t"

			"   VFMULAD.M1   VR24,VR63,VR27,VR27    \n\t"
			"|  VFMULAD.M2   VR25,VR63,VR28,VR28    \n\t"
			"|  VFMULAD.M3   VR26,VR63,VR29,VR29    \n\t"
			"		SNOP	5\n\t"
			"   VFMULAD.M1   VR27,VR63,VR28,VR28    \n\t"
			"		SNOP	5\n\t"
			"   VFMULAD.M1   VR28,VR63,VR29,VR57    \n\t"
			"		SNOP 5		\n\t" 
			"   VSTW        VR57,*+AR2[80]       \n\t"
			"		SNOP 3\n\t"
			:
			: [Vec]  "r" (vec),
			  [BE]	 "r" (BE_node),
			  [s0]   "r" (S0)
			: "R10","R11","R29","R30","R31","R32",\
			  "VR0","VR1","VR2","VR3","VR4","VR5","VR6","VR7","VR8",\
			  "VR9","VR10","VR11","VR12","VR13","VR14","VR15","VR16",\
			  "VR17","VR18","VR19","VR20","VR21","VR22","VR23","VR24",\
			  "VR25","VR26","VR27","VR28","VR29","VR30","VR31","VR32",\
			  "VR33","VR34","VR35","VR36","VR37","VR38","VR39","VR40",\
			  "VR41","VR42","VR43","VR44","VR45","VR46","VR47","VR48",\
			  "VR49","VR50","VR51","VR52","VR53","VR54","VR55","VR56",\
			  "VR57","VR58","VR59","VR60","VR61","VR62","VR63"
			);
}

