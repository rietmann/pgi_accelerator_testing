#include "omp.h"
#define t_max 1
/*
(w1_a[0][0]=((a[0][0][0][0][1]*(a[0][0][0][0][1]+1.0))*((a[0][0][0][0][1]+2.0)*0.16666666666666666)))
(w2_a[0][0]=(((a[0][0][0][0][1]-1.0)*(a[0][0][0][0][1]+1.0))*((a[0][0][0][0][1]+2.0)*-0.5)))
(w3_a[0][0]=(((a[0][0][0][0][1]-1.0)*a[0][0][0][0][1])*((a[0][0][0][0][1]+2.0)*0.5)))
(w4_a[0][0]=(((a[0][0][0][0][1]-1.0)*a[0][0][0][0][1])*((a[0][0][0][0][1]+1.0)*-0.16666666666666666)))
(w1_b[0][0]=((b[0][0][0][0][2]*(b[0][0][0][0][2]+1.0))*((b[0][0][0][0][2]+2.0)*0.16666666666666666)))
(w2_b[0][0]=(((b[0][0][0][0][2]-1.0)*(b[0][0][0][0][2]+1.0))*((b[0][0][0][0][2]+2.0)*-0.5)))
(w3_b[0][0]=(((b[0][0][0][0][2]-1.0)*b[0][0][0][0][2])*((b[0][0][0][0][2]+2.0)*0.5)))
(w4_b[0][0]=(((b[0][0][0][0][2]-1.0)*b[0][0][0][0][2])*((b[0][0][0][0][2]+1.0)*-0.16666666666666666)))
(w1_c[0][0]=((c[0][0][0][0][3]*(c[0][0][0][0][3]+1.0))*((c[0][0][0][0][3]+2.0)*0.16666666666666666)))
(w2_c[0][0]=(((c[0][0][0][0][3]-1.0)*(c[0][0][0][0][3]+1.0))*((c[0][0][0][0][3]+2.0)*-0.5)))
(w3_c[0][0]=(((c[0][0][0][0][3]-1.0)*c[0][0][0][0][3])*((c[0][0][0][0][3]+2.0)*0.5)))
(w4_c[0][0]=(((c[0][0][0][0][3]-1.0)*c[0][0][0][0][3])*((c[0][0][0][0][3]+1.0)*-0.16666666666666666)))
(u[0][0][0][1][0]=((((((((w1_a*w1_b)*(w1_c*u[-1][-1][-1][0][0]))+((w2_a*w1_b)*(w1_c*u[0][-1][-1][0][0])))+(((w3_a*w1_b)*(w1_c*u[1][-1][-1][0][0]))+((w4_a*w1_b)*(w1_c*u[2][-1][-1][0][0]))))+((((w1_a*w2_b)*(w1_c*u[-1][0][-1][0][0]))+((w2_a*w2_b)*(w1_c*u[0][0][-1][0][0])))+(((w3_a*w2_b)*(w1_c*u[1][0][-1][0][0]))+((w4_a*w2_b)*(w1_c*u[2][0][-1][0][0])))))+(((((w1_a*w3_b)*(w1_c*u[-1][1][-1][0][0]))+((w2_a*w3_b)*(w1_c*u[0][1][-1][0][0])))+(((w3_a*w3_b)*(w1_c*u[1][1][-1][0][0]))+((w4_a*w3_b)*(w1_c*u[2][1][-1][0][0]))))+((((w1_a*w4_b)*(w1_c*u[-1][2][-1][0][0]))+((w2_a*w4_b)*(w1_c*u[0][2][-1][0][0])))+(((w3_a*w4_b)*(w1_c*u[1][2][-1][0][0]))+((w4_a*w4_b)*(w1_c*u[2][2][-1][0][0]))))))+((((((w1_a*w1_b)*(w2_c*u[-1][-1][0][0][0]))+((w2_a*w1_b)*(w2_c*u[0][-1][0][0][0])))+(((w3_a*w1_b)*(w2_c*u[1][-1][0][0][0]))+((w4_a*w1_b)*(w2_c*u[2][-1][0][0][0]))))+((((w1_a*w2_b)*(w2_c*u[-1][0][0][0][0]))+((w2_a*w2_b)*(w2_c*u[0][0][0][0][0])))+(((w3_a*w2_b)*(w2_c*u[1][0][0][0][0]))+((w4_a*w2_b)*(w2_c*u[2][0][0][0][0])))))+(((((w1_a*w3_b)*(w2_c*u[-1][1][0][0][0]))+((w2_a*w3_b)*(w2_c*u[0][1][0][0][0])))+(((w3_a*w3_b)*(w2_c*u[1][1][0][0][0]))+((w4_a*w3_b)*(w2_c*u[2][1][0][0][0]))))+((((w1_a*w4_b)*(w2_c*u[-1][2][0][0][0]))+((w2_a*w4_b)*(w2_c*u[0][2][0][0][0])))+(((w3_a*w4_b)*(w2_c*u[1][2][0][0][0]))+((w4_a*w4_b)*(w2_c*u[2][2][0][0][0])))))))+(((((((w1_a*w1_b)*(w3_c*u[-1][-1][1][0][0]))+((w2_a*w1_b)*(w3_c*u[0][-1][1][0][0])))+(((w3_a*w1_b)*(w3_c*u[1][-1][1][0][0]))+((w4_a*w1_b)*(w3_c*u[2][-1][1][0][0]))))+((((w1_a*w2_b)*(w3_c*u[-1][0][1][0][0]))+((w2_a*w2_b)*(w3_c*u[0][0][1][0][0])))+(((w3_a*w2_b)*(w3_c*u[1][0][1][0][0]))+((w4_a*w2_b)*(w3_c*u[2][0][1][0][0])))))+(((((w1_a*w3_b)*(w3_c*u[-1][1][1][0][0]))+((w2_a*w3_b)*(w3_c*u[0][1][1][0][0])))+(((w3_a*w3_b)*(w3_c*u[1][1][1][0][0]))+((w4_a*w3_b)*(w3_c*u[2][1][1][0][0]))))+((((w1_a*w4_b)*(w3_c*u[-1][2][1][0][0]))+((w2_a*w4_b)*(w3_c*u[0][2][1][0][0])))+(((w3_a*w4_b)*(w3_c*u[1][2][1][0][0]))+((w4_a*w4_b)*(w3_c*u[2][2][1][0][0]))))))+((((((w1_a*w1_b)*(w4_c*u[-1][-1][2][0][0]))+((w2_a*w1_b)*(w4_c*u[0][-1][2][0][0])))+(((w3_a*w1_b)*(w4_c*u[1][-1][2][0][0]))+((w4_a*w1_b)*(w4_c*u[2][-1][2][0][0]))))+((((w1_a*w2_b)*(w4_c*u[-1][0][2][0][0]))+((w2_a*w2_b)*(w4_c*u[0][0][2][0][0])))+(((w3_a*w2_b)*(w4_c*u[1][0][2][0][0]))+((w4_a*w2_b)*(w4_c*u[2][0][2][0][0])))))+(((((w1_a*w3_b)*(w4_c*u[-1][1][2][0][0]))+((w2_a*w3_b)*(w4_c*u[0][1][2][0][0])))+(((w3_a*w3_b)*(w4_c*u[1][1][2][0][0]))+((w4_a*w3_b)*(w4_c*u[2][1][2][0][0]))))+((((w1_a*w4_b)*(w4_c*u[-1][2][2][0][0]))+((w2_a*w4_b)*(w4_c*u[0][2][2][0][0])))+(((w3_a*w4_b)*(w4_c*u[1][2][2][0][0]))+((w4_a*w4_b)*(w4_c*u[2][2][2][0][0])))))))))

*/
void tricubic_interpolation(double *  *  u_0_1_out, double *  u_0_0, double *  u_0_1, double *  a_1_0, double *  b_2_0, double *  c_3_0, int x_max, int y_max, int z_max)
{
	int _idx0;
	int _idx1;
	int _idx10;
	int _idx11;
	int _idx12;
	int _idx13;
	int _idx14;
	int _idx15;
	int _idx16;
	int _idx17;
	int _idx18;
	int _idx19;
	int _idx2;
	int _idx20;
	int _idx21;
	int _idx22;
	int _idx23;
	int _idx24;
	int _idx25;
	int _idx26;
	int _idx27;
	int _idx28;
	int _idx29;
	int _idx3;
	int _idx30;
	int _idx31;
	int _idx32;
	int _idx33;
	int _idx34;
	int _idx35;
	int _idx36;
	int _idx37;
	int _idx38;
	int _idx39;
	int _idx4;
	int _idx40;
	int _idx41;
	int _idx42;
	int _idx43;
	int _idx44;
	int _idx45;
	int _idx46;
	int _idx47;
	int _idx48;
	int _idx49;
	int _idx5;
	int _idx50;
	int _idx51;
	int _idx52;
	int _idx53;
	int _idx54;
	int _idx55;
	int _idx56;
	int _idx57;
	int _idx58;
	int _idx59;
	int _idx6;
	int _idx60;
	int _idx61;
	int _idx62;
	int _idx63;
	int _idx64;
	int _idx7;
	int _idx8;
	int _idx9;
	int p_idx_x;
	int p_idx_y;
	int p_idx_z;
	int t;
	const double *  __restrict__ const u__a_1[1] =  { a_1_0 } ;
	const double *  __restrict__ const u__b_2[1] =  { b_2_0 } ;
	const double *  __restrict__ const u__c_3[1] =  { c_3_0 } ;
	double *  __restrict__ const u__u_0[1] =  { u_0_0, u_0_1 } ;
	double w1_a;
	double w1_b;
	double w1_c;
	double w2_a;
	double w2_b;
	double w2_c;
	double w3_a;
	double w3_b;
	double w3_c;
	double w4_a;
	double w4_b;
	double w4_c;
	/*
	Initializations
	*/
	/*
	Implementation
	*/
	/*
	for t = 1..t_max by 1 parallel 1 <level 0> schedule  { ... }
	*/
	for (t=1; t<=t_max; t+=1)
	{
		/*
		for POINT p[t=t, s=(1, 1, 1)][0] of size [1, 1, 1] in u[t=t, s=(:, :, :)][0] parallel 1 <level 0> schedule default { ... }
		*/
		{
			/* Index bounds calculations for iterators in p[t=t, s=(1, 1, 1)][0] */
			for (p_idx_z=0; p_idx_z<z_max; p_idx_z+=1)
			{
				for (p_idx_y=0; p_idx_y<y_max; p_idx_y+=1)
				{
					for (p_idx_x=0; p_idx_x<x_max; p_idx_x+=1)
					{
						/* Index bounds calculations for iterators in p[t=t, s=(1, 1, 1)][0] */
						/*
						u[t=(t+1), s=p[t=?, s=?][0]][0]=stencil(u[t=t, s=p[t=?, s=?][0]][0])
						*/
						/* _idx0 = ((((p_idx_z*x_max)*y_max)+(p_idx_y*x_max))+p_idx_x) */
						{
							int __tmp0 = (p_idx_z*x_max);
							int __tmp1 = (__tmp0*y_max);
							int __tmp2 = (p_idx_y*x_max);
							int __tmp3 = (__tmp1+__tmp2);
							int __tmp4 = (__tmp3+p_idx_x);
							_idx0=__tmp4;
						}
						{
							double __tmp5 = (u__a_1[0][_idx0]+1.0);
							double __tmp6 = (u__a_1[0][_idx0]*__tmp5);
							double __tmp7 = (u__a_1[0][_idx0]+2.0);
							double __tmp8 = (__tmp7*0.16666666666666666);
							double __tmp9 = (__tmp6*__tmp8);
							w1_a=__tmp9;
						}
						{
							double __tmp10 = (u__a_1[0][_idx0]-1.0);
							double __tmp11 = (u__a_1[0][_idx0]+1.0);
							double __tmp12 = (__tmp10*__tmp11);
							double __tmp13 = (u__a_1[0][_idx0]+2.0);
							double __tmp14 = (__tmp13*-0.5);
							double __tmp15 = (__tmp12*__tmp14);
							w2_a=__tmp15;
						}
						{
							double __tmp16 = (u__a_1[0][_idx0]-1.0);
							double __tmp17 = (__tmp16*u__a_1[0][_idx0]);
							double __tmp18 = (u__a_1[0][_idx0]+2.0);
							double __tmp19 = (__tmp18*0.5);
							double __tmp20 = (__tmp17*__tmp19);
							w3_a=__tmp20;
						}
						{
							double __tmp21 = (u__a_1[0][_idx0]-1.0);
							double __tmp22 = (__tmp21*u__a_1[0][_idx0]);
							double __tmp23 = (u__a_1[0][_idx0]+1.0);
							double __tmp24 = (__tmp23*-0.16666666666666666);
							double __tmp25 = (__tmp22*__tmp24);
							w4_a=__tmp25;
						}
						{
							double __tmp26 = (u__b_2[0][_idx0]+1.0);
							double __tmp27 = (u__b_2[0][_idx0]*__tmp26);
							double __tmp28 = (u__b_2[0][_idx0]+2.0);
							double __tmp29 = (__tmp28*0.16666666666666666);
							double __tmp30 = (__tmp27*__tmp29);
							w1_b=__tmp30;
						}
						{
							double __tmp31 = (u__b_2[0][_idx0]-1.0);
							double __tmp32 = (u__b_2[0][_idx0]+1.0);
							double __tmp33 = (__tmp31*__tmp32);
							double __tmp34 = (u__b_2[0][_idx0]+2.0);
							double __tmp35 = (__tmp34*-0.5);
							double __tmp36 = (__tmp33*__tmp35);
							w2_b=__tmp36;
						}
						{
							double __tmp37 = (u__b_2[0][_idx0]-1.0);
							double __tmp38 = (__tmp37*u__b_2[0][_idx0]);
							double __tmp39 = (u__b_2[0][_idx0]+2.0);
							double __tmp40 = (__tmp39*0.5);
							double __tmp41 = (__tmp38*__tmp40);
							w3_b=__tmp41;
						}
						{
							double __tmp42 = (u__b_2[0][_idx0]-1.0);
							double __tmp43 = (__tmp42*u__b_2[0][_idx0]);
							double __tmp44 = (u__b_2[0][_idx0]+1.0);
							double __tmp45 = (__tmp44*-0.16666666666666666);
							double __tmp46 = (__tmp43*__tmp45);
							w4_b=__tmp46;
						}
						{
							double __tmp47 = (u__c_3[0][_idx0]+1.0);
							double __tmp48 = (u__c_3[0][_idx0]*__tmp47);
							double __tmp49 = (u__c_3[0][_idx0]+2.0);
							double __tmp50 = (__tmp49*0.16666666666666666);
							double __tmp51 = (__tmp48*__tmp50);
							w1_c=__tmp51;
						}
						{
							double __tmp52 = (u__c_3[0][_idx0]-1.0);
							double __tmp53 = (u__c_3[0][_idx0]+1.0);
							double __tmp54 = (__tmp52*__tmp53);
							double __tmp55 = (u__c_3[0][_idx0]+2.0);
							double __tmp56 = (__tmp55*-0.5);
							double __tmp57 = (__tmp54*__tmp56);
							w2_c=__tmp57;
						}
						{
							double __tmp58 = (u__c_3[0][_idx0]-1.0);
							double __tmp59 = (__tmp58*u__c_3[0][_idx0]);
							double __tmp60 = (u__c_3[0][_idx0]+2.0);
							double __tmp61 = (__tmp60*0.5);
							double __tmp62 = (__tmp59*__tmp61);
							w3_c=__tmp62;
						}
						{
							double __tmp63 = (u__c_3[0][_idx0]-1.0);
							double __tmp64 = (__tmp63*u__c_3[0][_idx0]);
							double __tmp65 = (u__c_3[0][_idx0]+1.0);
							double __tmp66 = (__tmp65*-0.16666666666666666);
							double __tmp67 = (__tmp64*__tmp66);
							w4_c=__tmp67;
						}
						/* _idx1 = (((((((p_idx_z*x_max)+(3*p_idx_z))*y_max)+(((3*p_idx_z)+p_idx_y)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x) */
						{
							int __tmp68 = (3*p_idx_z);
							int __tmp69 = (__tmp68*y_max);
							int __tmp70 = (_idx0+__tmp69);
							int __tmp71 = (3*p_idx_z);
							int __tmp72 = (__tmp71*x_max);
							int __tmp73 = (__tmp70+__tmp72);
							int __tmp74 = (9*p_idx_z);
							int __tmp75 = (__tmp73+__tmp74);
							int __tmp76 = (3*p_idx_y);
							int __tmp77 = (__tmp75+__tmp76);
							_idx1=__tmp77;
						}
						/* _idx2 = ((((((((p_idx_z*x_max)+(3*p_idx_z))*y_max)+(((3*p_idx_z)+p_idx_y)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+1) */
						{
							int __tmp78 = (_idx1+1);
							_idx2=__tmp78;
						}
						/* _idx3 = ((((((((p_idx_z*x_max)+(3*p_idx_z))*y_max)+(((3*p_idx_z)+p_idx_y)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+2) */
						{
							int __tmp79 = (_idx2+1);
							_idx3=__tmp79;
						}
						/* _idx4 = ((((((((p_idx_z*x_max)+(3*p_idx_z))*y_max)+(((3*p_idx_z)+p_idx_y)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+3) */
						{
							int __tmp80 = (_idx3+1);
							_idx4=__tmp80;
						}
						/* _idx5 = ((((((((p_idx_z*x_max)+(3*p_idx_z))*y_max)+((((3*p_idx_z)+p_idx_y)+1)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+3) */
						{
							int __tmp81 = (_idx3+x_max);
							int __tmp82 = (__tmp81+1);
							_idx5=__tmp82;
						}
						/* _idx6 = ((((((((p_idx_z*x_max)+(3*p_idx_z))*y_max)+((((3*p_idx_z)+p_idx_y)+1)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+4) */
						{
							int __tmp83 = (_idx3+x_max);
							int __tmp84 = (__tmp83+2);
							_idx6=__tmp84;
						}
						/* _idx7 = ((((((((p_idx_z*x_max)+(3*p_idx_z))*y_max)+((((3*p_idx_z)+p_idx_y)+1)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+5) */
						{
							int __tmp85 = (_idx3+x_max);
							int __tmp86 = (__tmp85+3);
							_idx7=__tmp86;
						}
						/* _idx8 = ((((((((p_idx_z*x_max)+(3*p_idx_z))*y_max)+((((3*p_idx_z)+p_idx_y)+1)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+6) */
						{
							int __tmp87 = (_idx3+x_max);
							int __tmp88 = (__tmp87+4);
							_idx8=__tmp88;
						}
						/* _idx9 = ((((((((p_idx_z*x_max)+(3*p_idx_z))*y_max)+((((3*p_idx_z)+p_idx_y)+2)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+6) */
						{
							int __tmp89 = (_idx6+x_max);
							int __tmp90 = (__tmp89+2);
							_idx9=__tmp90;
						}
						/* _idx10 = ((((((((p_idx_z*x_max)+(3*p_idx_z))*y_max)+((((3*p_idx_z)+p_idx_y)+2)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+7) */
						{
							int __tmp91 = (_idx6+x_max);
							int __tmp92 = (__tmp91+3);
							_idx10=__tmp92;
						}
						/* _idx11 = ((((((((p_idx_z*x_max)+(3*p_idx_z))*y_max)+((((3*p_idx_z)+p_idx_y)+2)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+8) */
						{
							int __tmp93 = (_idx10+1);
							_idx11=__tmp93;
						}
						/* _idx12 = ((((((((p_idx_z*x_max)+(3*p_idx_z))*y_max)+((((3*p_idx_z)+p_idx_y)+2)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+9) */
						{
							int __tmp94 = (_idx6+x_max);
							int __tmp95 = (__tmp94+5);
							_idx12=__tmp95;
						}
						/* _idx13 = ((((((((p_idx_z*x_max)+(3*p_idx_z))*y_max)+((((3*p_idx_z)+p_idx_y)+3)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+9) */
						{
							int __tmp96 = (_idx11+x_max);
							int __tmp97 = (__tmp96+1);
							_idx13=__tmp97;
						}
						/* _idx14 = ((((((((p_idx_z*x_max)+(3*p_idx_z))*y_max)+((((3*p_idx_z)+p_idx_y)+3)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+10) */
						{
							int __tmp98 = (_idx11+x_max);
							int __tmp99 = (__tmp98+2);
							_idx14=__tmp99;
						}
						/* _idx15 = ((((((((p_idx_z*x_max)+(3*p_idx_z))*y_max)+((((3*p_idx_z)+p_idx_y)+3)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+11) */
						{
							int __tmp100 = (_idx11+x_max);
							int __tmp101 = (__tmp100+3);
							_idx15=__tmp101;
						}
						/* _idx16 = ((((((((p_idx_z*x_max)+(3*p_idx_z))*y_max)+((((3*p_idx_z)+p_idx_y)+3)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+12) */
						{
							int __tmp102 = (_idx11+x_max);
							int __tmp103 = (__tmp102+4);
							_idx16=__tmp103;
						}
						/* _idx17 = ((((((((((p_idx_z+1)*x_max)+(3*p_idx_z))+3)*y_max)+((((3*p_idx_z)+p_idx_y)+3)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+9) */
						{
							int __tmp104 = (x_max+3);
							int __tmp105 = (__tmp104*y_max);
							int __tmp106 = (_idx13+__tmp105);
							_idx17=__tmp106;
						}
						/* _idx18 = ((((((((((p_idx_z+1)*x_max)+(3*p_idx_z))+3)*y_max)+((((3*p_idx_z)+p_idx_y)+3)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+10) */
						{
							int __tmp107 = (_idx17+1);
							_idx18=__tmp107;
						}
						/* _idx19 = ((((((((((p_idx_z+1)*x_max)+(3*p_idx_z))+3)*y_max)+((((3*p_idx_z)+p_idx_y)+3)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+11) */
						{
							int __tmp108 = (_idx18+1);
							_idx19=__tmp108;
						}
						/* _idx20 = ((((((((((p_idx_z+1)*x_max)+(3*p_idx_z))+3)*y_max)+((((3*p_idx_z)+p_idx_y)+3)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+12) */
						{
							int __tmp109 = (_idx18+2);
							_idx20=__tmp109;
						}
						/* _idx21 = ((((((((((p_idx_z+1)*x_max)+(3*p_idx_z))+3)*y_max)+((((3*p_idx_z)+p_idx_y)+4)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+12) */
						{
							int __tmp110 = (_idx18+x_max);
							int __tmp111 = (__tmp110+2);
							_idx21=__tmp111;
						}
						/* _idx22 = ((((((((((p_idx_z+1)*x_max)+(3*p_idx_z))+3)*y_max)+((((3*p_idx_z)+p_idx_y)+4)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+13) */
						{
							int __tmp112 = (_idx18+x_max);
							int __tmp113 = (__tmp112+3);
							_idx22=__tmp113;
						}
						/* _idx23 = ((((((((((p_idx_z+1)*x_max)+(3*p_idx_z))+3)*y_max)+((((3*p_idx_z)+p_idx_y)+4)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+14) */
						{
							int __tmp114 = (_idx18+x_max);
							int __tmp115 = (__tmp114+4);
							_idx23=__tmp115;
						}
						/* _idx24 = ((((((((((p_idx_z+1)*x_max)+(3*p_idx_z))+3)*y_max)+((((3*p_idx_z)+p_idx_y)+4)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+15) */
						{
							int __tmp116 = (_idx18+x_max);
							int __tmp117 = (__tmp116+5);
							_idx24=__tmp117;
						}
						/* _idx25 = ((((((((((p_idx_z+1)*x_max)+(3*p_idx_z))+3)*y_max)+((((3*p_idx_z)+p_idx_y)+5)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+15) */
						{
							int __tmp118 = (_idx22+x_max);
							int __tmp119 = (__tmp118+2);
							_idx25=__tmp119;
						}
						/* _idx26 = ((((((((((p_idx_z+1)*x_max)+(3*p_idx_z))+3)*y_max)+((((3*p_idx_z)+p_idx_y)+5)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+16) */
						{
							int __tmp120 = (_idx25+1);
							_idx26=__tmp120;
						}
						/* _idx27 = ((((((((((p_idx_z+1)*x_max)+(3*p_idx_z))+3)*y_max)+((((3*p_idx_z)+p_idx_y)+5)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+17) */
						{
							int __tmp121 = (_idx25+2);
							_idx27=__tmp121;
						}
						/* _idx28 = ((((((((((p_idx_z+1)*x_max)+(3*p_idx_z))+3)*y_max)+((((3*p_idx_z)+p_idx_y)+5)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+18) */
						{
							int __tmp122 = (_idx27+1);
							_idx28=__tmp122;
						}
						/* _idx29 = ((((((((((p_idx_z+1)*x_max)+(3*p_idx_z))+3)*y_max)+((((3*p_idx_z)+p_idx_y)+6)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+18) */
						{
							int __tmp123 = (_idx27+x_max);
							int __tmp124 = (__tmp123+1);
							_idx29=__tmp124;
						}
						/* _idx30 = ((((((((((p_idx_z+1)*x_max)+(3*p_idx_z))+3)*y_max)+((((3*p_idx_z)+p_idx_y)+6)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+19) */
						{
							int __tmp125 = (_idx27+x_max);
							int __tmp126 = (__tmp125+2);
							_idx30=__tmp126;
						}
						/* _idx31 = ((((((((((p_idx_z+1)*x_max)+(3*p_idx_z))+3)*y_max)+((((3*p_idx_z)+p_idx_y)+6)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+20) */
						{
							int __tmp127 = (_idx27+x_max);
							int __tmp128 = (__tmp127+3);
							_idx31=__tmp128;
						}
						/* _idx32 = ((((((((((p_idx_z+1)*x_max)+(3*p_idx_z))+3)*y_max)+((((3*p_idx_z)+p_idx_y)+6)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+21) */
						{
							int __tmp129 = (_idx27+x_max);
							int __tmp130 = (__tmp129+4);
							_idx32=__tmp130;
						}
						/* _idx33 = ((((((((((p_idx_z+2)*x_max)+(3*p_idx_z))+6)*y_max)+((((3*p_idx_z)+p_idx_y)+6)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+18) */
						{
							int __tmp131 = (x_max+3);
							int __tmp132 = (__tmp131*y_max);
							int __tmp133 = (_idx29+__tmp132);
							_idx33=__tmp133;
						}
						/* _idx34 = ((((((((((p_idx_z+2)*x_max)+(3*p_idx_z))+6)*y_max)+((((3*p_idx_z)+p_idx_y)+6)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+19) */
						{
							int __tmp134 = (_idx33+1);
							_idx34=__tmp134;
						}
						/* _idx35 = ((((((((((p_idx_z+2)*x_max)+(3*p_idx_z))+6)*y_max)+((((3*p_idx_z)+p_idx_y)+6)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+20) */
						{
							int __tmp135 = (_idx34+1);
							_idx35=__tmp135;
						}
						/* _idx36 = ((((((((((p_idx_z+2)*x_max)+(3*p_idx_z))+6)*y_max)+((((3*p_idx_z)+p_idx_y)+6)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+21) */
						{
							int __tmp136 = (_idx35+1);
							_idx36=__tmp136;
						}
						/* _idx37 = ((((((((((p_idx_z+2)*x_max)+(3*p_idx_z))+6)*y_max)+((((3*p_idx_z)+p_idx_y)+7)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+21) */
						{
							int __tmp137 = (_idx35+x_max);
							int __tmp138 = (__tmp137+1);
							_idx37=__tmp138;
						}
						/* _idx38 = ((((((((((p_idx_z+2)*x_max)+(3*p_idx_z))+6)*y_max)+((((3*p_idx_z)+p_idx_y)+7)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+22) */
						{
							int __tmp139 = (_idx35+x_max);
							int __tmp140 = (__tmp139+2);
							_idx38=__tmp140;
						}
						/* _idx39 = ((((((((((p_idx_z+2)*x_max)+(3*p_idx_z))+6)*y_max)+((((3*p_idx_z)+p_idx_y)+7)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+23) */
						{
							int __tmp141 = (_idx35+x_max);
							int __tmp142 = (__tmp141+3);
							_idx39=__tmp142;
						}
						/* _idx40 = ((((((((((p_idx_z+2)*x_max)+(3*p_idx_z))+6)*y_max)+((((3*p_idx_z)+p_idx_y)+7)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+24) */
						{
							int __tmp143 = (_idx35+x_max);
							int __tmp144 = (__tmp143+4);
							_idx40=__tmp144;
						}
						/* _idx41 = ((((((((((p_idx_z+2)*x_max)+(3*p_idx_z))+6)*y_max)+((((3*p_idx_z)+p_idx_y)+8)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+24) */
						{
							int __tmp145 = (_idx39+x_max);
							int __tmp146 = (__tmp145+1);
							_idx41=__tmp146;
						}
						/* _idx42 = ((((((((((p_idx_z+2)*x_max)+(3*p_idx_z))+6)*y_max)+((((3*p_idx_z)+p_idx_y)+8)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+25) */
						{
							int __tmp147 = (_idx39+x_max);
							int __tmp148 = (__tmp147+2);
							_idx42=__tmp148;
						}
						/* _idx43 = ((((((((((p_idx_z+2)*x_max)+(3*p_idx_z))+6)*y_max)+((((3*p_idx_z)+p_idx_y)+8)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+26) */
						{
							int __tmp149 = (_idx42+1);
							_idx43=__tmp149;
						}
						/* _idx44 = ((((((((((p_idx_z+2)*x_max)+(3*p_idx_z))+6)*y_max)+((((3*p_idx_z)+p_idx_y)+8)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+27) */
						{
							int __tmp150 = (_idx42+2);
							_idx44=__tmp150;
						}
						/* _idx45 = ((((((((((p_idx_z+2)*x_max)+(3*p_idx_z))+6)*y_max)+((((3*p_idx_z)+p_idx_y)+9)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+27) */
						{
							int __tmp151 = (_idx42+x_max);
							int __tmp152 = (__tmp151+2);
							_idx45=__tmp152;
						}
						/* _idx46 = ((((((((((p_idx_z+2)*x_max)+(3*p_idx_z))+6)*y_max)+((((3*p_idx_z)+p_idx_y)+9)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+28) */
						{
							int __tmp153 = (_idx42+x_max);
							int __tmp154 = (__tmp153+3);
							_idx46=__tmp154;
						}
						/* _idx47 = ((((((((((p_idx_z+2)*x_max)+(3*p_idx_z))+6)*y_max)+((((3*p_idx_z)+p_idx_y)+9)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+29) */
						{
							int __tmp155 = (_idx42+x_max);
							int __tmp156 = (__tmp155+4);
							_idx47=__tmp156;
						}
						/* _idx48 = ((((((((((p_idx_z+2)*x_max)+(3*p_idx_z))+6)*y_max)+((((3*p_idx_z)+p_idx_y)+9)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+30) */
						{
							int __tmp157 = (_idx47+1);
							_idx48=__tmp157;
						}
						/* _idx49 = ((((((((((p_idx_z+3)*x_max)+(3*p_idx_z))+9)*y_max)+((((3*p_idx_z)+p_idx_y)+9)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+27) */
						{
							int __tmp158 = (x_max+3);
							int __tmp159 = (__tmp158*y_max);
							int __tmp160 = (_idx45+__tmp159);
							_idx49=__tmp160;
						}
						/* _idx50 = ((((((((((p_idx_z+3)*x_max)+(3*p_idx_z))+9)*y_max)+((((3*p_idx_z)+p_idx_y)+9)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+28) */
						{
							int __tmp161 = (_idx49+1);
							_idx50=__tmp161;
						}
						/* _idx51 = ((((((((((p_idx_z+3)*x_max)+(3*p_idx_z))+9)*y_max)+((((3*p_idx_z)+p_idx_y)+9)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+29) */
						{
							int __tmp162 = (_idx49+2);
							_idx51=__tmp162;
						}
						/* _idx52 = ((((((((((p_idx_z+3)*x_max)+(3*p_idx_z))+9)*y_max)+((((3*p_idx_z)+p_idx_y)+9)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+30) */
						{
							int __tmp163 = (_idx51+1);
							_idx52=__tmp163;
						}
						/* _idx53 = ((((((((((p_idx_z+3)*x_max)+(3*p_idx_z))+9)*y_max)+((((3*p_idx_z)+p_idx_y)+10)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+30) */
						{
							int __tmp164 = (_idx51+x_max);
							int __tmp165 = (__tmp164+1);
							_idx53=__tmp165;
						}
						/* _idx54 = ((((((((((p_idx_z+3)*x_max)+(3*p_idx_z))+9)*y_max)+((((3*p_idx_z)+p_idx_y)+10)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+31) */
						{
							int __tmp166 = (_idx51+x_max);
							int __tmp167 = (__tmp166+2);
							_idx54=__tmp167;
						}
						/* _idx55 = ((((((((((p_idx_z+3)*x_max)+(3*p_idx_z))+9)*y_max)+((((3*p_idx_z)+p_idx_y)+10)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+32) */
						{
							int __tmp168 = (_idx51+x_max);
							int __tmp169 = (__tmp168+3);
							_idx55=__tmp169;
						}
						/* _idx56 = ((((((((((p_idx_z+3)*x_max)+(3*p_idx_z))+9)*y_max)+((((3*p_idx_z)+p_idx_y)+10)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+33) */
						{
							int __tmp170 = (_idx55+1);
							_idx56=__tmp170;
						}
						/* _idx57 = ((((((((((p_idx_z+3)*x_max)+(3*p_idx_z))+9)*y_max)+((((3*p_idx_z)+p_idx_y)+11)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+33) */
						{
							int __tmp171 = (_idx55+x_max);
							int __tmp172 = (__tmp171+1);
							_idx57=__tmp172;
						}
						/* _idx58 = ((((((((((p_idx_z+3)*x_max)+(3*p_idx_z))+9)*y_max)+((((3*p_idx_z)+p_idx_y)+11)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+34) */
						{
							int __tmp173 = (_idx55+x_max);
							int __tmp174 = (__tmp173+2);
							_idx58=__tmp174;
						}
						/* _idx59 = ((((((((((p_idx_z+3)*x_max)+(3*p_idx_z))+9)*y_max)+((((3*p_idx_z)+p_idx_y)+11)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+35) */
						{
							int __tmp175 = (_idx58+1);
							_idx59=__tmp175;
						}
						/* _idx60 = ((((((((((p_idx_z+3)*x_max)+(3*p_idx_z))+9)*y_max)+((((3*p_idx_z)+p_idx_y)+11)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+36) */
						{
							int __tmp176 = (_idx58+2);
							_idx60=__tmp176;
						}
						/* _idx61 = ((((((((((p_idx_z+3)*x_max)+(3*p_idx_z))+9)*y_max)+((((3*p_idx_z)+p_idx_y)+12)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+36) */
						{
							int __tmp177 = (_idx58+x_max);
							int __tmp178 = (__tmp177+2);
							_idx61=__tmp178;
						}
						/* _idx62 = ((((((((((p_idx_z+3)*x_max)+(3*p_idx_z))+9)*y_max)+((((3*p_idx_z)+p_idx_y)+12)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+37) */
						{
							int __tmp179 = (_idx58+x_max);
							int __tmp180 = (__tmp179+3);
							_idx62=__tmp180;
						}
						/* _idx63 = ((((((((((p_idx_z+3)*x_max)+(3*p_idx_z))+9)*y_max)+((((3*p_idx_z)+p_idx_y)+12)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+38) */
						{
							int __tmp181 = (_idx58+x_max);
							int __tmp182 = (__tmp181+4);
							_idx63=__tmp182;
						}
						/* _idx64 = ((((((((((p_idx_z+3)*x_max)+(3*p_idx_z))+9)*y_max)+((((3*p_idx_z)+p_idx_y)+12)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+39) */
						{
							int __tmp183 = (_idx58+x_max);
							int __tmp184 = (__tmp183+5);
							_idx64=__tmp184;
						}
						u__u_0[1][_idx22]=((((((((w1_a*w1_b)*(w1_c*u__u_0[0][_idx1]))+((w2_a*w1_b)*(w1_c*u__u_0[0][_idx2])))+(((w3_a*w1_b)*(w1_c*u__u_0[0][_idx3]))+((w4_a*w1_b)*(w1_c*u__u_0[0][_idx4]))))+((((w1_a*w2_b)*(w1_c*u__u_0[0][_idx5]))+((w2_a*w2_b)*(w1_c*u__u_0[0][_idx6])))+(((w3_a*w2_b)*(w1_c*u__u_0[0][_idx7]))+((w4_a*w2_b)*(w1_c*u__u_0[0][_idx8])))))+(((((w1_a*w3_b)*(w1_c*u__u_0[0][_idx9]))+((w2_a*w3_b)*(w1_c*u__u_0[0][_idx10])))+(((w3_a*w3_b)*(w1_c*u__u_0[0][_idx11]))+((w4_a*w3_b)*(w1_c*u__u_0[0][_idx12]))))+((((w1_a*w4_b)*(w1_c*u__u_0[0][_idx13]))+((w2_a*w4_b)*(w1_c*u__u_0[0][_idx14])))+(((w3_a*w4_b)*(w1_c*u__u_0[0][_idx15]))+((w4_a*w4_b)*(w1_c*u__u_0[0][_idx16]))))))+((((((w1_a*w1_b)*(w2_c*u__u_0[0][_idx17]))+((w2_a*w1_b)*(w2_c*u__u_0[0][_idx18])))+(((w3_a*w1_b)*(w2_c*u__u_0[0][_idx19]))+((w4_a*w1_b)*(w2_c*u__u_0[0][_idx20]))))+((((w1_a*w2_b)*(w2_c*u__u_0[0][_idx21]))+((w2_a*w2_b)*(w2_c*u__u_0[0][_idx22])))+(((w3_a*w2_b)*(w2_c*u__u_0[0][_idx23]))+((w4_a*w2_b)*(w2_c*u__u_0[0][_idx24])))))+(((((w1_a*w3_b)*(w2_c*u__u_0[0][_idx25]))+((w2_a*w3_b)*(w2_c*u__u_0[0][_idx26])))+(((w3_a*w3_b)*(w2_c*u__u_0[0][_idx27]))+((w4_a*w3_b)*(w2_c*u__u_0[0][_idx28]))))+((((w1_a*w4_b)*(w2_c*u__u_0[0][_idx29]))+((w2_a*w4_b)*(w2_c*u__u_0[0][_idx30])))+(((w3_a*w4_b)*(w2_c*u__u_0[0][_idx31]))+((w4_a*w4_b)*(w2_c*u__u_0[0][_idx32])))))))+(((((((w1_a*w1_b)*(w3_c*u__u_0[0][_idx33]))+((w2_a*w1_b)*(w3_c*u__u_0[0][_idx34])))+(((w3_a*w1_b)*(w3_c*u__u_0[0][_idx35]))+((w4_a*w1_b)*(w3_c*u__u_0[0][_idx36]))))+((((w1_a*w2_b)*(w3_c*u__u_0[0][_idx37]))+((w2_a*w2_b)*(w3_c*u__u_0[0][_idx38])))+(((w3_a*w2_b)*(w3_c*u__u_0[0][_idx39]))+((w4_a*w2_b)*(w3_c*u__u_0[0][_idx40])))))+(((((w1_a*w3_b)*(w3_c*u__u_0[0][_idx41]))+((w2_a*w3_b)*(w3_c*u__u_0[0][_idx42])))+(((w3_a*w3_b)*(w3_c*u__u_0[0][_idx43]))+((w4_a*w3_b)*(w3_c*u__u_0[0][_idx44]))))+((((w1_a*w4_b)*(w3_c*u__u_0[0][_idx45]))+((w2_a*w4_b)*(w3_c*u__u_0[0][_idx46])))+(((w3_a*w4_b)*(w3_c*u__u_0[0][_idx47]))+((w4_a*w4_b)*(w3_c*u__u_0[0][_idx48]))))))+((((((w1_a*w1_b)*(w4_c*u__u_0[0][_idx49]))+((w2_a*w1_b)*(w4_c*u__u_0[0][_idx50])))+(((w3_a*w1_b)*(w4_c*u__u_0[0][_idx51]))+((w4_a*w1_b)*(w4_c*u__u_0[0][_idx52]))))+((((w1_a*w2_b)*(w4_c*u__u_0[0][_idx53]))+((w2_a*w2_b)*(w4_c*u__u_0[0][_idx54])))+(((w3_a*w2_b)*(w4_c*u__u_0[0][_idx55]))+((w4_a*w2_b)*(w4_c*u__u_0[0][_idx56])))))+(((((w1_a*w3_b)*(w4_c*u__u_0[0][_idx57]))+((w2_a*w3_b)*(w4_c*u__u_0[0][_idx58])))+(((w3_a*w3_b)*(w4_c*u__u_0[0][_idx59]))+((w4_a*w3_b)*(w4_c*u__u_0[0][_idx60]))))+((((w1_a*w4_b)*(w4_c*u__u_0[0][_idx61]))+((w2_a*w4_b)*(w4_c*u__u_0[0][_idx62])))+(((w3_a*w4_b)*(w4_c*u__u_0[0][_idx63]))+((w4_a*w4_b)*(w4_c*u__u_0[0][_idx64]))))))));
					}
				}
			}
		}
	}
}

void initialize(double *  u_0_0, double *  u_0_1, double *  a_1_0, double *  b_2_0, double *  c_3_0, int x_max, int y_max, int z_max)
{
	int _idx0;
	int _idx1;
	int _idx10;
	int _idx11;
	int _idx12;
	int _idx13;
	int _idx14;
	int _idx15;
	int _idx16;
	int _idx17;
	int _idx18;
	int _idx19;
	int _idx2;
	int _idx20;
	int _idx21;
	int _idx22;
	int _idx23;
	int _idx24;
	int _idx25;
	int _idx26;
	int _idx27;
	int _idx28;
	int _idx29;
	int _idx3;
	int _idx30;
	int _idx31;
	int _idx32;
	int _idx33;
	int _idx34;
	int _idx35;
	int _idx36;
	int _idx37;
	int _idx38;
	int _idx39;
	int _idx4;
	int _idx40;
	int _idx41;
	int _idx42;
	int _idx43;
	int _idx44;
	int _idx45;
	int _idx46;
	int _idx47;
	int _idx48;
	int _idx49;
	int _idx5;
	int _idx50;
	int _idx51;
	int _idx52;
	int _idx53;
	int _idx54;
	int _idx55;
	int _idx56;
	int _idx57;
	int _idx58;
	int _idx59;
	int _idx6;
	int _idx60;
	int _idx61;
	int _idx62;
	int _idx63;
	int _idx64;
	int _idx7;
	int _idx8;
	int _idx9;
	int p_idx_x;
	int p_idx_y;
	int p_idx_z;
	int t;
	double *  __restrict__ const u__a_1[1] =  { a_1_0 } ;
	double *  __restrict__ const u__b_2[1] =  { b_2_0 } ;
	double *  __restrict__ const u__c_3[1] =  { c_3_0 } ;
	double *  __restrict__ const u__u_0[1] =  { u_0_0, u_0_1 } ;
	double w1_a;
	double w1_b;
	double w1_c;
	double w2_a;
	double w2_b;
	double w2_c;
	double w3_a;
	double w3_b;
	double w3_c;
	double w4_a;
	double w4_b;
	double w4_c;
	/*
	Initializations
	*/
	/*
	Implementation
	*/
	/*
	for t = 1..t_max by 1 parallel 1 <level 0> schedule  { ... }
	*/
	for (t=1; t<=t_max; t+=1)
	{
		/*
		for POINT p[t=t, s=(1, 1, 1)][0] of size [1, 1, 1] in u[t=t, s=(:, :, :)][0] parallel 1 <level 0> schedule default { ... }
		*/
		{
			/* Index bounds calculations for iterators in p[t=t, s=(1, 1, 1)][0] */
			for (p_idx_z=0; p_idx_z<z_max; p_idx_z+=1)
			{
				for (p_idx_y=0; p_idx_y<y_max; p_idx_y+=1)
				{
					for (p_idx_x=0; p_idx_x<x_max; p_idx_x+=1)
					{
						/* Index bounds calculations for iterators in p[t=t, s=(1, 1, 1)][0] */
						/*
						u[t=(t+1), s=p[t=?, s=?][0]][0]=stencil(u[t=t, s=p[t=?, s=?][0]][0])
						*/
						/* _idx0 = ((((p_idx_z*x_max)*y_max)+(p_idx_y*x_max))+p_idx_x) */
						_idx0=((((p_idx_z*x_max)*y_max)+(p_idx_y*x_max))+p_idx_x);
						u__a_1[0][_idx0]=0.2;
						w1_a=0.1;
						u__a_1[0][_idx0]=0.2;
						w2_a=0.1;
						u__a_1[0][_idx0]=0.2;
						w3_a=0.1;
						u__a_1[0][_idx0]=0.2;
						w4_a=0.1;
						u__b_2[0][_idx0]=0.30000000000000004;
						w1_b=0.1;
						u__b_2[0][_idx0]=0.30000000000000004;
						w2_b=0.1;
						u__b_2[0][_idx0]=0.30000000000000004;
						w3_b=0.1;
						u__b_2[0][_idx0]=0.30000000000000004;
						w4_b=0.1;
						u__c_3[0][_idx0]=0.4;
						w1_c=0.1;
						u__c_3[0][_idx0]=0.4;
						w2_c=0.1;
						u__c_3[0][_idx0]=0.4;
						w3_c=0.1;
						u__c_3[0][_idx0]=0.4;
						w4_c=0.1;
						/* _idx1 = (((((((p_idx_z*x_max)+(3*p_idx_z))*y_max)+(((3*p_idx_z)+p_idx_y)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x) */
						_idx1=((((_idx0+((3*p_idx_z)*y_max))+((3*p_idx_z)*x_max))+(9*p_idx_z))+(3*p_idx_y));
						u__u_0[0][_idx1]=0.1;
						/* _idx2 = ((((((((((p_idx_z+1)*x_max)+(3*p_idx_z))+3)*y_max)+((((3*p_idx_z)+p_idx_y)+3)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+9) */
						_idx2=(((_idx1+((x_max+3)*y_max))+(3*x_max))+9);
						u__u_0[0][_idx2]=0.1;
						/* _idx3 = ((((((((((p_idx_z+2)*x_max)+(3*p_idx_z))+6)*y_max)+((((3*p_idx_z)+p_idx_y)+6)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+18) */
						_idx3=(((_idx2+((x_max+3)*y_max))+(3*x_max))+9);
						u__u_0[0][_idx3]=0.1;
						/* _idx4 = ((((((((((p_idx_z+3)*x_max)+(3*p_idx_z))+9)*y_max)+((((3*p_idx_z)+p_idx_y)+9)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+27) */
						_idx4=(((_idx3+((x_max+3)*y_max))+(3*x_max))+9);
						u__u_0[0][_idx4]=0.1;
						/* _idx5 = ((((((((p_idx_z*x_max)+(3*p_idx_z))*y_max)+((((3*p_idx_z)+p_idx_y)+1)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+3) */
						_idx5=((_idx1+x_max)+3);
						u__u_0[0][_idx5]=0.1;
						/* _idx6 = ((((((((((p_idx_z+1)*x_max)+(3*p_idx_z))+3)*y_max)+((((3*p_idx_z)+p_idx_y)+4)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+12) */
						_idx6=((_idx2+x_max)+3);
						u__u_0[0][_idx6]=0.1;
						/* _idx7 = ((((((((((p_idx_z+2)*x_max)+(3*p_idx_z))+6)*y_max)+((((3*p_idx_z)+p_idx_y)+7)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+21) */
						_idx7=((_idx3+x_max)+3);
						u__u_0[0][_idx7]=0.1;
						/* _idx8 = ((((((((((p_idx_z+3)*x_max)+(3*p_idx_z))+9)*y_max)+((((3*p_idx_z)+p_idx_y)+10)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+30) */
						_idx8=((_idx4+x_max)+3);
						u__u_0[0][_idx8]=0.1;
						/* _idx9 = ((((((((p_idx_z*x_max)+(3*p_idx_z))*y_max)+((((3*p_idx_z)+p_idx_y)+2)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+6) */
						_idx9=((_idx5+x_max)+3);
						u__u_0[0][_idx9]=0.1;
						/* _idx10 = ((((((((((p_idx_z+1)*x_max)+(3*p_idx_z))+3)*y_max)+((((3*p_idx_z)+p_idx_y)+5)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+15) */
						_idx10=((_idx6+x_max)+3);
						u__u_0[0][_idx10]=0.1;
						/* _idx11 = ((((((((((p_idx_z+2)*x_max)+(3*p_idx_z))+6)*y_max)+((((3*p_idx_z)+p_idx_y)+8)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+24) */
						_idx11=((_idx7+x_max)+3);
						u__u_0[0][_idx11]=0.1;
						/* _idx12 = ((((((((((p_idx_z+3)*x_max)+(3*p_idx_z))+9)*y_max)+((((3*p_idx_z)+p_idx_y)+11)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+33) */
						_idx12=((_idx8+x_max)+3);
						u__u_0[0][_idx12]=0.1;
						/* _idx13 = ((((((((p_idx_z*x_max)+(3*p_idx_z))*y_max)+((((3*p_idx_z)+p_idx_y)+3)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+9) */
						_idx13=((_idx9+x_max)+3);
						u__u_0[0][_idx13]=0.1;
						/* _idx14 = ((((((((((p_idx_z+1)*x_max)+(3*p_idx_z))+3)*y_max)+((((3*p_idx_z)+p_idx_y)+6)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+18) */
						_idx14=((_idx10+x_max)+3);
						u__u_0[0][_idx14]=0.1;
						/* _idx15 = ((((((((((p_idx_z+2)*x_max)+(3*p_idx_z))+6)*y_max)+((((3*p_idx_z)+p_idx_y)+9)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+27) */
						_idx15=((_idx11+x_max)+3);
						u__u_0[0][_idx15]=0.1;
						/* _idx16 = ((((((((((p_idx_z+3)*x_max)+(3*p_idx_z))+9)*y_max)+((((3*p_idx_z)+p_idx_y)+12)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+36) */
						_idx16=((_idx12+x_max)+3);
						u__u_0[0][_idx16]=0.1;
						/* _idx17 = ((((((((p_idx_z*x_max)+(3*p_idx_z))*y_max)+(((3*p_idx_z)+p_idx_y)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+1) */
						_idx17=(_idx1+1);
						u__u_0[0][_idx17]=0.1;
						/* _idx18 = ((((((((((p_idx_z+1)*x_max)+(3*p_idx_z))+3)*y_max)+((((3*p_idx_z)+p_idx_y)+3)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+10) */
						_idx18=(_idx2+1);
						u__u_0[0][_idx18]=0.1;
						/* _idx19 = ((((((((((p_idx_z+2)*x_max)+(3*p_idx_z))+6)*y_max)+((((3*p_idx_z)+p_idx_y)+6)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+19) */
						_idx19=(_idx3+1);
						u__u_0[0][_idx19]=0.1;
						/* _idx20 = ((((((((((p_idx_z+3)*x_max)+(3*p_idx_z))+9)*y_max)+((((3*p_idx_z)+p_idx_y)+9)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+28) */
						_idx20=(_idx4+1);
						u__u_0[0][_idx20]=0.1;
						/* _idx21 = ((((((((p_idx_z*x_max)+(3*p_idx_z))*y_max)+((((3*p_idx_z)+p_idx_y)+1)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+4) */
						_idx21=((_idx17+x_max)+3);
						u__u_0[0][_idx21]=0.1;
						/* _idx22 = ((((((((((p_idx_z+1)*x_max)+(3*p_idx_z))+3)*y_max)+((((3*p_idx_z)+p_idx_y)+4)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+13) */
						_idx22=((_idx18+x_max)+3);
						u__u_0[0][_idx22]=0.1;
						/* _idx23 = ((((((((((p_idx_z+2)*x_max)+(3*p_idx_z))+6)*y_max)+((((3*p_idx_z)+p_idx_y)+7)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+22) */
						_idx23=((_idx19+x_max)+3);
						u__u_0[0][_idx23]=0.1;
						/* _idx24 = ((((((((((p_idx_z+3)*x_max)+(3*p_idx_z))+9)*y_max)+((((3*p_idx_z)+p_idx_y)+10)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+31) */
						_idx24=(_idx8+1);
						u__u_0[0][_idx24]=0.1;
						/* _idx25 = ((((((((p_idx_z*x_max)+(3*p_idx_z))*y_max)+((((3*p_idx_z)+p_idx_y)+2)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+7) */
						_idx25=((_idx5+x_max)+4);
						u__u_0[0][_idx25]=0.1;
						/* _idx26 = ((((((((((p_idx_z+1)*x_max)+(3*p_idx_z))+3)*y_max)+((((3*p_idx_z)+p_idx_y)+5)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+16) */
						_idx26=((_idx22+x_max)+3);
						u__u_0[0][_idx26]=0.1;
						/* _idx27 = ((((((((((p_idx_z+2)*x_max)+(3*p_idx_z))+6)*y_max)+((((3*p_idx_z)+p_idx_y)+8)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+25) */
						_idx27=((_idx23+x_max)+3);
						u__u_0[0][_idx27]=0.1;
						/* _idx28 = ((((((((((p_idx_z+3)*x_max)+(3*p_idx_z))+9)*y_max)+((((3*p_idx_z)+p_idx_y)+11)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+34) */
						_idx28=((_idx8+x_max)+4);
						u__u_0[0][_idx28]=0.1;
						/* _idx29 = ((((((((p_idx_z*x_max)+(3*p_idx_z))*y_max)+((((3*p_idx_z)+p_idx_y)+3)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+10) */
						_idx29=((_idx25+x_max)+3);
						u__u_0[0][_idx29]=0.1;
						/* _idx30 = ((((((((((p_idx_z+1)*x_max)+(3*p_idx_z))+3)*y_max)+((((3*p_idx_z)+p_idx_y)+6)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+19) */
						_idx30=((_idx10+x_max)+4);
						u__u_0[0][_idx30]=0.1;
						/* _idx31 = ((((((((((p_idx_z+2)*x_max)+(3*p_idx_z))+6)*y_max)+((((3*p_idx_z)+p_idx_y)+9)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+28) */
						_idx31=((_idx27+x_max)+3);
						u__u_0[0][_idx31]=0.1;
						/* _idx32 = ((((((((((p_idx_z+3)*x_max)+(3*p_idx_z))+9)*y_max)+((((3*p_idx_z)+p_idx_y)+12)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+37) */
						_idx32=((_idx28+x_max)+3);
						u__u_0[0][_idx32]=0.1;
						/* _idx33 = ((((((((p_idx_z*x_max)+(3*p_idx_z))*y_max)+(((3*p_idx_z)+p_idx_y)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+2) */
						_idx33=(_idx17+1);
						u__u_0[0][_idx33]=0.1;
						/* _idx34 = ((((((((((p_idx_z+1)*x_max)+(3*p_idx_z))+3)*y_max)+((((3*p_idx_z)+p_idx_y)+3)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+11) */
						_idx34=(_idx18+1);
						u__u_0[0][_idx34]=0.1;
						/* _idx35 = ((((((((((p_idx_z+2)*x_max)+(3*p_idx_z))+6)*y_max)+((((3*p_idx_z)+p_idx_y)+6)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+20) */
						_idx35=(_idx19+1);
						u__u_0[0][_idx35]=0.1;
						/* _idx36 = ((((((((((p_idx_z+3)*x_max)+(3*p_idx_z))+9)*y_max)+((((3*p_idx_z)+p_idx_y)+9)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+29) */
						_idx36=(_idx4+2);
						u__u_0[0][_idx36]=0.1;
						/* _idx37 = ((((((((p_idx_z*x_max)+(3*p_idx_z))*y_max)+((((3*p_idx_z)+p_idx_y)+1)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+5) */
						_idx37=((_idx17+x_max)+4);
						u__u_0[0][_idx37]=0.1;
						/* _idx38 = ((((((((((p_idx_z+1)*x_max)+(3*p_idx_z))+3)*y_max)+((((3*p_idx_z)+p_idx_y)+4)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+14) */
						_idx38=((_idx18+x_max)+4);
						u__u_0[0][_idx38]=0.1;
						/* _idx39 = ((((((((((p_idx_z+2)*x_max)+(3*p_idx_z))+6)*y_max)+((((3*p_idx_z)+p_idx_y)+7)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+23) */
						_idx39=((_idx19+x_max)+4);
						u__u_0[0][_idx39]=0.1;
						/* _idx40 = ((((((((((p_idx_z+3)*x_max)+(3*p_idx_z))+9)*y_max)+((((3*p_idx_z)+p_idx_y)+10)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+32) */
						_idx40=((_idx36+x_max)+3);
						u__u_0[0][_idx40]=0.1;
						/* _idx41 = ((((((((p_idx_z*x_max)+(3*p_idx_z))*y_max)+((((3*p_idx_z)+p_idx_y)+2)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+8) */
						_idx41=((_idx37+x_max)+3);
						u__u_0[0][_idx41]=0.1;
						/* _idx42 = ((((((((((p_idx_z+1)*x_max)+(3*p_idx_z))+3)*y_max)+((((3*p_idx_z)+p_idx_y)+5)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+17) */
						_idx42=((_idx22+x_max)+4);
						u__u_0[0][_idx42]=0.1;
						/* _idx43 = ((((((((((p_idx_z+2)*x_max)+(3*p_idx_z))+6)*y_max)+((((3*p_idx_z)+p_idx_y)+8)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+26) */
						_idx43=((_idx39+x_max)+3);
						u__u_0[0][_idx43]=0.1;
						/* _idx44 = ((((((((((p_idx_z+3)*x_max)+(3*p_idx_z))+9)*y_max)+((((3*p_idx_z)+p_idx_y)+11)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+35) */
						_idx44=(_idx28+1);
						u__u_0[0][_idx44]=0.1;
						/* _idx45 = ((((((((p_idx_z*x_max)+(3*p_idx_z))*y_max)+((((3*p_idx_z)+p_idx_y)+3)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+11) */
						_idx45=((_idx41+x_max)+3);
						u__u_0[0][_idx45]=0.1;
						/* _idx46 = ((((((((((p_idx_z+1)*x_max)+(3*p_idx_z))+3)*y_max)+((((3*p_idx_z)+p_idx_y)+6)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+20) */
						_idx46=((_idx42+x_max)+3);
						u__u_0[0][_idx46]=0.1;
						/* _idx47 = ((((((((((p_idx_z+2)*x_max)+(3*p_idx_z))+6)*y_max)+((((3*p_idx_z)+p_idx_y)+9)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+29) */
						_idx47=(_idx31+1);
						u__u_0[0][_idx47]=0.1;
						/* _idx48 = ((((((((((p_idx_z+3)*x_max)+(3*p_idx_z))+9)*y_max)+((((3*p_idx_z)+p_idx_y)+12)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+38) */
						_idx48=((_idx28+x_max)+4);
						u__u_0[0][_idx48]=0.1;
						/* _idx49 = ((((((((p_idx_z*x_max)+(3*p_idx_z))*y_max)+(((3*p_idx_z)+p_idx_y)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+3) */
						_idx49=(_idx17+2);
						u__u_0[0][_idx49]=0.1;
						/* _idx50 = ((((((((((p_idx_z+1)*x_max)+(3*p_idx_z))+3)*y_max)+((((3*p_idx_z)+p_idx_y)+3)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+12) */
						_idx50=(_idx18+2);
						u__u_0[0][_idx50]=0.1;
						/* _idx51 = ((((((((((p_idx_z+2)*x_max)+(3*p_idx_z))+6)*y_max)+((((3*p_idx_z)+p_idx_y)+6)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+21) */
						_idx51=(_idx19+2);
						u__u_0[0][_idx51]=0.1;
						/* _idx52 = ((((((((((p_idx_z+3)*x_max)+(3*p_idx_z))+9)*y_max)+((((3*p_idx_z)+p_idx_y)+9)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+30) */
						_idx52=(_idx36+1);
						u__u_0[0][_idx52]=0.1;
						/* _idx53 = ((((((((p_idx_z*x_max)+(3*p_idx_z))*y_max)+((((3*p_idx_z)+p_idx_y)+1)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+6) */
						_idx53=(_idx37+1);
						u__u_0[0][_idx53]=0.1;
						/* _idx54 = ((((((((((p_idx_z+1)*x_max)+(3*p_idx_z))+3)*y_max)+((((3*p_idx_z)+p_idx_y)+4)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+15) */
						_idx54=((_idx18+x_max)+5);
						u__u_0[0][_idx54]=0.1;
						/* _idx55 = ((((((((((p_idx_z+2)*x_max)+(3*p_idx_z))+6)*y_max)+((((3*p_idx_z)+p_idx_y)+7)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+24) */
						_idx55=((_idx19+x_max)+5);
						u__u_0[0][_idx55]=0.1;
						/* _idx56 = ((((((((((p_idx_z+3)*x_max)+(3*p_idx_z))+9)*y_max)+((((3*p_idx_z)+p_idx_y)+10)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+33) */
						_idx56=(_idx40+1);
						u__u_0[0][_idx56]=0.1;
						/* _idx57 = ((((((((p_idx_z*x_max)+(3*p_idx_z))*y_max)+((((3*p_idx_z)+p_idx_y)+2)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+9) */
						_idx57=(_idx41+1);
						u__u_0[0][_idx57]=0.1;
						/* _idx58 = ((((((((((p_idx_z+1)*x_max)+(3*p_idx_z))+3)*y_max)+((((3*p_idx_z)+p_idx_y)+5)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+18) */
						_idx58=(_idx42+1);
						u__u_0[0][_idx58]=0.1;
						/* _idx59 = ((((((((((p_idx_z+2)*x_max)+(3*p_idx_z))+6)*y_max)+((((3*p_idx_z)+p_idx_y)+8)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+27) */
						_idx59=((_idx39+x_max)+4);
						u__u_0[0][_idx59]=0.1;
						/* _idx60 = ((((((((((p_idx_z+3)*x_max)+(3*p_idx_z))+9)*y_max)+((((3*p_idx_z)+p_idx_y)+11)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+36) */
						_idx60=(_idx28+2);
						u__u_0[0][_idx60]=0.1;
						/* _idx61 = ((((((((p_idx_z*x_max)+(3*p_idx_z))*y_max)+((((3*p_idx_z)+p_idx_y)+3)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+12) */
						_idx61=((_idx41+x_max)+4);
						u__u_0[0][_idx61]=0.1;
						/* _idx62 = ((((((((((p_idx_z+1)*x_max)+(3*p_idx_z))+3)*y_max)+((((3*p_idx_z)+p_idx_y)+6)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+21) */
						_idx62=((_idx42+x_max)+4);
						u__u_0[0][_idx62]=0.1;
						/* _idx63 = ((((((((((p_idx_z+2)*x_max)+(3*p_idx_z))+6)*y_max)+((((3*p_idx_z)+p_idx_y)+9)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+30) */
						_idx63=(_idx47+1);
						u__u_0[0][_idx63]=0.1;
						/* _idx64 = ((((((((((p_idx_z+3)*x_max)+(3*p_idx_z))+9)*y_max)+((((3*p_idx_z)+p_idx_y)+12)*x_max))+(9*p_idx_z))+(3*p_idx_y))+p_idx_x)+39) */
						_idx64=((_idx28+x_max)+5);
						u__u_0[0][_idx64]=0.1;
						u__u_0[1][_idx22]=1.1;
					}
				}
			}
		}
	}
}

