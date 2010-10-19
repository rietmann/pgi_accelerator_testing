#include "omp.h"
#define t_max 1
/*
(ux[0][0][0][0][1]=(alpha*(u[1][0][0][0][0]+u[-1][0][0][0][0])))
(uy[0][0][0][0][2]=(beta*(u[0][1][0][0][0]+u[0][-1][0][0][0])))
(uz[0][0][0][0][3]=(gamma*(u[0][0][1][0][0]+u[0][0][-1][0][0])))

*/
void gradient(float *  *  ux_1_0_out, float *  *  uy_2_0_out, float *  *  uz_3_0_out, float *  u_0_0, float *  ux_1_0, float *  uy_2_0, float *  uz_3_0, float alpha, float beta, float gamma, int x_max, int y_max, int z_max)
{
	int _idx0;
	int _idx1;
	int _idx2;
	int _idx3;
	int _idx4;
	int _idx5;
	int _idx6;
	int p_idx_x;
	int p_idx_y;
	int p_idx_z;
	int t;
	float *  __restrict__ const u__u_0[1] =  { u_0_0 } ;
	float *  __restrict__ const u__ux_1[1] =  { ux_1_0 } ;
	float *  __restrict__ const u__uy_2[1] =  { uy_2_0 } ;
	float *  __restrict__ const u__uz_3[1] =  { uz_3_0 } ;
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
						/* _idx0 = ((((((((((p_idx_z+1)*x_max)+(2*p_idx_z))+2)*y_max)+((((2*p_idx_z)+p_idx_y)+3)*x_max))+(4*p_idx_z))+(2*p_idx_y))+p_idx_x)+8) */
						{
							int __tmp0 = (p_idx_z+1);
							int __tmp1 = (__tmp0*x_max);
							int __tmp2 = (2*p_idx_z);
							int __tmp3 = (__tmp1+__tmp2);
							int __tmp4 = (__tmp3+2);
							int __tmp5 = (__tmp4*y_max);
							int __tmp6 = (2*p_idx_z);
							int __tmp7 = (__tmp6+p_idx_y);
							int __tmp8 = (__tmp7+3);
							int __tmp9 = (__tmp8*x_max);
							int __tmp10 = (__tmp5+__tmp9);
							int __tmp11 = (4*p_idx_z);
							int __tmp12 = (__tmp10+__tmp11);
							int __tmp13 = (2*p_idx_y);
							int __tmp14 = (__tmp12+__tmp13);
							int __tmp15 = (__tmp14+p_idx_x);
							int __tmp16 = (__tmp15+8);
							_idx0=__tmp16;
						}
						/* _idx1 = ((((((((((p_idx_z+1)*x_max)+(2*p_idx_z))+2)*y_max)+((((2*p_idx_z)+p_idx_y)+3)*x_max))+(4*p_idx_z))+(2*p_idx_y))+p_idx_x)+6) */
						{
							int __tmp17 = (_idx0-2);
							_idx1=__tmp17;
						}
						/* _idx2 = ((((p_idx_z*x_max)*y_max)+(p_idx_y*x_max))+p_idx_x) */
						{
							int __tmp18 = (p_idx_z*y_max);
							int __tmp19 = (__tmp18+p_idx_y);
							int __tmp20 = (x_max*__tmp19);
							int __tmp21 = (__tmp20+p_idx_x);
							_idx2=__tmp21;
						}
						u__ux_1[0][_idx2]=(alpha*(u__u_0[0][_idx0]+u__u_0[0][_idx1]));
						/* _idx3 = ((((((((((p_idx_z+1)*x_max)+(2*p_idx_z))+2)*y_max)+((((2*p_idx_z)+p_idx_y)+4)*x_max))+(4*p_idx_z))+(2*p_idx_y))+p_idx_x)+9) */
						{
							int __tmp22 = (_idx1+x_max);
							int __tmp23 = (__tmp22+3);
							_idx3=__tmp23;
						}
						/* _idx4 = ((((((((((p_idx_z+1)*x_max)+(2*p_idx_z))+2)*y_max)+((((2*p_idx_z)+p_idx_y)+2)*x_max))+(4*p_idx_z))+(2*p_idx_y))+p_idx_x)+5) */
						{
							int __tmp24 = (_idx1-x_max);
							int __tmp25 = (__tmp24-1);
							_idx4=__tmp25;
						}
						u__uy_2[0][_idx2]=(beta*(u__u_0[0][_idx3]+u__u_0[0][_idx4]));
						/* _idx5 = ((((((((((p_idx_z+2)*x_max)+(2*p_idx_z))+4)*y_max)+((((2*p_idx_z)+p_idx_y)+5)*x_max))+(4*p_idx_z))+(2*p_idx_y))+p_idx_x)+11) */
						{
							int __tmp26 = (x_max+2);
							int __tmp27 = (__tmp26*y_max);
							int __tmp28 = (_idx4+__tmp27);
							int __tmp29 = (3*x_max);
							int __tmp30 = (__tmp28+__tmp29);
							int __tmp31 = (__tmp30+6);
							_idx5=__tmp31;
						}
						/* _idx6 = ((((((((p_idx_z*x_max)+(2*p_idx_z))*y_max)+((((2*p_idx_z)+p_idx_y)+1)*x_max))+(4*p_idx_z))+(2*p_idx_y))+p_idx_x)+3) */
						{
							int __tmp32 = (( - x_max)-2);
							int __tmp33 = (__tmp32*y_max);
							int __tmp34 = (_idx4+__tmp33);
							int __tmp35 = (__tmp34-x_max);
							int __tmp36 = (__tmp35-2);
							_idx6=__tmp36;
						}
						u__uz_3[0][_idx2]=(gamma*(u__u_0[0][_idx5]+u__u_0[0][_idx6]));
					}
				}
			}
		}
	}
}

void initialize(float *  u_0_0, float *  ux_1_0, float *  uy_2_0, float *  uz_3_0, float alpha, float beta, float gamma, int x_max, int y_max, int z_max)
{
	int _idx0;
	int _idx1;
	int _idx2;
	int _idx3;
	int _idx4;
	int _idx5;
	int _idx6;
	int p_idx_x;
	int p_idx_y;
	int p_idx_z;
	int t;
	float *  __restrict__ const u__u_0[1] =  { u_0_0 } ;
	float *  __restrict__ const u__ux_1[1] =  { ux_1_0 } ;
	float *  __restrict__ const u__uy_2[1] =  { uy_2_0 } ;
	float *  __restrict__ const u__uz_3[1] =  { uz_3_0 } ;
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
						/* _idx0 = ((((((((((p_idx_z+1)*x_max)+(2*p_idx_z))+2)*y_max)+((((2*p_idx_z)+p_idx_y)+3)*x_max))+(4*p_idx_z))+(2*p_idx_y))+p_idx_x)+6) */
						_idx0=((((((((((p_idx_z+1)*x_max)+(2*p_idx_z))+2)*y_max)+((((2*p_idx_z)+p_idx_y)+3)*x_max))+(4*p_idx_z))+(2*p_idx_y))+p_idx_x)+6);
						u__u_0[0][_idx0]=0.1;
						/* _idx1 = ((((((((((p_idx_z+1)*x_max)+(2*p_idx_z))+2)*y_max)+((((2*p_idx_z)+p_idx_y)+3)*x_max))+(4*p_idx_z))+(2*p_idx_y))+p_idx_x)+8) */
						_idx1=(_idx0+2);
						u__u_0[0][_idx1]=0.1;
						/* _idx2 = ((((p_idx_z*x_max)*y_max)+(p_idx_y*x_max))+p_idx_x) */
						_idx2=((x_max*((p_idx_z*y_max)+p_idx_y))+p_idx_x);
						u__ux_1[0][_idx2]=0.2;
						/* _idx3 = ((((((((((p_idx_z+1)*x_max)+(2*p_idx_z))+2)*y_max)+((((2*p_idx_z)+p_idx_y)+2)*x_max))+(4*p_idx_z))+(2*p_idx_y))+p_idx_x)+5) */
						_idx3=((_idx0-x_max)-1);
						u__u_0[0][_idx3]=0.1;
						/* _idx4 = ((((((((((p_idx_z+1)*x_max)+(2*p_idx_z))+2)*y_max)+((((2*p_idx_z)+p_idx_y)+4)*x_max))+(4*p_idx_z))+(2*p_idx_y))+p_idx_x)+9) */
						_idx4=((_idx0+x_max)+3);
						u__u_0[0][_idx4]=0.1;
						u__uy_2[0][_idx2]=0.30000000000000004;
						/* _idx5 = ((((((((p_idx_z*x_max)+(2*p_idx_z))*y_max)+((((2*p_idx_z)+p_idx_y)+1)*x_max))+(4*p_idx_z))+(2*p_idx_y))+p_idx_x)+3) */
						_idx5=(((_idx4+((( - x_max)-2)*y_max))-(3*x_max))-6);
						u__u_0[0][_idx5]=0.1;
						/* _idx6 = ((((((((((p_idx_z+2)*x_max)+(2*p_idx_z))+4)*y_max)+((((2*p_idx_z)+p_idx_y)+5)*x_max))+(4*p_idx_z))+(2*p_idx_y))+p_idx_x)+11) */
						_idx6=(((_idx4+((x_max+2)*y_max))+x_max)+2);
						u__u_0[0][_idx6]=0.1;
						u__uz_3[0][_idx2]=0.4;
					}
				}
			}
		}
	}
}

