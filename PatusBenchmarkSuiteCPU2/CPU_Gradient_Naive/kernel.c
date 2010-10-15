#include "omp.h"

/*
(ux[0][0][0][0][1]=(alpha*(u[1][0][0][0][0]+u[-1][0][0][0][0])))
(uy[0][0][0][0][2]=(beta*(u[0][1][0][0][0]+u[0][-1][0][0][0])))
(uz[0][0][0][0][3]=(gamma*(u[0][0][1][0][0]+u[0][0][-1][0][0])))

*/
void gradient(float *  *  ux_1_0_out, float *  *  uy_2_0_out, float *  *  uz_3_0_out, float *  u_0_0, float *  ux_1_0, float *  uy_2_0, float *  uz_3_0, float alpha, float beta, float gamma, int x_max, int y_max, int z_max, int t_max)
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
	
#pragma acc region copyin(u_0_0[0:z_max*y_max*x_max]) copyout(ux_1_0[0:z_max*y_max*x_max],uy_2_0[0:z_max*y_max*x_max],uz_3_0[0:z_max*y_max*x_max])
	{
#pragma acc for seq 
	for (t=1; t<=t_max; t++)	
		{
			/* Index bounds calculations for iterators in p[t=t, s=(1, 1, 1)][0] */
#pragma acc for independent		  
		  for (p_idx_z=0; p_idx_z<z_max; p_idx_z+=1)
			{
#pragma acc for independent
			  for (p_idx_y=0; p_idx_y<y_max; p_idx_y+=1)
				{
#pragma acc for independent
				  for (p_idx_x=0; p_idx_x<x_max; p_idx_x+=1)
					{
					
						_idx0=((((((((((p_idx_z+1)*x_max)+(2*p_idx_z))+2)*y_max)+((((2*p_idx_z)+p_idx_y)+3)*x_max))+(4*p_idx_z))+(2*p_idx_y))+p_idx_x)+8);
						_idx1=(_idx0-2);
						_idx2=((x_max*((p_idx_z*y_max)+p_idx_y))+p_idx_x);
						ux_1_0[_idx2]=(alpha*(u_0_0[_idx0]+u_0_0[_idx1]));
						_idx3=((_idx1+x_max)+3);
						_idx4=((_idx1-x_max)-1);
						uy_2_0[_idx2]=(beta*(u_0_0[_idx3]+u_0_0[_idx4]));
						_idx5=(((_idx4+((x_max+2)*y_max))+(3*x_max))+6);
					
						_idx6=(((_idx4+((( - x_max)-2)*y_max))-x_max)-2);
						uz_3_0[_idx2]=(gamma*(u_0_0[_idx5]+u_0_0[_idx6]));
					}
				}
			}
		}
	}
	*ux_1_0_out = ux_1_0;
	*uy_2_0_out = uy_2_0;
	*uz_3_0_out = uz_3_0;
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

void gradient_cpu(float *  *  ux_1_0_out, float *  *  uy_2_0_out, float *  *  uz_3_0_out, float *  u_0_0, float *  ux_1_0, float *  uy_2_0, float *  uz_3_0, float alpha, float beta, float gamma, int x_max, int y_max, int z_max, int t_max)
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
	
	for (t=1; t<=t_max; t++)	
		{
			/* Index bounds calculations for iterators in p[t=t, s=(1, 1, 1)][0] */
		  for (p_idx_z=0; p_idx_z<z_max; p_idx_z+=1)
			{
			  for (p_idx_y=0; p_idx_y<y_max; p_idx_y+=1)
				{
				  for (p_idx_x=0; p_idx_x<x_max; p_idx_x+=1)
					{
					
						_idx0=((((((((((p_idx_z+1)*x_max)+(2*p_idx_z))+2)*y_max)+((((2*p_idx_z)+p_idx_y)+3)*x_max))+(4*p_idx_z))+(2*p_idx_y))+p_idx_x)+8);
						_idx1=(_idx0-2);
						_idx2=((x_max*((p_idx_z*y_max)+p_idx_y))+p_idx_x);
						ux_1_0[_idx2]=(alpha*(u_0_0[_idx0]+u_0_0[_idx1]));
						_idx3=((_idx1+x_max)+3);
						_idx4=((_idx1-x_max)-1);
						uy_2_0[_idx2]=(beta*(u_0_0[_idx3]+u_0_0[_idx4]));
						_idx5=(((_idx4+((x_max+2)*y_max))+(3*x_max))+6);
					
						_idx6=(((_idx4+((( - x_max)-2)*y_max))-x_max)-2);
						uz_3_0[_idx2]=(gamma*(u_0_0[_idx5]+u_0_0[_idx6]));
					}
				}
			}
		}

	*ux_1_0_out = ux_1_0;
	*uy_2_0_out = uy_2_0;
	*uz_3_0_out = uz_3_0;
}
