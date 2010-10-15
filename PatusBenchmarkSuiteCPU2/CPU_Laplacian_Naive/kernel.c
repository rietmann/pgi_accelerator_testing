#include "omp.h"

void laplacian(float *  *  u_0_1_out, float *  u_0_0, float *  u_0_1, int x_max, int y_max, int z_max, int t_max)
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

	/*
	Initializations
	*/
	/*
	Implementation
	*/
	/*
	for t = 1..t_max by 1 parallel 1 <level 0> schedule  { ... }
	*/
#pragma acc region copyin(u_0_0[0:x_max*y_max*z_max]) copyout(u_0_1[0:x_max*y_max*z_max])
	{
#pragma acc for seq
	  for (t=1; t<=t_max; t+=1)
	    {
#pragma acc for independent
			for (p_idx_z=0; p_idx_z<z_max; p_idx_z+=1)
			{
#pragma acc for independent
			  for (p_idx_y=0; p_idx_y<y_max; p_idx_y+=1)
				{
#pragma acc for independent 
				  for (p_idx_x=0; p_idx_x<x_max; p_idx_x+=1)
					{
						/* Index bounds calculations for iterators in p[t=t, s=(1, 1, 1)][0] */
						/*
						u[t=(t+1), s=p[t=?, s=?][0]][0]=stencil(u[t=t, s=p[t=?, s=?][0]][0])
						*/
						/* _idx0 = ((((((((((p_idx_z+1)*x_max)+(2*p_idx_z))+2)*y_max)+((((2*p_idx_z)+p_idx_y)+3)*x_max))+(4*p_idx_z))+(2*p_idx_y))+p_idx_x)+8) */
						_idx0=((((((((((p_idx_z+1)*x_max)+(2*p_idx_z))+2)*y_max)+((((2*p_idx_z)+p_idx_y)+3)*x_max))+(4*p_idx_z))+(2*p_idx_y))+p_idx_x)+8);
						/* _idx1 = ((((((((((p_idx_z+1)*x_max)+(2*p_idx_z))+2)*y_max)+((((2*p_idx_z)+p_idx_y)+3)*x_max))+(4*p_idx_z))+(2*p_idx_y))+p_idx_x)+6) */
						_idx1=(_idx0-2);
						/* _idx2 = ((((((((((p_idx_z+1)*x_max)+(2*p_idx_z))+2)*y_max)+((((2*p_idx_z)+p_idx_y)+4)*x_max))+(4*p_idx_z))+(2*p_idx_y))+p_idx_x)+9) */
						_idx2=((_idx1+x_max)+3);
						/* _idx3 = ((((((((((p_idx_z+1)*x_max)+(2*p_idx_z))+2)*y_max)+((((2*p_idx_z)+p_idx_y)+2)*x_max))+(4*p_idx_z))+(2*p_idx_y))+p_idx_x)+5) */
						_idx3=((_idx1-x_max)-1);
						/* _idx4 = ((((((((((p_idx_z+2)*x_max)+(2*p_idx_z))+4)*y_max)+((((2*p_idx_z)+p_idx_y)+5)*x_max))+(4*p_idx_z))+(2*p_idx_y))+p_idx_x)+11) */
						_idx4=(((_idx3+((x_max+2)*y_max))+(3*x_max))+6);
						/* _idx5 = ((((((((p_idx_z*x_max)+(2*p_idx_z))*y_max)+((((2*p_idx_z)+p_idx_y)+1)*x_max))+(4*p_idx_z))+(2*p_idx_y))+p_idx_x)+3) */
						_idx5=(((_idx3+((( - x_max)-2)*y_max))-x_max)-2);
						/* _idx6 = ((((((((((p_idx_z+1)*x_max)+(2*p_idx_z))+2)*y_max)+((((2*p_idx_z)+p_idx_y)+3)*x_max))+(4*p_idx_z))+(2*p_idx_y))+p_idx_x)+7) */
						_idx6=((_idx3+x_max)+2);
						u_0_1[_idx6]=((((u_0_0[_idx0]+(u_0_0[_idx1]+u_0_0[_idx2]))+(u_0_0[_idx3]+(u_0_0[_idx4]+u_0_0[_idx5])))*0.25)-u_0_0[_idx6]);
					}
				}
			}
		}
	}
	*u_0_1_out = u_0_1;
}

void initialize(float *  u_0_0, float *  u_0_1, int x_max, int y_max, int z_max)
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
	float *  __restrict__ const u__u_0[2] =  { u_0_0, u_0_1 } ;
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
						/* _idx1 = ((((((((((p_idx_z+1)*x_max)+(2*p_idx_z))+2)*y_max)+((((2*p_idx_z)+p_idx_y)+2)*x_max))+(4*p_idx_z))+(2*p_idx_y))+p_idx_x)+5) */
						_idx1=((_idx0-x_max)-1);
						u__u_0[0][_idx1]=0.1;
						/* _idx2 = ((((((((p_idx_z*x_max)+(2*p_idx_z))*y_max)+((((2*p_idx_z)+p_idx_y)+1)*x_max))+(4*p_idx_z))+(2*p_idx_y))+p_idx_x)+3) */
						_idx2=(((_idx1+((( - x_max)-2)*y_max))-x_max)-2);
						u__u_0[0][_idx2]=0.1;
						/* _idx3 = ((((((((((p_idx_z+1)*x_max)+(2*p_idx_z))+2)*y_max)+((((2*p_idx_z)+p_idx_y)+3)*x_max))+(4*p_idx_z))+(2*p_idx_y))+p_idx_x)+7) */
						_idx3=((_idx1+x_max)+2);
						u__u_0[0][_idx3]=0.1;
						/* _idx4 = ((((((((((p_idx_z+2)*x_max)+(2*p_idx_z))+4)*y_max)+((((2*p_idx_z)+p_idx_y)+5)*x_max))+(4*p_idx_z))+(2*p_idx_y))+p_idx_x)+11) */
						_idx4=(((_idx1+((x_max+2)*y_max))+(3*x_max))+6);
						u__u_0[0][_idx4]=0.1;
						/* _idx5 = ((((((((((p_idx_z+1)*x_max)+(2*p_idx_z))+2)*y_max)+((((2*p_idx_z)+p_idx_y)+4)*x_max))+(4*p_idx_z))+(2*p_idx_y))+p_idx_x)+9) */
						_idx5=((_idx0+x_max)+3);
						u__u_0[0][_idx5]=0.1;
						/* _idx6 = ((((((((((p_idx_z+1)*x_max)+(2*p_idx_z))+2)*y_max)+((((2*p_idx_z)+p_idx_y)+3)*x_max))+(4*p_idx_z))+(2*p_idx_y))+p_idx_x)+8) */
						_idx6=((_idx1+x_max)+3);
						u__u_0[0][_idx6]=0.1;
						u__u_0[1][_idx3]=1.1;
					}
				}
			}
		}
	}
}

void laplacian_cpu(float *  *  u_0_1_out, float *  u_0_0, float *  u_0_1, int x_max, int y_max, int z_max, int t_max)
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
	
	  for (t=1; t<=t_max; t+=1)
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
						_idx0=((((((((((p_idx_z+1)*x_max)+(2*p_idx_z))+2)*y_max)+((((2*p_idx_z)+p_idx_y)+3)*x_max))+(4*p_idx_z))+(2*p_idx_y))+p_idx_x)+8);
						/* _idx1 = ((((((((((p_idx_z+1)*x_max)+(2*p_idx_z))+2)*y_max)+((((2*p_idx_z)+p_idx_y)+3)*x_max))+(4*p_idx_z))+(2*p_idx_y))+p_idx_x)+6) */
						_idx1=(_idx0-2);
						/* _idx2 = ((((((((((p_idx_z+1)*x_max)+(2*p_idx_z))+2)*y_max)+((((2*p_idx_z)+p_idx_y)+4)*x_max))+(4*p_idx_z))+(2*p_idx_y))+p_idx_x)+9) */
						_idx2=((_idx1+x_max)+3);
						/* _idx3 = ((((((((((p_idx_z+1)*x_max)+(2*p_idx_z))+2)*y_max)+((((2*p_idx_z)+p_idx_y)+2)*x_max))+(4*p_idx_z))+(2*p_idx_y))+p_idx_x)+5) */
						_idx3=((_idx1-x_max)-1);
						/* _idx4 = ((((((((((p_idx_z+2)*x_max)+(2*p_idx_z))+4)*y_max)+((((2*p_idx_z)+p_idx_y)+5)*x_max))+(4*p_idx_z))+(2*p_idx_y))+p_idx_x)+11) */
						_idx4=(((_idx3+((x_max+2)*y_max))+(3*x_max))+6);
						/* _idx5 = ((((((((p_idx_z*x_max)+(2*p_idx_z))*y_max)+((((2*p_idx_z)+p_idx_y)+1)*x_max))+(4*p_idx_z))+(2*p_idx_y))+p_idx_x)+3) */
						_idx5=(((_idx3+((( - x_max)-2)*y_max))-x_max)-2);
						/* _idx6 = ((((((((((p_idx_z+1)*x_max)+(2*p_idx_z))+2)*y_max)+((((2*p_idx_z)+p_idx_y)+3)*x_max))+(4*p_idx_z))+(2*p_idx_y))+p_idx_x)+7) */
						_idx6=((_idx3+x_max)+2);
						u_0_1[_idx6]=((((u_0_0[_idx0]+(u_0_0[_idx1]+u_0_0[_idx2]))+(u_0_0[_idx3]+(u_0_0[_idx4]+u_0_0[_idx5])))*0.25)-u_0_0[_idx6]);
					}
				}
			}
		}
	}
	*u_0_1_out = u_0_1;
}
