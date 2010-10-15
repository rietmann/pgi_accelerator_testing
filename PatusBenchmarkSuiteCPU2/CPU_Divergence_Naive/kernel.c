#include "omp.h"
/* #define t_max 1 */
/*
(u[0][0][0][0][0]=((alpha*(ux[1][0][0][0][1]-ux[-1][0][0][0][1]))+((beta*(uy[0][1][0][0][2]-uy[0][-1][0][0][2]))+(gamma*(uz[0][0][1][0][3]-uz[0][0][-1][0][3])))))

*/
void divergence(float *  *  u_0_0_out, float *  u_0_0, float *  ux_1_0, float *  uy_2_0, float *  uz_3_0, float alpha, float beta, float gamma, int x_max, int y_max, int z_max, int t_max)
{
  int _idx0;
  int _idx1;
  int _idx2;
  int _idx3;
  int _idx4;
  int _idx5;
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
	
#pragma acc region copyin(ux_1_0[0:(x_max*y_max*z_max)],uy_2_0[0:(x_max*y_max*z_max)],uz_3_0[0:(x_max*y_max*z_max)]) copyout(u_0_0[0:(x_max*y_max*z_max)])
  {
#pragma acc for seq 
    for (t=1; t<=t_max; t++)
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
		    _idx0=(((((((p_idx_z*x_max)+(2*p_idx_z))*y_max)+(p_idx_y*x_max))+(2*p_idx_y))+p_idx_x)+2);
		    _idx1=(_idx0-2);
		    _idx2=((x_max*(((p_idx_z*(y_max+2))+p_idx_y)+2))+p_idx_x);
		    _idx3=(_idx2-(2*x_max));
		    _idx4=((_idx3+((2*x_max)*y_max))-((2*p_idx_z)*x_max));
		    _idx5=(_idx3-((2*p_idx_z)*x_max));
		    u_0_0[_idx5]=((alpha*(ux_1_0[_idx0]-ux_1_0[_idx1]))+((beta*(uy_2_0[_idx2]-uy_2_0[_idx3]))+(gamma*(uz_3_0[_idx4]-uz_3_0[_idx5]))));
		  }
	      }
	  }
      }
  }
  *u_0_0_out = u_0_0;
}

void initialize(float *  u_0_0, float *  ux_1_0, float *  uy_2_0, float *  uz_3_0, float alpha, float beta, float gamma, int x_max, int y_max, int z_max)
{
	int _idx0;
	int _idx1;
	int _idx2;
	int _idx3;
	int _idx4;
	int _idx5;
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
						/* _idx0 = ((((((p_idx_z*x_max)+(2*p_idx_z))*y_max)+(p_idx_y*x_max))+(2*p_idx_y))+p_idx_x) */
						_idx0=((((((p_idx_z*x_max)+(2*p_idx_z))*y_max)+(p_idx_y*x_max))+(2*p_idx_y))+p_idx_x);
						u__ux_1[0][_idx0]=0.2;
						/* _idx1 = ((((p_idx_z*x_max)*y_max)+(((2*p_idx_z)+p_idx_y)*x_max))+p_idx_x) */
						_idx1=((x_max*((p_idx_z*(y_max+2))+p_idx_y))+p_idx_x);
						u__uy_2[0][_idx1]=0.30000000000000004;
						/* _idx2 = ((((p_idx_z*x_max)*y_max)+(p_idx_y*x_max))+p_idx_x) */
						_idx2=(_idx1-((2*p_idx_z)*x_max));
						u__uz_3[0][_idx2]=0.4;
						/* _idx3 = (((((p_idx_z+2)*x_max)*y_max)+(p_idx_y*x_max))+p_idx_x) */
						_idx3=(_idx2+((2*x_max)*y_max));
						u__uz_3[0][_idx3]=0.4;
						/* _idx4 = ((((p_idx_z*x_max)*y_max)+((((2*p_idx_z)+p_idx_y)+2)*x_max))+p_idx_x) */
						_idx4=(_idx1+(2*x_max));
						u__uy_2[0][_idx4]=0.30000000000000004;
						/* _idx5 = (((((((p_idx_z*x_max)+(2*p_idx_z))*y_max)+(p_idx_y*x_max))+(2*p_idx_y))+p_idx_x)+2) */
						_idx5=(_idx0+2);
						u__ux_1[0][_idx5]=0.2;
						u__u_0[0][_idx2]=0.1;
					}
				}
			}
		}
	}
}

void divergence_cpu(float *  *  u_0_0_out, float *  u_0_0, float *  ux_1_0, float *  uy_2_0, float *  uz_3_0, float alpha, float beta, float gamma, int x_max, int y_max, int z_max,int t_max)
{
	int _idx0;
	int _idx1;
	int _idx2;
	int _idx3;
	int _idx4;
	int _idx5;
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
	
	  for (t=1; t<=t_max; t++)
	    {
	      for (p_idx_z=0; p_idx_z<z_max; p_idx_z+=1)
		{
		  for (p_idx_y=0; p_idx_y<y_max; p_idx_y+=1)
		    {
		      for (p_idx_x=0; p_idx_x<x_max; p_idx_x+=1)
			{
			  _idx0=(((((((p_idx_z*x_max)+(2*p_idx_z))*y_max)+(p_idx_y*x_max))+(2*p_idx_y))+p_idx_x)+2);
			  _idx1=(_idx0-2);
			  _idx2=((x_max*(((p_idx_z*(y_max+2))+p_idx_y)+2))+p_idx_x);
			  _idx3=(_idx2-(2*x_max));
			  _idx4=((_idx3+((2*x_max)*y_max))-((2*p_idx_z)*x_max));
			  _idx5=(_idx3-((2*p_idx_z)*x_max));
			  u_0_0[_idx5]=((alpha*(ux_1_0[_idx0]-ux_1_0[_idx1]))+((beta*(uy_2_0[_idx2]-uy_2_0[_idx3]))+(gamma*(uz_3_0[_idx4]-uz_3_0[_idx5]))));
			}
		    }
		}
	    }
	
	*u_0_0_out = u_0_0;
}
