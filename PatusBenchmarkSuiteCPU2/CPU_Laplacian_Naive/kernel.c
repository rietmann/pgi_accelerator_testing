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
  
#pragma acc region copyin(u_0_0[0:(x_max+2)*(y_max+2)*(z_max+2)-1]) copyout(u_0_1[0:(x_max+2)*(y_max+2)*(z_max+2)-1])
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
		    int temp,temp2;
		    /* Index bounds calculations for iterators in p[t=t, s=(1, 1, 1)][0] */

		    /* _idx0=((((((((((p_idx_z+1)*x_max)+(2*p_idx_z))+2)*y_max)+((((2*p_idx_z)+p_idx_y)+3)*x_max))+(4*p_idx_z))+(2*p_idx_y))+p_idx_x)+8); */
		    
		    /* _idx0=(((((((((())())))+(((()))))+())+())+)+8); */

		    /* _idx0 = p_idx_z+1; */
		    /* _idx0 = _idx0*x_max; */
		    /* temp = 2*p_idx_z; */
		    /* _idx0 += + temp; */
		    /* _idx0 += 2; */
		    /* _idx0 = _idx0*y_max; */
		    /* temp = 2*p_idx_z; */
		    /* temp += p_idx_y; */
		    /* temp += 3; */
		    /* temp = temp*x_max; */
		    /* _idx0 += temp; */
		    /* temp = 4*p_idx_z; */
		    /* _idx0 += temp; */

		    /* /\* PGCC-ACC BUG: if replace this line *\/ */
		    /* /\* _idx0 = _idx0 + 2*p_idx_y + p_idx_x + 8; *\/ */
		    /* /\* for the following 4 lines, we get incorrect results *\/ */
		    
		    /* temp = 2*p_idx_y; */
		    /* _idx0 += temp; */
		    /* _idx0 += p_idx_x; */
		    /* _idx0 += 8; */
		    
		    
		    
		    /* _idx1=_idx0-2; */
		    
		    /* /\* _idx2=((_idx1+x_max)+3); *\/ */
		    /* _idx2 = _idx1+x_max; */
		    /* _idx2 += 3;		     */
		    
		    /* /\* _idx3=((_idx1-x_max)-1); *\/ */

		    /* _idx3 = _idx1-x_max; */
		    /* _idx3 -= 1; */
		    
		    /* /\* _idx4=(((_idx3+((x_max+2)*y_max))+(3*x_max))+6); *\/ */
		    /* _idx4 = x_max+2; */
		    /* _idx4 = _idx4*y_max; */
		    /* _idx4 += _idx3; */
		    /* temp = 3*x_max; */
		    /* _idx4 += temp; */
		    /* _idx4 += 6; */
		    
		    /* /\* _idx5=(((_idx3+((( - x_max)-2)*y_max))-x_max)-2); *\/ */

		    /* _idx5 = -x_max-2; */
		    /* _idx5 *= y_max; */
		    /* _idx5 += _idx3; */
		    /* _idx5 += -x_max; */
		    /* _idx5 += -2; */
		    
		    /* /\* _idx6=((_idx3+x_max)+2);		     *\/ */

		    /* _idx6 = _idx3 + x_max; */
		    /* _idx6 += 2; */

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
		    /* _idx2 = ((((((((((p_idx_z+1)*x_max)+(2*p_idx_z))+2)*y_max)+((((2*p_idx_z)+p_idx_y)+4)*x_max))+(4*p_idx_z))+(2*p_idx_y))+p_idx_x)+9) */
		    {
		      int __tmp18 = (_idx1+x_max);
		      int __tmp19 = (__tmp18+3);
		      _idx2=__tmp19;
		    }
		    /* _idx3 = ((((((((((p_idx_z+1)*x_max)+(2*p_idx_z))+2)*y_max)+((((2*p_idx_z)+p_idx_y)+2)*x_max))+(4*p_idx_z))+(2*p_idx_y))+p_idx_x)+5) */
		    {
		      int __tmp20 = (_idx1-x_max);
		      int __tmp21 = (__tmp20-1);
		      _idx3=__tmp21;
		    }
		    /* _idx4 = ((((((((((p_idx_z+2)*x_max)+(2*p_idx_z))+4)*y_max)+((((2*p_idx_z)+p_idx_y)+5)*x_max))+(4*p_idx_z))+(2*p_idx_y))+p_idx_x)+11) */
		    {
		      int __tmp22 = (x_max+2);
		      int __tmp23 = (__tmp22*y_max);
		      int __tmp24 = (_idx3+__tmp23);
		      int __tmp25 = (3*x_max);
		      int __tmp26 = (__tmp24+__tmp25);
		      int __tmp27 = (__tmp26+6);
		      _idx4=__tmp27;
		    }
		    /* _idx5 = ((((((((p_idx_z*x_max)+(2*p_idx_z))*y_max)+((((2*p_idx_z)+p_idx_y)+1)*x_max))+(4*p_idx_z))+(2*p_idx_y))+p_idx_x)+3) */
		    {
		      int __tmp28 = (( - x_max)-2);
		      int __tmp29 = (__tmp28*y_max);
		      int __tmp30 = (_idx3+__tmp29);
		      int __tmp31 = (__tmp30-x_max);
		      int __tmp32 = (__tmp31-2);
		      _idx5=__tmp32;
		    }
		    /* _idx6 = ((((((((((p_idx_z+1)*x_max)+(2*p_idx_z))+2)*y_max)+((((2*p_idx_z)+p_idx_y)+3)*x_max))+(4*p_idx_z))+(2*p_idx_y))+p_idx_x)+7) */
		    {
		      int __tmp33 = (_idx3+x_max);
		      int __tmp34 = (__tmp33+2);
		      _idx6=__tmp34;
		    }

		    
		    u_0_1[_idx6]=((((u_0_0[_idx0]+(u_0_0[_idx1]+u_0_0[_idx2]))+(u_0_0[_idx3]+(u_0_0[_idx4]+u_0_0[_idx5])))*0.25)-u_0_0[_idx6]);

		    /* test _idx0 */
		    /* u_0_1[p_idx_x + p_idx_y*x_max + p_idx_z*x_max*y_max] = _idx3; */
		    
		    /* for cut and pasting 

		    _idx0=((((((((((p_idx_z+1)*x_max)+(2*p_idx_z))+2)*y_max)+((((2*p_idx_z)+p_idx_y)+3)*x_max))+(4*p_idx_z))+(2*p_idx_y))+p_idx_x)+8);
		    _idx1 = ( ((((((((((p_idx_z+1)*x_max)+(2*p_idx_z))+2)*y_max)+((((2*p_idx_z)+p_idx_y)+3)*x_max))+(4*p_idx_z))+(2*p_idx_y))+p_idx_x)+8)-2 );
		    _idx2 = ( ( ( ((((((((((p_idx_z+1)*x_max)+(2*p_idx_z))+2)*y_max)+((((2*p_idx_z)+p_idx_y)+3)*x_max))+(4*p_idx_z))+(2*p_idx_y))+p_idx_x)+8)-2 ) + x_max ) + 3);
		    _idx3 = ( (( ((((((((((p_idx_z+1)*x_max)+(2*p_idx_z))+2)*y_max)+((((2*p_idx_z)+p_idx_y)+3)*x_max))+(4*p_idx_z))+(2*p_idx_y))+p_idx_x)+8)-2 ) - x_max) - 1 );
		    _idx4 = ((( ( (( ((((((((((p_idx_z+1)*x_max)+(2*p_idx_z))+2)*y_max)+((((2*p_idx_z)+p_idx_y)+3)*x_max))+(4*p_idx_z))+(2*p_idx_y))+p_idx_x)+8)-2 ) - x_max) - 1 ) + ((x_max+2)*y_max))+(3*x_max))+6);
		    _idx5 = ((( ( (( ((((((((((p_idx_z+1)*x_max)+(2*p_idx_z))+2)*y_max)+((((2*p_idx_z)+p_idx_y)+3)*x_max))+(4*p_idx_z))+(2*p_idx_y))+p_idx_x)+8)-2 ) - x_max) - 1 ) + ((( - x_max)-2)*y_max))-x_max)-2);
		    _idx6 = (( ( (( ((((((((((p_idx_z+1)*x_max)+(2*p_idx_z))+2)*y_max)+((((2*p_idx_z)+p_idx_y)+3)*x_max))+(4*p_idx_z))+(2*p_idx_y))+p_idx_x)+8)-2 ) - x_max) - 1 ) + x_max)+2);

		    end for cut and pasting */

		    /* didn't work -- error: unknown shared array tile size
		    /* u_0_1[(( ( (( ((((((((((p_idx_z+1)*x_max)+(2*p_idx_z))+2)*y_max)+((((2*p_idx_z)+p_idx_y)+3)*x_max))+(4*p_idx_z))+(2*p_idx_y))+p_idx_x)+8)-2 ) - x_max) - 1 ) + x_max)+2)] = ((((u_0_0[((((((((((p_idx_z+1)*x_max)+(2*p_idx_z))+2)*y_max)+((((2*p_idx_z)+p_idx_y)+3)*x_max))+(4*p_idx_z))+(2*p_idx_y))+p_idx_x)+8)]+(u_0_0[( ((((((((((p_idx_z+1)*x_max)+(2*p_idx_z))+2)*y_max)+((((2*p_idx_z)+p_idx_y)+3)*x_max))+(4*p_idx_z))+(2*p_idx_y))+p_idx_x)+8)-2 )]+u_0_0[( ( ( ((((((((((p_idx_z+1)*x_max)+(2*p_idx_z))+2)*y_max)+((((2*p_idx_z)+p_idx_y)+3)*x_max))+(4*p_idx_z))+(2*p_idx_y))+p_idx_x)+8)-2 ) + x_max ) + 3)]))+(u_0_0[( (( ((((((((((p_idx_z+1)*x_max)+(2*p_idx_z))+2)*y_max)+((((2*p_idx_z)+p_idx_y)+3)*x_max))+(4*p_idx_z))+(2*p_idx_y))+p_idx_x)+8)-2 ) - x_max) - 1 )]+(u_0_0[((( ( (( ((((((((((p_idx_z+1)*x_max)+(2*p_idx_z))+2)*y_max)+((((2*p_idx_z)+p_idx_y)+3)*x_max))+(4*p_idx_z))+(2*p_idx_y))+p_idx_x)+8)-2 ) - x_max) - 1 ) + ((x_max+2)*y_max))+(3*x_max))+6)]+u_0_0[((( ( (( ((((((((((p_idx_z+1)*x_max)+(2*p_idx_z))+2)*y_max)+((((2*p_idx_z)+p_idx_y)+3)*x_max))+(4*p_idx_z))+(2*p_idx_y))+p_idx_x)+8)-2 ) - x_max) - 1 ) + ((( - x_max)-2)*y_max))-x_max)-2)])))*0.25)-u_0_0[(( ( (( ((((((((((p_idx_z+1)*x_max)+(2*p_idx_z))+2)*y_max)+((((2*p_idx_z)+p_idx_y)+3)*x_max))+(4*p_idx_z))+(2*p_idx_y))+p_idx_x)+8)-2 ) - x_max) - 1 ) + x_max)+2)]); */
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
	  int count=0;
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
						
						
						/* debugging */
						/* u_0_0[count] = p_idx_x + p_idx_y*x_max + p_idx_z*x_max*y_max; */
						/* count++; */
						
						/* u_0_1[p_idx_x + p_idx_y*x_max + p_idx_z*x_max*y_max] = _idx3; */
						/* u_0_1[p_idx_x + p_idx_y*x_max + p_idx_z*x_max] = p_idx_x + p_idx_y*x_max + p_idx_z*x_max; */
						/* end debugging */
						
						/* printf("u_0_0[_idx0]=%f,u_0_0[_idx2]=%f,u_0_0[_idx1]=%f\n",u_0_0[_idx0],u_0_0[_idx2],u_0_0[_idx1]); */
						
					}
				}
			}
		}
	}
	*u_0_1_out = u_0_1;
}
