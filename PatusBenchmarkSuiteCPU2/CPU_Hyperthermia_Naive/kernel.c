#include "omp.h"

/*
(T[0][0][0][1][0]=((((T[0][0][0][0][0]*((c[0][0][0][0][1]*T[0][0][0][0][0])+c[0][0][0][0][2]))+c[0][0][0][0][3])+((c[0][0][0][0][4]*T[-1][0][0][0][0])+(c[0][0][0][0][5]*T[1][0][0][0][0])))+(((c[0][0][0][0][6]*T[0][-1][0][0][0])+(c[0][0][0][0][7]*T[0][1][0][0][0]))+((c[0][0][0][0][8]*T[0][0][-1][0][0])+(c[0][0][0][0][9]*T[0][0][1][0][0])))))
*/

void hyperthermia(float *  *  T_0_1_out, float *  T_0_0, float *  T_0_1, float *  c_1_0, float *  c_2_0, float *  c_3_0, float *  c_4_0, float *  c_5_0, float *  c_6_0, float *  c_7_0, float *  c_8_0, float *  c_9_0, int x_max, int y_max, int z_max, int t_max)
{
	int _idx0;
	int _idx1;
	int _idx2;
	int _idx3;
	int _idx4;
	int _idx5;
	int _idx6;
	int _idx7;
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
#pragma acc region copyin(T_0_0[0:(x_max+2)*(y_max+2)*(z_max+2)],c_1_0[0:x_max*y_max*z_max],c_2_0[0:x_max*y_max*z_max],c_3_0[0:x_max*y_max*z_max],c_4_0[0:x_max*y_max*z_max],c_5_0[0:x_max*y_max*z_max],c_6_0[0:x_max*y_max*z_max],c_7_0[0:x_max*y_max*z_max],c_8_0[0:x_max*y_max*z_max],c_9_0[0:x_max*y_max*z_max]) copyout(T_0_1[0:(x_max+2)*(y_max+2)*(z_max+2)])
	{
#pragma acc for seq
	  for (t=1; t<=t_max; t+=1)
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
						/* Index bounds calculations for iterators in p[t=t, s=(1, 1, 1)][0] */
						/*
						u[t=(t+1), s=p[t=?, s=?][0]][0]=stencil(u[t=t, s=p[t=?, s=?][0]][0])
						*/
						/* _idx0 = ((((((((((p_idx_z+1)*x_max)+(2*p_idx_z))+2)*y_max)+((((2*p_idx_z)+p_idx_y)+3)*x_max))+(4*p_idx_z))+(2*p_idx_y))+p_idx_x)+7) */
						/* _idx0=((((((((((p_idx_z+1)*x_max)+(2*p_idx_z))+2)*y_max)+((((2*p_idx_z)+p_idx_y)+3)*x_max))+(4*p_idx_z))+(2*p_idx_y))+p_idx_x)+7); */
						/* /\* _idx1 = ((((p_idx_z*x_max)*y_max)+(p_idx_y*x_max))+p_idx_x) *\/ */
						/* _idx1=((x_max*((p_idx_z*y_max)+p_idx_y))+p_idx_x); */
						/* /\* _idx2 = ((((((((((p_idx_z+1)*x_max)+(2*p_idx_z))+2)*y_max)+((((2*p_idx_z)+p_idx_y)+3)*x_max))+(4*p_idx_z))+(2*p_idx_y))+p_idx_x)+6) *\/ */
						/* _idx2=(_idx0-1); */
						/* /\* _idx3 = ((((((((((p_idx_z+1)*x_max)+(2*p_idx_z))+2)*y_max)+((((2*p_idx_z)+p_idx_y)+3)*x_max))+(4*p_idx_z))+(2*p_idx_y))+p_idx_x)+8) *\/ */
						/* _idx3=(_idx2+2); */
						/* /\* _idx4 = ((((((((((p_idx_z+1)*x_max)+(2*p_idx_z))+2)*y_max)+((((2*p_idx_z)+p_idx_y)+2)*x_max))+(4*p_idx_z))+(2*p_idx_y))+p_idx_x)+5) *\/ */
						/* _idx4=((_idx2-x_max)-1); */
						/* /\* _idx5 = ((((((((((p_idx_z+1)*x_max)+(2*p_idx_z))+2)*y_max)+((((2*p_idx_z)+p_idx_y)+4)*x_max))+(4*p_idx_z))+(2*p_idx_y))+p_idx_x)+9) *\/ */
						/* _idx5=((_idx2+x_max)+3); */
						/* /\* _idx6 = ((((((((p_idx_z*x_max)+(2*p_idx_z))*y_max)+((((2*p_idx_z)+p_idx_y)+1)*x_max))+(4*p_idx_z))+(2*p_idx_y))+p_idx_x)+3) *\/ */
						/* _idx6=(((_idx5+((( - x_max)-2)*y_max))-(3*x_max))-6); */
						/* /\* _idx7 = ((((((((((p_idx_z+2)*x_max)+(2*p_idx_z))+4)*y_max)+((((2*p_idx_z)+p_idx_y)+5)*x_max))+(4*p_idx_z))+(2*p_idx_y))+p_idx_x)+11) *\/ */
						/* _idx7=(((_idx5+((x_max+2)*y_max))+x_max)+2); */

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
						  int __tmp16 = (__tmp15+7);
						  _idx0=__tmp16;
						}
						/* _idx1 = ((((p_idx_z*x_max)*y_max)+(p_idx_y*x_max))+p_idx_x) */
						{
						  int __tmp17 = (p_idx_z*y_max);
						  int __tmp18 = (__tmp17+p_idx_y);
						  int __tmp19 = (x_max*__tmp18);
						  int __tmp20 = (__tmp19+p_idx_x);
						  _idx1=__tmp20;
						}
						/* _idx2 = ((((((((((p_idx_z+1)*x_max)+(2*p_idx_z))+2)*y_max)+((((2*p_idx_z)+p_idx_y)+3)*x_max))+(4*p_idx_z))+(2*p_idx_y))+p_idx_x)+6) */
						{
						  int __tmp21 = (_idx0-1);
						  _idx2=__tmp21;
						}
						/* _idx3 = ((((((((((p_idx_z+1)*x_max)+(2*p_idx_z))+2)*y_max)+((((2*p_idx_z)+p_idx_y)+3)*x_max))+(4*p_idx_z))+(2*p_idx_y))+p_idx_x)+8) */
						{
						  int __tmp22 = (_idx2+2);
						  _idx3=__tmp22;
						}
						/* _idx4 = ((((((((((p_idx_z+1)*x_max)+(2*p_idx_z))+2)*y_max)+((((2*p_idx_z)+p_idx_y)+2)*x_max))+(4*p_idx_z))+(2*p_idx_y))+p_idx_x)+5) */
						{
						  int __tmp23 = (_idx2-x_max);
						  int __tmp24 = (__tmp23-1);
						  _idx4=__tmp24;
						}
						/* _idx5 = ((((((((((p_idx_z+1)*x_max)+(2*p_idx_z))+2)*y_max)+((((2*p_idx_z)+p_idx_y)+4)*x_max))+(4*p_idx_z))+(2*p_idx_y))+p_idx_x)+9) */
						{
						  int __tmp25 = (_idx2+x_max);
						  int __tmp26 = (__tmp25+3);
						  _idx5=__tmp26;
						}
						/* _idx6 = ((((((((p_idx_z*x_max)+(2*p_idx_z))*y_max)+((((2*p_idx_z)+p_idx_y)+1)*x_max))+(4*p_idx_z))+(2*p_idx_y))+p_idx_x)+3) */
						{
						  int __tmp27 = (( - x_max)-2);
						  int __tmp28 = (__tmp27*y_max);
						  int __tmp29 = (_idx5+__tmp28);
						  int __tmp30 = (3*x_max);
						  int __tmp31 = (__tmp29-__tmp30);
						  int __tmp32 = (__tmp31-6);
						  _idx6=__tmp32;
						}
						/* _idx7 = ((((((((((p_idx_z+2)*x_max)+(2*p_idx_z))+4)*y_max)+((((2*p_idx_z)+p_idx_y)+5)*x_max))+(4*p_idx_z))+(2*p_idx_y))+p_idx_x)+11) */
						{
						  int __tmp33 = (x_max+2);
						  int __tmp34 = (__tmp33*y_max);
						  int __tmp35 = (_idx5+__tmp34);
						  int __tmp36 = (__tmp35+x_max);
						  int __tmp37 = (__tmp36+2);
						  _idx7=__tmp37;
						}


						T_0_1[_idx0]=((((T_0_0[_idx0]*((c_1_0[_idx1]*T_0_0[_idx0])+c_2_0[_idx1]))+c_3_0[_idx1])+((c_4_0[_idx1]*T_0_0[_idx2])+(c_5_0[_idx1]*T_0_0[_idx3])))+(((c_6_0[_idx1]*T_0_0[_idx4])+(c_7_0[_idx1]*T_0_0[_idx5]))+((c_8_0[_idx1]*T_0_0[_idx6])+(c_9_0[_idx1]*T_0_0[_idx7]))));
					}
				}
			}
		}
	}
	*T_0_1_out = T_0_1;
}

void initialize(float *  T_0_0, float *  T_0_1, float *  c_1_0, float *  c_2_0, float *  c_3_0, float *  c_4_0, float *  c_5_0, float *  c_6_0, float *  c_7_0, float *  c_8_0, float *  c_9_0, int x_max, int y_max, int z_max)
{
	int _idx0;
	int _idx1;
	int _idx2;
	int _idx3;
	int _idx4;
	int _idx5;
	int _idx6;
	int _idx7;
	int p_idx_x;
	int p_idx_y;
	int p_idx_z;
	int t;
	float *  __restrict__ const u__T_0[2] =  { T_0_0, T_0_1 } ;
	float *  __restrict__ const u__c_1[1] =  { c_1_0 } ;
	float *  __restrict__ const u__c_2[1] =  { c_2_0 } ;
	float *  __restrict__ const u__c_3[1] =  { c_3_0 } ;
	float *  __restrict__ const u__c_4[1] =  { c_4_0 } ;
	float *  __restrict__ const u__c_5[1] =  { c_5_0 } ;
	float *  __restrict__ const u__c_6[1] =  { c_6_0 } ;
	float *  __restrict__ const u__c_7[1] =  { c_7_0 } ;
	float *  __restrict__ const u__c_8[1] =  { c_8_0 } ;
	float *  __restrict__ const u__c_9[1] =  { c_9_0 } ;
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
						u__T_0[0][_idx0]=0.1;
						/* _idx1 = ((((((((((p_idx_z+1)*x_max)+(2*p_idx_z))+2)*y_max)+((((2*p_idx_z)+p_idx_y)+2)*x_max))+(4*p_idx_z))+(2*p_idx_y))+p_idx_x)+5) */
						_idx1=((_idx0-x_max)-1);
						u__T_0[0][_idx1]=0.1;
						/* _idx2 = ((((((((p_idx_z*x_max)+(2*p_idx_z))*y_max)+((((2*p_idx_z)+p_idx_y)+1)*x_max))+(4*p_idx_z))+(2*p_idx_y))+p_idx_x)+3) */
						_idx2=(((_idx1+((( - x_max)-2)*y_max))-x_max)-2);
						u__T_0[0][_idx2]=0.1;
						/* _idx3 = ((((((((((p_idx_z+1)*x_max)+(2*p_idx_z))+2)*y_max)+((((2*p_idx_z)+p_idx_y)+3)*x_max))+(4*p_idx_z))+(2*p_idx_y))+p_idx_x)+7) */
						_idx3=((_idx1+x_max)+2);
						u__T_0[0][_idx3]=0.1;
						/* _idx4 = ((((p_idx_z*x_max)*y_max)+(p_idx_y*x_max))+p_idx_x) */
						_idx4=((x_max*((p_idx_z*y_max)+p_idx_y))+p_idx_x);
						u__c_1[0][_idx4]=0.2;
						u__c_2[0][_idx4]=0.30000000000000004;
						u__c_3[0][_idx4]=0.4;
						u__c_4[0][_idx4]=0.5;
						u__c_5[0][_idx4]=0.6000000000000001;
						u__c_6[0][_idx4]=0.7000000000000001;
						u__c_7[0][_idx4]=0.8;
						u__c_8[0][_idx4]=0.9;
						u__c_9[0][_idx4]=1.0;
						/* _idx5 = ((((((((((p_idx_z+2)*x_max)+(2*p_idx_z))+4)*y_max)+((((2*p_idx_z)+p_idx_y)+5)*x_max))+(4*p_idx_z))+(2*p_idx_y))+p_idx_x)+11) */
						_idx5=(((_idx1+((x_max+2)*y_max))+(3*x_max))+6);
						u__T_0[0][_idx5]=0.1;
						/* _idx6 = ((((((((((p_idx_z+1)*x_max)+(2*p_idx_z))+2)*y_max)+((((2*p_idx_z)+p_idx_y)+4)*x_max))+(4*p_idx_z))+(2*p_idx_y))+p_idx_x)+9) */
						_idx6=((_idx0+x_max)+3);
						u__T_0[0][_idx6]=0.1;
						/* _idx7 = ((((((((((p_idx_z+1)*x_max)+(2*p_idx_z))+2)*y_max)+((((2*p_idx_z)+p_idx_y)+3)*x_max))+(4*p_idx_z))+(2*p_idx_y))+p_idx_x)+8) */
						_idx7=((_idx1+x_max)+3);
						u__T_0[0][_idx7]=0.1;
						u__T_0[1][_idx3]=1.1;
					}
				}
			}
		}
	}
}

void hyperthermia_cpu(float *  *  T_0_1_out, float *  T_0_0, float *  T_0_1, float *  c_1_0, float *  c_2_0, float *  c_3_0, float *  c_4_0, float *  c_5_0, float *  c_6_0, float *  c_7_0, float *  c_8_0, float *  c_9_0, int x_max, int y_max, int z_max,int t_max)
{
  int _idx0;
  int _idx1;
  int _idx2;
  int _idx3;
  int _idx4;
  int _idx5;
  int _idx6;
  int _idx7;
  int p_idx_x;
  int p_idx_y;
  int p_idx_z;
  int t;

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
		  /* _idx0 = ((((((((((p_idx_z+1)*x_max)+(2*p_idx_z))+2)*y_max)+((((2*p_idx_z)+p_idx_y)+3)*x_max))+(4*p_idx_z))+(2*p_idx_y))+p_idx_x)+7) */
		  _idx0=((((((((((p_idx_z+1)*x_max)+(2*p_idx_z))+2)*y_max)+((((2*p_idx_z)+p_idx_y)+3)*x_max))+(4*p_idx_z))+(2*p_idx_y))+p_idx_x)+7);
		  /* _idx1 = ((((p_idx_z*x_max)*y_max)+(p_idx_y*x_max))+p_idx_x) */
		  _idx1=((x_max*((p_idx_z*y_max)+p_idx_y))+p_idx_x);
		  /* _idx2 = ((((((((((p_idx_z+1)*x_max)+(2*p_idx_z))+2)*y_max)+((((2*p_idx_z)+p_idx_y)+3)*x_max))+(4*p_idx_z))+(2*p_idx_y))+p_idx_x)+6) */
		  _idx2=(_idx0-1);
		  /* _idx3 = ((((((((((p_idx_z+1)*x_max)+(2*p_idx_z))+2)*y_max)+((((2*p_idx_z)+p_idx_y)+3)*x_max))+(4*p_idx_z))+(2*p_idx_y))+p_idx_x)+8) */
		  _idx3=(_idx2+2);
		  /* _idx4 = ((((((((((p_idx_z+1)*x_max)+(2*p_idx_z))+2)*y_max)+((((2*p_idx_z)+p_idx_y)+2)*x_max))+(4*p_idx_z))+(2*p_idx_y))+p_idx_x)+5) */
		  _idx4=((_idx2-x_max)-1);
		  /* _idx5 = ((((((((((p_idx_z+1)*x_max)+(2*p_idx_z))+2)*y_max)+((((2*p_idx_z)+p_idx_y)+4)*x_max))+(4*p_idx_z))+(2*p_idx_y))+p_idx_x)+9) */
		  _idx5=((_idx2+x_max)+3);
		  /* _idx6 = ((((((((p_idx_z*x_max)+(2*p_idx_z))*y_max)+((((2*p_idx_z)+p_idx_y)+1)*x_max))+(4*p_idx_z))+(2*p_idx_y))+p_idx_x)+3) */
		  _idx6=(((_idx5+((( - x_max)-2)*y_max))-(3*x_max))-6);
		  /* _idx7 = ((((((((((p_idx_z+2)*x_max)+(2*p_idx_z))+4)*y_max)+((((2*p_idx_z)+p_idx_y)+5)*x_max))+(4*p_idx_z))+(2*p_idx_y))+p_idx_x)+11) */
		  _idx7=(((_idx5+((x_max+2)*y_max))+x_max)+2);
		  T_0_1[_idx0]=((((T_0_0[_idx0]*((c_1_0[_idx1]*T_0_0[_idx0])+c_2_0[_idx1]))+c_3_0[_idx1])+((c_4_0[_idx1]*T_0_0[_idx2])+(c_5_0[_idx1]*T_0_0[_idx3])))+(((c_6_0[_idx1]*T_0_0[_idx4])+(c_7_0[_idx1]*T_0_0[_idx5]))+((c_8_0[_idx1]*T_0_0[_idx6])+(c_9_0[_idx1]*T_0_0[_idx7]))));
		}
	    }
	}
    }
  *T_0_1_out = T_0_1;
}
