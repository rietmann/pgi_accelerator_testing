#include "omp.h"


void upstream_5_3d(double *  *  u_0_1_out, double *  u_0_0, double *  u_0_1, double a, int x_max, int y_max, int z_max, int t_max)
{
	int _idx0;
	int _idx1;
	int _idx10;
	int _idx11;
	int _idx12;
	int _idx13;
	int _idx14;
	int _idx15;
	int _idx2;
	int _idx3;
	int _idx4;
	int _idx5;
	int _idx6;
	int _idx7;
	int _idx8;
	int _idx9;
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
#pragma acc region copyin(u_0_0[0:(x_max+5)*(y_max+5)*(z_max+5)]) copy(u_0_1[0:(x_max+5)*(y_max+5)*(z_max+5)])
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
						/* _idx0 = ((((((((((p_idx_z+3)*x_max)+(5*p_idx_z))+15)*y_max)+((((5*p_idx_z)+p_idx_y)+18)*x_max))+(25*p_idx_z))+(5*p_idx_y))+p_idx_x)+90) */
						/* _idx0=((((((((((p_idx_z+3)*x_max)+(5*p_idx_z))+15)*y_max)+((((5*p_idx_z)+p_idx_y)+18)*x_max))+(25*p_idx_z))+(5*p_idx_y))+p_idx_x)+90); */
						/* /\* _idx1 = ((((((((((p_idx_z+3)*x_max)+(5*p_idx_z))+15)*y_max)+((((5*p_idx_z)+p_idx_y)+15)*x_max))+(25*p_idx_z))+(5*p_idx_y))+p_idx_x)+78) *\/ */
						/* _idx1=((_idx0-(3*x_max))-12); */
						/* /\* _idx2 = ((((((((p_idx_z*x_max)+(5*p_idx_z))*y_max)+((((5*p_idx_z)+p_idx_y)+3)*x_max))+(25*p_idx_z))+(5*p_idx_y))+p_idx_x)+18) *\/ */
						/* _idx2=(((_idx1+(((-3*x_max)-15)*y_max))-(12*x_max))-60); */
						/* /\* _idx3 = ((((((((((p_idx_z+3)*x_max)+(5*p_idx_z))+15)*y_max)+((((5*p_idx_z)+p_idx_y)+18)*x_max))+(25*p_idx_z))+(5*p_idx_y))+p_idx_x)+91) *\/ */
						/* _idx3=(_idx0+1); */
						/* /\* _idx4 = ((((((((((p_idx_z+3)*x_max)+(5*p_idx_z))+15)*y_max)+((((5*p_idx_z)+p_idx_y)+16)*x_max))+(25*p_idx_z))+(5*p_idx_y))+p_idx_x)+83) *\/ */
						/* _idx4=((_idx1+x_max)+5); */
						/* /\* _idx5 = ((((((((((p_idx_z+1)*x_max)+(5*p_idx_z))+5)*y_max)+((((5*p_idx_z)+p_idx_y)+8)*x_max))+(25*p_idx_z))+(5*p_idx_y))+p_idx_x)+43) *\/ */
						/* _idx5=(((_idx2+((x_max+5)*y_max))+(5*x_max))+25); */
						/* /\* _idx6 = ((((((((((p_idx_z+3)*x_max)+(5*p_idx_z))+15)*y_max)+((((5*p_idx_z)+p_idx_y)+18)*x_max))+(25*p_idx_z))+(5*p_idx_y))+p_idx_x)+92) *\/ */
						/* _idx6=(_idx3+1); */
						/* /\* _idx7 = ((((((((((p_idx_z+3)*x_max)+(5*p_idx_z))+15)*y_max)+((((5*p_idx_z)+p_idx_y)+17)*x_max))+(25*p_idx_z))+(5*p_idx_y))+p_idx_x)+88) *\/ */
						/* _idx7=((_idx4+x_max)+5); */
						/* /\* _idx8 = ((((((((((p_idx_z+2)*x_max)+(5*p_idx_z))+10)*y_max)+((((5*p_idx_z)+p_idx_y)+13)*x_max))+(25*p_idx_z))+(5*p_idx_y))+p_idx_x)+68) *\/ */
						/* _idx8=(((_idx5+((x_max+5)*y_max))+(5*x_max))+25); */
						/* /\* _idx9 = ((((((((((p_idx_z+3)*x_max)+(5*p_idx_z))+15)*y_max)+((((5*p_idx_z)+p_idx_y)+18)*x_max))+(25*p_idx_z))+(5*p_idx_y))+p_idx_x)+93) *\/ */
						/* _idx9=((_idx7+x_max)+5); */
						/* /\* _idx10 = ((((((((((p_idx_z+3)*x_max)+(5*p_idx_z))+15)*y_max)+((((5*p_idx_z)+p_idx_y)+18)*x_max))+(25*p_idx_z))+(5*p_idx_y))+p_idx_x)+94) *\/ */
						/* _idx10=((_idx7+x_max)+6); */
						/* /\* _idx11 = ((((((((((p_idx_z+3)*x_max)+(5*p_idx_z))+15)*y_max)+((((5*p_idx_z)+p_idx_y)+19)*x_max))+(25*p_idx_z))+(5*p_idx_y))+p_idx_x)+98) *\/ */
						/* _idx11=((_idx10+x_max)+4); */
						/* /\* _idx12 = ((((((((((p_idx_z+4)*x_max)+(5*p_idx_z))+20)*y_max)+((((5*p_idx_z)+p_idx_y)+23)*x_max))+(25*p_idx_z))+(5*p_idx_y))+p_idx_x)+118) *\/ */
						/* _idx12=(((_idx1+((x_max+5)*y_max))+(8*x_max))+40); */
						/* /\* _idx13 = ((((((((((p_idx_z+3)*x_max)+(5*p_idx_z))+15)*y_max)+((((5*p_idx_z)+p_idx_y)+18)*x_max))+(25*p_idx_z))+(5*p_idx_y))+p_idx_x)+95) *\/ */
						/* _idx13=((_idx7+x_max)+7); */
						/* /\* _idx14 = ((((((((((p_idx_z+3)*x_max)+(5*p_idx_z))+15)*y_max)+((((5*p_idx_z)+p_idx_y)+20)*x_max))+(25*p_idx_z))+(5*p_idx_y))+p_idx_x)+103) *\/ */
						/* _idx14=((_idx11+x_max)+5); */
						/* /\* _idx15 = ((((((((((p_idx_z+5)*x_max)+(5*p_idx_z))+25)*y_max)+((((5*p_idx_z)+p_idx_y)+28)*x_max))+(25*p_idx_z))+(5*p_idx_y))+p_idx_x)+143) *\/ */
						/* _idx15=(((_idx12+((x_max+5)*y_max))+(5*x_max))+25); */
					  
						{
						  int __tmp0 = (p_idx_z+3);
						  int __tmp1 = (__tmp0*x_max);
						  int __tmp2 = (5*p_idx_z);
						  int __tmp3 = (__tmp1+__tmp2);
						  int __tmp4 = (__tmp3+15);
						  int __tmp5 = (__tmp4*y_max);
						  int __tmp6 = (5*p_idx_z);
						  int __tmp7 = (__tmp6+p_idx_y);
						  int __tmp8 = (__tmp7+18);
						  int __tmp9 = (__tmp8*x_max);
						  int __tmp10 = (__tmp5+__tmp9);
						  int __tmp11 = (25*p_idx_z);
						  int __tmp12 = (__tmp10+__tmp11);
						  int __tmp13 = (5*p_idx_y);
						  int __tmp14 = (__tmp12+__tmp13);
						  int __tmp15 = (__tmp14+p_idx_x);
						  int __tmp16 = (__tmp15+90);
						  _idx0=__tmp16;
						}
						/* _idx1 = ((((((((((p_idx_z+3)*x_max)+(5*p_idx_z))+15)*y_max)+((((5*p_idx_z)+p_idx_y)+15)*x_max))+(25*p_idx_z))+(5*p_idx_y))+p_idx_x)+78) */
						{
						  int __tmp17 = (3*x_max);
						  int __tmp18 = (_idx0-__tmp17);
						  int __tmp19 = (__tmp18-12);
						  _idx1=__tmp19;
						}
						/* _idx2 = ((((((((p_idx_z*x_max)+(5*p_idx_z))*y_max)+((((5*p_idx_z)+p_idx_y)+3)*x_max))+(25*p_idx_z))+(5*p_idx_y))+p_idx_x)+18) */
						{
						  int __tmp20 = (-3*x_max);
						  int __tmp21 = (__tmp20-15);
						  int __tmp22 = (__tmp21*y_max);
						  int __tmp23 = (_idx1+__tmp22);
						  int __tmp24 = (12*x_max);
						  int __tmp25 = (__tmp23-__tmp24);
						  int __tmp26 = (__tmp25-60);
						  _idx2=__tmp26;
						}
						/* _idx3 = ((((((((((p_idx_z+3)*x_max)+(5*p_idx_z))+15)*y_max)+((((5*p_idx_z)+p_idx_y)+18)*x_max))+(25*p_idx_z))+(5*p_idx_y))+p_idx_x)+91) */
						{
						  int __tmp27 = (_idx0+1);
						  _idx3=__tmp27;
						}
						/* _idx4 = ((((((((((p_idx_z+3)*x_max)+(5*p_idx_z))+15)*y_max)+((((5*p_idx_z)+p_idx_y)+16)*x_max))+(25*p_idx_z))+(5*p_idx_y))+p_idx_x)+83) */
						{
						  int __tmp28 = (_idx1+x_max);
						  int __tmp29 = (__tmp28+5);
						  _idx4=__tmp29;
						}
						/* _idx5 = ((((((((((p_idx_z+1)*x_max)+(5*p_idx_z))+5)*y_max)+((((5*p_idx_z)+p_idx_y)+8)*x_max))+(25*p_idx_z))+(5*p_idx_y))+p_idx_x)+43) */
						{
						  int __tmp30 = (x_max+5);
						  int __tmp31 = (__tmp30*y_max);
						  int __tmp32 = (_idx2+__tmp31);
						  int __tmp33 = (5*x_max);
						  int __tmp34 = (__tmp32+__tmp33);
						  int __tmp35 = (__tmp34+25);
						  _idx5=__tmp35;
						}
						/* _idx6 = ((((((((((p_idx_z+3)*x_max)+(5*p_idx_z))+15)*y_max)+((((5*p_idx_z)+p_idx_y)+18)*x_max))+(25*p_idx_z))+(5*p_idx_y))+p_idx_x)+92) */
						{
						  int __tmp36 = (_idx3+1);
						  _idx6=__tmp36;
						}
						/* _idx7 = ((((((((((p_idx_z+3)*x_max)+(5*p_idx_z))+15)*y_max)+((((5*p_idx_z)+p_idx_y)+17)*x_max))+(25*p_idx_z))+(5*p_idx_y))+p_idx_x)+88) */
						{
						  int __tmp37 = (_idx4+x_max);
						  int __tmp38 = (__tmp37+5);
						  _idx7=__tmp38;
						}
						/* _idx8 = ((((((((((p_idx_z+2)*x_max)+(5*p_idx_z))+10)*y_max)+((((5*p_idx_z)+p_idx_y)+13)*x_max))+(25*p_idx_z))+(5*p_idx_y))+p_idx_x)+68) */
						{
						  int __tmp39 = (x_max+5);
						  int __tmp40 = (__tmp39*y_max);
						  int __tmp41 = (_idx5+__tmp40);
						  int __tmp42 = (5*x_max);
						  int __tmp43 = (__tmp41+__tmp42);
						  int __tmp44 = (__tmp43+25);
						  _idx8=__tmp44;
						}
						/* _idx9 = ((((((((((p_idx_z+3)*x_max)+(5*p_idx_z))+15)*y_max)+((((5*p_idx_z)+p_idx_y)+18)*x_max))+(25*p_idx_z))+(5*p_idx_y))+p_idx_x)+93) */
						{
						  int __tmp45 = (_idx7+x_max);
						  int __tmp46 = (__tmp45+5);
						  _idx9=__tmp46;
						}
						/* _idx10 = ((((((((((p_idx_z+3)*x_max)+(5*p_idx_z))+15)*y_max)+((((5*p_idx_z)+p_idx_y)+18)*x_max))+(25*p_idx_z))+(5*p_idx_y))+p_idx_x)+94) */
						{
						  int __tmp47 = (_idx7+x_max);
						  int __tmp48 = (__tmp47+6);
						  _idx10=__tmp48;
						}
						/* _idx11 = ((((((((((p_idx_z+3)*x_max)+(5*p_idx_z))+15)*y_max)+((((5*p_idx_z)+p_idx_y)+19)*x_max))+(25*p_idx_z))+(5*p_idx_y))+p_idx_x)+98) */
						{
						  int __tmp49 = (_idx10+x_max);
						  int __tmp50 = (__tmp49+4);
						  _idx11=__tmp50;
						}
						/* _idx12 = ((((((((((p_idx_z+4)*x_max)+(5*p_idx_z))+20)*y_max)+((((5*p_idx_z)+p_idx_y)+23)*x_max))+(25*p_idx_z))+(5*p_idx_y))+p_idx_x)+118) */
						{
						  int __tmp51 = (x_max+5);
						  int __tmp52 = (__tmp51*y_max);
						  int __tmp53 = (_idx1+__tmp52);
						  int __tmp54 = (8*x_max);
						  int __tmp55 = (__tmp53+__tmp54);
						  int __tmp56 = (__tmp55+40);
						  _idx12=__tmp56;
						}
						/* _idx13 = ((((((((((p_idx_z+3)*x_max)+(5*p_idx_z))+15)*y_max)+((((5*p_idx_z)+p_idx_y)+18)*x_max))+(25*p_idx_z))+(5*p_idx_y))+p_idx_x)+95) */
						{
						  int __tmp57 = (_idx7+x_max);
						  int __tmp58 = (__tmp57+7);
						  _idx13=__tmp58;
						}
						/* _idx14 = ((((((((((p_idx_z+3)*x_max)+(5*p_idx_z))+15)*y_max)+((((5*p_idx_z)+p_idx_y)+20)*x_max))+(25*p_idx_z))+(5*p_idx_y))+p_idx_x)+103) */
						{
						  int __tmp59 = (_idx11+x_max);
						  int __tmp60 = (__tmp59+5);
						  _idx14=__tmp60;
						}
						/* _idx15 = ((((((((((p_idx_z+5)*x_max)+(5*p_idx_z))+25)*y_max)+((((5*p_idx_z)+p_idx_y)+28)*x_max))+(25*p_idx_z))+(5*p_idx_y))+p_idx_x)+143) */
						{
						  int __tmp61 = (x_max+5);
						  int __tmp62 = (__tmp61*y_max);
						  int __tmp63 = (_idx12+__tmp62);
						  int __tmp64 = (5*x_max);
						  int __tmp65 = (__tmp63+__tmp64);
						  int __tmp66 = (__tmp65+25);
						  _idx15=__tmp66;
						}
						
						u_0_1[_idx9]=(a*((((u_0_0[_idx0]+(u_0_0[_idx1]+u_0_0[_idx2]))*-2.0)+(((u_0_0[_idx3]+(u_0_0[_idx4]+u_0_0[_idx5]))*15.0)+((u_0_0[_idx6]+(u_0_0[_idx7]+u_0_0[_idx8]))*-60.0)))+((u_0_0[_idx9]*20.0)+(((u_0_0[_idx10]+(u_0_0[_idx11]+u_0_0[_idx12]))*30.0)+((u_0_0[_idx13]+(u_0_0[_idx14]+u_0_0[_idx15]))*-3.0)))));
					}
				}
			}
		}
	}
	*u_0_1_out = u_0_1;
}

void initialize(double *  u_0_0, double *  u_0_1, double a, int x_max, int y_max, int z_max)
{
	int _idx0;
	int _idx1;
	int _idx10;
	int _idx11;
	int _idx12;
	int _idx13;
	int _idx14;
	int _idx15;
	int _idx2;
	int _idx3;
	int _idx4;
	int _idx5;
	int _idx6;
	int _idx7;
	int _idx8;
	int _idx9;
	int p_idx_x;
	int p_idx_y;
	int p_idx_z;
	int t;
	double *  __restrict__ const u__u_0[2] =  { u_0_0, u_0_1 } ;
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
						/* _idx0 = ((((((((((p_idx_z+3)*x_max)+(5*p_idx_z))+15)*y_max)+((((5*p_idx_z)+p_idx_y)+18)*x_max))+(25*p_idx_z))+(5*p_idx_y))+p_idx_x)+90) */
						_idx0=((((((((((p_idx_z+3)*x_max)+(5*p_idx_z))+15)*y_max)+((((5*p_idx_z)+p_idx_y)+18)*x_max))+(25*p_idx_z))+(5*p_idx_y))+p_idx_x)+90);
						u__u_0[0][_idx0]=0.1;
						/* _idx1 = ((((((((((p_idx_z+3)*x_max)+(5*p_idx_z))+15)*y_max)+((((5*p_idx_z)+p_idx_y)+18)*x_max))+(25*p_idx_z))+(5*p_idx_y))+p_idx_x)+91) */
						_idx1=(_idx0+1);
						u__u_0[0][_idx1]=0.1;
						/* _idx2 = ((((((((((p_idx_z+3)*x_max)+(5*p_idx_z))+15)*y_max)+((((5*p_idx_z)+p_idx_y)+18)*x_max))+(25*p_idx_z))+(5*p_idx_y))+p_idx_x)+92) */
						_idx2=(_idx1+1);
						u__u_0[0][_idx2]=0.1;
						/* _idx3 = ((((((((((p_idx_z+3)*x_max)+(5*p_idx_z))+15)*y_max)+((((5*p_idx_z)+p_idx_y)+15)*x_max))+(25*p_idx_z))+(5*p_idx_y))+p_idx_x)+78) */
						_idx3=((_idx1-(3*x_max))-13);
						u__u_0[0][_idx3]=0.1;
						/* _idx4 = ((((((((((p_idx_z+3)*x_max)+(5*p_idx_z))+15)*y_max)+((((5*p_idx_z)+p_idx_y)+16)*x_max))+(25*p_idx_z))+(5*p_idx_y))+p_idx_x)+83) */
						_idx4=((_idx3+x_max)+5);
						u__u_0[0][_idx4]=0.1;
						/* _idx5 = ((((((((((p_idx_z+3)*x_max)+(5*p_idx_z))+15)*y_max)+((((5*p_idx_z)+p_idx_y)+17)*x_max))+(25*p_idx_z))+(5*p_idx_y))+p_idx_x)+88) */
						_idx5=((_idx4+x_max)+5);
						u__u_0[0][_idx5]=0.1;
						/* _idx6 = ((((((((p_idx_z*x_max)+(5*p_idx_z))*y_max)+((((5*p_idx_z)+p_idx_y)+3)*x_max))+(25*p_idx_z))+(5*p_idx_y))+p_idx_x)+18) */
						_idx6=(((_idx5+(((-3*x_max)-15)*y_max))-(14*x_max))-70);
						u__u_0[0][_idx6]=0.1;
						/* _idx7 = ((((((((((p_idx_z+1)*x_max)+(5*p_idx_z))+5)*y_max)+((((5*p_idx_z)+p_idx_y)+8)*x_max))+(25*p_idx_z))+(5*p_idx_y))+p_idx_x)+43) */
						_idx7=(((_idx6+((x_max+5)*y_max))+(5*x_max))+25);
						u__u_0[0][_idx7]=0.1;
						/* _idx8 = ((((((((((p_idx_z+2)*x_max)+(5*p_idx_z))+10)*y_max)+((((5*p_idx_z)+p_idx_y)+13)*x_max))+(25*p_idx_z))+(5*p_idx_y))+p_idx_x)+68) */
						_idx8=(((_idx7+((x_max+5)*y_max))+(5*x_max))+25);
						u__u_0[0][_idx8]=0.1;
						/* _idx9 = ((((((((((p_idx_z+3)*x_max)+(5*p_idx_z))+15)*y_max)+((((5*p_idx_z)+p_idx_y)+18)*x_max))+(25*p_idx_z))+(5*p_idx_y))+p_idx_x)+93) */
						_idx9=((_idx5+x_max)+5);
						u__u_0[0][_idx9]=0.1;
						/* _idx10 = ((((((((((p_idx_z+4)*x_max)+(5*p_idx_z))+20)*y_max)+((((5*p_idx_z)+p_idx_y)+23)*x_max))+(25*p_idx_z))+(5*p_idx_y))+p_idx_x)+118) */
						_idx10=(((_idx5+((x_max+5)*y_max))+(6*x_max))+30);
						u__u_0[0][_idx10]=0.1;
						/* _idx11 = ((((((((((p_idx_z+5)*x_max)+(5*p_idx_z))+25)*y_max)+((((5*p_idx_z)+p_idx_y)+28)*x_max))+(25*p_idx_z))+(5*p_idx_y))+p_idx_x)+143) */
						_idx11=(((_idx10+((x_max+5)*y_max))+(5*x_max))+25);
						u__u_0[0][_idx11]=0.1;
						/* _idx12 = ((((((((((p_idx_z+3)*x_max)+(5*p_idx_z))+15)*y_max)+((((5*p_idx_z)+p_idx_y)+19)*x_max))+(25*p_idx_z))+(5*p_idx_y))+p_idx_x)+98) */
						_idx12=((_idx9+x_max)+5);
						u__u_0[0][_idx12]=0.1;
						/* _idx13 = ((((((((((p_idx_z+3)*x_max)+(5*p_idx_z))+15)*y_max)+((((5*p_idx_z)+p_idx_y)+20)*x_max))+(25*p_idx_z))+(5*p_idx_y))+p_idx_x)+103) */
						_idx13=((_idx12+x_max)+5);
						u__u_0[0][_idx13]=0.1;
						/* _idx14 = ((((((((((p_idx_z+3)*x_max)+(5*p_idx_z))+15)*y_max)+((((5*p_idx_z)+p_idx_y)+18)*x_max))+(25*p_idx_z))+(5*p_idx_y))+p_idx_x)+94) */
						_idx14=((_idx5+x_max)+6);
						u__u_0[0][_idx14]=0.1;
						/* _idx15 = ((((((((((p_idx_z+3)*x_max)+(5*p_idx_z))+15)*y_max)+((((5*p_idx_z)+p_idx_y)+18)*x_max))+(25*p_idx_z))+(5*p_idx_y))+p_idx_x)+95) */
						_idx15=((_idx5+x_max)+7);
						u__u_0[0][_idx15]=0.1;
						u__u_0[1][_idx9]=1.1;
					}
				}
			}
		}
	}
}

void upstream_5_3d_cpu(double *  *  u_0_1_out, double *  u_0_0, double *  u_0_1, double a, int x_max, int y_max, int z_max, int t_max)
{
	int _idx0;
	int _idx1;
	int _idx10;
	int _idx11;
	int _idx12;
	int _idx13;
	int _idx14;
	int _idx15;
	int _idx2;
	int _idx3;
	int _idx4;
	int _idx5;
	int _idx6;
	int _idx7;
	int _idx8;
	int _idx9;
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
		/*
		for POINT p[t=t, s=(1, 1, 1)][0] of size [1, 1, 1] in u[t=t, s=(:, :, :)][0] parallel 1 <level 0> schedule default { ... }
		*/
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
						/* _idx0 = ((((((((((p_idx_z+3)*x_max)+(5*p_idx_z))+15)*y_max)+((((5*p_idx_z)+p_idx_y)+18)*x_max))+(25*p_idx_z))+(5*p_idx_y))+p_idx_x)+90) */
						_idx0=((((((((((p_idx_z+3)*x_max)+(5*p_idx_z))+15)*y_max)+((((5*p_idx_z)+p_idx_y)+18)*x_max))+(25*p_idx_z))+(5*p_idx_y))+p_idx_x)+90);
						/* _idx1 = ((((((((((p_idx_z+3)*x_max)+(5*p_idx_z))+15)*y_max)+((((5*p_idx_z)+p_idx_y)+15)*x_max))+(25*p_idx_z))+(5*p_idx_y))+p_idx_x)+78) */
						_idx1=((_idx0-(3*x_max))-12);
						/* _idx2 = ((((((((p_idx_z*x_max)+(5*p_idx_z))*y_max)+((((5*p_idx_z)+p_idx_y)+3)*x_max))+(25*p_idx_z))+(5*p_idx_y))+p_idx_x)+18) */
						_idx2=(((_idx1+(((-3*x_max)-15)*y_max))-(12*x_max))-60);
						/* _idx3 = ((((((((((p_idx_z+3)*x_max)+(5*p_idx_z))+15)*y_max)+((((5*p_idx_z)+p_idx_y)+18)*x_max))+(25*p_idx_z))+(5*p_idx_y))+p_idx_x)+91) */
						_idx3=(_idx0+1);
						/* _idx4 = ((((((((((p_idx_z+3)*x_max)+(5*p_idx_z))+15)*y_max)+((((5*p_idx_z)+p_idx_y)+16)*x_max))+(25*p_idx_z))+(5*p_idx_y))+p_idx_x)+83) */
						_idx4=((_idx1+x_max)+5);
						/* _idx5 = ((((((((((p_idx_z+1)*x_max)+(5*p_idx_z))+5)*y_max)+((((5*p_idx_z)+p_idx_y)+8)*x_max))+(25*p_idx_z))+(5*p_idx_y))+p_idx_x)+43) */
						_idx5=(((_idx2+((x_max+5)*y_max))+(5*x_max))+25);
						/* _idx6 = ((((((((((p_idx_z+3)*x_max)+(5*p_idx_z))+15)*y_max)+((((5*p_idx_z)+p_idx_y)+18)*x_max))+(25*p_idx_z))+(5*p_idx_y))+p_idx_x)+92) */
						_idx6=(_idx3+1);
						/* _idx7 = ((((((((((p_idx_z+3)*x_max)+(5*p_idx_z))+15)*y_max)+((((5*p_idx_z)+p_idx_y)+17)*x_max))+(25*p_idx_z))+(5*p_idx_y))+p_idx_x)+88) */
						_idx7=((_idx4+x_max)+5);
						/* _idx8 = ((((((((((p_idx_z+2)*x_max)+(5*p_idx_z))+10)*y_max)+((((5*p_idx_z)+p_idx_y)+13)*x_max))+(25*p_idx_z))+(5*p_idx_y))+p_idx_x)+68) */
						_idx8=(((_idx5+((x_max+5)*y_max))+(5*x_max))+25);
						/* _idx9 = ((((((((((p_idx_z+3)*x_max)+(5*p_idx_z))+15)*y_max)+((((5*p_idx_z)+p_idx_y)+18)*x_max))+(25*p_idx_z))+(5*p_idx_y))+p_idx_x)+93) */
						_idx9=((_idx7+x_max)+5);
						/* _idx10 = ((((((((((p_idx_z+3)*x_max)+(5*p_idx_z))+15)*y_max)+((((5*p_idx_z)+p_idx_y)+18)*x_max))+(25*p_idx_z))+(5*p_idx_y))+p_idx_x)+94) */
						_idx10=((_idx7+x_max)+6);
						/* _idx11 = ((((((((((p_idx_z+3)*x_max)+(5*p_idx_z))+15)*y_max)+((((5*p_idx_z)+p_idx_y)+19)*x_max))+(25*p_idx_z))+(5*p_idx_y))+p_idx_x)+98) */
						_idx11=((_idx10+x_max)+4);
						/* _idx12 = ((((((((((p_idx_z+4)*x_max)+(5*p_idx_z))+20)*y_max)+((((5*p_idx_z)+p_idx_y)+23)*x_max))+(25*p_idx_z))+(5*p_idx_y))+p_idx_x)+118) */
						_idx12=(((_idx1+((x_max+5)*y_max))+(8*x_max))+40);
						/* _idx13 = ((((((((((p_idx_z+3)*x_max)+(5*p_idx_z))+15)*y_max)+((((5*p_idx_z)+p_idx_y)+18)*x_max))+(25*p_idx_z))+(5*p_idx_y))+p_idx_x)+95) */
						_idx13=((_idx7+x_max)+7);
						/* _idx14 = ((((((((((p_idx_z+3)*x_max)+(5*p_idx_z))+15)*y_max)+((((5*p_idx_z)+p_idx_y)+20)*x_max))+(25*p_idx_z))+(5*p_idx_y))+p_idx_x)+103) */
						_idx14=((_idx11+x_max)+5);
						/* _idx15 = ((((((((((p_idx_z+5)*x_max)+(5*p_idx_z))+25)*y_max)+((((5*p_idx_z)+p_idx_y)+28)*x_max))+(25*p_idx_z))+(5*p_idx_y))+p_idx_x)+143) */
						_idx15=(((_idx12+((x_max+5)*y_max))+(5*x_max))+25);
						u_0_1[_idx9]=(a*((((u_0_0[_idx0]+(u_0_0[_idx1]+u_0_0[_idx2]))*-2.0)+(((u_0_0[_idx3]+(u_0_0[_idx4]+u_0_0[_idx5]))*15.0)+((u_0_0[_idx6]+(u_0_0[_idx7]+u_0_0[_idx8]))*-60.0)))+((u_0_0[_idx9]*20.0)+(((u_0_0[_idx10]+(u_0_0[_idx11]+u_0_0[_idx12]))*30.0)+((u_0_0[_idx13]+(u_0_0[_idx14]+u_0_0[_idx15]))*-3.0)))));
					}
				}
			}
		}
	}
	*u_0_1_out = u_0_1;
}
