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
#pragma acc region copyin(u_0_0[0:(x_max+5)*(y_max+5)*(z_max+5)]) copyout(u_0_1[0:(x_max+5)*(y_max+5)*(z_max+5)])
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
