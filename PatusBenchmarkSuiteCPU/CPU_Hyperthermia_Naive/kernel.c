#include "omp.h"
#define t_max 1
/*
(T[0][0][0][1][0]=((((T[0][0][0][0][0]*((c[0][0][0][0][1]*T[0][0][0][0][0])+c[0][0][0][0][2]))+c[0][0][0][0][3])+((c[0][0][0][0][4]*T[-1][0][0][0][0])+(c[0][0][0][0][5]*T[1][0][0][0][0])))+(((c[0][0][0][0][6]*T[0][-1][0][0][0])+(c[0][0][0][0][7]*T[0][1][0][0][0]))+((c[0][0][0][0][8]*T[0][0][-1][0][0])+(c[0][0][0][0][9]*T[0][0][1][0][0])))))

*/
void hyperthermia(float *  *  T_0_1_out, float *  T_0_0, float *  T_0_1, float *  c_1_0, float *  c_2_0, float *  c_3_0, float *  c_4_0, float *  c_5_0, float *  c_6_0, float *  c_7_0, float *  c_8_0, float *  c_9_0, int x_max, int y_max, int z_max)
{
	int _idx0;
	int _idx1;
	int _idx2;
	int _idx3;
	int _idx4;
	int _idx5;
	int _idx6;
	int _idx7;
	int dimidx0 = omp_get_thread_num();
	int dimsize0 = omp_get_num_threads();
	int p_idx;
	int p_idx_x;
	int p_idx_x_max;
	int p_idx_y;
	int p_idx_y_max;
	int p_idx_z;
	int p_idx_z_max;
	int p_numblocks;
	int t;
	float *  __restrict__ const u__T_0[16] =  { T_0_0, T_0_1 } ;
	const float *  __restrict__ const u__c_1[16] =  { c_1_0 } ;
	const float *  __restrict__ const u__c_2[16] =  { c_2_0 } ;
	const float *  __restrict__ const u__c_3[16] =  { c_3_0 } ;
	const float *  __restrict__ const u__c_4[16] =  { c_4_0 } ;
	const float *  __restrict__ const u__c_5[16] =  { c_5_0 } ;
	const float *  __restrict__ const u__c_6[16] =  { c_6_0 } ;
	const float *  __restrict__ const u__c_7[16] =  { c_7_0 } ;
	const float *  __restrict__ const u__c_8[16] =  { c_8_0 } ;
	const float *  __restrict__ const u__c_9[16] =  { c_9_0 } ;
	/*
	Initializations
	*/
	p_numblocks=((((int)x_max)*((int)y_max))*((int)z_max));
	/*
	Implementation
	*/
	/*
	for t = 1..t_max by 1 parallel 1 <level 0> schedule  { ... }
	*/
	for (t=1; t<=t_max; t+=1)
	{
		/*
		for p_idx = dimidx0..(p_numblocks-1) by dimsize0 parallel 1 <level 1> schedule  { ... }
		*/
		for (p_idx=dimidx0; p_idx<=(p_numblocks-1); p_idx+=dimsize0)
		{
			/* Index bounds calculations for iterators in p[t=t, s=(1, 1, 1)][0] */
			p_idx_x=(0+((p_idx%((int)x_max))*1));
			p_idx_x_max=(p_idx_x+1);
			p_idx_y=(0+(((p_idx/((int)x_max))%((int)y_max))*1));
			p_idx_y_max=(p_idx_y+1);
			p_idx_z=(0+(((p_idx/(((int)x_max)*((int)y_max)))%((int)z_max))*1));
			p_idx_z_max=(p_idx_z+1);
			/*
			u[t=(t+1), s=p[t=?, s=?][0]][0]=stencil(u[t=t, s=p[t=?, s=?][0]][0])
			*/
			/* _idx0 = (((((((((p_idx_z+1)*x_max)+(((2*p_idx_z)+2)*t))*y_max)+((((((2*p_idx_z)+2)*t)+p_idx_y)+1)*x_max))+(((4*p_idx_z)+4)*(t*t)))+(((2*p_idx_y)+2)*t))+p_idx_x)+1) */
			_idx0=(((((((((p_idx_z+1)*x_max)+(((2*p_idx_z)+2)*t))*y_max)+((((((2*p_idx_z)+2)*t)+p_idx_y)+1)*x_max))+(((4*p_idx_z)+4)*(t*t)))+(((2*p_idx_y)+2)*t))+p_idx_x)+1);
			/* _idx1 = ((((p_idx_z*x_max)*y_max)+(p_idx_y*x_max))+p_idx_x) */
			_idx1=(((((_idx0+(((((-2*p_idx_z)-2)*t)-x_max)*y_max))+(((((-2*p_idx_z)-2)*t)-1)*x_max))+(((-4*p_idx_z)-4)*(t*t)))+(((-2*p_idx_y)-2)*t))-1);
			/* _idx2 = ((((((((p_idx_z+1)*x_max)+(((2*p_idx_z)+2)*t))*y_max)+((((((2*p_idx_z)+2)*t)+p_idx_y)+1)*x_max))+(((4*p_idx_z)+4)*(t*t)))+(((2*p_idx_y)+2)*t))+p_idx_x) */
			_idx2=(_idx0-1);
			/* _idx3 = (((((((((p_idx_z+1)*x_max)+(((2*p_idx_z)+2)*t))*y_max)+((((((2*p_idx_z)+2)*t)+p_idx_y)+1)*x_max))+(((4*p_idx_z)+4)*(t*t)))+(((2*p_idx_y)+2)*t))+p_idx_x)+2) */
			_idx3=(_idx2+2);
			/* _idx4 = (((((((((p_idx_z+1)*x_max)+(((2*p_idx_z)+2)*t))*y_max)+(((((2*p_idx_z)+2)*t)+p_idx_y)*x_max))+(((4*p_idx_z)+4)*(t*t)))+((2*p_idx_y)*t))+p_idx_x)+1) */
			_idx4=((_idx0-x_max)-(2*t));
			/* _idx5 = (((((((((p_idx_z+1)*x_max)+(((2*p_idx_z)+2)*t))*y_max)+((((((2*p_idx_z)+2)*t)+p_idx_y)+2)*x_max))+(((4*p_idx_z)+4)*(t*t)))+(((2*p_idx_y)+4)*t))+p_idx_x)+1) */
			_idx5=((_idx0+x_max)+(2*t));
			/* _idx6 = ((((((((p_idx_z*x_max)+((2*p_idx_z)*t))*y_max)+(((((2*p_idx_z)*t)+p_idx_y)+1)*x_max))+((4*p_idx_z)*(t*t)))+(((2*p_idx_y)+2)*t))+p_idx_x)+1) */
			_idx6=((((_idx2+((( - x_max)-(2*t))*y_max))-((2*t)*x_max))-(4*(t*t)))+1);
			/* _idx7 = (((((((((p_idx_z+2)*x_max)+(((2*p_idx_z)+4)*t))*y_max)+((((((2*p_idx_z)+4)*t)+p_idx_y)+1)*x_max))+(((4*p_idx_z)+8)*(t*t)))+(((2*p_idx_y)+2)*t))+p_idx_x)+1) */
			_idx7=(((_idx0+((x_max+(2*t))*y_max))+((2*t)*x_max))+(4*(t*t)));
			u__T_0[t][_idx0]=((((u__T_0[(t-1)][_idx0]*((u__c_1[(t-1)][_idx1]*u__T_0[(t-1)][_idx0])+u__c_2[(t-1)][_idx1]))+u__c_3[(t-1)][_idx1])+((u__c_4[(t-1)][_idx1]*u__T_0[(t-1)][_idx2])+(u__c_5[(t-1)][_idx1]*u__T_0[(t-1)][_idx3])))+(((u__c_6[(t-1)][_idx1]*u__T_0[(t-1)][_idx4])+(u__c_7[(t-1)][_idx1]*u__T_0[(t-1)][_idx5]))+((u__c_8[(t-1)][_idx1]*u__T_0[(t-1)][_idx6])+(u__c_9[(t-1)][_idx1]*u__T_0[(t-1)][_idx7]))));
		}
#pragma omp barrier 
	}
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
	int dimidx0 = omp_get_thread_num();
	int dimsize0 = omp_get_num_threads();
	int p_idx;
	int p_idx_x;
	int p_idx_x_max;
	int p_idx_y;
	int p_idx_y_max;
	int p_idx_z;
	int p_idx_z_max;
	int p_numblocks;
	int t;
	float *  __restrict__ const u__T_0[16] =  { T_0_0, T_0_1 } ;
	float *  __restrict__ const u__c_1[16] =  { c_1_0 } ;
	float *  __restrict__ const u__c_2[16] =  { c_2_0 } ;
	float *  __restrict__ const u__c_3[16] =  { c_3_0 } ;
	float *  __restrict__ const u__c_4[16] =  { c_4_0 } ;
	float *  __restrict__ const u__c_5[16] =  { c_5_0 } ;
	float *  __restrict__ const u__c_6[16] =  { c_6_0 } ;
	float *  __restrict__ const u__c_7[16] =  { c_7_0 } ;
	float *  __restrict__ const u__c_8[16] =  { c_8_0 } ;
	float *  __restrict__ const u__c_9[16] =  { c_9_0 } ;
	/*
	Initializations
	*/
	p_numblocks=((((int)x_max)*((int)y_max))*((int)z_max));
	/*
	Implementation
	*/
	/*
	for t = 1..t_max by 1 parallel 1 <level 0> schedule  { ... }
	*/
	for (t=1; t<=t_max; t+=1)
	{
		/*
		for p_idx = dimidx0..(p_numblocks-1) by dimsize0 parallel 1 <level 1> schedule  { ... }
		*/
		for (p_idx=dimidx0; p_idx<=(p_numblocks-1); p_idx+=dimsize0)
		{
			/* Index bounds calculations for iterators in p[t=t, s=(1, 1, 1)][0] */
			p_idx_x=(0+((p_idx%((int)x_max))*1));
			p_idx_x_max=(p_idx_x+1);
			p_idx_y=(0+(((p_idx/((int)x_max))%((int)y_max))*1));
			p_idx_y_max=(p_idx_y+1);
			p_idx_z=(0+(((p_idx/(((int)x_max)*((int)y_max)))%((int)z_max))*1));
			p_idx_z_max=(p_idx_z+1);
			/*
			u[t=(t+1), s=p[t=?, s=?][0]][0]=stencil(u[t=t, s=p[t=?, s=?][0]][0])
			*/
			/* _idx0 = ((((((((p_idx_z+1)*x_max)+(((2*p_idx_z)+2)*t))*y_max)+((((((2*p_idx_z)+2)*t)+p_idx_y)+1)*x_max))+(((4*p_idx_z)+4)*(t*t)))+(((2*p_idx_y)+2)*t))+p_idx_x) */
			_idx0=((((((((p_idx_z+1)*x_max)+(((2*p_idx_z)+2)*t))*y_max)+((((((2*p_idx_z)+2)*t)+p_idx_y)+1)*x_max))+(((4*p_idx_z)+4)*(t*t)))+(((2*p_idx_y)+2)*t))+p_idx_x);
			u__T_0[(t-1)][_idx0]=0.1;
			/* _idx1 = (((((((((p_idx_z+1)*x_max)+(((2*p_idx_z)+2)*t))*y_max)+(((((2*p_idx_z)+2)*t)+p_idx_y)*x_max))+(((4*p_idx_z)+4)*(t*t)))+((2*p_idx_y)*t))+p_idx_x)+1) */
			_idx1=(((_idx0-x_max)-(2*t))+1);
			u__T_0[(t-1)][_idx1]=0.1;
			/* _idx2 = ((((((((p_idx_z*x_max)+((2*p_idx_z)*t))*y_max)+(((((2*p_idx_z)*t)+p_idx_y)+1)*x_max))+((4*p_idx_z)*(t*t)))+(((2*p_idx_y)+2)*t))+p_idx_x)+1) */
			_idx2=((((_idx0+((( - x_max)-(2*t))*y_max))-((2*t)*x_max))-(4*(t*t)))+1);
			u__T_0[(t-1)][_idx2]=0.1;
			/* _idx3 = (((((((((p_idx_z+1)*x_max)+(((2*p_idx_z)+2)*t))*y_max)+((((((2*p_idx_z)+2)*t)+p_idx_y)+1)*x_max))+(((4*p_idx_z)+4)*(t*t)))+(((2*p_idx_y)+2)*t))+p_idx_x)+1) */
			_idx3=(_idx0+1);
			u__T_0[(t-1)][_idx3]=0.1;
			/* _idx4 = ((((p_idx_z*x_max)*y_max)+(p_idx_y*x_max))+p_idx_x) */
			_idx4=(((((_idx2-(((2*p_idx_z)*t)*y_max))+((((-2*p_idx_z)*t)-1)*x_max))-((4*p_idx_z)*(t*t)))+(((-2*p_idx_y)-2)*t))-1);
			u__c_1[(t-1)][_idx4]=0.2;
			u__c_2[(t-1)][_idx4]=0.30000000000000004;
			u__c_3[(t-1)][_idx4]=0.4;
			u__c_4[(t-1)][_idx4]=0.5;
			u__c_5[(t-1)][_idx4]=0.6000000000000001;
			u__c_6[(t-1)][_idx4]=0.7000000000000001;
			u__c_7[(t-1)][_idx4]=0.8;
			u__c_8[(t-1)][_idx4]=0.9;
			u__c_9[(t-1)][_idx4]=1.0;
			/* _idx5 = (((((((((p_idx_z+2)*x_max)+(((2*p_idx_z)+4)*t))*y_max)+((((((2*p_idx_z)+4)*t)+p_idx_y)+1)*x_max))+(((4*p_idx_z)+8)*(t*t)))+(((2*p_idx_y)+2)*t))+p_idx_x)+1) */
			_idx5=(((_idx3+((x_max+(2*t))*y_max))+((2*t)*x_max))+(4*(t*t)));
			u__T_0[(t-1)][_idx5]=0.1;
			/* _idx6 = (((((((((p_idx_z+1)*x_max)+(((2*p_idx_z)+2)*t))*y_max)+((((((2*p_idx_z)+2)*t)+p_idx_y)+2)*x_max))+(((4*p_idx_z)+4)*(t*t)))+(((2*p_idx_y)+4)*t))+p_idx_x)+1) */
			_idx6=((_idx3+x_max)+(2*t));
			u__T_0[(t-1)][_idx6]=0.1;
			/* _idx7 = (((((((((p_idx_z+1)*x_max)+(((2*p_idx_z)+2)*t))*y_max)+((((((2*p_idx_z)+2)*t)+p_idx_y)+1)*x_max))+(((4*p_idx_z)+4)*(t*t)))+(((2*p_idx_y)+2)*t))+p_idx_x)+2) */
			_idx7=(_idx0+2);
			u__T_0[(t-1)][_idx7]=0.1;
			u__T_0[t][_idx3]=1.1;
		}
#pragma omp barrier 
	}
}

