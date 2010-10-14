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
	const double *  __restrict__ const u__a_1[16] =  { a_1_0 } ;
	const double *  __restrict__ const u__b_2[16] =  { b_2_0 } ;
	const double *  __restrict__ const u__c_3[16] =  { c_3_0 } ;
	double *  __restrict__ const u__u_0[16] =  { u_0_0, u_0_1 } ;
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
			/* _idx0 = ((((p_idx_z*x_max)*y_max)+(p_idx_y*x_max))+p_idx_x) */
			_idx0=((((p_idx_z*x_max)*y_max)+(p_idx_y*x_max))+p_idx_x);
			w1_a=((u__a_1[(t-1)][_idx0]*(u__a_1[(t-1)][_idx0]+1.0))*((u__a_1[(t-1)][_idx0]+2.0)*0.16666666666666666));
			w2_a=(((u__a_1[(t-1)][_idx0]-1.0)*(u__a_1[(t-1)][_idx0]+1.0))*((u__a_1[(t-1)][_idx0]+2.0)*-0.5));
			w3_a=(((u__a_1[(t-1)][_idx0]-1.0)*u__a_1[(t-1)][_idx0])*((u__a_1[(t-1)][_idx0]+2.0)*0.5));
			w4_a=(((u__a_1[(t-1)][_idx0]-1.0)*u__a_1[(t-1)][_idx0])*((u__a_1[(t-1)][_idx0]+1.0)*-0.16666666666666666));
			w1_b=((u__b_2[(t-1)][_idx0]*(u__b_2[(t-1)][_idx0]+1.0))*((u__b_2[(t-1)][_idx0]+2.0)*0.16666666666666666));
			w2_b=(((u__b_2[(t-1)][_idx0]-1.0)*(u__b_2[(t-1)][_idx0]+1.0))*((u__b_2[(t-1)][_idx0]+2.0)*-0.5));
			w3_b=(((u__b_2[(t-1)][_idx0]-1.0)*u__b_2[(t-1)][_idx0])*((u__b_2[(t-1)][_idx0]+2.0)*0.5));
			w4_b=(((u__b_2[(t-1)][_idx0]-1.0)*u__b_2[(t-1)][_idx0])*((u__b_2[(t-1)][_idx0]+1.0)*-0.16666666666666666));
			w1_c=((u__c_3[(t-1)][_idx0]*(u__c_3[(t-1)][_idx0]+1.0))*((u__c_3[(t-1)][_idx0]+2.0)*0.16666666666666666));
			w2_c=(((u__c_3[(t-1)][_idx0]-1.0)*(u__c_3[(t-1)][_idx0]+1.0))*((u__c_3[(t-1)][_idx0]+2.0)*-0.5));
			w3_c=(((u__c_3[(t-1)][_idx0]-1.0)*u__c_3[(t-1)][_idx0])*((u__c_3[(t-1)][_idx0]+2.0)*0.5));
			w4_c=(((u__c_3[(t-1)][_idx0]-1.0)*u__c_3[(t-1)][_idx0])*((u__c_3[(t-1)][_idx0]+1.0)*-0.16666666666666666));
			/* _idx1 = (((((((p_idx_z*x_max)+((3*p_idx_z)*t))*y_max)+((((3*p_idx_z)*t)+p_idx_y)*x_max))+((9*p_idx_z)*(t*t)))+((3*p_idx_y)*t))+p_idx_x) */
			_idx1=((((_idx0+(((3*p_idx_z)*t)*y_max))+(((3*p_idx_z)*t)*x_max))+((9*p_idx_z)*(t*t)))+((3*p_idx_y)*t));
			/* _idx2 = ((((((((p_idx_z*x_max)+((3*p_idx_z)*t))*y_max)+((((3*p_idx_z)*t)+p_idx_y)*x_max))+((9*p_idx_z)*(t*t)))+((3*p_idx_y)*t))+p_idx_x)+1) */
			_idx2=(_idx1+1);
			/* _idx3 = ((((((((p_idx_z*x_max)+((3*p_idx_z)*t))*y_max)+((((3*p_idx_z)*t)+p_idx_y)*x_max))+((9*p_idx_z)*(t*t)))+((3*p_idx_y)*t))+p_idx_x)+2) */
			_idx3=(_idx2+1);
			/* _idx4 = ((((((((p_idx_z*x_max)+((3*p_idx_z)*t))*y_max)+((((3*p_idx_z)*t)+p_idx_y)*x_max))+((9*p_idx_z)*(t*t)))+((3*p_idx_y)*t))+p_idx_x)+3) */
			_idx4=(_idx2+2);
			/* _idx5 = (((((((p_idx_z*x_max)+((3*p_idx_z)*t))*y_max)+(((((3*p_idx_z)*t)+p_idx_y)+1)*x_max))+((9*p_idx_z)*(t*t)))+(((3*p_idx_y)+3)*t))+p_idx_x) */
			_idx5=((_idx1+x_max)+(3*t));
			/* _idx6 = ((((((((p_idx_z*x_max)+((3*p_idx_z)*t))*y_max)+(((((3*p_idx_z)*t)+p_idx_y)+1)*x_max))+((9*p_idx_z)*(t*t)))+(((3*p_idx_y)+3)*t))+p_idx_x)+1) */
			_idx6=(_idx5+1);
			/* _idx7 = ((((((((p_idx_z*x_max)+((3*p_idx_z)*t))*y_max)+(((((3*p_idx_z)*t)+p_idx_y)+1)*x_max))+((9*p_idx_z)*(t*t)))+(((3*p_idx_y)+3)*t))+p_idx_x)+2) */
			_idx7=(_idx5+2);
			/* _idx8 = ((((((((p_idx_z*x_max)+((3*p_idx_z)*t))*y_max)+(((((3*p_idx_z)*t)+p_idx_y)+1)*x_max))+((9*p_idx_z)*(t*t)))+(((3*p_idx_y)+3)*t))+p_idx_x)+3) */
			_idx8=(_idx7+1);
			/* _idx9 = (((((((p_idx_z*x_max)+((3*p_idx_z)*t))*y_max)+(((((3*p_idx_z)*t)+p_idx_y)+2)*x_max))+((9*p_idx_z)*(t*t)))+(((3*p_idx_y)+6)*t))+p_idx_x) */
			_idx9=((_idx5+x_max)+(3*t));
			/* _idx10 = ((((((((p_idx_z*x_max)+((3*p_idx_z)*t))*y_max)+(((((3*p_idx_z)*t)+p_idx_y)+2)*x_max))+((9*p_idx_z)*(t*t)))+(((3*p_idx_y)+6)*t))+p_idx_x)+1) */
			_idx10=(_idx9+1);
			/* _idx11 = ((((((((p_idx_z*x_max)+((3*p_idx_z)*t))*y_max)+(((((3*p_idx_z)*t)+p_idx_y)+2)*x_max))+((9*p_idx_z)*(t*t)))+(((3*p_idx_y)+6)*t))+p_idx_x)+2) */
			_idx11=(_idx10+1);
			/* _idx12 = ((((((((p_idx_z*x_max)+((3*p_idx_z)*t))*y_max)+(((((3*p_idx_z)*t)+p_idx_y)+2)*x_max))+((9*p_idx_z)*(t*t)))+(((3*p_idx_y)+6)*t))+p_idx_x)+3) */
			_idx12=(_idx11+1);
			/* _idx13 = (((((((p_idx_z*x_max)+((3*p_idx_z)*t))*y_max)+(((((3*p_idx_z)*t)+p_idx_y)+3)*x_max))+((9*p_idx_z)*(t*t)))+(((3*p_idx_y)+9)*t))+p_idx_x) */
			_idx13=((_idx9+x_max)+(3*t));
			/* _idx14 = ((((((((p_idx_z*x_max)+((3*p_idx_z)*t))*y_max)+(((((3*p_idx_z)*t)+p_idx_y)+3)*x_max))+((9*p_idx_z)*(t*t)))+(((3*p_idx_y)+9)*t))+p_idx_x)+1) */
			_idx14=(_idx13+1);
			/* _idx15 = ((((((((p_idx_z*x_max)+((3*p_idx_z)*t))*y_max)+(((((3*p_idx_z)*t)+p_idx_y)+3)*x_max))+((9*p_idx_z)*(t*t)))+(((3*p_idx_y)+9)*t))+p_idx_x)+2) */
			_idx15=(_idx14+1);
			/* _idx16 = ((((((((p_idx_z*x_max)+((3*p_idx_z)*t))*y_max)+(((((3*p_idx_z)*t)+p_idx_y)+3)*x_max))+((9*p_idx_z)*(t*t)))+(((3*p_idx_y)+9)*t))+p_idx_x)+3) */
			_idx16=(_idx14+2);
			/* _idx17 = ((((((((p_idx_z+1)*x_max)+(((3*p_idx_z)+3)*t))*y_max)+(((((3*p_idx_z)+3)*t)+p_idx_y)*x_max))+(((9*p_idx_z)+9)*(t*t)))+((3*p_idx_y)*t))+p_idx_x) */
			_idx17=(((_idx1+((x_max+(3*t))*y_max))+((3*t)*x_max))+(9*(t*t)));
			/* _idx18 = (((((((((p_idx_z+1)*x_max)+(((3*p_idx_z)+3)*t))*y_max)+(((((3*p_idx_z)+3)*t)+p_idx_y)*x_max))+(((9*p_idx_z)+9)*(t*t)))+((3*p_idx_y)*t))+p_idx_x)+1) */
			_idx18=(_idx17+1);
			/* _idx19 = (((((((((p_idx_z+1)*x_max)+(((3*p_idx_z)+3)*t))*y_max)+(((((3*p_idx_z)+3)*t)+p_idx_y)*x_max))+(((9*p_idx_z)+9)*(t*t)))+((3*p_idx_y)*t))+p_idx_x)+2) */
			_idx19=(_idx18+1);
			/* _idx20 = (((((((((p_idx_z+1)*x_max)+(((3*p_idx_z)+3)*t))*y_max)+(((((3*p_idx_z)+3)*t)+p_idx_y)*x_max))+(((9*p_idx_z)+9)*(t*t)))+((3*p_idx_y)*t))+p_idx_x)+3) */
			_idx20=(_idx19+1);
			/* _idx21 = ((((((((p_idx_z+1)*x_max)+(((3*p_idx_z)+3)*t))*y_max)+((((((3*p_idx_z)+3)*t)+p_idx_y)+1)*x_max))+(((9*p_idx_z)+9)*(t*t)))+(((3*p_idx_y)+3)*t))+p_idx_x) */
			_idx21=((_idx17+x_max)+(3*t));
			/* _idx22 = (((((((((p_idx_z+1)*x_max)+(((3*p_idx_z)+3)*t))*y_max)+((((((3*p_idx_z)+3)*t)+p_idx_y)+1)*x_max))+(((9*p_idx_z)+9)*(t*t)))+(((3*p_idx_y)+3)*t))+p_idx_x)+1) */
			_idx22=(_idx21+1);
			/* _idx23 = (((((((((p_idx_z+1)*x_max)+(((3*p_idx_z)+3)*t))*y_max)+((((((3*p_idx_z)+3)*t)+p_idx_y)+1)*x_max))+(((9*p_idx_z)+9)*(t*t)))+(((3*p_idx_y)+3)*t))+p_idx_x)+2) */
			_idx23=(_idx22+1);
			/* _idx24 = (((((((((p_idx_z+1)*x_max)+(((3*p_idx_z)+3)*t))*y_max)+((((((3*p_idx_z)+3)*t)+p_idx_y)+1)*x_max))+(((9*p_idx_z)+9)*(t*t)))+(((3*p_idx_y)+3)*t))+p_idx_x)+3) */
			_idx24=(_idx23+1);
			/* _idx25 = ((((((((p_idx_z+1)*x_max)+(((3*p_idx_z)+3)*t))*y_max)+((((((3*p_idx_z)+3)*t)+p_idx_y)+2)*x_max))+(((9*p_idx_z)+9)*(t*t)))+(((3*p_idx_y)+6)*t))+p_idx_x) */
			_idx25=((_idx21+x_max)+(3*t));
			/* _idx26 = (((((((((p_idx_z+1)*x_max)+(((3*p_idx_z)+3)*t))*y_max)+((((((3*p_idx_z)+3)*t)+p_idx_y)+2)*x_max))+(((9*p_idx_z)+9)*(t*t)))+(((3*p_idx_y)+6)*t))+p_idx_x)+1) */
			_idx26=(_idx25+1);
			/* _idx27 = (((((((((p_idx_z+1)*x_max)+(((3*p_idx_z)+3)*t))*y_max)+((((((3*p_idx_z)+3)*t)+p_idx_y)+2)*x_max))+(((9*p_idx_z)+9)*(t*t)))+(((3*p_idx_y)+6)*t))+p_idx_x)+2) */
			_idx27=(_idx25+2);
			/* _idx28 = (((((((((p_idx_z+1)*x_max)+(((3*p_idx_z)+3)*t))*y_max)+((((((3*p_idx_z)+3)*t)+p_idx_y)+2)*x_max))+(((9*p_idx_z)+9)*(t*t)))+(((3*p_idx_y)+6)*t))+p_idx_x)+3) */
			_idx28=(_idx27+1);
			/* _idx29 = ((((((((p_idx_z+1)*x_max)+(((3*p_idx_z)+3)*t))*y_max)+((((((3*p_idx_z)+3)*t)+p_idx_y)+3)*x_max))+(((9*p_idx_z)+9)*(t*t)))+(((3*p_idx_y)+9)*t))+p_idx_x) */
			_idx29=((_idx25+x_max)+(3*t));
			/* _idx30 = (((((((((p_idx_z+1)*x_max)+(((3*p_idx_z)+3)*t))*y_max)+((((((3*p_idx_z)+3)*t)+p_idx_y)+3)*x_max))+(((9*p_idx_z)+9)*(t*t)))+(((3*p_idx_y)+9)*t))+p_idx_x)+1) */
			_idx30=(_idx29+1);
			/* _idx31 = (((((((((p_idx_z+1)*x_max)+(((3*p_idx_z)+3)*t))*y_max)+((((((3*p_idx_z)+3)*t)+p_idx_y)+3)*x_max))+(((9*p_idx_z)+9)*(t*t)))+(((3*p_idx_y)+9)*t))+p_idx_x)+2) */
			_idx31=(_idx29+2);
			/* _idx32 = (((((((((p_idx_z+1)*x_max)+(((3*p_idx_z)+3)*t))*y_max)+((((((3*p_idx_z)+3)*t)+p_idx_y)+3)*x_max))+(((9*p_idx_z)+9)*(t*t)))+(((3*p_idx_y)+9)*t))+p_idx_x)+3) */
			_idx32=(_idx29+3);
			/* _idx33 = ((((((((p_idx_z+2)*x_max)+(((3*p_idx_z)+6)*t))*y_max)+(((((3*p_idx_z)+6)*t)+p_idx_y)*x_max))+(((9*p_idx_z)+18)*(t*t)))+((3*p_idx_y)*t))+p_idx_x) */
			_idx33=(((_idx17+((x_max+(3*t))*y_max))+((3*t)*x_max))+(9*(t*t)));
			/* _idx34 = (((((((((p_idx_z+2)*x_max)+(((3*p_idx_z)+6)*t))*y_max)+(((((3*p_idx_z)+6)*t)+p_idx_y)*x_max))+(((9*p_idx_z)+18)*(t*t)))+((3*p_idx_y)*t))+p_idx_x)+1) */
			_idx34=(_idx33+1);
			/* _idx35 = (((((((((p_idx_z+2)*x_max)+(((3*p_idx_z)+6)*t))*y_max)+(((((3*p_idx_z)+6)*t)+p_idx_y)*x_max))+(((9*p_idx_z)+18)*(t*t)))+((3*p_idx_y)*t))+p_idx_x)+2) */
			_idx35=(_idx34+1);
			/* _idx36 = (((((((((p_idx_z+2)*x_max)+(((3*p_idx_z)+6)*t))*y_max)+(((((3*p_idx_z)+6)*t)+p_idx_y)*x_max))+(((9*p_idx_z)+18)*(t*t)))+((3*p_idx_y)*t))+p_idx_x)+3) */
			_idx36=(_idx34+2);
			/* _idx37 = ((((((((p_idx_z+2)*x_max)+(((3*p_idx_z)+6)*t))*y_max)+((((((3*p_idx_z)+6)*t)+p_idx_y)+1)*x_max))+(((9*p_idx_z)+18)*(t*t)))+(((3*p_idx_y)+3)*t))+p_idx_x) */
			_idx37=((_idx33+x_max)+(3*t));
			/* _idx38 = (((((((((p_idx_z+2)*x_max)+(((3*p_idx_z)+6)*t))*y_max)+((((((3*p_idx_z)+6)*t)+p_idx_y)+1)*x_max))+(((9*p_idx_z)+18)*(t*t)))+(((3*p_idx_y)+3)*t))+p_idx_x)+1) */
			_idx38=(_idx37+1);
			/* _idx39 = (((((((((p_idx_z+2)*x_max)+(((3*p_idx_z)+6)*t))*y_max)+((((((3*p_idx_z)+6)*t)+p_idx_y)+1)*x_max))+(((9*p_idx_z)+18)*(t*t)))+(((3*p_idx_y)+3)*t))+p_idx_x)+2) */
			_idx39=(_idx37+2);
			/* _idx40 = (((((((((p_idx_z+2)*x_max)+(((3*p_idx_z)+6)*t))*y_max)+((((((3*p_idx_z)+6)*t)+p_idx_y)+1)*x_max))+(((9*p_idx_z)+18)*(t*t)))+(((3*p_idx_y)+3)*t))+p_idx_x)+3) */
			_idx40=(_idx39+1);
			/* _idx41 = ((((((((p_idx_z+2)*x_max)+(((3*p_idx_z)+6)*t))*y_max)+((((((3*p_idx_z)+6)*t)+p_idx_y)+2)*x_max))+(((9*p_idx_z)+18)*(t*t)))+(((3*p_idx_y)+6)*t))+p_idx_x) */
			_idx41=((_idx37+x_max)+(3*t));
			/* _idx42 = (((((((((p_idx_z+2)*x_max)+(((3*p_idx_z)+6)*t))*y_max)+((((((3*p_idx_z)+6)*t)+p_idx_y)+2)*x_max))+(((9*p_idx_z)+18)*(t*t)))+(((3*p_idx_y)+6)*t))+p_idx_x)+1) */
			_idx42=(_idx41+1);
			/* _idx43 = (((((((((p_idx_z+2)*x_max)+(((3*p_idx_z)+6)*t))*y_max)+((((((3*p_idx_z)+6)*t)+p_idx_y)+2)*x_max))+(((9*p_idx_z)+18)*(t*t)))+(((3*p_idx_y)+6)*t))+p_idx_x)+2) */
			_idx43=(_idx41+2);
			/* _idx44 = (((((((((p_idx_z+2)*x_max)+(((3*p_idx_z)+6)*t))*y_max)+((((((3*p_idx_z)+6)*t)+p_idx_y)+2)*x_max))+(((9*p_idx_z)+18)*(t*t)))+(((3*p_idx_y)+6)*t))+p_idx_x)+3) */
			_idx44=(_idx41+3);
			/* _idx45 = ((((((((p_idx_z+2)*x_max)+(((3*p_idx_z)+6)*t))*y_max)+((((((3*p_idx_z)+6)*t)+p_idx_y)+3)*x_max))+(((9*p_idx_z)+18)*(t*t)))+(((3*p_idx_y)+9)*t))+p_idx_x) */
			_idx45=((_idx41+x_max)+(3*t));
			/* _idx46 = (((((((((p_idx_z+2)*x_max)+(((3*p_idx_z)+6)*t))*y_max)+((((((3*p_idx_z)+6)*t)+p_idx_y)+3)*x_max))+(((9*p_idx_z)+18)*(t*t)))+(((3*p_idx_y)+9)*t))+p_idx_x)+1) */
			_idx46=(_idx45+1);
			/* _idx47 = (((((((((p_idx_z+2)*x_max)+(((3*p_idx_z)+6)*t))*y_max)+((((((3*p_idx_z)+6)*t)+p_idx_y)+3)*x_max))+(((9*p_idx_z)+18)*(t*t)))+(((3*p_idx_y)+9)*t))+p_idx_x)+2) */
			_idx47=(_idx45+2);
			/* _idx48 = (((((((((p_idx_z+2)*x_max)+(((3*p_idx_z)+6)*t))*y_max)+((((((3*p_idx_z)+6)*t)+p_idx_y)+3)*x_max))+(((9*p_idx_z)+18)*(t*t)))+(((3*p_idx_y)+9)*t))+p_idx_x)+3) */
			_idx48=(_idx45+3);
			/* _idx49 = ((((((((p_idx_z+3)*x_max)+(((3*p_idx_z)+9)*t))*y_max)+(((((3*p_idx_z)+9)*t)+p_idx_y)*x_max))+(((9*p_idx_z)+27)*(t*t)))+((3*p_idx_y)*t))+p_idx_x) */
			_idx49=(((_idx33+((x_max+(3*t))*y_max))+((3*t)*x_max))+(9*(t*t)));
			/* _idx50 = (((((((((p_idx_z+3)*x_max)+(((3*p_idx_z)+9)*t))*y_max)+(((((3*p_idx_z)+9)*t)+p_idx_y)*x_max))+(((9*p_idx_z)+27)*(t*t)))+((3*p_idx_y)*t))+p_idx_x)+1) */
			_idx50=(_idx49+1);
			/* _idx51 = (((((((((p_idx_z+3)*x_max)+(((3*p_idx_z)+9)*t))*y_max)+(((((3*p_idx_z)+9)*t)+p_idx_y)*x_max))+(((9*p_idx_z)+27)*(t*t)))+((3*p_idx_y)*t))+p_idx_x)+2) */
			_idx51=(_idx50+1);
			/* _idx52 = (((((((((p_idx_z+3)*x_max)+(((3*p_idx_z)+9)*t))*y_max)+(((((3*p_idx_z)+9)*t)+p_idx_y)*x_max))+(((9*p_idx_z)+27)*(t*t)))+((3*p_idx_y)*t))+p_idx_x)+3) */
			_idx52=(_idx50+2);
			/* _idx53 = ((((((((p_idx_z+3)*x_max)+(((3*p_idx_z)+9)*t))*y_max)+((((((3*p_idx_z)+9)*t)+p_idx_y)+1)*x_max))+(((9*p_idx_z)+27)*(t*t)))+(((3*p_idx_y)+3)*t))+p_idx_x) */
			_idx53=((_idx49+x_max)+(3*t));
			/* _idx54 = (((((((((p_idx_z+3)*x_max)+(((3*p_idx_z)+9)*t))*y_max)+((((((3*p_idx_z)+9)*t)+p_idx_y)+1)*x_max))+(((9*p_idx_z)+27)*(t*t)))+(((3*p_idx_y)+3)*t))+p_idx_x)+1) */
			_idx54=(_idx53+1);
			/* _idx55 = (((((((((p_idx_z+3)*x_max)+(((3*p_idx_z)+9)*t))*y_max)+((((((3*p_idx_z)+9)*t)+p_idx_y)+1)*x_max))+(((9*p_idx_z)+27)*(t*t)))+(((3*p_idx_y)+3)*t))+p_idx_x)+2) */
			_idx55=(_idx53+2);
			/* _idx56 = (((((((((p_idx_z+3)*x_max)+(((3*p_idx_z)+9)*t))*y_max)+((((((3*p_idx_z)+9)*t)+p_idx_y)+1)*x_max))+(((9*p_idx_z)+27)*(t*t)))+(((3*p_idx_y)+3)*t))+p_idx_x)+3) */
			_idx56=(_idx53+3);
			/* _idx57 = ((((((((p_idx_z+3)*x_max)+(((3*p_idx_z)+9)*t))*y_max)+((((((3*p_idx_z)+9)*t)+p_idx_y)+2)*x_max))+(((9*p_idx_z)+27)*(t*t)))+(((3*p_idx_y)+6)*t))+p_idx_x) */
			_idx57=((_idx53+x_max)+(3*t));
			/* _idx58 = (((((((((p_idx_z+3)*x_max)+(((3*p_idx_z)+9)*t))*y_max)+((((((3*p_idx_z)+9)*t)+p_idx_y)+2)*x_max))+(((9*p_idx_z)+27)*(t*t)))+(((3*p_idx_y)+6)*t))+p_idx_x)+1) */
			_idx58=(_idx57+1);
			/* _idx59 = (((((((((p_idx_z+3)*x_max)+(((3*p_idx_z)+9)*t))*y_max)+((((((3*p_idx_z)+9)*t)+p_idx_y)+2)*x_max))+(((9*p_idx_z)+27)*(t*t)))+(((3*p_idx_y)+6)*t))+p_idx_x)+2) */
			_idx59=(_idx57+2);
			/* _idx60 = (((((((((p_idx_z+3)*x_max)+(((3*p_idx_z)+9)*t))*y_max)+((((((3*p_idx_z)+9)*t)+p_idx_y)+2)*x_max))+(((9*p_idx_z)+27)*(t*t)))+(((3*p_idx_y)+6)*t))+p_idx_x)+3) */
			_idx60=(_idx57+3);
			/* _idx61 = ((((((((p_idx_z+3)*x_max)+(((3*p_idx_z)+9)*t))*y_max)+((((((3*p_idx_z)+9)*t)+p_idx_y)+3)*x_max))+(((9*p_idx_z)+27)*(t*t)))+(((3*p_idx_y)+9)*t))+p_idx_x) */
			_idx61=((_idx57+x_max)+(3*t));
			/* _idx62 = (((((((((p_idx_z+3)*x_max)+(((3*p_idx_z)+9)*t))*y_max)+((((((3*p_idx_z)+9)*t)+p_idx_y)+3)*x_max))+(((9*p_idx_z)+27)*(t*t)))+(((3*p_idx_y)+9)*t))+p_idx_x)+1) */
			_idx62=(_idx61+1);
			/* _idx63 = (((((((((p_idx_z+3)*x_max)+(((3*p_idx_z)+9)*t))*y_max)+((((((3*p_idx_z)+9)*t)+p_idx_y)+3)*x_max))+(((9*p_idx_z)+27)*(t*t)))+(((3*p_idx_y)+9)*t))+p_idx_x)+2) */
			_idx63=(_idx61+2);
			/* _idx64 = (((((((((p_idx_z+3)*x_max)+(((3*p_idx_z)+9)*t))*y_max)+((((((3*p_idx_z)+9)*t)+p_idx_y)+3)*x_max))+(((9*p_idx_z)+27)*(t*t)))+(((3*p_idx_y)+9)*t))+p_idx_x)+3) */
			_idx64=(_idx63+1);
			u__u_0[t][_idx22]=((((((((w1_a*w1_b)*(w1_c*u__u_0[(t-1)][_idx1]))+((w2_a*w1_b)*(w1_c*u__u_0[(t-1)][_idx2])))+(((w3_a*w1_b)*(w1_c*u__u_0[(t-1)][_idx3]))+((w4_a*w1_b)*(w1_c*u__u_0[(t-1)][_idx4]))))+((((w1_a*w2_b)*(w1_c*u__u_0[(t-1)][_idx5]))+((w2_a*w2_b)*(w1_c*u__u_0[(t-1)][_idx6])))+(((w3_a*w2_b)*(w1_c*u__u_0[(t-1)][_idx7]))+((w4_a*w2_b)*(w1_c*u__u_0[(t-1)][_idx8])))))+(((((w1_a*w3_b)*(w1_c*u__u_0[(t-1)][_idx9]))+((w2_a*w3_b)*(w1_c*u__u_0[(t-1)][_idx10])))+(((w3_a*w3_b)*(w1_c*u__u_0[(t-1)][_idx11]))+((w4_a*w3_b)*(w1_c*u__u_0[(t-1)][_idx12]))))+((((w1_a*w4_b)*(w1_c*u__u_0[(t-1)][_idx13]))+((w2_a*w4_b)*(w1_c*u__u_0[(t-1)][_idx14])))+(((w3_a*w4_b)*(w1_c*u__u_0[(t-1)][_idx15]))+((w4_a*w4_b)*(w1_c*u__u_0[(t-1)][_idx16]))))))+((((((w1_a*w1_b)*(w2_c*u__u_0[(t-1)][_idx17]))+((w2_a*w1_b)*(w2_c*u__u_0[(t-1)][_idx18])))+(((w3_a*w1_b)*(w2_c*u__u_0[(t-1)][_idx19]))+((w4_a*w1_b)*(w2_c*u__u_0[(t-1)][_idx20]))))+((((w1_a*w2_b)*(w2_c*u__u_0[(t-1)][_idx21]))+((w2_a*w2_b)*(w2_c*u__u_0[(t-1)][_idx22])))+(((w3_a*w2_b)*(w2_c*u__u_0[(t-1)][_idx23]))+((w4_a*w2_b)*(w2_c*u__u_0[(t-1)][_idx24])))))+(((((w1_a*w3_b)*(w2_c*u__u_0[(t-1)][_idx25]))+((w2_a*w3_b)*(w2_c*u__u_0[(t-1)][_idx26])))+(((w3_a*w3_b)*(w2_c*u__u_0[(t-1)][_idx27]))+((w4_a*w3_b)*(w2_c*u__u_0[(t-1)][_idx28]))))+((((w1_a*w4_b)*(w2_c*u__u_0[(t-1)][_idx29]))+((w2_a*w4_b)*(w2_c*u__u_0[(t-1)][_idx30])))+(((w3_a*w4_b)*(w2_c*u__u_0[(t-1)][_idx31]))+((w4_a*w4_b)*(w2_c*u__u_0[(t-1)][_idx32])))))))+(((((((w1_a*w1_b)*(w3_c*u__u_0[(t-1)][_idx33]))+((w2_a*w1_b)*(w3_c*u__u_0[(t-1)][_idx34])))+(((w3_a*w1_b)*(w3_c*u__u_0[(t-1)][_idx35]))+((w4_a*w1_b)*(w3_c*u__u_0[(t-1)][_idx36]))))+((((w1_a*w2_b)*(w3_c*u__u_0[(t-1)][_idx37]))+((w2_a*w2_b)*(w3_c*u__u_0[(t-1)][_idx38])))+(((w3_a*w2_b)*(w3_c*u__u_0[(t-1)][_idx39]))+((w4_a*w2_b)*(w3_c*u__u_0[(t-1)][_idx40])))))+(((((w1_a*w3_b)*(w3_c*u__u_0[(t-1)][_idx41]))+((w2_a*w3_b)*(w3_c*u__u_0[(t-1)][_idx42])))+(((w3_a*w3_b)*(w3_c*u__u_0[(t-1)][_idx43]))+((w4_a*w3_b)*(w3_c*u__u_0[(t-1)][_idx44]))))+((((w1_a*w4_b)*(w3_c*u__u_0[(t-1)][_idx45]))+((w2_a*w4_b)*(w3_c*u__u_0[(t-1)][_idx46])))+(((w3_a*w4_b)*(w3_c*u__u_0[(t-1)][_idx47]))+((w4_a*w4_b)*(w3_c*u__u_0[(t-1)][_idx48]))))))+((((((w1_a*w1_b)*(w4_c*u__u_0[(t-1)][_idx49]))+((w2_a*w1_b)*(w4_c*u__u_0[(t-1)][_idx50])))+(((w3_a*w1_b)*(w4_c*u__u_0[(t-1)][_idx51]))+((w4_a*w1_b)*(w4_c*u__u_0[(t-1)][_idx52]))))+((((w1_a*w2_b)*(w4_c*u__u_0[(t-1)][_idx53]))+((w2_a*w2_b)*(w4_c*u__u_0[(t-1)][_idx54])))+(((w3_a*w2_b)*(w4_c*u__u_0[(t-1)][_idx55]))+((w4_a*w2_b)*(w4_c*u__u_0[(t-1)][_idx56])))))+(((((w1_a*w3_b)*(w4_c*u__u_0[(t-1)][_idx57]))+((w2_a*w3_b)*(w4_c*u__u_0[(t-1)][_idx58])))+(((w3_a*w3_b)*(w4_c*u__u_0[(t-1)][_idx59]))+((w4_a*w3_b)*(w4_c*u__u_0[(t-1)][_idx60]))))+((((w1_a*w4_b)*(w4_c*u__u_0[(t-1)][_idx61]))+((w2_a*w4_b)*(w4_c*u__u_0[(t-1)][_idx62])))+(((w3_a*w4_b)*(w4_c*u__u_0[(t-1)][_idx63]))+((w4_a*w4_b)*(w4_c*u__u_0[(t-1)][_idx64]))))))));
		}
#pragma omp barrier 
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
	double *  __restrict__ const u__a_1[16] =  { a_1_0 } ;
	double *  __restrict__ const u__b_2[16] =  { b_2_0 } ;
	double *  __restrict__ const u__c_3[16] =  { c_3_0 } ;
	double *  __restrict__ const u__u_0[16] =  { u_0_0, u_0_1 } ;
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
			/* _idx0 = ((((p_idx_z*x_max)*y_max)+(p_idx_y*x_max))+p_idx_x) */
			_idx0=((((p_idx_z*x_max)*y_max)+(p_idx_y*x_max))+p_idx_x);
			u__a_1[(t-1)][_idx0]=0.2;
			w1_a=0.1;
			u__a_1[(t-1)][_idx0]=0.2;
			w2_a=0.1;
			u__a_1[(t-1)][_idx0]=0.2;
			w3_a=0.1;
			u__a_1[(t-1)][_idx0]=0.2;
			w4_a=0.1;
			u__b_2[(t-1)][_idx0]=0.30000000000000004;
			w1_b=0.1;
			u__b_2[(t-1)][_idx0]=0.30000000000000004;
			w2_b=0.1;
			u__b_2[(t-1)][_idx0]=0.30000000000000004;
			w3_b=0.1;
			u__b_2[(t-1)][_idx0]=0.30000000000000004;
			w4_b=0.1;
			u__c_3[(t-1)][_idx0]=0.4;
			w1_c=0.1;
			u__c_3[(t-1)][_idx0]=0.4;
			w2_c=0.1;
			u__c_3[(t-1)][_idx0]=0.4;
			w3_c=0.1;
			u__c_3[(t-1)][_idx0]=0.4;
			w4_c=0.1;
			/* _idx1 = (((((((p_idx_z*x_max)+((3*p_idx_z)*t))*y_max)+((((3*p_idx_z)*t)+p_idx_y)*x_max))+((9*p_idx_z)*(t*t)))+((3*p_idx_y)*t))+p_idx_x) */
			_idx1=((((_idx0+(((3*p_idx_z)*t)*y_max))+(((3*p_idx_z)*t)*x_max))+((9*p_idx_z)*(t*t)))+((3*p_idx_y)*t));
			u__u_0[(t-1)][_idx1]=0.1;
			/* _idx2 = ((((((((p_idx_z+1)*x_max)+(((3*p_idx_z)+3)*t))*y_max)+(((((3*p_idx_z)+3)*t)+p_idx_y)*x_max))+(((9*p_idx_z)+9)*(t*t)))+((3*p_idx_y)*t))+p_idx_x) */
			_idx2=(((_idx1+((x_max+(3*t))*y_max))+((3*t)*x_max))+(9*(t*t)));
			u__u_0[(t-1)][_idx2]=0.1;
			/* _idx3 = ((((((((p_idx_z+2)*x_max)+(((3*p_idx_z)+6)*t))*y_max)+(((((3*p_idx_z)+6)*t)+p_idx_y)*x_max))+(((9*p_idx_z)+18)*(t*t)))+((3*p_idx_y)*t))+p_idx_x) */
			_idx3=(((_idx2+((x_max+(3*t))*y_max))+((3*t)*x_max))+(9*(t*t)));
			u__u_0[(t-1)][_idx3]=0.1;
			/* _idx4 = ((((((((p_idx_z+3)*x_max)+(((3*p_idx_z)+9)*t))*y_max)+(((((3*p_idx_z)+9)*t)+p_idx_y)*x_max))+(((9*p_idx_z)+27)*(t*t)))+((3*p_idx_y)*t))+p_idx_x) */
			_idx4=(((_idx3+((x_max+(3*t))*y_max))+((3*t)*x_max))+(9*(t*t)));
			u__u_0[(t-1)][_idx4]=0.1;
			/* _idx5 = (((((((p_idx_z*x_max)+((3*p_idx_z)*t))*y_max)+(((((3*p_idx_z)*t)+p_idx_y)+1)*x_max))+((9*p_idx_z)*(t*t)))+(((3*p_idx_y)+3)*t))+p_idx_x) */
			_idx5=((_idx1+x_max)+(3*t));
			u__u_0[(t-1)][_idx5]=0.1;
			/* _idx6 = ((((((((p_idx_z+1)*x_max)+(((3*p_idx_z)+3)*t))*y_max)+((((((3*p_idx_z)+3)*t)+p_idx_y)+1)*x_max))+(((9*p_idx_z)+9)*(t*t)))+(((3*p_idx_y)+3)*t))+p_idx_x) */
			_idx6=((_idx2+x_max)+(3*t));
			u__u_0[(t-1)][_idx6]=0.1;
			/* _idx7 = ((((((((p_idx_z+2)*x_max)+(((3*p_idx_z)+6)*t))*y_max)+((((((3*p_idx_z)+6)*t)+p_idx_y)+1)*x_max))+(((9*p_idx_z)+18)*(t*t)))+(((3*p_idx_y)+3)*t))+p_idx_x) */
			_idx7=((_idx3+x_max)+(3*t));
			u__u_0[(t-1)][_idx7]=0.1;
			/* _idx8 = ((((((((p_idx_z+3)*x_max)+(((3*p_idx_z)+9)*t))*y_max)+((((((3*p_idx_z)+9)*t)+p_idx_y)+1)*x_max))+(((9*p_idx_z)+27)*(t*t)))+(((3*p_idx_y)+3)*t))+p_idx_x) */
			_idx8=((_idx4+x_max)+(3*t));
			u__u_0[(t-1)][_idx8]=0.1;
			/* _idx9 = (((((((p_idx_z*x_max)+((3*p_idx_z)*t))*y_max)+(((((3*p_idx_z)*t)+p_idx_y)+2)*x_max))+((9*p_idx_z)*(t*t)))+(((3*p_idx_y)+6)*t))+p_idx_x) */
			_idx9=((_idx5+x_max)+(3*t));
			u__u_0[(t-1)][_idx9]=0.1;
			/* _idx10 = ((((((((p_idx_z+1)*x_max)+(((3*p_idx_z)+3)*t))*y_max)+((((((3*p_idx_z)+3)*t)+p_idx_y)+2)*x_max))+(((9*p_idx_z)+9)*(t*t)))+(((3*p_idx_y)+6)*t))+p_idx_x) */
			_idx10=((_idx6+x_max)+(3*t));
			u__u_0[(t-1)][_idx10]=0.1;
			/* _idx11 = ((((((((p_idx_z+2)*x_max)+(((3*p_idx_z)+6)*t))*y_max)+((((((3*p_idx_z)+6)*t)+p_idx_y)+2)*x_max))+(((9*p_idx_z)+18)*(t*t)))+(((3*p_idx_y)+6)*t))+p_idx_x) */
			_idx11=((_idx7+x_max)+(3*t));
			u__u_0[(t-1)][_idx11]=0.1;
			/* _idx12 = ((((((((p_idx_z+3)*x_max)+(((3*p_idx_z)+9)*t))*y_max)+((((((3*p_idx_z)+9)*t)+p_idx_y)+2)*x_max))+(((9*p_idx_z)+27)*(t*t)))+(((3*p_idx_y)+6)*t))+p_idx_x) */
			_idx12=((_idx8+x_max)+(3*t));
			u__u_0[(t-1)][_idx12]=0.1;
			/* _idx13 = (((((((p_idx_z*x_max)+((3*p_idx_z)*t))*y_max)+(((((3*p_idx_z)*t)+p_idx_y)+3)*x_max))+((9*p_idx_z)*(t*t)))+(((3*p_idx_y)+9)*t))+p_idx_x) */
			_idx13=((_idx9+x_max)+(3*t));
			u__u_0[(t-1)][_idx13]=0.1;
			/* _idx14 = ((((((((p_idx_z+1)*x_max)+(((3*p_idx_z)+3)*t))*y_max)+((((((3*p_idx_z)+3)*t)+p_idx_y)+3)*x_max))+(((9*p_idx_z)+9)*(t*t)))+(((3*p_idx_y)+9)*t))+p_idx_x) */
			_idx14=((_idx10+x_max)+(3*t));
			u__u_0[(t-1)][_idx14]=0.1;
			/* _idx15 = ((((((((p_idx_z+2)*x_max)+(((3*p_idx_z)+6)*t))*y_max)+((((((3*p_idx_z)+6)*t)+p_idx_y)+3)*x_max))+(((9*p_idx_z)+18)*(t*t)))+(((3*p_idx_y)+9)*t))+p_idx_x) */
			_idx15=((_idx11+x_max)+(3*t));
			u__u_0[(t-1)][_idx15]=0.1;
			/* _idx16 = ((((((((p_idx_z+3)*x_max)+(((3*p_idx_z)+9)*t))*y_max)+((((((3*p_idx_z)+9)*t)+p_idx_y)+3)*x_max))+(((9*p_idx_z)+27)*(t*t)))+(((3*p_idx_y)+9)*t))+p_idx_x) */
			_idx16=((_idx12+x_max)+(3*t));
			u__u_0[(t-1)][_idx16]=0.1;
			/* _idx17 = ((((((((p_idx_z*x_max)+((3*p_idx_z)*t))*y_max)+((((3*p_idx_z)*t)+p_idx_y)*x_max))+((9*p_idx_z)*(t*t)))+((3*p_idx_y)*t))+p_idx_x)+1) */
			_idx17=(_idx1+1);
			u__u_0[(t-1)][_idx17]=0.1;
			/* _idx18 = (((((((((p_idx_z+1)*x_max)+(((3*p_idx_z)+3)*t))*y_max)+(((((3*p_idx_z)+3)*t)+p_idx_y)*x_max))+(((9*p_idx_z)+9)*(t*t)))+((3*p_idx_y)*t))+p_idx_x)+1) */
			_idx18=(_idx2+1);
			u__u_0[(t-1)][_idx18]=0.1;
			/* _idx19 = (((((((((p_idx_z+2)*x_max)+(((3*p_idx_z)+6)*t))*y_max)+(((((3*p_idx_z)+6)*t)+p_idx_y)*x_max))+(((9*p_idx_z)+18)*(t*t)))+((3*p_idx_y)*t))+p_idx_x)+1) */
			_idx19=(_idx3+1);
			u__u_0[(t-1)][_idx19]=0.1;
			/* _idx20 = (((((((((p_idx_z+3)*x_max)+(((3*p_idx_z)+9)*t))*y_max)+(((((3*p_idx_z)+9)*t)+p_idx_y)*x_max))+(((9*p_idx_z)+27)*(t*t)))+((3*p_idx_y)*t))+p_idx_x)+1) */
			_idx20=(_idx4+1);
			u__u_0[(t-1)][_idx20]=0.1;
			/* _idx21 = ((((((((p_idx_z*x_max)+((3*p_idx_z)*t))*y_max)+(((((3*p_idx_z)*t)+p_idx_y)+1)*x_max))+((9*p_idx_z)*(t*t)))+(((3*p_idx_y)+3)*t))+p_idx_x)+1) */
			_idx21=(_idx5+1);
			u__u_0[(t-1)][_idx21]=0.1;
			/* _idx22 = (((((((((p_idx_z+1)*x_max)+(((3*p_idx_z)+3)*t))*y_max)+((((((3*p_idx_z)+3)*t)+p_idx_y)+1)*x_max))+(((9*p_idx_z)+9)*(t*t)))+(((3*p_idx_y)+3)*t))+p_idx_x)+1) */
			_idx22=(_idx6+1);
			u__u_0[(t-1)][_idx22]=0.1;
			/* _idx23 = (((((((((p_idx_z+2)*x_max)+(((3*p_idx_z)+6)*t))*y_max)+((((((3*p_idx_z)+6)*t)+p_idx_y)+1)*x_max))+(((9*p_idx_z)+18)*(t*t)))+(((3*p_idx_y)+3)*t))+p_idx_x)+1) */
			_idx23=(_idx7+1);
			u__u_0[(t-1)][_idx23]=0.1;
			/* _idx24 = (((((((((p_idx_z+3)*x_max)+(((3*p_idx_z)+9)*t))*y_max)+((((((3*p_idx_z)+9)*t)+p_idx_y)+1)*x_max))+(((9*p_idx_z)+27)*(t*t)))+(((3*p_idx_y)+3)*t))+p_idx_x)+1) */
			_idx24=(_idx8+1);
			u__u_0[(t-1)][_idx24]=0.1;
			/* _idx25 = ((((((((p_idx_z*x_max)+((3*p_idx_z)*t))*y_max)+(((((3*p_idx_z)*t)+p_idx_y)+2)*x_max))+((9*p_idx_z)*(t*t)))+(((3*p_idx_y)+6)*t))+p_idx_x)+1) */
			_idx25=(_idx9+1);
			u__u_0[(t-1)][_idx25]=0.1;
			/* _idx26 = (((((((((p_idx_z+1)*x_max)+(((3*p_idx_z)+3)*t))*y_max)+((((((3*p_idx_z)+3)*t)+p_idx_y)+2)*x_max))+(((9*p_idx_z)+9)*(t*t)))+(((3*p_idx_y)+6)*t))+p_idx_x)+1) */
			_idx26=(_idx10+1);
			u__u_0[(t-1)][_idx26]=0.1;
			/* _idx27 = (((((((((p_idx_z+2)*x_max)+(((3*p_idx_z)+6)*t))*y_max)+((((((3*p_idx_z)+6)*t)+p_idx_y)+2)*x_max))+(((9*p_idx_z)+18)*(t*t)))+(((3*p_idx_y)+6)*t))+p_idx_x)+1) */
			_idx27=(_idx11+1);
			u__u_0[(t-1)][_idx27]=0.1;
			/* _idx28 = (((((((((p_idx_z+3)*x_max)+(((3*p_idx_z)+9)*t))*y_max)+((((((3*p_idx_z)+9)*t)+p_idx_y)+2)*x_max))+(((9*p_idx_z)+27)*(t*t)))+(((3*p_idx_y)+6)*t))+p_idx_x)+1) */
			_idx28=(_idx12+1);
			u__u_0[(t-1)][_idx28]=0.1;
			/* _idx29 = ((((((((p_idx_z*x_max)+((3*p_idx_z)*t))*y_max)+(((((3*p_idx_z)*t)+p_idx_y)+3)*x_max))+((9*p_idx_z)*(t*t)))+(((3*p_idx_y)+9)*t))+p_idx_x)+1) */
			_idx29=(_idx13+1);
			u__u_0[(t-1)][_idx29]=0.1;
			/* _idx30 = (((((((((p_idx_z+1)*x_max)+(((3*p_idx_z)+3)*t))*y_max)+((((((3*p_idx_z)+3)*t)+p_idx_y)+3)*x_max))+(((9*p_idx_z)+9)*(t*t)))+(((3*p_idx_y)+9)*t))+p_idx_x)+1) */
			_idx30=(_idx14+1);
			u__u_0[(t-1)][_idx30]=0.1;
			/* _idx31 = (((((((((p_idx_z+2)*x_max)+(((3*p_idx_z)+6)*t))*y_max)+((((((3*p_idx_z)+6)*t)+p_idx_y)+3)*x_max))+(((9*p_idx_z)+18)*(t*t)))+(((3*p_idx_y)+9)*t))+p_idx_x)+1) */
			_idx31=(_idx15+1);
			u__u_0[(t-1)][_idx31]=0.1;
			/* _idx32 = (((((((((p_idx_z+3)*x_max)+(((3*p_idx_z)+9)*t))*y_max)+((((((3*p_idx_z)+9)*t)+p_idx_y)+3)*x_max))+(((9*p_idx_z)+27)*(t*t)))+(((3*p_idx_y)+9)*t))+p_idx_x)+1) */
			_idx32=(_idx16+1);
			u__u_0[(t-1)][_idx32]=0.1;
			/* _idx33 = ((((((((p_idx_z*x_max)+((3*p_idx_z)*t))*y_max)+((((3*p_idx_z)*t)+p_idx_y)*x_max))+((9*p_idx_z)*(t*t)))+((3*p_idx_y)*t))+p_idx_x)+2) */
			_idx33=(_idx17+1);
			u__u_0[(t-1)][_idx33]=0.1;
			/* _idx34 = (((((((((p_idx_z+1)*x_max)+(((3*p_idx_z)+3)*t))*y_max)+(((((3*p_idx_z)+3)*t)+p_idx_y)*x_max))+(((9*p_idx_z)+9)*(t*t)))+((3*p_idx_y)*t))+p_idx_x)+2) */
			_idx34=(_idx2+2);
			u__u_0[(t-1)][_idx34]=0.1;
			/* _idx35 = (((((((((p_idx_z+2)*x_max)+(((3*p_idx_z)+6)*t))*y_max)+(((((3*p_idx_z)+6)*t)+p_idx_y)*x_max))+(((9*p_idx_z)+18)*(t*t)))+((3*p_idx_y)*t))+p_idx_x)+2) */
			_idx35=(_idx19+1);
			u__u_0[(t-1)][_idx35]=0.1;
			/* _idx36 = (((((((((p_idx_z+3)*x_max)+(((3*p_idx_z)+9)*t))*y_max)+(((((3*p_idx_z)+9)*t)+p_idx_y)*x_max))+(((9*p_idx_z)+27)*(t*t)))+((3*p_idx_y)*t))+p_idx_x)+2) */
			_idx36=(_idx20+1);
			u__u_0[(t-1)][_idx36]=0.1;
			/* _idx37 = ((((((((p_idx_z*x_max)+((3*p_idx_z)*t))*y_max)+(((((3*p_idx_z)*t)+p_idx_y)+1)*x_max))+((9*p_idx_z)*(t*t)))+(((3*p_idx_y)+3)*t))+p_idx_x)+2) */
			_idx37=(_idx21+1);
			u__u_0[(t-1)][_idx37]=0.1;
			/* _idx38 = (((((((((p_idx_z+1)*x_max)+(((3*p_idx_z)+3)*t))*y_max)+((((((3*p_idx_z)+3)*t)+p_idx_y)+1)*x_max))+(((9*p_idx_z)+9)*(t*t)))+(((3*p_idx_y)+3)*t))+p_idx_x)+2) */
			_idx38=(_idx6+2);
			u__u_0[(t-1)][_idx38]=0.1;
			/* _idx39 = (((((((((p_idx_z+2)*x_max)+(((3*p_idx_z)+6)*t))*y_max)+((((((3*p_idx_z)+6)*t)+p_idx_y)+1)*x_max))+(((9*p_idx_z)+18)*(t*t)))+(((3*p_idx_y)+3)*t))+p_idx_x)+2) */
			_idx39=(_idx7+2);
			u__u_0[(t-1)][_idx39]=0.1;
			/* _idx40 = (((((((((p_idx_z+3)*x_max)+(((3*p_idx_z)+9)*t))*y_max)+((((((3*p_idx_z)+9)*t)+p_idx_y)+1)*x_max))+(((9*p_idx_z)+27)*(t*t)))+(((3*p_idx_y)+3)*t))+p_idx_x)+2) */
			_idx40=(_idx8+2);
			u__u_0[(t-1)][_idx40]=0.1;
			/* _idx41 = ((((((((p_idx_z*x_max)+((3*p_idx_z)*t))*y_max)+(((((3*p_idx_z)*t)+p_idx_y)+2)*x_max))+((9*p_idx_z)*(t*t)))+(((3*p_idx_y)+6)*t))+p_idx_x)+2) */
			_idx41=(_idx25+1);
			u__u_0[(t-1)][_idx41]=0.1;
			/* _idx42 = (((((((((p_idx_z+1)*x_max)+(((3*p_idx_z)+3)*t))*y_max)+((((((3*p_idx_z)+3)*t)+p_idx_y)+2)*x_max))+(((9*p_idx_z)+9)*(t*t)))+(((3*p_idx_y)+6)*t))+p_idx_x)+2) */
			_idx42=(_idx10+2);
			u__u_0[(t-1)][_idx42]=0.1;
			/* _idx43 = (((((((((p_idx_z+2)*x_max)+(((3*p_idx_z)+6)*t))*y_max)+((((((3*p_idx_z)+6)*t)+p_idx_y)+2)*x_max))+(((9*p_idx_z)+18)*(t*t)))+(((3*p_idx_y)+6)*t))+p_idx_x)+2) */
			_idx43=(_idx27+1);
			u__u_0[(t-1)][_idx43]=0.1;
			/* _idx44 = (((((((((p_idx_z+3)*x_max)+(((3*p_idx_z)+9)*t))*y_max)+((((((3*p_idx_z)+9)*t)+p_idx_y)+2)*x_max))+(((9*p_idx_z)+27)*(t*t)))+(((3*p_idx_y)+6)*t))+p_idx_x)+2) */
			_idx44=(_idx12+2);
			u__u_0[(t-1)][_idx44]=0.1;
			/* _idx45 = ((((((((p_idx_z*x_max)+((3*p_idx_z)*t))*y_max)+(((((3*p_idx_z)*t)+p_idx_y)+3)*x_max))+((9*p_idx_z)*(t*t)))+(((3*p_idx_y)+9)*t))+p_idx_x)+2) */
			_idx45=(_idx29+1);
			u__u_0[(t-1)][_idx45]=0.1;
			/* _idx46 = (((((((((p_idx_z+1)*x_max)+(((3*p_idx_z)+3)*t))*y_max)+((((((3*p_idx_z)+3)*t)+p_idx_y)+3)*x_max))+(((9*p_idx_z)+9)*(t*t)))+(((3*p_idx_y)+9)*t))+p_idx_x)+2) */
			_idx46=(_idx30+1);
			u__u_0[(t-1)][_idx46]=0.1;
			/* _idx47 = (((((((((p_idx_z+2)*x_max)+(((3*p_idx_z)+6)*t))*y_max)+((((((3*p_idx_z)+6)*t)+p_idx_y)+3)*x_max))+(((9*p_idx_z)+18)*(t*t)))+(((3*p_idx_y)+9)*t))+p_idx_x)+2) */
			_idx47=(_idx15+2);
			u__u_0[(t-1)][_idx47]=0.1;
			/* _idx48 = (((((((((p_idx_z+3)*x_max)+(((3*p_idx_z)+9)*t))*y_max)+((((((3*p_idx_z)+9)*t)+p_idx_y)+3)*x_max))+(((9*p_idx_z)+27)*(t*t)))+(((3*p_idx_y)+9)*t))+p_idx_x)+2) */
			_idx48=(_idx16+2);
			u__u_0[(t-1)][_idx48]=0.1;
			/* _idx49 = ((((((((p_idx_z*x_max)+((3*p_idx_z)*t))*y_max)+((((3*p_idx_z)*t)+p_idx_y)*x_max))+((9*p_idx_z)*(t*t)))+((3*p_idx_y)*t))+p_idx_x)+3) */
			_idx49=(_idx17+2);
			u__u_0[(t-1)][_idx49]=0.1;
			/* _idx50 = (((((((((p_idx_z+1)*x_max)+(((3*p_idx_z)+3)*t))*y_max)+(((((3*p_idx_z)+3)*t)+p_idx_y)*x_max))+(((9*p_idx_z)+9)*(t*t)))+((3*p_idx_y)*t))+p_idx_x)+3) */
			_idx50=(_idx2+3);
			u__u_0[(t-1)][_idx50]=0.1;
			/* _idx51 = (((((((((p_idx_z+2)*x_max)+(((3*p_idx_z)+6)*t))*y_max)+(((((3*p_idx_z)+6)*t)+p_idx_y)*x_max))+(((9*p_idx_z)+18)*(t*t)))+((3*p_idx_y)*t))+p_idx_x)+3) */
			_idx51=(_idx19+2);
			u__u_0[(t-1)][_idx51]=0.1;
			/* _idx52 = (((((((((p_idx_z+3)*x_max)+(((3*p_idx_z)+9)*t))*y_max)+(((((3*p_idx_z)+9)*t)+p_idx_y)*x_max))+(((9*p_idx_z)+27)*(t*t)))+((3*p_idx_y)*t))+p_idx_x)+3) */
			_idx52=(_idx20+2);
			u__u_0[(t-1)][_idx52]=0.1;
			/* _idx53 = ((((((((p_idx_z*x_max)+((3*p_idx_z)*t))*y_max)+(((((3*p_idx_z)*t)+p_idx_y)+1)*x_max))+((9*p_idx_z)*(t*t)))+(((3*p_idx_y)+3)*t))+p_idx_x)+3) */
			_idx53=(_idx21+2);
			u__u_0[(t-1)][_idx53]=0.1;
			/* _idx54 = (((((((((p_idx_z+1)*x_max)+(((3*p_idx_z)+3)*t))*y_max)+((((((3*p_idx_z)+3)*t)+p_idx_y)+1)*x_max))+(((9*p_idx_z)+9)*(t*t)))+(((3*p_idx_y)+3)*t))+p_idx_x)+3) */
			_idx54=(_idx6+3);
			u__u_0[(t-1)][_idx54]=0.1;
			/* _idx55 = (((((((((p_idx_z+2)*x_max)+(((3*p_idx_z)+6)*t))*y_max)+((((((3*p_idx_z)+6)*t)+p_idx_y)+1)*x_max))+(((9*p_idx_z)+18)*(t*t)))+(((3*p_idx_y)+3)*t))+p_idx_x)+3) */
			_idx55=(_idx7+3);
			u__u_0[(t-1)][_idx55]=0.1;
			/* _idx56 = (((((((((p_idx_z+3)*x_max)+(((3*p_idx_z)+9)*t))*y_max)+((((((3*p_idx_z)+9)*t)+p_idx_y)+1)*x_max))+(((9*p_idx_z)+27)*(t*t)))+(((3*p_idx_y)+3)*t))+p_idx_x)+3) */
			_idx56=(_idx8+3);
			u__u_0[(t-1)][_idx56]=0.1;
			/* _idx57 = ((((((((p_idx_z*x_max)+((3*p_idx_z)*t))*y_max)+(((((3*p_idx_z)*t)+p_idx_y)+2)*x_max))+((9*p_idx_z)*(t*t)))+(((3*p_idx_y)+6)*t))+p_idx_x)+3) */
			_idx57=(_idx25+2);
			u__u_0[(t-1)][_idx57]=0.1;
			/* _idx58 = (((((((((p_idx_z+1)*x_max)+(((3*p_idx_z)+3)*t))*y_max)+((((((3*p_idx_z)+3)*t)+p_idx_y)+2)*x_max))+(((9*p_idx_z)+9)*(t*t)))+(((3*p_idx_y)+6)*t))+p_idx_x)+3) */
			_idx58=(_idx42+1);
			u__u_0[(t-1)][_idx58]=0.1;
			/* _idx59 = (((((((((p_idx_z+2)*x_max)+(((3*p_idx_z)+6)*t))*y_max)+((((((3*p_idx_z)+6)*t)+p_idx_y)+2)*x_max))+(((9*p_idx_z)+18)*(t*t)))+(((3*p_idx_y)+6)*t))+p_idx_x)+3) */
			_idx59=(_idx43+1);
			u__u_0[(t-1)][_idx59]=0.1;
			/* _idx60 = (((((((((p_idx_z+3)*x_max)+(((3*p_idx_z)+9)*t))*y_max)+((((((3*p_idx_z)+9)*t)+p_idx_y)+2)*x_max))+(((9*p_idx_z)+27)*(t*t)))+(((3*p_idx_y)+6)*t))+p_idx_x)+3) */
			_idx60=(_idx12+3);
			u__u_0[(t-1)][_idx60]=0.1;
			/* _idx61 = ((((((((p_idx_z*x_max)+((3*p_idx_z)*t))*y_max)+(((((3*p_idx_z)*t)+p_idx_y)+3)*x_max))+((9*p_idx_z)*(t*t)))+(((3*p_idx_y)+9)*t))+p_idx_x)+3) */
			_idx61=(_idx29+2);
			u__u_0[(t-1)][_idx61]=0.1;
			/* _idx62 = (((((((((p_idx_z+1)*x_max)+(((3*p_idx_z)+3)*t))*y_max)+((((((3*p_idx_z)+3)*t)+p_idx_y)+3)*x_max))+(((9*p_idx_z)+9)*(t*t)))+(((3*p_idx_y)+9)*t))+p_idx_x)+3) */
			_idx62=(_idx46+1);
			u__u_0[(t-1)][_idx62]=0.1;
			/* _idx63 = (((((((((p_idx_z+2)*x_max)+(((3*p_idx_z)+6)*t))*y_max)+((((((3*p_idx_z)+6)*t)+p_idx_y)+3)*x_max))+(((9*p_idx_z)+18)*(t*t)))+(((3*p_idx_y)+9)*t))+p_idx_x)+3) */
			_idx63=(_idx15+3);
			u__u_0[(t-1)][_idx63]=0.1;
			/* _idx64 = (((((((((p_idx_z+3)*x_max)+(((3*p_idx_z)+9)*t))*y_max)+((((((3*p_idx_z)+9)*t)+p_idx_y)+3)*x_max))+(((9*p_idx_z)+27)*(t*t)))+(((3*p_idx_y)+9)*t))+p_idx_x)+3) */
			_idx64=(_idx48+1);
			u__u_0[(t-1)][_idx64]=0.1;
			u__u_0[t][_idx22]=1.1;
		}
#pragma omp barrier 
	}
}

