#define t_max 1
#define t 1 

/*
(T[0][0][0][1][0]=((((T[0][0][0][0][0]*((c[0][0][0][0][1]*T[0][0][0][0][0])+c[0][0][0][0][2]))+c[0][0][0][0][3])+((c[0][0][0][0][4]*T[-1][0][0][0][0])+(c[0][0][0][0][5]*T[1][0][0][0][0])))+(((c[0][0][0][0][6]*T[0][-1][0][0][0])+(c[0][0][0][0][7]*T[0][1][0][0][0]))+((c[0][0][0][0][8]*T[0][0][-1][0][0])+(c[0][0][0][0][9]*T[0][0][1][0][0])))))

*/
__global__ void hyperthermia(float *  *  T_0_1_out, float *  T_0_0, float *  T_0_1, float *  c_1_0, float *  c_2_0, float *  c_3_0, float *  c_4_0, float *  c_5_0, float *  c_6_0, float *  c_7_0, float *  c_8_0, float *  c_9_0, int x_max, int y_max, int z_max, int tbx, int tby, int tbz, int c)
{
/*
	const float *  const u__c_1[16] =  { c_1_0 } ;
	const float *  const u__c_2[16] =  { c_2_0 } ;
	const float *  const u__c_3[16] =  { c_3_0 } ;
	const float *  const u__c_4[16] =  { c_4_0 } ;
	const float *  const u__c_5[16] =  { c_5_0 } ;
	const float *  const u__c_6[16] =  { c_6_0 } ;
	const float *  const u__c_7[16] =  { c_7_0 } ;
	const float *  const u__c_8[16] =  { c_8_0 } ;
	const float *  const u__c_9[16] =  { c_9_0 } ;
	float *  const u__T_0[16] =  { T_0_0, T_0_1 } ;
*/
	int _idx0;
	int _idx1;
	int _idx2;
	int _idx3;
	int _idx4;
	int _idx5;
	int _idx6;
	int _idx7;
	int chunk_idx_x;
	int chunk_idx_x_max;
	int chunk_idx_y;
	int chunk_idx_y_max;
	int chunk_idx_z;
	int chunk_idx_z_max;
	int idx_1_2;
	int size_1_1;
	int size_1_2;
	//int t;
	int thd_idx_x;
	int thd_idx_y;
	int thd_idx_z;
	int thdblks_idx_x;
	int thdblks_idx_x_max;
	int thdblks_idx_y;
	int thdblks_idx_y_max;
	int thdblks_idx_z;
	int thdblks_idx_z_max;
	int tmp;
	/*
	Initializations
	*/
	size_1_1=(y_max/blockDim.y);
	size_1_2=(z_max/blockDim.z);
	idx_1_2=(blockIdx.y/size_1_2);
	tmp=(blockIdx.y-(idx_1_2*size_1_2));
	chunk_idx_x=(c*(threadIdx.x+(blockDim.x*blockIdx.x)));
	chunk_idx_x_max=(chunk_idx_x+c);
	chunk_idx_y=(threadIdx.y+(tmp*blockDim.y));
	chunk_idx_y_max=(chunk_idx_y+1);
	chunk_idx_z=(threadIdx.z+(idx_1_2*blockDim.z));
	chunk_idx_z_max=(chunk_idx_z+1);
	thdblks_idx_x=(tbx*(threadIdx.x+(blockDim.x*blockIdx.x)));
	thdblks_idx_x_max=(thdblks_idx_x+tbx);
	thdblks_idx_y=(tby*(threadIdx.y+(tmp*blockDim.y)));
	thdblks_idx_y_max=(thdblks_idx_y+tby);
	thdblks_idx_z=(tbz*(threadIdx.z+(idx_1_2*blockDim.z)));
	thdblks_idx_z_max=(thdblks_idx_z+tbz);
	/*
	Implementation
	*/
	/*
	for t = 1..t_max by 1 parallel 1 <level 0> schedule  { ... }
	*/
	//for (t=1; t<=t_max; t+=1)
	{
		/* Index bounds calculations for iterators in thdblks[t=t, s=(tbx, tby, tbz)][0] */
		/* Index bounds calculations for iterators in chunk[t=t, s=(c, 1, 1)][0] */
		/*
		for POINT thd[t=t, s=(1, 1, 1)][0] of size [1, 1, 1] in chunk[t=t, s=(:, :, :)][0] parallel 1 <level 2> schedule default { ... }
		*/
		{
			/* Index bounds calculations for iterators in thd[t=t, s=(1, 1, 1)][0] */
			thd_idx_z=chunk_idx_z;
			thd_idx_y=chunk_idx_y;
			for (thd_idx_x=chunk_idx_x; thd_idx_x<(chunk_idx_x_max-0); thd_idx_x+=1)
			{
				/* Index bounds calculations for iterators in thd[t=t, s=(1, 1, 1)][0] */
				/*
				u[t=(t+1), s=thd[t=?, s=?][0]][0]=stencil(u[t=t, s=thd[t=?, s=?][0]][0])
				*/
				/* _idx0 = ((((((((((((thd_idx_z+1)*x_max)+((2*t)*thd_idx_z))+(2*t))*y_max)+((((((2*t)*thd_idx_z)+thd_idx_y)+(2*t))+1)*x_max))+((4*(t*t))*thd_idx_z))+((2*t)*thd_idx_y))+thd_idx_x)+(4*(t*t)))+(2*t))+1) */
				_idx0=((((((((((((thd_idx_z+1)*x_max)+((2*t)*thd_idx_z))+(2*t))*y_max)+((((((2*t)*thd_idx_z)+thd_idx_y)+(2*t))+1)*x_max))+((4*(t*t))*thd_idx_z))+((2*t)*thd_idx_y))+thd_idx_x)+(4*(t*t)))+(2*t))+1);
				/* _idx1 = ((((thd_idx_z*x_max)*y_max)+(thd_idx_y*x_max))+thd_idx_x) */
				_idx1=(((((((_idx0+(((( - x_max)-((2*t)*thd_idx_z))-(2*t))*y_max))+(((((-2*t)*thd_idx_z)-(2*t))-1)*x_max))-((4*(t*t))*thd_idx_z))-((2*t)*thd_idx_y))-(4*(t*t)))-(2*t))-1);
				/* _idx2 = (((((((((((thd_idx_z+1)*x_max)+((2*t)*thd_idx_z))+(2*t))*y_max)+((((((2*t)*thd_idx_z)+thd_idx_y)+(2*t))+1)*x_max))+((4*(t*t))*thd_idx_z))+((2*t)*thd_idx_y))+thd_idx_x)+(4*(t*t)))+(2*t)) */
				_idx2=(_idx0-1);
				/* _idx3 = ((((((((((((thd_idx_z+1)*x_max)+((2*t)*thd_idx_z))+(2*t))*y_max)+((((((2*t)*thd_idx_z)+thd_idx_y)+(2*t))+1)*x_max))+((4*(t*t))*thd_idx_z))+((2*t)*thd_idx_y))+thd_idx_x)+(4*(t*t)))+(2*t))+2) */
				_idx3=(_idx2+2);
				/* _idx4 = (((((((((((thd_idx_z+1)*x_max)+((2*t)*thd_idx_z))+(2*t))*y_max)+(((((2*t)*thd_idx_z)+thd_idx_y)+(2*t))*x_max))+((4*(t*t))*thd_idx_z))+((2*t)*thd_idx_y))+thd_idx_x)+(4*(t*t)))+1) */
				_idx4=((_idx0-x_max)-(2*t));
				/* _idx5 = ((((((((((((thd_idx_z+1)*x_max)+((2*t)*thd_idx_z))+(2*t))*y_max)+((((((2*t)*thd_idx_z)+thd_idx_y)+(2*t))+2)*x_max))+((4*(t*t))*thd_idx_z))+((2*t)*thd_idx_y))+thd_idx_x)+(4*(t*t)))+(4*t))+1) */
				_idx5=((_idx0+x_max)+(2*t));
				/* _idx6 = (((((((((thd_idx_z*x_max)+((2*t)*thd_idx_z))*y_max)+(((((2*t)*thd_idx_z)+thd_idx_y)+1)*x_max))+((4*(t*t))*thd_idx_z))+((2*t)*thd_idx_y))+thd_idx_x)+(2*t))+1) */
				_idx6=((((_idx2+((( - x_max)-(2*t))*y_max))-((2*t)*x_max))-(4*(t*t)))+1);
				/* _idx7 = ((((((((((((thd_idx_z+2)*x_max)+((2*t)*thd_idx_z))+(4*t))*y_max)+((((((2*t)*thd_idx_z)+thd_idx_y)+(4*t))+1)*x_max))+((4*(t*t))*thd_idx_z))+((2*t)*thd_idx_y))+thd_idx_x)+(8*(t*t)))+(2*t))+1) */
				_idx7=(((_idx0+((x_max+(2*t))*y_max))+((2*t)*x_max))+(4*(t*t)));
//				u__T_0[t][_idx0]=((((u__T_0[(t-1)][_idx0]*((u__c_1[(t-1)][_idx1]*u__T_0[(t-1)][_idx0])+u__c_2[(t-1)][_idx1]))+u__c_3[(t-1)][_idx1])+((u__c_4[(t-1)][_idx1]*u__T_0[(t-1)][_idx2])+(u__c_5[(t-1)][_idx1]*u__T_0[(t-1)][_idx3])))+(((u__c_6[(t-1)][_idx1]*u__T_0[(t-1)][_idx4])+(u__c_7[(t-1)][_idx1]*u__T_0[(t-1)][_idx5]))+((u__c_8[(t-1)][_idx1]*u__T_0[(t-1)][_idx6])+(u__c_9[(t-1)][_idx1]*u__T_0[(t-1)][_idx7]))));
				T_0_1[_idx0]=((((T_0_0[_idx0]*((c_1_0[_idx1]*T_0_0[_idx0])+c_2_0[_idx1]))+c_3_0[_idx1])+((c_4_0[_idx1]*T_0_0[_idx2])+(c_5_0[_idx1]*T_0_0[_idx3])))+(((c_6_0[_idx1]*T_0_0[_idx4])+(c_7_0[_idx1]*T_0_0[_idx5]))+((c_8_0[_idx1]*T_0_0[_idx6])+(c_9_0[_idx1]*T_0_0[_idx7]))));
			}
		}
	}
}

__global__ void initialize(float *  T_0_0, float *  T_0_1, float *  c_1_0, float *  c_2_0, float *  c_3_0, float *  c_4_0, float *  c_5_0, float *  c_6_0, float *  c_7_0, float *  c_8_0, float *  c_9_0, int x_max, int y_max, int z_max, int tbx, int tby, int tbz, int c)
{
	 float *  const u__c_1[16] =  { c_1_0 } ;
	 float *  const u__c_2[16] =  { c_2_0 } ;
	 float *  const u__c_3[16] =  { c_3_0 } ;
	 float *  const u__c_4[16] =  { c_4_0 } ;
	 float *  const u__c_5[16] =  { c_5_0 } ;
	 float *  const u__c_6[16] =  { c_6_0 } ;
	 float *  const u__c_7[16] =  { c_7_0 } ;
	 float *  const u__c_8[16] =  { c_8_0 } ;
	 float *  const u__c_9[16] =  { c_9_0 } ;
	float *  const u__T_0[16] =  { T_0_0, T_0_1 } ;
	int _idx0;
	int _idx1;
	int _idx2;
	int _idx3;
	int _idx4;
	int _idx5;
	int _idx6;
	int _idx7;
	int chunk_idx_x;
	int chunk_idx_x_max;
	int chunk_idx_y;
	int chunk_idx_y_max;
	int chunk_idx_z;
	int chunk_idx_z_max;
	int idx_1_2;
	int size_1_1;
	int size_1_2;
	//int t;
	int thd_idx_x;
	int thd_idx_y;
	int thd_idx_z;
	int thdblks_idx_x;
	int thdblks_idx_x_max;
	int thdblks_idx_y;
	int thdblks_idx_y_max;
	int thdblks_idx_z;
	int thdblks_idx_z_max;
	int tmp;
	/*
	Initializations
	*/
	size_1_1=(y_max/blockDim.y);
	size_1_2=(z_max/blockDim.z);
	idx_1_2=(blockIdx.y/size_1_2);
	tmp=(blockIdx.y-(idx_1_2*size_1_2));
	chunk_idx_x=(c*(threadIdx.x+(blockDim.x*blockIdx.x)));
	chunk_idx_x_max=(chunk_idx_x+c);
	chunk_idx_y=(threadIdx.y+(tmp*blockDim.y));
	chunk_idx_y_max=(chunk_idx_y+1);
	chunk_idx_z=(threadIdx.z+(idx_1_2*blockDim.z));
	chunk_idx_z_max=(chunk_idx_z+1);
	thdblks_idx_x=(tbx*(threadIdx.x+(blockDim.x*blockIdx.x)));
	thdblks_idx_x_max=(thdblks_idx_x+tbx);
	thdblks_idx_y=(tby*(threadIdx.y+(tmp*blockDim.y)));
	thdblks_idx_y_max=(thdblks_idx_y+tby);
	thdblks_idx_z=(tbz*(threadIdx.z+(idx_1_2*blockDim.z)));
	thdblks_idx_z_max=(thdblks_idx_z+tbz);
	/*
	Implementation
	*/
	/*
	for t = 1..t_max by 1 parallel 1 <level 0> schedule  { ... }
	*/
	//for (t=1; t<=t_max; t+=1)
	{
		/* Index bounds calculations for iterators in thdblks[t=t, s=(tbx, tby, tbz)][0] */
		/* Index bounds calculations for iterators in chunk[t=t, s=(c, 1, 1)][0] */
		/*
		for POINT thd[t=t, s=(1, 1, 1)][0] of size [1, 1, 1] in chunk[t=t, s=(:, :, :)][0] parallel 1 <level 2> schedule default { ... }
		*/
		{
			/* Index bounds calculations for iterators in thd[t=t, s=(1, 1, 1)][0] */
			thd_idx_z=chunk_idx_z;
			thd_idx_y=chunk_idx_y;
			for (thd_idx_x=chunk_idx_x; thd_idx_x<(chunk_idx_x_max-0); thd_idx_x+=1)
			{
				/* Index bounds calculations for iterators in thd[t=t, s=(1, 1, 1)][0] */
				/*
				u[t=(t+1), s=thd[t=?, s=?][0]][0]=stencil(u[t=t, s=thd[t=?, s=?][0]][0])
				*/
				/* _idx0 = (((((((((((thd_idx_z+1)*x_max)+((2*t)*thd_idx_z))+(2*t))*y_max)+((((((2*t)*thd_idx_z)+thd_idx_y)+(2*t))+1)*x_max))+((4*(t*t))*thd_idx_z))+((2*t)*thd_idx_y))+thd_idx_x)+(4*(t*t)))+(2*t)) */
				_idx0=(((((((((((thd_idx_z+1)*x_max)+((2*t)*thd_idx_z))+(2*t))*y_max)+((((((2*t)*thd_idx_z)+thd_idx_y)+(2*t))+1)*x_max))+((4*(t*t))*thd_idx_z))+((2*t)*thd_idx_y))+thd_idx_x)+(4*(t*t)))+(2*t));
				u__T_0[(t-1)][_idx0]=0.1;
				/* _idx1 = (((((((((((thd_idx_z+1)*x_max)+((2*t)*thd_idx_z))+(2*t))*y_max)+(((((2*t)*thd_idx_z)+thd_idx_y)+(2*t))*x_max))+((4*(t*t))*thd_idx_z))+((2*t)*thd_idx_y))+thd_idx_x)+(4*(t*t)))+1) */
				_idx1=(((_idx0-x_max)-(2*t))+1);
				u__T_0[(t-1)][_idx1]=0.1;
				/* _idx2 = (((((((((thd_idx_z*x_max)+((2*t)*thd_idx_z))*y_max)+(((((2*t)*thd_idx_z)+thd_idx_y)+1)*x_max))+((4*(t*t))*thd_idx_z))+((2*t)*thd_idx_y))+thd_idx_x)+(2*t))+1) */
				_idx2=((((_idx0+((( - x_max)-(2*t))*y_max))-((2*t)*x_max))-(4*(t*t)))+1);
				u__T_0[(t-1)][_idx2]=0.1;
				/* _idx3 = ((((((((((((thd_idx_z+1)*x_max)+((2*t)*thd_idx_z))+(2*t))*y_max)+((((((2*t)*thd_idx_z)+thd_idx_y)+(2*t))+1)*x_max))+((4*(t*t))*thd_idx_z))+((2*t)*thd_idx_y))+thd_idx_x)+(4*(t*t)))+(2*t))+1) */
				_idx3=(_idx0+1);
				u__T_0[(t-1)][_idx3]=0.1;
				/* _idx4 = ((((thd_idx_z*x_max)*y_max)+(thd_idx_y*x_max))+thd_idx_x) */
				_idx4=((((((_idx2-(((2*t)*thd_idx_z)*y_max))+((((-2*t)*thd_idx_z)-1)*x_max))-((4*(t*t))*thd_idx_z))-((2*t)*thd_idx_y))-(2*t))-1);
				u__c_1[(t-1)][_idx4]=0.2;
				u__c_2[(t-1)][_idx4]=0.30000000000000004;
				u__c_3[(t-1)][_idx4]=0.4;
				u__c_4[(t-1)][_idx4]=0.5;
				u__c_5[(t-1)][_idx4]=0.6000000000000001;
				u__c_6[(t-1)][_idx4]=0.7000000000000001;
				u__c_7[(t-1)][_idx4]=0.8;
				u__c_8[(t-1)][_idx4]=0.9;
				u__c_9[(t-1)][_idx4]=1.0;
				/* _idx5 = ((((((((((((thd_idx_z+2)*x_max)+((2*t)*thd_idx_z))+(4*t))*y_max)+((((((2*t)*thd_idx_z)+thd_idx_y)+(4*t))+1)*x_max))+((4*(t*t))*thd_idx_z))+((2*t)*thd_idx_y))+thd_idx_x)+(8*(t*t)))+(2*t))+1) */
				_idx5=(((_idx3+((x_max+(2*t))*y_max))+((2*t)*x_max))+(4*(t*t)));
				u__T_0[(t-1)][_idx5]=0.1;
				/* _idx6 = ((((((((((((thd_idx_z+1)*x_max)+((2*t)*thd_idx_z))+(2*t))*y_max)+((((((2*t)*thd_idx_z)+thd_idx_y)+(2*t))+2)*x_max))+((4*(t*t))*thd_idx_z))+((2*t)*thd_idx_y))+thd_idx_x)+(4*(t*t)))+(4*t))+1) */
				_idx6=((_idx3+x_max)+(2*t));
				u__T_0[(t-1)][_idx6]=0.1;
				/* _idx7 = ((((((((((((thd_idx_z+1)*x_max)+((2*t)*thd_idx_z))+(2*t))*y_max)+((((((2*t)*thd_idx_z)+thd_idx_y)+(2*t))+1)*x_max))+((4*(t*t))*thd_idx_z))+((2*t)*thd_idx_y))+thd_idx_x)+(4*(t*t)))+(2*t))+2) */
				_idx7=(_idx0+2);
				u__T_0[(t-1)][_idx7]=0.1;
				u__T_0[t][_idx3]=1.1;
			}
		}
	}
}

