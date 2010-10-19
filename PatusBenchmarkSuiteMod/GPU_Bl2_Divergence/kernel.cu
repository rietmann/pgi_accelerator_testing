#define t_max 1
#define t 1

/*
(u[0][0][0][0][0]=((alpha*(ux[1][0][0][0][1]-ux[-1][0][0][0][1]))+((beta*(uy[0][1][0][0][2]-uy[0][-1][0][0][2]))+(gamma*(uz[0][0][1][0][3]-uz[0][0][-1][0][3])))))

*/
__global__ void divergence(float *  *  u_0_0_out, float *  u_0_0, float *  ux_1_0, float *  uy_2_0, float *  uz_3_0, float alpha, float beta, float gamma, int x_max, int y_max, int z_max, int tbx, int tby, int tbz, int c)
{
/*
	float *  const u__u_0[16] =  { u_0_0 } ;
	float *  const u__ux_1[16] =  { ux_1_0 } ;
	float *  const u__uy_2[16] =  { uy_2_0 } ;
	float *  const u__uz_3[16] =  { uz_3_0 } ;
	*/
	
	int _idx0;
	int _idx1;
	int _idx2;
	int _idx3;
	int _idx4;
	int _idx5;
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
				/* _idx0 = (((((((thd_idx_z*x_max)+((2*t)*thd_idx_z))*y_max)+(thd_idx_y*x_max))+((2*t)*thd_idx_y))+thd_idx_x)+2) */
				_idx0=(((((((thd_idx_z*x_max)+((2*t)*thd_idx_z))*y_max)+(thd_idx_y*x_max))+((2*t)*thd_idx_y))+thd_idx_x)+2);
				/* _idx1 = ((((((thd_idx_z*x_max)+((2*t)*thd_idx_z))*y_max)+(thd_idx_y*x_max))+((2*t)*thd_idx_y))+thd_idx_x) */
				_idx1=(_idx0-2);
				/* _idx2 = ((((thd_idx_z*x_max)*y_max)+(((((2*t)*thd_idx_z)+thd_idx_y)+2)*x_max))+thd_idx_x) */
				_idx2=(((_idx1-(((2*t)*thd_idx_z)*y_max))+((((2*t)*thd_idx_z)+2)*x_max))-((2*t)*thd_idx_y));
				/* _idx3 = ((((thd_idx_z*x_max)*y_max)+((((2*t)*thd_idx_z)+thd_idx_y)*x_max))+thd_idx_x) */
				_idx3=(_idx2-(2*x_max));
				/* _idx4 = (((((thd_idx_z+2)*x_max)*y_max)+(thd_idx_y*x_max))+thd_idx_x) */
				_idx4=((_idx3+((2*x_max)*y_max))-(((2*t)*thd_idx_z)*x_max));
				/* _idx5 = ((((thd_idx_z*x_max)*y_max)+(thd_idx_y*x_max))+thd_idx_x) */
				_idx5=(_idx4-((2*x_max)*y_max));
				u_0_0[_idx5]=((alpha*(ux_1_0[_idx0]-ux_1_0[_idx1]))+((beta*(uy_2_0[_idx2]-uy_2_0[_idx3]))+(gamma*(uz_3_0[_idx4]-uz_3_0[_idx5]))));
			}
		}
	}
}

__global__ void initialize(float *  u_0_0, float *  ux_1_0, float *  uy_2_0, float *  uz_3_0, float alpha, float beta, float gamma, int x_max, int y_max, int z_max, int tbx, int tby, int tbz, int c)
{
	float *  const u__u_0[16] =  { u_0_0 } ;
	float *  const u__ux_1[16] =  { ux_1_0 } ;
	float *  const u__uy_2[16] =  { uy_2_0 } ;
	float *  const u__uz_3[16] =  { uz_3_0 } ;
	int _idx0;
	int _idx1;
	int _idx2;
	int _idx3;
	int _idx4;
	int _idx5;
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
				/* _idx0 = ((((((thd_idx_z*x_max)+((2*t)*thd_idx_z))*y_max)+(thd_idx_y*x_max))+((2*t)*thd_idx_y))+thd_idx_x) */
				_idx0=((((((thd_idx_z*x_max)+((2*t)*thd_idx_z))*y_max)+(thd_idx_y*x_max))+((2*t)*thd_idx_y))+thd_idx_x);
				u__ux_1[(t-1)][_idx0]=0.2;
				/* _idx1 = ((((thd_idx_z*x_max)*y_max)+((((2*t)*thd_idx_z)+thd_idx_y)*x_max))+thd_idx_x) */
				_idx1=(((_idx0-(((2*t)*thd_idx_z)*y_max))+(((2*t)*thd_idx_z)*x_max))-((2*t)*thd_idx_y));
				u__uy_2[(t-1)][_idx1]=0.30000000000000004;
				/* _idx2 = ((((thd_idx_z*x_max)*y_max)+(thd_idx_y*x_max))+thd_idx_x) */
				_idx2=(_idx1-(((2*t)*thd_idx_z)*x_max));
				u__uz_3[(t-1)][_idx2]=0.4;
				/* _idx3 = (((((thd_idx_z+2)*x_max)*y_max)+(thd_idx_y*x_max))+thd_idx_x) */
				_idx3=(_idx2+((2*x_max)*y_max));
				u__uz_3[(t-1)][_idx3]=0.4;
				/* _idx4 = ((((thd_idx_z*x_max)*y_max)+(((((2*t)*thd_idx_z)+thd_idx_y)+2)*x_max))+thd_idx_x) */
				_idx4=(_idx1+(2*x_max));
				u__uy_2[(t-1)][_idx4]=0.30000000000000004;
				/* _idx5 = (((((((thd_idx_z*x_max)+((2*t)*thd_idx_z))*y_max)+(thd_idx_y*x_max))+((2*t)*thd_idx_y))+thd_idx_x)+2) */
				_idx5=(_idx0+2);
				u__ux_1[(t-1)][_idx5]=0.2;
				u__u_0[(t-1)][_idx2]=0.1;
			}
		}
	}
}

