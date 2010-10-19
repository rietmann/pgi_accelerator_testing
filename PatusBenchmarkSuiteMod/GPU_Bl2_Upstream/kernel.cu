#define t_max 1
#define t 1

/*
(u[0][0][0][1][0]=(a*((((u[-3][0][0][0][0]+(u[0][-3][0][0][0]+u[0][0][-3][0][0]))*-2.0)+(((u[-2][0][0][0][0]+(u[0][-2][0][0][0]+u[0][0][-2][0][0]))*15.0)+((u[-1][0][0][0][0]+(u[0][-1][0][0][0]+u[0][0][-1][0][0]))*-60.0)))+((u[0][0][0][0][0]*20.0)+(((u[1][0][0][0][0]+(u[0][1][0][0][0]+u[0][0][1][0][0]))*30.0)+((u[2][0][0][0][0]+(u[0][2][0][0][0]+u[0][0][2][0][0]))*-3.0))))))

*/
__global__ void upstream_5_3d(double *  *  u_0_1_out, double *  u_0_0, double *  u_0_1, double a, int x_max, int y_max, int z_max, int tbx, int tby, int tbz, int c)
{
	//double *  const u__u_0[16] =  { u_0_0, u_0_1 } ;
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
				/* _idx0 = (((((((((((thd_idx_z+3)*x_max)+((5*t)*thd_idx_z))+(15*t))*y_max)+((((((5*t)*thd_idx_z)+thd_idx_y)+(15*t))+3)*x_max))+((25*(t*t))*thd_idx_z))+((5*t)*thd_idx_y))+thd_idx_x)+(75*(t*t)))+(15*t)) */
				_idx0=(((((((((((thd_idx_z+3)*x_max)+((5*t)*thd_idx_z))+(15*t))*y_max)+((((((5*t)*thd_idx_z)+thd_idx_y)+(15*t))+3)*x_max))+((25*(t*t))*thd_idx_z))+((5*t)*thd_idx_y))+thd_idx_x)+(75*(t*t)))+(15*t));
				/* _idx1 = (((((((((((thd_idx_z+3)*x_max)+((5*t)*thd_idx_z))+(15*t))*y_max)+(((((5*t)*thd_idx_z)+thd_idx_y)+(15*t))*x_max))+((25*(t*t))*thd_idx_z))+((5*t)*thd_idx_y))+thd_idx_x)+(75*(t*t)))+3) */
				_idx1=(((_idx0-(3*x_max))-(15*t))+3);
				/* _idx2 = (((((((((thd_idx_z*x_max)+((5*t)*thd_idx_z))*y_max)+(((((5*t)*thd_idx_z)+thd_idx_y)+3)*x_max))+((25*(t*t))*thd_idx_z))+((5*t)*thd_idx_y))+thd_idx_x)+(15*t))+3) */
				_idx2=((((_idx0+(((-3*x_max)-(15*t))*y_max))-((15*t)*x_max))-(75*(t*t)))+3);
				/* _idx3 = ((((((((((((thd_idx_z+3)*x_max)+((5*t)*thd_idx_z))+(15*t))*y_max)+((((((5*t)*thd_idx_z)+thd_idx_y)+(15*t))+3)*x_max))+((25*(t*t))*thd_idx_z))+((5*t)*thd_idx_y))+thd_idx_x)+(75*(t*t)))+(15*t))+1) */
				_idx3=(_idx0+1);
				/* _idx4 = ((((((((((((thd_idx_z+3)*x_max)+((5*t)*thd_idx_z))+(15*t))*y_max)+((((((5*t)*thd_idx_z)+thd_idx_y)+(15*t))+1)*x_max))+((25*(t*t))*thd_idx_z))+((5*t)*thd_idx_y))+thd_idx_x)+(75*(t*t)))+(5*t))+3) */
				_idx4=((_idx1+x_max)+(5*t));
				/* _idx5 = ((((((((((((thd_idx_z+1)*x_max)+((5*t)*thd_idx_z))+(5*t))*y_max)+((((((5*t)*thd_idx_z)+thd_idx_y)+(5*t))+3)*x_max))+((25*(t*t))*thd_idx_z))+((5*t)*thd_idx_y))+thd_idx_x)+(25*(t*t)))+(15*t))+3) */
				_idx5=(((_idx2+((x_max+(5*t))*y_max))+((5*t)*x_max))+(25*(t*t)));
				/* _idx6 = ((((((((((((thd_idx_z+3)*x_max)+((5*t)*thd_idx_z))+(15*t))*y_max)+((((((5*t)*thd_idx_z)+thd_idx_y)+(15*t))+3)*x_max))+((25*(t*t))*thd_idx_z))+((5*t)*thd_idx_y))+thd_idx_x)+(75*(t*t)))+(15*t))+2) */
				_idx6=(_idx3+1);
				/* _idx7 = ((((((((((((thd_idx_z+3)*x_max)+((5*t)*thd_idx_z))+(15*t))*y_max)+((((((5*t)*thd_idx_z)+thd_idx_y)+(15*t))+2)*x_max))+((25*(t*t))*thd_idx_z))+((5*t)*thd_idx_y))+thd_idx_x)+(75*(t*t)))+(10*t))+3) */
				_idx7=((_idx4+x_max)+(5*t));
				/* _idx8 = ((((((((((((thd_idx_z+2)*x_max)+((5*t)*thd_idx_z))+(10*t))*y_max)+((((((5*t)*thd_idx_z)+thd_idx_y)+(10*t))+3)*x_max))+((25*(t*t))*thd_idx_z))+((5*t)*thd_idx_y))+thd_idx_x)+(50*(t*t)))+(15*t))+3) */
				_idx8=(((_idx5+((x_max+(5*t))*y_max))+((5*t)*x_max))+(25*(t*t)));
				/* _idx9 = ((((((((((((thd_idx_z+3)*x_max)+((5*t)*thd_idx_z))+(15*t))*y_max)+((((((5*t)*thd_idx_z)+thd_idx_y)+(15*t))+3)*x_max))+((25*(t*t))*thd_idx_z))+((5*t)*thd_idx_y))+thd_idx_x)+(75*(t*t)))+(15*t))+3) */
				_idx9=(_idx3+2);
				/* _idx10 = ((((((((((((thd_idx_z+3)*x_max)+((5*t)*thd_idx_z))+(15*t))*y_max)+((((((5*t)*thd_idx_z)+thd_idx_y)+(15*t))+3)*x_max))+((25*(t*t))*thd_idx_z))+((5*t)*thd_idx_y))+thd_idx_x)+(75*(t*t)))+(15*t))+4) */
				_idx10=(_idx3+3);
				/* _idx11 = ((((((((((((thd_idx_z+3)*x_max)+((5*t)*thd_idx_z))+(15*t))*y_max)+((((((5*t)*thd_idx_z)+thd_idx_y)+(15*t))+4)*x_max))+((25*(t*t))*thd_idx_z))+((5*t)*thd_idx_y))+thd_idx_x)+(75*(t*t)))+(20*t))+3) */
				_idx11=((_idx9+x_max)+(5*t));
				/* _idx12 = ((((((((((((thd_idx_z+4)*x_max)+((5*t)*thd_idx_z))+(20*t))*y_max)+((((((5*t)*thd_idx_z)+thd_idx_y)+(20*t))+3)*x_max))+((25*(t*t))*thd_idx_z))+((5*t)*thd_idx_y))+thd_idx_x)+(100*(t*t)))+(15*t))+3) */
				_idx12=(((_idx9+((x_max+(5*t))*y_max))+((5*t)*x_max))+(25*(t*t)));
				/* _idx13 = ((((((((((((thd_idx_z+3)*x_max)+((5*t)*thd_idx_z))+(15*t))*y_max)+((((((5*t)*thd_idx_z)+thd_idx_y)+(15*t))+3)*x_max))+((25*(t*t))*thd_idx_z))+((5*t)*thd_idx_y))+thd_idx_x)+(75*(t*t)))+(15*t))+5) */
				_idx13=(_idx3+4);
				/* _idx14 = ((((((((((((thd_idx_z+3)*x_max)+((5*t)*thd_idx_z))+(15*t))*y_max)+((((((5*t)*thd_idx_z)+thd_idx_y)+(15*t))+5)*x_max))+((25*(t*t))*thd_idx_z))+((5*t)*thd_idx_y))+thd_idx_x)+(75*(t*t)))+(25*t))+3) */
				_idx14=((_idx11+x_max)+(5*t));
				/* _idx15 = ((((((((((((thd_idx_z+5)*x_max)+((5*t)*thd_idx_z))+(25*t))*y_max)+((((((5*t)*thd_idx_z)+thd_idx_y)+(25*t))+3)*x_max))+((25*(t*t))*thd_idx_z))+((5*t)*thd_idx_y))+thd_idx_x)+(125*(t*t)))+(15*t))+3) */
				_idx15=(((_idx12+((x_max+(5*t))*y_max))+((5*t)*x_max))+(25*(t*t)));
				u_0_1[_idx9]=(a*((((u_0_0[_idx0]+(u_0_0[_idx1]+u_0_0[_idx2]))*-2.0)+(((u_0_0[_idx3]+(u_0_0[_idx4]+u_0_0[_idx5]))*15.0)+((u_0_0[_idx6]+(u_0_0[_idx7]+u_0_0[_idx8]))*-60.0)))+((u_0_0[_idx9]*20.0)+(((u_0_0[_idx10]+(u_0_0[_idx11]+u_0_0[_idx12]))*30.0)+((u_0_0[_idx13]+(u_0_0[_idx14]+u_0_0[_idx15]))*-3.0)))));
			}
		}
	}
}

__global__ void initialize(double *  u_0_0, double *  u_0_1, double a, int x_max, int y_max, int z_max, int tbx, int tby, int tbz, int c)
{
	double *  const u__u_0[16] =  { u_0_0, u_0_1 } ;
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
				/* _idx0 = (((((((((((thd_idx_z+3)*x_max)+((5*t)*thd_idx_z))+(15*t))*y_max)+((((((5*t)*thd_idx_z)+thd_idx_y)+(15*t))+3)*x_max))+((25*(t*t))*thd_idx_z))+((5*t)*thd_idx_y))+thd_idx_x)+(75*(t*t)))+(15*t)) */
				_idx0=(((((((((((thd_idx_z+3)*x_max)+((5*t)*thd_idx_z))+(15*t))*y_max)+((((((5*t)*thd_idx_z)+thd_idx_y)+(15*t))+3)*x_max))+((25*(t*t))*thd_idx_z))+((5*t)*thd_idx_y))+thd_idx_x)+(75*(t*t)))+(15*t));
				u__u_0[(t-1)][_idx0]=0.1;
				/* _idx1 = ((((((((((((thd_idx_z+3)*x_max)+((5*t)*thd_idx_z))+(15*t))*y_max)+((((((5*t)*thd_idx_z)+thd_idx_y)+(15*t))+3)*x_max))+((25*(t*t))*thd_idx_z))+((5*t)*thd_idx_y))+thd_idx_x)+(75*(t*t)))+(15*t))+1) */
				_idx1=(_idx0+1);
				u__u_0[(t-1)][_idx1]=0.1;
				/* _idx2 = ((((((((((((thd_idx_z+3)*x_max)+((5*t)*thd_idx_z))+(15*t))*y_max)+((((((5*t)*thd_idx_z)+thd_idx_y)+(15*t))+3)*x_max))+((25*(t*t))*thd_idx_z))+((5*t)*thd_idx_y))+thd_idx_x)+(75*(t*t)))+(15*t))+2) */
				_idx2=(_idx1+1);
				u__u_0[(t-1)][_idx2]=0.1;
				/* _idx3 = (((((((((((thd_idx_z+3)*x_max)+((5*t)*thd_idx_z))+(15*t))*y_max)+(((((5*t)*thd_idx_z)+thd_idx_y)+(15*t))*x_max))+((25*(t*t))*thd_idx_z))+((5*t)*thd_idx_y))+thd_idx_x)+(75*(t*t)))+3) */
				_idx3=(((_idx1-(3*x_max))-(15*t))+2);
				u__u_0[(t-1)][_idx3]=0.1;
				/* _idx4 = ((((((((((((thd_idx_z+3)*x_max)+((5*t)*thd_idx_z))+(15*t))*y_max)+((((((5*t)*thd_idx_z)+thd_idx_y)+(15*t))+1)*x_max))+((25*(t*t))*thd_idx_z))+((5*t)*thd_idx_y))+thd_idx_x)+(75*(t*t)))+(5*t))+3) */
				_idx4=((_idx3+x_max)+(5*t));
				u__u_0[(t-1)][_idx4]=0.1;
				/* _idx5 = ((((((((((((thd_idx_z+3)*x_max)+((5*t)*thd_idx_z))+(15*t))*y_max)+((((((5*t)*thd_idx_z)+thd_idx_y)+(15*t))+2)*x_max))+((25*(t*t))*thd_idx_z))+((5*t)*thd_idx_y))+thd_idx_x)+(75*(t*t)))+(10*t))+3) */
				_idx5=((_idx4+x_max)+(5*t));
				u__u_0[(t-1)][_idx5]=0.1;
				/* _idx6 = (((((((((thd_idx_z*x_max)+((5*t)*thd_idx_z))*y_max)+(((((5*t)*thd_idx_z)+thd_idx_y)+3)*x_max))+((25*(t*t))*thd_idx_z))+((5*t)*thd_idx_y))+thd_idx_x)+(15*t))+3) */
				_idx6=((((_idx1+(((-3*x_max)-(15*t))*y_max))-((15*t)*x_max))-(75*(t*t)))+2);
				u__u_0[(t-1)][_idx6]=0.1;
				/* _idx7 = ((((((((((((thd_idx_z+1)*x_max)+((5*t)*thd_idx_z))+(5*t))*y_max)+((((((5*t)*thd_idx_z)+thd_idx_y)+(5*t))+3)*x_max))+((25*(t*t))*thd_idx_z))+((5*t)*thd_idx_y))+thd_idx_x)+(25*(t*t)))+(15*t))+3) */
				_idx7=(((_idx6+((x_max+(5*t))*y_max))+((5*t)*x_max))+(25*(t*t)));
				u__u_0[(t-1)][_idx7]=0.1;
				/* _idx8 = ((((((((((((thd_idx_z+2)*x_max)+((5*t)*thd_idx_z))+(10*t))*y_max)+((((((5*t)*thd_idx_z)+thd_idx_y)+(10*t))+3)*x_max))+((25*(t*t))*thd_idx_z))+((5*t)*thd_idx_y))+thd_idx_x)+(50*(t*t)))+(15*t))+3) */
				_idx8=(((_idx7+((x_max+(5*t))*y_max))+((5*t)*x_max))+(25*(t*t)));
				u__u_0[(t-1)][_idx8]=0.1;
				/* _idx9 = ((((((((((((thd_idx_z+3)*x_max)+((5*t)*thd_idx_z))+(15*t))*y_max)+((((((5*t)*thd_idx_z)+thd_idx_y)+(15*t))+3)*x_max))+((25*(t*t))*thd_idx_z))+((5*t)*thd_idx_y))+thd_idx_x)+(75*(t*t)))+(15*t))+3) */
				_idx9=(_idx1+2);
				u__u_0[(t-1)][_idx9]=0.1;
				/* _idx10 = ((((((((((((thd_idx_z+4)*x_max)+((5*t)*thd_idx_z))+(20*t))*y_max)+((((((5*t)*thd_idx_z)+thd_idx_y)+(20*t))+3)*x_max))+((25*(t*t))*thd_idx_z))+((5*t)*thd_idx_y))+thd_idx_x)+(100*(t*t)))+(15*t))+3) */
				_idx10=(((_idx9+((x_max+(5*t))*y_max))+((5*t)*x_max))+(25*(t*t)));
				u__u_0[(t-1)][_idx10]=0.1;
				/* _idx11 = ((((((((((((thd_idx_z+5)*x_max)+((5*t)*thd_idx_z))+(25*t))*y_max)+((((((5*t)*thd_idx_z)+thd_idx_y)+(25*t))+3)*x_max))+((25*(t*t))*thd_idx_z))+((5*t)*thd_idx_y))+thd_idx_x)+(125*(t*t)))+(15*t))+3) */
				_idx11=(((_idx10+((x_max+(5*t))*y_max))+((5*t)*x_max))+(25*(t*t)));
				u__u_0[(t-1)][_idx11]=0.1;
				/* _idx12 = ((((((((((((thd_idx_z+3)*x_max)+((5*t)*thd_idx_z))+(15*t))*y_max)+((((((5*t)*thd_idx_z)+thd_idx_y)+(15*t))+4)*x_max))+((25*(t*t))*thd_idx_z))+((5*t)*thd_idx_y))+thd_idx_x)+(75*(t*t)))+(20*t))+3) */
				_idx12=((_idx9+x_max)+(5*t));
				u__u_0[(t-1)][_idx12]=0.1;
				/* _idx13 = ((((((((((((thd_idx_z+3)*x_max)+((5*t)*thd_idx_z))+(15*t))*y_max)+((((((5*t)*thd_idx_z)+thd_idx_y)+(15*t))+5)*x_max))+((25*(t*t))*thd_idx_z))+((5*t)*thd_idx_y))+thd_idx_x)+(75*(t*t)))+(25*t))+3) */
				_idx13=((_idx12+x_max)+(5*t));
				u__u_0[(t-1)][_idx13]=0.1;
				/* _idx14 = ((((((((((((thd_idx_z+3)*x_max)+((5*t)*thd_idx_z))+(15*t))*y_max)+((((((5*t)*thd_idx_z)+thd_idx_y)+(15*t))+3)*x_max))+((25*(t*t))*thd_idx_z))+((5*t)*thd_idx_y))+thd_idx_x)+(75*(t*t)))+(15*t))+4) */
				_idx14=(_idx1+3);
				u__u_0[(t-1)][_idx14]=0.1;
				/* _idx15 = ((((((((((((thd_idx_z+3)*x_max)+((5*t)*thd_idx_z))+(15*t))*y_max)+((((((5*t)*thd_idx_z)+thd_idx_y)+(15*t))+3)*x_max))+((25*(t*t))*thd_idx_z))+((5*t)*thd_idx_y))+thd_idx_x)+(75*(t*t)))+(15*t))+5) */
				_idx15=(_idx1+4);
				u__u_0[(t-1)][_idx15]=0.1;
				u__u_0[t][_idx9]=1.1;
			}
		}
	}
}

