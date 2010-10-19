#define t_max 1
#define t 1

/*
(ux[0][0][0][0][1]=(alpha*(u[1][0][0][0][0]+u[-1][0][0][0][0])))
(uy[0][0][0][0][2]=(beta*(u[0][1][0][0][0]+u[0][-1][0][0][0])))
(uz[0][0][0][0][3]=(gamma*(u[0][0][1][0][0]+u[0][0][-1][0][0])))

*/
__global__ void gradient(float *  *  ux_1_0_out, float *  *  uy_2_0_out, float *  *  uz_3_0_out, float *  u_0_0, float *  ux_1_0, float *  uy_2_0, float *  uz_3_0, float alpha, float beta, float gamma, int x_max, int y_max, int z_max)
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
	int _idx6;
	int idx_1_2;
	int p_idx_x;
	int p_idx_x_max;
	int p_idx_y;
	int p_idx_y_max;
	int p_idx_z;
	int p_idx_z_max;
	int size_1_1;
	int size_1_2;
	//int t;
	int tmp;
	/*
	Initializations
	*/
	size_1_1=(y_max/blockDim.y);
	size_1_2=(z_max/blockDim.z);
	idx_1_2=(blockIdx.y/size_1_2);
	tmp=(blockIdx.y-(idx_1_2*size_1_2));
	p_idx_x=(threadIdx.x+(blockDim.x*blockIdx.x));
	p_idx_x_max=(p_idx_x+1);
	p_idx_y=(threadIdx.y+(tmp*blockDim.y));
	p_idx_y_max=(p_idx_y+1);
	p_idx_z=(threadIdx.z+(idx_1_2*blockDim.z));
	p_idx_z_max=(p_idx_z+1);
	/*
	Implementation
	*/
	/*
	for t = 1..t_max by 1 parallel 1 <level 0> schedule  { ... }
	*/
	//for (t=1; t<=t_max; t+=1)
	{
		/* Index bounds calculations for iterators in p[t=t, s=(1, 1, 1)][0] */
		/*
		u[t=(t+1), s=p[t=?, s=?][0]][0]=stencil(u[t=t, s=p[t=?, s=?][0]][0])
		*/
		/* _idx0 = (((((((((p_idx_z+1)*x_max)+(((2*p_idx_z)+2)*t))*y_max)+((((((2*p_idx_z)+2)*t)+p_idx_y)+1)*x_max))+(((4*p_idx_z)+4)*(t*t)))+(((2*p_idx_y)+2)*t))+p_idx_x)+2) */
		_idx0=(((((((((p_idx_z+1)*x_max)+(((2*p_idx_z)+2)*t))*y_max)+((((((2*p_idx_z)+2)*t)+p_idx_y)+1)*x_max))+(((4*p_idx_z)+4)*(t*t)))+(((2*p_idx_y)+2)*t))+p_idx_x)+2);
		/* _idx1 = ((((((((p_idx_z+1)*x_max)+(((2*p_idx_z)+2)*t))*y_max)+((((((2*p_idx_z)+2)*t)+p_idx_y)+1)*x_max))+(((4*p_idx_z)+4)*(t*t)))+(((2*p_idx_y)+2)*t))+p_idx_x) */
		_idx1=(_idx0-2);
		/* _idx2 = ((((p_idx_z*x_max)*y_max)+(p_idx_y*x_max))+p_idx_x) */
		_idx2=((((_idx1+(((((-2*p_idx_z)-2)*t)-x_max)*y_max))+(((((-2*p_idx_z)-2)*t)-1)*x_max))+(((-4*p_idx_z)-4)*(t*t)))+(((-2*p_idx_y)-2)*t));
		ux_1_0[_idx2]=(alpha*(u_0_0[_idx0]+u_0_0[_idx1]));
		/* _idx3 = (((((((((p_idx_z+1)*x_max)+(((2*p_idx_z)+2)*t))*y_max)+((((((2*p_idx_z)+2)*t)+p_idx_y)+2)*x_max))+(((4*p_idx_z)+4)*(t*t)))+(((2*p_idx_y)+4)*t))+p_idx_x)+1) */
		_idx3=(((_idx1+x_max)+(2*t))+1);
		/* _idx4 = (((((((((p_idx_z+1)*x_max)+(((2*p_idx_z)+2)*t))*y_max)+(((((2*p_idx_z)+2)*t)+p_idx_y)*x_max))+(((4*p_idx_z)+4)*(t*t)))+((2*p_idx_y)*t))+p_idx_x)+1) */
		_idx4=(((_idx1-x_max)-(2*t))+1);
		uy_2_0[_idx2]=(beta*(u_0_0[_idx3]+u_0_0[_idx4]));
		/* _idx5 = (((((((((p_idx_z+2)*x_max)+(((2*p_idx_z)+4)*t))*y_max)+((((((2*p_idx_z)+4)*t)+p_idx_y)+1)*x_max))+(((4*p_idx_z)+8)*(t*t)))+(((2*p_idx_y)+2)*t))+p_idx_x)+1) */
		_idx5=((((_idx4+((x_max+(2*t))*y_max))+(((2*t)+1)*x_max))+(4*(t*t)))+(2*t));
		/* _idx6 = ((((((((p_idx_z*x_max)+((2*p_idx_z)*t))*y_max)+(((((2*p_idx_z)*t)+p_idx_y)+1)*x_max))+((4*p_idx_z)*(t*t)))+(((2*p_idx_y)+2)*t))+p_idx_x)+1) */
		_idx6=((((_idx1+((( - x_max)-(2*t))*y_max))-((2*t)*x_max))-(4*(t*t)))+1);
		uz_3_0[_idx2]=(gamma*(u_0_0[_idx5]+u_0_0[_idx6]));
	}
}

__global__ void initialize(float *  u_0_0, float *  ux_1_0, float *  uy_2_0, float *  uz_3_0, float alpha, float beta, float gamma, int x_max, int y_max, int z_max)
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
	int _idx6;
	int idx_1_2;
	int p_idx_x;
	int p_idx_x_max;
	int p_idx_y;
	int p_idx_y_max;
	int p_idx_z;
	int p_idx_z_max;
	int size_1_1;
	int size_1_2;
	//int t;
	int tmp;
	/*
	Initializations
	*/
	size_1_1=(y_max/blockDim.y);
	size_1_2=(z_max/blockDim.z);
	idx_1_2=(blockIdx.y/size_1_2);
	tmp=(blockIdx.y-(idx_1_2*size_1_2));
	p_idx_x=(threadIdx.x+(blockDim.x*blockIdx.x));
	p_idx_x_max=(p_idx_x+1);
	p_idx_y=(threadIdx.y+(tmp*blockDim.y));
	p_idx_y_max=(p_idx_y+1);
	p_idx_z=(threadIdx.z+(idx_1_2*blockDim.z));
	p_idx_z_max=(p_idx_z+1);
	/*
	Implementation
	*/
	/*
	for t = 1..t_max by 1 parallel 1 <level 0> schedule  { ... }
	*/
	//for (t=1; t<=t_max; t+=1)
	{
		/* Index bounds calculations for iterators in p[t=t, s=(1, 1, 1)][0] */
		/*
		u[t=(t+1), s=p[t=?, s=?][0]][0]=stencil(u[t=t, s=p[t=?, s=?][0]][0])
		*/
		/* _idx0 = ((((((((p_idx_z+1)*x_max)+(((2*p_idx_z)+2)*t))*y_max)+((((((2*p_idx_z)+2)*t)+p_idx_y)+1)*x_max))+(((4*p_idx_z)+4)*(t*t)))+(((2*p_idx_y)+2)*t))+p_idx_x) */
		_idx0=((((((((p_idx_z+1)*x_max)+(((2*p_idx_z)+2)*t))*y_max)+((((((2*p_idx_z)+2)*t)+p_idx_y)+1)*x_max))+(((4*p_idx_z)+4)*(t*t)))+(((2*p_idx_y)+2)*t))+p_idx_x);
		u_0_0[_idx0]=0.1;
		/* _idx1 = (((((((((p_idx_z+1)*x_max)+(((2*p_idx_z)+2)*t))*y_max)+((((((2*p_idx_z)+2)*t)+p_idx_y)+1)*x_max))+(((4*p_idx_z)+4)*(t*t)))+(((2*p_idx_y)+2)*t))+p_idx_x)+2) */
		_idx1=(_idx0+2);
		u_0_0[_idx1]=0.1;
		/* _idx2 = ((((p_idx_z*x_max)*y_max)+(p_idx_y*x_max))+p_idx_x) */
		_idx2=((((_idx0+(((((-2*p_idx_z)-2)*t)-x_max)*y_max))+(((((-2*p_idx_z)-2)*t)-1)*x_max))+(((-4*p_idx_z)-4)*(t*t)))+(((-2*p_idx_y)-2)*t));
		u__ux_1[(t-1)][_idx2]=0.2;
		/* _idx3 = (((((((((p_idx_z+1)*x_max)+(((2*p_idx_z)+2)*t))*y_max)+(((((2*p_idx_z)+2)*t)+p_idx_y)*x_max))+(((4*p_idx_z)+4)*(t*t)))+((2*p_idx_y)*t))+p_idx_x)+1) */
		_idx3=(((_idx0-x_max)-(2*t))+1);
		u_0_0[_idx3]=0.1;
		/* _idx4 = (((((((((p_idx_z+1)*x_max)+(((2*p_idx_z)+2)*t))*y_max)+((((((2*p_idx_z)+2)*t)+p_idx_y)+2)*x_max))+(((4*p_idx_z)+4)*(t*t)))+(((2*p_idx_y)+4)*t))+p_idx_x)+1) */
		_idx4=((_idx3+(2*x_max))+(4*t));
		u_0_0[_idx4]=0.1;
		u__uy_2[(t-1)][_idx2]=0.30000000000000004;
		/* _idx5 = ((((((((p_idx_z*x_max)+((2*p_idx_z)*t))*y_max)+(((((2*p_idx_z)*t)+p_idx_y)+1)*x_max))+((4*p_idx_z)*(t*t)))+(((2*p_idx_y)+2)*t))+p_idx_x)+1) */
		_idx5=((((_idx0+((( - x_max)-(2*t))*y_max))-((2*t)*x_max))-(4*(t*t)))+1);
		u_0_0[_idx5]=0.1;
		/* _idx6 = (((((((((p_idx_z+2)*x_max)+(((2*p_idx_z)+4)*t))*y_max)+((((((2*p_idx_z)+4)*t)+p_idx_y)+1)*x_max))+(((4*p_idx_z)+8)*(t*t)))+(((2*p_idx_y)+2)*t))+p_idx_x)+1) */
		_idx6=((((_idx3+((x_max+(2*t))*y_max))+(((2*t)+1)*x_max))+(4*(t*t)))+(2*t));
		u_0_0[_idx6]=0.1;
		u__uz_3[(t-1)][_idx2]=0.4;
	}
}

