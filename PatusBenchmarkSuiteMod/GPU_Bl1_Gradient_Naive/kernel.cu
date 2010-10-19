#define t_max 1
#define t 1

/*
(ux[0][0][0][0][1]=(alpha*(u[1][0][0][0][0]+u[-1][0][0][0][0])))
(uy[0][0][0][0][2]=(beta*(u[0][1][0][0][0]+u[0][-1][0][0][0])))
(uz[0][0][0][0][3]=(gamma*(u[0][0][1][0][0]+u[0][0][-1][0][0])))

*/
__global__ void gradient(float *  *  ux_1_0_out, float *  *  uy_2_0_out, float *  *  uz_3_0_out, float *  u_0_0, float *  ux_1_0, float *  uy_2_0, float *  uz_3_0, float alpha, float beta, float gamma, int x_max, int y_max, int z_max, int cbx)
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
	int pt_idx_x;
	int pt_idx_y;
	int pt_idx_z;
	int size_1_1;
	int size_1_2;
	//int t;
	int tmp;
	int v_idx_x;
	int v_idx_x_max;
	int v_idx_y;
	int v_idx_y_max;
	int v_idx_z;
	int v_idx_z_max;
	/*
	Initializations
	*/
	size_1_1=(y_max/blockDim.y);
	size_1_2=(z_max/blockDim.z);
	idx_1_2=(blockIdx.y/size_1_2);
	tmp=(blockIdx.y-(idx_1_2*size_1_2));
	v_idx_x=(cbx*(threadIdx.x+(blockDim.x*blockIdx.x)));
	v_idx_x_max=(v_idx_x+cbx);
	v_idx_y=(threadIdx.y+(tmp*blockDim.y));
	v_idx_y_max=(v_idx_y+1);
	v_idx_z=(threadIdx.z+(idx_1_2*blockDim.z));
	v_idx_z_max=(v_idx_z+1);
	/*
	Implementation
	*/
	/*
	for t = 1..t_max by 1 parallel 1 <level 0> schedule  { ... }
	*/
	//for (t=1; t<=t_max; t+=1)
	{
		/* Index bounds calculations for iterators in v[t=t, s=(cbx, 1, 1)][0] */
		/*
		for POINT pt[t=t, s=(1, 1, 1)][0] of size [1, 1, 1] in v[t=t, s=(:, :, :)][0] parallel 1 <level 1> schedule default { ... }
		*/
		{
			/* Index bounds calculations for iterators in pt[t=t, s=(1, 1, 1)][0] */
			pt_idx_z=v_idx_z;
			pt_idx_y=v_idx_y;
			for (pt_idx_x=v_idx_x; pt_idx_x<(v_idx_x_max-0); pt_idx_x+=1)
			{
				/* Index bounds calculations for iterators in pt[t=t, s=(1, 1, 1)][0] */
				/*
				v[t=(t+1), s=pt[t=?, s=?][0]][0]=stencil(v[t=t, s=pt[t=?, s=?][0]][0])
				*/
				/* _idx0 = (((((((((pt_idx_z+1)*x_max)+(((2*pt_idx_z)+2)*t))*y_max)+((((((2*pt_idx_z)+2)*t)+pt_idx_y)+1)*x_max))+(((4*pt_idx_z)+4)*(t*t)))+(((2*pt_idx_y)+2)*t))+pt_idx_x)+2) */
				_idx0=(((((((((pt_idx_z+1)*x_max)+(((2*pt_idx_z)+2)*t))*y_max)+((((((2*pt_idx_z)+2)*t)+pt_idx_y)+1)*x_max))+(((4*pt_idx_z)+4)*(t*t)))+(((2*pt_idx_y)+2)*t))+pt_idx_x)+2);
				/* _idx1 = ((((((((pt_idx_z+1)*x_max)+(((2*pt_idx_z)+2)*t))*y_max)+((((((2*pt_idx_z)+2)*t)+pt_idx_y)+1)*x_max))+(((4*pt_idx_z)+4)*(t*t)))+(((2*pt_idx_y)+2)*t))+pt_idx_x) */
				_idx1=(_idx0-2);
				/* _idx2 = ((((pt_idx_z*x_max)*y_max)+(pt_idx_y*x_max))+pt_idx_x) */
				_idx2=((((_idx1+(((((-2*pt_idx_z)-2)*t)-x_max)*y_max))+(((((-2*pt_idx_z)-2)*t)-1)*x_max))+(((-4*pt_idx_z)-4)*(t*t)))+(((-2*pt_idx_y)-2)*t));
				ux_1_0[_idx2]=(alpha*(u_0_0[_idx0]+u_0_0[_idx1]));
				/* _idx3 = (((((((((pt_idx_z+1)*x_max)+(((2*pt_idx_z)+2)*t))*y_max)+((((((2*pt_idx_z)+2)*t)+pt_idx_y)+2)*x_max))+(((4*pt_idx_z)+4)*(t*t)))+(((2*pt_idx_y)+4)*t))+pt_idx_x)+1) */
				_idx3=(((_idx1+x_max)+(2*t))+1);
				/* _idx4 = (((((((((pt_idx_z+1)*x_max)+(((2*pt_idx_z)+2)*t))*y_max)+(((((2*pt_idx_z)+2)*t)+pt_idx_y)*x_max))+(((4*pt_idx_z)+4)*(t*t)))+((2*pt_idx_y)*t))+pt_idx_x)+1) */
				_idx4=((_idx3-(2*x_max))-(4*t));
				uy_2_0[_idx2]=(beta*(u_0_0[_idx3]+u_0_0[_idx4]));
				/* _idx5 = (((((((((pt_idx_z+2)*x_max)+(((2*pt_idx_z)+4)*t))*y_max)+((((((2*pt_idx_z)+4)*t)+pt_idx_y)+1)*x_max))+(((4*pt_idx_z)+8)*(t*t)))+(((2*pt_idx_y)+2)*t))+pt_idx_x)+1) */
				_idx5=((((_idx3+((x_max+(2*t))*y_max))+(((2*t)-1)*x_max))+(4*(t*t)))-(2*t));
				/* _idx6 = ((((((((pt_idx_z*x_max)+((2*pt_idx_z)*t))*y_max)+(((((2*pt_idx_z)*t)+pt_idx_y)+1)*x_max))+((4*pt_idx_z)*(t*t)))+(((2*pt_idx_y)+2)*t))+pt_idx_x)+1) */
				_idx6=((((_idx1+((( - x_max)-(2*t))*y_max))-((2*t)*x_max))-(4*(t*t)))+1);
				uz_3_0[_idx2]=(gamma*(u_0_0[_idx5]+u_0_0[_idx6]));
			}
		}
	}
}

__global__ void initialize(float *  u_0_0, float *  ux_1_0, float *  uy_2_0, float *  uz_3_0, float alpha, float beta, float gamma, int x_max, int y_max, int z_max, int cbx)
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
	int pt_idx_x;
	int pt_idx_y;
	int pt_idx_z;
	int size_1_1;
	int size_1_2;
	//int t;
	int tmp;
	int v_idx_x;
	int v_idx_x_max;
	int v_idx_y;
	int v_idx_y_max;
	int v_idx_z;
	int v_idx_z_max;
	/*
	Initializations
	*/
	size_1_1=(y_max/blockDim.y);
	size_1_2=(z_max/blockDim.z);
	idx_1_2=(blockIdx.y/size_1_2);
	tmp=(blockIdx.y-(idx_1_2*size_1_2));
	v_idx_x=(cbx*(threadIdx.x+(blockDim.x*blockIdx.x)));
	v_idx_x_max=(v_idx_x+cbx);
	v_idx_y=(threadIdx.y+(tmp*blockDim.y));
	v_idx_y_max=(v_idx_y+1);
	v_idx_z=(threadIdx.z+(idx_1_2*blockDim.z));
	v_idx_z_max=(v_idx_z+1);
	/*
	Implementation
	*/
	/*
	for t = 1..t_max by 1 parallel 1 <level 0> schedule  { ... }
	*/
	//for (t=1; t<=t_max; t+=1)
	{
		/* Index bounds calculations for iterators in v[t=t, s=(cbx, 1, 1)][0] */
		/*
		for POINT pt[t=t, s=(1, 1, 1)][0] of size [1, 1, 1] in v[t=t, s=(:, :, :)][0] parallel 1 <level 1> schedule default { ... }
		*/
		{
			/* Index bounds calculations for iterators in pt[t=t, s=(1, 1, 1)][0] */
			pt_idx_z=v_idx_z;
			pt_idx_y=v_idx_y;
			for (pt_idx_x=v_idx_x; pt_idx_x<(v_idx_x_max-0); pt_idx_x+=1)
			{
				/* Index bounds calculations for iterators in pt[t=t, s=(1, 1, 1)][0] */
				/*
				v[t=(t+1), s=pt[t=?, s=?][0]][0]=stencil(v[t=t, s=pt[t=?, s=?][0]][0])
				*/
				/* _idx0 = ((((((((pt_idx_z+1)*x_max)+(((2*pt_idx_z)+2)*t))*y_max)+((((((2*pt_idx_z)+2)*t)+pt_idx_y)+1)*x_max))+(((4*pt_idx_z)+4)*(t*t)))+(((2*pt_idx_y)+2)*t))+pt_idx_x) */
				_idx0=((((((((pt_idx_z+1)*x_max)+(((2*pt_idx_z)+2)*t))*y_max)+((((((2*pt_idx_z)+2)*t)+pt_idx_y)+1)*x_max))+(((4*pt_idx_z)+4)*(t*t)))+(((2*pt_idx_y)+2)*t))+pt_idx_x);
				u_0_0[_idx0]=0.1;
				/* _idx1 = (((((((((pt_idx_z+1)*x_max)+(((2*pt_idx_z)+2)*t))*y_max)+((((((2*pt_idx_z)+2)*t)+pt_idx_y)+1)*x_max))+(((4*pt_idx_z)+4)*(t*t)))+(((2*pt_idx_y)+2)*t))+pt_idx_x)+2) */
				_idx1=(_idx0+2);
				u_0_0[_idx1]=0.1;
				/* _idx2 = ((((pt_idx_z*x_max)*y_max)+(pt_idx_y*x_max))+pt_idx_x) */
				_idx2=((((_idx0+(((((-2*pt_idx_z)-2)*t)-x_max)*y_max))+(((((-2*pt_idx_z)-2)*t)-1)*x_max))+(((-4*pt_idx_z)-4)*(t*t)))+(((-2*pt_idx_y)-2)*t));
				ux_1_0[_idx2]=0.2;
				/* _idx3 = (((((((((pt_idx_z+1)*x_max)+(((2*pt_idx_z)+2)*t))*y_max)+(((((2*pt_idx_z)+2)*t)+pt_idx_y)*x_max))+(((4*pt_idx_z)+4)*(t*t)))+((2*pt_idx_y)*t))+pt_idx_x)+1) */
				_idx3=(((_idx0-x_max)-(2*t))+1);
				u_0_0[_idx3]=0.1;
				/* _idx4 = (((((((((pt_idx_z+1)*x_max)+(((2*pt_idx_z)+2)*t))*y_max)+((((((2*pt_idx_z)+2)*t)+pt_idx_y)+2)*x_max))+(((4*pt_idx_z)+4)*(t*t)))+(((2*pt_idx_y)+4)*t))+pt_idx_x)+1) */
				_idx4=((_idx3+(2*x_max))+(4*t));
				u_0_0[_idx4]=0.1;
				uy_2_0[_idx2]=0.30000000000000004;
				/* _idx5 = ((((((((pt_idx_z*x_max)+((2*pt_idx_z)*t))*y_max)+(((((2*pt_idx_z)*t)+pt_idx_y)+1)*x_max))+((4*pt_idx_z)*(t*t)))+(((2*pt_idx_y)+2)*t))+pt_idx_x)+1) */
				_idx5=((((_idx0+((( - x_max)-(2*t))*y_max))-((2*t)*x_max))-(4*(t*t)))+1);
				u_0_0[_idx5]=0.1;
				/* _idx6 = (((((((((pt_idx_z+2)*x_max)+(((2*pt_idx_z)+4)*t))*y_max)+((((((2*pt_idx_z)+4)*t)+pt_idx_y)+1)*x_max))+(((4*pt_idx_z)+8)*(t*t)))+(((2*pt_idx_y)+2)*t))+pt_idx_x)+1) */
				_idx6=((((_idx4+((x_max+(2*t))*y_max))+(((2*t)-1)*x_max))+(4*(t*t)))-(2*t));
				u_0_0[_idx6]=0.1;
				uz_3_0[_idx2]=0.4;
			}
		}
	}
}

