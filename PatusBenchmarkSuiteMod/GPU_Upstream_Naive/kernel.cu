#define t_max 1
#define t 1

/*
(u[0][0][0][1][0]=(a*((((u[-3][0][0][0][0]+(u[0][-3][0][0][0]+u[0][0][-3][0][0]))*-2.0)+(((u[-2][0][0][0][0]+(u[0][-2][0][0][0]+u[0][0][-2][0][0]))*15.0)+((u[-1][0][0][0][0]+(u[0][-1][0][0][0]+u[0][0][-1][0][0]))*-60.0)))+((u[0][0][0][0][0]*20.0)+(((u[1][0][0][0][0]+(u[0][1][0][0][0]+u[0][0][1][0][0]))*30.0)+((u[2][0][0][0][0]+(u[0][2][0][0][0]+u[0][0][2][0][0]))*-3.0))))))

*/
__global__ void upstream_5_3d(double *  *  u_0_1_out, double *  u_0_0, double *  u_0_1, double a, int x_max, int y_max, int z_max)
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
		/* _idx0 = ((((((((p_idx_z+3)*x_max)+(((5*p_idx_z)+15)*t))*y_max)+((((((5*p_idx_z)+15)*t)+p_idx_y)+3)*x_max))+(((25*p_idx_z)+75)*(t*t)))+(((5*p_idx_y)+15)*t))+p_idx_x) */
		_idx0=((((((((p_idx_z+3)*x_max)+(((5*p_idx_z)+15)*t))*y_max)+((((((5*p_idx_z)+15)*t)+p_idx_y)+3)*x_max))+(((25*p_idx_z)+75)*(t*t)))+(((5*p_idx_y)+15)*t))+p_idx_x);
		/* _idx1 = (((((((((p_idx_z+3)*x_max)+(((5*p_idx_z)+15)*t))*y_max)+(((((5*p_idx_z)+15)*t)+p_idx_y)*x_max))+(((25*p_idx_z)+75)*(t*t)))+((5*p_idx_y)*t))+p_idx_x)+3) */
		_idx1=(((_idx0-(3*x_max))-(15*t))+3);
		/* _idx2 = ((((((((p_idx_z*x_max)+((5*p_idx_z)*t))*y_max)+(((((5*p_idx_z)*t)+p_idx_y)+3)*x_max))+((25*p_idx_z)*(t*t)))+(((5*p_idx_y)+15)*t))+p_idx_x)+3) */
		_idx2=((((_idx0+(((-3*x_max)-(15*t))*y_max))-((15*t)*x_max))-(75*(t*t)))+3);
		/* _idx3 = (((((((((p_idx_z+3)*x_max)+(((5*p_idx_z)+15)*t))*y_max)+((((((5*p_idx_z)+15)*t)+p_idx_y)+3)*x_max))+(((25*p_idx_z)+75)*(t*t)))+(((5*p_idx_y)+15)*t))+p_idx_x)+1) */
		_idx3=(_idx0+1);
		/* _idx4 = (((((((((p_idx_z+3)*x_max)+(((5*p_idx_z)+15)*t))*y_max)+((((((5*p_idx_z)+15)*t)+p_idx_y)+1)*x_max))+(((25*p_idx_z)+75)*(t*t)))+(((5*p_idx_y)+5)*t))+p_idx_x)+3) */
		_idx4=((_idx1+x_max)+(5*t));
		/* _idx5 = (((((((((p_idx_z+1)*x_max)+(((5*p_idx_z)+5)*t))*y_max)+((((((5*p_idx_z)+5)*t)+p_idx_y)+3)*x_max))+(((25*p_idx_z)+25)*(t*t)))+(((5*p_idx_y)+15)*t))+p_idx_x)+3) */
		_idx5=(((_idx2+((x_max+(5*t))*y_max))+((5*t)*x_max))+(25*(t*t)));
		/* _idx6 = (((((((((p_idx_z+3)*x_max)+(((5*p_idx_z)+15)*t))*y_max)+((((((5*p_idx_z)+15)*t)+p_idx_y)+3)*x_max))+(((25*p_idx_z)+75)*(t*t)))+(((5*p_idx_y)+15)*t))+p_idx_x)+2) */
		_idx6=(_idx3+1);
		/* _idx7 = (((((((((p_idx_z+3)*x_max)+(((5*p_idx_z)+15)*t))*y_max)+((((((5*p_idx_z)+15)*t)+p_idx_y)+2)*x_max))+(((25*p_idx_z)+75)*(t*t)))+(((5*p_idx_y)+10)*t))+p_idx_x)+3) */
		_idx7=((_idx4+x_max)+(5*t));
		/* _idx8 = (((((((((p_idx_z+2)*x_max)+(((5*p_idx_z)+10)*t))*y_max)+((((((5*p_idx_z)+10)*t)+p_idx_y)+3)*x_max))+(((25*p_idx_z)+50)*(t*t)))+(((5*p_idx_y)+15)*t))+p_idx_x)+3) */
		_idx8=(((_idx5+((x_max+(5*t))*y_max))+((5*t)*x_max))+(25*(t*t)));
		/* _idx9 = (((((((((p_idx_z+3)*x_max)+(((5*p_idx_z)+15)*t))*y_max)+((((((5*p_idx_z)+15)*t)+p_idx_y)+3)*x_max))+(((25*p_idx_z)+75)*(t*t)))+(((5*p_idx_y)+15)*t))+p_idx_x)+3) */
		_idx9=(_idx3+2);
		/* _idx10 = (((((((((p_idx_z+3)*x_max)+(((5*p_idx_z)+15)*t))*y_max)+((((((5*p_idx_z)+15)*t)+p_idx_y)+3)*x_max))+(((25*p_idx_z)+75)*(t*t)))+(((5*p_idx_y)+15)*t))+p_idx_x)+4) */
		_idx10=(_idx3+3);
		/* _idx11 = (((((((((p_idx_z+3)*x_max)+(((5*p_idx_z)+15)*t))*y_max)+((((((5*p_idx_z)+15)*t)+p_idx_y)+4)*x_max))+(((25*p_idx_z)+75)*(t*t)))+(((5*p_idx_y)+20)*t))+p_idx_x)+3) */
		_idx11=((_idx9+x_max)+(5*t));
		/* _idx12 = (((((((((p_idx_z+4)*x_max)+(((5*p_idx_z)+20)*t))*y_max)+((((((5*p_idx_z)+20)*t)+p_idx_y)+3)*x_max))+(((25*p_idx_z)+100)*(t*t)))+(((5*p_idx_y)+15)*t))+p_idx_x)+3) */
		_idx12=(((_idx9+((x_max+(5*t))*y_max))+((5*t)*x_max))+(25*(t*t)));
		/* _idx13 = (((((((((p_idx_z+3)*x_max)+(((5*p_idx_z)+15)*t))*y_max)+((((((5*p_idx_z)+15)*t)+p_idx_y)+3)*x_max))+(((25*p_idx_z)+75)*(t*t)))+(((5*p_idx_y)+15)*t))+p_idx_x)+5) */
		_idx13=(_idx3+4);
		/* _idx14 = (((((((((p_idx_z+3)*x_max)+(((5*p_idx_z)+15)*t))*y_max)+((((((5*p_idx_z)+15)*t)+p_idx_y)+5)*x_max))+(((25*p_idx_z)+75)*(t*t)))+(((5*p_idx_y)+25)*t))+p_idx_x)+3) */
		_idx14=((_idx11+x_max)+(5*t));
		/* _idx15 = (((((((((p_idx_z+5)*x_max)+(((5*p_idx_z)+25)*t))*y_max)+((((((5*p_idx_z)+25)*t)+p_idx_y)+3)*x_max))+(((25*p_idx_z)+125)*(t*t)))+(((5*p_idx_y)+15)*t))+p_idx_x)+3) */
		_idx15=(((_idx12+((x_max+(5*t))*y_max))+((5*t)*x_max))+(25*(t*t)));
		u_0_1[_idx9]=(a*((((u_0_0[_idx0]+(u_0_0[_idx1]+u_0_0[_idx2]))*-2.0)+(((u_0_0[_idx3]+(u_0_0[_idx4]+u_0_0[_idx5]))*15.0)+((u_0_0[_idx6]+(u_0_0[_idx7]+u_0_0[_idx8]))*-60.0)))+((u_0_0[_idx9]*20.0)+(((u_0_0[_idx10]+(u_0_0[_idx11]+u_0_0[_idx12]))*30.0)+((u_0_0[_idx13]+(u_0_0[_idx14]+u_0_0[_idx15]))*-3.0)))));
	}
}

__global__ void initialize(double *  u_0_0, double *  u_0_1, double a, int x_max, int y_max, int z_max)
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
		/* _idx0 = ((((((((p_idx_z+3)*x_max)+(((5*p_idx_z)+15)*t))*y_max)+((((((5*p_idx_z)+15)*t)+p_idx_y)+3)*x_max))+(((25*p_idx_z)+75)*(t*t)))+(((5*p_idx_y)+15)*t))+p_idx_x) */
		_idx0=((((((((p_idx_z+3)*x_max)+(((5*p_idx_z)+15)*t))*y_max)+((((((5*p_idx_z)+15)*t)+p_idx_y)+3)*x_max))+(((25*p_idx_z)+75)*(t*t)))+(((5*p_idx_y)+15)*t))+p_idx_x);
		u_0_0[_idx0]=0.1;
		/* _idx1 = (((((((((p_idx_z+3)*x_max)+(((5*p_idx_z)+15)*t))*y_max)+((((((5*p_idx_z)+15)*t)+p_idx_y)+3)*x_max))+(((25*p_idx_z)+75)*(t*t)))+(((5*p_idx_y)+15)*t))+p_idx_x)+1) */
		_idx1=(_idx0+1);
		u_0_0[_idx1]=0.1;
		/* _idx2 = (((((((((p_idx_z+3)*x_max)+(((5*p_idx_z)+15)*t))*y_max)+((((((5*p_idx_z)+15)*t)+p_idx_y)+3)*x_max))+(((25*p_idx_z)+75)*(t*t)))+(((5*p_idx_y)+15)*t))+p_idx_x)+2) */
		_idx2=(_idx1+1);
		u_0_0[_idx2]=0.1;
		/* _idx3 = (((((((((p_idx_z+3)*x_max)+(((5*p_idx_z)+15)*t))*y_max)+(((((5*p_idx_z)+15)*t)+p_idx_y)*x_max))+(((25*p_idx_z)+75)*(t*t)))+((5*p_idx_y)*t))+p_idx_x)+3) */
		_idx3=(((_idx1-(3*x_max))-(15*t))+2);
		u_0_0[_idx3]=0.1;
		/* _idx4 = (((((((((p_idx_z+3)*x_max)+(((5*p_idx_z)+15)*t))*y_max)+((((((5*p_idx_z)+15)*t)+p_idx_y)+1)*x_max))+(((25*p_idx_z)+75)*(t*t)))+(((5*p_idx_y)+5)*t))+p_idx_x)+3) */
		_idx4=((_idx3+x_max)+(5*t));
		u_0_0[_idx4]=0.1;
		/* _idx5 = (((((((((p_idx_z+3)*x_max)+(((5*p_idx_z)+15)*t))*y_max)+((((((5*p_idx_z)+15)*t)+p_idx_y)+2)*x_max))+(((25*p_idx_z)+75)*(t*t)))+(((5*p_idx_y)+10)*t))+p_idx_x)+3) */
		_idx5=((_idx4+x_max)+(5*t));
		u_0_0[_idx5]=0.1;
		/* _idx6 = ((((((((p_idx_z*x_max)+((5*p_idx_z)*t))*y_max)+(((((5*p_idx_z)*t)+p_idx_y)+3)*x_max))+((25*p_idx_z)*(t*t)))+(((5*p_idx_y)+15)*t))+p_idx_x)+3) */
		_idx6=((((_idx1+(((-3*x_max)-(15*t))*y_max))-((15*t)*x_max))-(75*(t*t)))+2);
		u_0_0[_idx6]=0.1;
		/* _idx7 = (((((((((p_idx_z+1)*x_max)+(((5*p_idx_z)+5)*t))*y_max)+((((((5*p_idx_z)+5)*t)+p_idx_y)+3)*x_max))+(((25*p_idx_z)+25)*(t*t)))+(((5*p_idx_y)+15)*t))+p_idx_x)+3) */
		_idx7=(((_idx6+((x_max+(5*t))*y_max))+((5*t)*x_max))+(25*(t*t)));
		u_0_0[_idx7]=0.1;
		/* _idx8 = (((((((((p_idx_z+2)*x_max)+(((5*p_idx_z)+10)*t))*y_max)+((((((5*p_idx_z)+10)*t)+p_idx_y)+3)*x_max))+(((25*p_idx_z)+50)*(t*t)))+(((5*p_idx_y)+15)*t))+p_idx_x)+3) */
		_idx8=(((_idx7+((x_max+(5*t))*y_max))+((5*t)*x_max))+(25*(t*t)));
		u_0_0[_idx8]=0.1;
		/* _idx9 = (((((((((p_idx_z+3)*x_max)+(((5*p_idx_z)+15)*t))*y_max)+((((((5*p_idx_z)+15)*t)+p_idx_y)+3)*x_max))+(((25*p_idx_z)+75)*(t*t)))+(((5*p_idx_y)+15)*t))+p_idx_x)+3) */
		_idx9=(_idx1+2);
		u_0_0[_idx9]=0.1;
		/* _idx10 = (((((((((p_idx_z+4)*x_max)+(((5*p_idx_z)+20)*t))*y_max)+((((((5*p_idx_z)+20)*t)+p_idx_y)+3)*x_max))+(((25*p_idx_z)+100)*(t*t)))+(((5*p_idx_y)+15)*t))+p_idx_x)+3) */
		_idx10=(((_idx9+((x_max+(5*t))*y_max))+((5*t)*x_max))+(25*(t*t)));
		u_0_0[_idx10]=0.1;
		/* _idx11 = (((((((((p_idx_z+5)*x_max)+(((5*p_idx_z)+25)*t))*y_max)+((((((5*p_idx_z)+25)*t)+p_idx_y)+3)*x_max))+(((25*p_idx_z)+125)*(t*t)))+(((5*p_idx_y)+15)*t))+p_idx_x)+3) */
		_idx11=(((_idx10+((x_max+(5*t))*y_max))+((5*t)*x_max))+(25*(t*t)));
		u_0_0[_idx11]=0.1;
		/* _idx12 = (((((((((p_idx_z+3)*x_max)+(((5*p_idx_z)+15)*t))*y_max)+((((((5*p_idx_z)+15)*t)+p_idx_y)+4)*x_max))+(((25*p_idx_z)+75)*(t*t)))+(((5*p_idx_y)+20)*t))+p_idx_x)+3) */
		_idx12=((_idx9+x_max)+(5*t));
		u_0_0[_idx12]=0.1;
		/* _idx13 = (((((((((p_idx_z+3)*x_max)+(((5*p_idx_z)+15)*t))*y_max)+((((((5*p_idx_z)+15)*t)+p_idx_y)+5)*x_max))+(((25*p_idx_z)+75)*(t*t)))+(((5*p_idx_y)+25)*t))+p_idx_x)+3) */
		_idx13=((_idx12+x_max)+(5*t));
		u_0_0[_idx13]=0.1;
		/* _idx14 = (((((((((p_idx_z+3)*x_max)+(((5*p_idx_z)+15)*t))*y_max)+((((((5*p_idx_z)+15)*t)+p_idx_y)+3)*x_max))+(((25*p_idx_z)+75)*(t*t)))+(((5*p_idx_y)+15)*t))+p_idx_x)+4) */
		_idx14=(_idx1+3);
		u_0_0[_idx14]=0.1;
		/* _idx15 = (((((((((p_idx_z+3)*x_max)+(((5*p_idx_z)+15)*t))*y_max)+((((((5*p_idx_z)+15)*t)+p_idx_y)+3)*x_max))+(((25*p_idx_z)+75)*(t*t)))+(((5*p_idx_y)+15)*t))+p_idx_x)+5) */
		_idx15=(_idx1+4);
		u_0_0[_idx15]=0.1;
		u__u_0[t][_idx9]=1.1;
	}
}

