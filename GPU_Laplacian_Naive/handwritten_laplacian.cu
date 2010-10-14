/* CUDA_SAFE_CALL (cudaMalloc ((void**) &u0_device, Nx*sizeof (cuComplex))); */
/* nls_stencil<<< blocks, threads >>> (u_device,u0_device); */
/* CUT_DEVICE_INIT (argc, argv); */
/* CUDA_SAFE_CALL (cudaMemcpy (u,u_device, N*NUMBER_OF_DOMAINS*sizeof (cuComplex), */
/* 				  cudaMemcpyDeviceToHost)); */

__global__ void laplacian_gpu(float *  *  u_0_1_out, float * u_0_0, float * u_0_1, int x_max, int y_max, int z_max)
{
  int _idx0;
  int _idx1;
  int _idx2;
  int _idx3;
  int _idx4;
  int _idx5;
  int _idx6;
  int p_idx_x;
  int p_idx_y;
  int p_idx_z;
  //	int t;
  //float *  __restrict__ const u__u_0[16] =  { u_0_0, u_0_1 } ;
  int u_idx_x;
  int u_idx_x_max;
  int u_idx_y;
  int u_idx_y_max;
  int u_idx_z;
  int u_idx_z_max;

  __shared__ float us[BLOCKx*BLOCKy*BLOCKz];

  id = 
  
  us[

  
      /* Index bounds calculations for iterators in p[t=t, s=(1, 1, 1)][0] */
      /* for (p_idx_z=0; p_idx_z<z_max; p_idx_z+=1) */
      /* 	{ */
      /* 	  for (p_idx_y=0; p_idx_y<y_max; p_idx_y+=1) */
      /* 	    { */
      /* 	      for (p_idx_x=0; p_idx_x<x_max; p_idx_x+=1) */
      /* 		{ */
		  
  _idx0=(((((((((p_idx_z+1)*x_max)+(((2*p_idx_z)+2)))*y_max)+((((((2*p_idx_z)+2)*t)+p_idx_y)+1)*x_max))+(((4*p_idx_z)+4)))+(((2*p_idx_y)+2)))+p_idx_x)+2);
  _idx1=_idx0-2;
  _idx2=(((_idx1+x_max)+2)+1);
  _idx3=(((_idx1-x_max)-2)+1);
  _idx4=((((_idx3+((x_max+2)*y_max))+((2+1)*x_max))+4)+2);
  _idx5=((((_idx1+((( -x_max)-2)*y_max))-(2*x_max))-4)+1);
  _idx6=(_idx1+1);

  u_0_1[_idx6]=((((u_0_0[_idx0]+(u_0_0[_idx1]+u_0_0[_idx2]))+(u_0_0[_idx3]+(u_0_0[_idx4]+u_0_0[_idx5])))*0.25)-u_0_0[_idx6]);
}
}
}
    
  
  *u_0_1_out = u_0_1;
}

