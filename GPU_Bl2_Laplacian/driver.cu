#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include <cutil.h>

#include "patusrt.h"

// forward_decls -->
__global__ void initialize(float *  u_0_0, float *  u_0_1, int x_max, int y_max, int z_max, int tbx, int tby, int tbz, int c);
__global__ void laplacian(float *  *  u_0_1_out, float *  u_0_0, float *  u_0_1, int x_max, int y_max, int z_max, int tbx, int tby, int tbz, int c);
void initialize_cpu(float *  u_0_0, float *  u_0_1, int x_max, int y_max, int z_max);
void laplacian_cpu(float *  *  u_0_1_out, float *  u_0_0, float *  u_0_1, int x_max, int y_max, int z_max, int t_max);
// <--


int main (int argc, char** argv)
{
	int i;
	cudaError_t res;
	
	// prepare grids
	// declare_grids -->
	float *  u_0_1_out;
	float *  u_0_0;
	float *  u_0_1;
	float * u_0_1_out_host;
	if ((argc!=8))
	{
		printf("Wrong number of parameters. Syntax:\n%s <x_max> <y_max> <z_max> <tbx> <tby> <tbz> <c>\n", argv[0]);
		exit(-1);
	}
	int x_max = atoi(argv[1]);
	int y_max = atoi(argv[2]);
	int z_max = atoi(argv[3]);
	int tbx = atoi(argv[4]);
	int tby = atoi(argv[5]);
	int tbz = atoi(argv[6]);
	int c = atoi(argv[7]);
	// <--
	
	// allocate_grids -->
	u_0_0=((float * )malloc(((((x_max+4)*(y_max+4))*(z_max+4))*sizeof (float))));
	u_0_1=((float * )malloc(((((x_max+4)*(y_max+4))*(z_max+4))*sizeof (float))));
	u_0_1_out_host = ((float * )malloc(((((x_max+4)*(y_max+4))*(z_max+4))*sizeof (float))));
	
	// <--
	
	
	// declare_GPU_grids -->
	float *  u_0_1_out_gpu;
	float *  u_0_0_gpu;
	float *  u_0_1_gpu;

	dim3 thds(tbx, tby, tbz);
	dim3 blks((x_max/tbx), ((y_max*z_max)/(tby*tbz)), 1);
	// <--
	
	// allocate_GPU_grids -->
	cudaMalloc(((void *  * )( & u_0_1_out_gpu)), ((((x_max+4)*(y_max+4))*(z_max+4))*sizeof (float * )));
	cudaMalloc(((void *  * )( & u_0_0_gpu)), ((((x_max+4)*(y_max+4))*(z_max+4))*sizeof (float)));
	cudaMalloc(((void *  * )( & u_0_1_gpu)), ((((x_max+4)*(y_max+4))*(z_max+4))*sizeof (float)));
	// <--
	
	
	
	// initialize_grids on host-->
	initialize_cpu(u_0_0, u_0_1, x_max, y_max, z_max);
	/* initialize<<<blks, thds>>>(u_0_0_gpu, u_0_1_gpu, x_max, y_max, z_max, tbx, tby, tbz, c); */
	// <--

	
	// copy_grids_to_GPU -->
	cudaMemcpy(((void * )u_0_0_gpu), ((void * )u_0_0), ((((x_max+4)*(y_max+4))*(z_max+4))*sizeof (float)), cudaMemcpyHostToDevice);
	cudaMemcpy(((void * )u_0_1_gpu), ((void * )u_0_1), ((((x_max+4)*(y_max+4))*(z_max+4))*sizeof (float)), cudaMemcpyHostToDevice);
	// <--
	
	cudaThreadSynchronize ();
	res = cudaGetLastError ();
	if (res != cudaSuccess)
	{
		printf ("CUDA Error [Initialization]: %s.\n", cudaGetErrorString (res));
	}
	
	long nFlopsPerStencil = 7;
	int iterations = 100;
	long nGridPointsCount = iterations * ((x_max*y_max)*z_max);
	long nBytesTransferred = iterations * (((((x_max+4)*(y_max+4))*(z_max+4))*sizeof (float))+(((x_max*y_max)*z_max)*sizeof (float)));
	
//	cudaFuncSetCacheConfig (laplacian, cudaFuncCachePreferShared);
	cudaFuncSetCacheConfig (laplacian, cudaFuncCachePreferL1);
	
	// warm up
	// compute_stencil -->
	laplacian<<<blks, thds>>>(( & u_0_1_out_gpu), u_0_0_gpu, u_0_1_gpu, x_max, y_max, z_max, tbx, tby, tbz, c);
	// <--
	
	cudaThreadSynchronize ();
	res = cudaGetLastError ();
	if (res != cudaSuccess)
	{
		printf ("CUDA Error [Stencil]: %s.\n", cudaGetErrorString (res));
	}
	
	// run the benchmark
	tic ();
	for (i = 0; i < iterations; i++)
	{
		// compute_stencil -->
		laplacian<<<blks, thds>>>(( & u_0_1_out_gpu), u_0_0_gpu, u_0_1_gpu, x_max, y_max, z_max, tbx, tby, tbz, c);
		// <--
		
		cudaThreadSynchronize ();
	}
	toc (nFlopsPerStencil, nGridPointsCount, nBytesTransferred);

	cudaMemcpy(((void * )u_0_1_out_host), ((void * )u_0_1_out_gpu), ((((x_max+4)*(y_max+4))*(z_max+4))*sizeof (float)), cudaMemcpyDeviceToHost);


	//warm up
	laplacian_cpu(( & u_0_1_out), u_0_0, u_0_1, x_max, y_max, z_max,iterations);
	tic ();
	// compute_stencil -->
	laplacian_cpu(( & u_0_1_out), u_0_0, u_0_1, x_max, y_max, z_max,iterations);
	// <--
	toc (nFlopsPerStencil, nGridPointsCount, nBytesTransferred);

	int error_count=0;
	int x,y,z;
	for(y=1;y<x_max+1;y++) {
	  for(x=1;x<x_max+1;x++) {
	    for(z=1;z<y_max+1;z++) {
	      i = x + (x_max+2)*y + (x_max+2)*(y_max+2)*z;
	      if(fabs(u_0_1_out_host[i] - u_0_1[i])>0.001) {
		error_count++;
		printf("%dth error encountered at u[%d]: |%f-%f|=%5.16f\n",error_count,i,u_0_1_out_host[i],u_0_1[i],fabs(u_0_1_out_host[i] - u_0_1[i]));
		if(error_count>30) {
		  printf("too many errors\n"); exit(1);
		}
	      }
	    }
	  }
	}
	
	
	printf("Error Check Successful. No errors encountered.\n");	  
	/* saving results */
	FILE* fp;
	char filename[100];
	sprintf(filename,"/users/rietmann/tmp/stencil_laplacian_cuda_%.3d_%.3d_%.3d.dat",x_max,y_max,z_max);
	fp = fopen(filename,"w");
	if(!fp) {printf("couldn't open file: %s\n",filename); exit(1);}
	fprintf(fp,"%d\n",x_max);
	fprintf(fp,"%d\n",y_max);
	fprintf(fp,"%d\n",z_max);	  
	for(i=0;i<(x_max+4)*(y_max+4)*(z_max+4);i++) {
	  fprintf(fp,"%2.6f\n", u_0_1_out_host[i]);
	}
	fclose(fp);
	
	// free memory
	// deallocate_grids -->
	cudaFree(((void * )u_0_1_out_gpu));
	cudaFree(((void * )u_0_0_gpu));
	cudaFree(((void * )u_0_1_gpu));
	free(u_0_0);
	free(u_0_1);
	// <--
	
	
	cudaThreadExit ();
	return EXIT_SUCCESS;
}

void initialize_cpu(float *  u_0_0, float *  u_0_1, int x_max, int y_max, int z_max)
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
	int t;
	float *  __restrict__ const u__u_0[2] =  { u_0_0, u_0_1 } ;
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
						/* _idx0 = ((((((((((p_idx_z+1)*x_max)+(2*p_idx_z))+2)*y_max)+((((2*p_idx_z)+p_idx_y)+3)*x_max))+(4*p_idx_z))+(2*p_idx_y))+p_idx_x)+6) */
						_idx0=((((((((((p_idx_z+1)*x_max)+(2*p_idx_z))+2)*y_max)+((((2*p_idx_z)+p_idx_y)+3)*x_max))+(4*p_idx_z))+(2*p_idx_y))+p_idx_x)+6);
						u__u_0[0][_idx0]=0.1;
						/* _idx1 = ((((((((((p_idx_z+1)*x_max)+(2*p_idx_z))+2)*y_max)+((((2*p_idx_z)+p_idx_y)+2)*x_max))+(4*p_idx_z))+(2*p_idx_y))+p_idx_x)+5) */
						_idx1=((_idx0-x_max)-1);
						u__u_0[0][_idx1]=0.1;
						/* _idx2 = ((((((((p_idx_z*x_max)+(2*p_idx_z))*y_max)+((((2*p_idx_z)+p_idx_y)+1)*x_max))+(4*p_idx_z))+(2*p_idx_y))+p_idx_x)+3) */
						_idx2=(((_idx1+((( - x_max)-2)*y_max))-x_max)-2);
						u__u_0[0][_idx2]=0.1;
						/* _idx3 = ((((((((((p_idx_z+1)*x_max)+(2*p_idx_z))+2)*y_max)+((((2*p_idx_z)+p_idx_y)+3)*x_max))+(4*p_idx_z))+(2*p_idx_y))+p_idx_x)+7) */
						_idx3=((_idx1+x_max)+2);
						u__u_0[0][_idx3]=0.1;
						/* _idx4 = ((((((((((p_idx_z+2)*x_max)+(2*p_idx_z))+4)*y_max)+((((2*p_idx_z)+p_idx_y)+5)*x_max))+(4*p_idx_z))+(2*p_idx_y))+p_idx_x)+11) */
						_idx4=(((_idx1+((x_max+2)*y_max))+(3*x_max))+6);
						u__u_0[0][_idx4]=0.1;
						/* _idx5 = ((((((((((p_idx_z+1)*x_max)+(2*p_idx_z))+2)*y_max)+((((2*p_idx_z)+p_idx_y)+4)*x_max))+(4*p_idx_z))+(2*p_idx_y))+p_idx_x)+9) */
						_idx5=((_idx0+x_max)+3);
						u__u_0[0][_idx5]=0.1;
						/* _idx6 = ((((((((((p_idx_z+1)*x_max)+(2*p_idx_z))+2)*y_max)+((((2*p_idx_z)+p_idx_y)+3)*x_max))+(4*p_idx_z))+(2*p_idx_y))+p_idx_x)+8) */
						_idx6=((_idx1+x_max)+3);
						u__u_0[0][_idx6]=0.1;
						u__u_0[1][_idx3]=1.1;
					}
				}
			}
		}
	}
}

void laplacian_cpu(float *  *  u_0_1_out, float *  u_0_0, float *  u_0_1, int x_max, int y_max, int z_max, int t_max)
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
	  int count=0;
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
						/* _idx0 = ((((((((((p_idx_z+1)*x_max)+(2*p_idx_z))+2)*y_max)+((((2*p_idx_z)+p_idx_y)+3)*x_max))+(4*p_idx_z))+(2*p_idx_y))+p_idx_x)+8) */
						_idx0=((((((((((p_idx_z+1)*x_max)+(2*p_idx_z))+2)*y_max)+((((2*p_idx_z)+p_idx_y)+3)*x_max))+(4*p_idx_z))+(2*p_idx_y))+p_idx_x)+8);
						/* _idx1 = ((((((((((p_idx_z+1)*x_max)+(2*p_idx_z))+2)*y_max)+((((2*p_idx_z)+p_idx_y)+3)*x_max))+(4*p_idx_z))+(2*p_idx_y))+p_idx_x)+6) */
						_idx1=(_idx0-2);
						/* _idx2 = ((((((((((p_idx_z+1)*x_max)+(2*p_idx_z))+2)*y_max)+((((2*p_idx_z)+p_idx_y)+4)*x_max))+(4*p_idx_z))+(2*p_idx_y))+p_idx_x)+9) */
						_idx2=((_idx1+x_max)+3);
						/* _idx3 = ((((((((((p_idx_z+1)*x_max)+(2*p_idx_z))+2)*y_max)+((((2*p_idx_z)+p_idx_y)+2)*x_max))+(4*p_idx_z))+(2*p_idx_y))+p_idx_x)+5) */
						_idx3=((_idx1-x_max)-1);
						/* _idx4 = ((((((((((p_idx_z+2)*x_max)+(2*p_idx_z))+4)*y_max)+((((2*p_idx_z)+p_idx_y)+5)*x_max))+(4*p_idx_z))+(2*p_idx_y))+p_idx_x)+11) */
						_idx4=(((_idx3+((x_max+2)*y_max))+(3*x_max))+6);
						/* _idx5 = ((((((((p_idx_z*x_max)+(2*p_idx_z))*y_max)+((((2*p_idx_z)+p_idx_y)+1)*x_max))+(4*p_idx_z))+(2*p_idx_y))+p_idx_x)+3) */
						_idx5=(((_idx3+((( - x_max)-2)*y_max))-x_max)-2);
						/* _idx6 = ((((((((((p_idx_z+1)*x_max)+(2*p_idx_z))+2)*y_max)+((((2*p_idx_z)+p_idx_y)+3)*x_max))+(4*p_idx_z))+(2*p_idx_y))+p_idx_x)+7) */
						_idx6=((_idx3+x_max)+2);
						u_0_1[_idx6]=((((u_0_0[_idx0]+(u_0_0[_idx1]+u_0_0[_idx2]))+(u_0_0[_idx3]+(u_0_0[_idx4]+u_0_0[_idx5])))*0.25)-u_0_0[_idx6]);
						
						
						/* debugging */
						/* u_0_0[count] = p_idx_x + p_idx_y*x_max + p_idx_z*x_max*y_max; */
						/* count++; */
						
						/* u_0_1[p_idx_x + p_idx_y*x_max + p_idx_z*x_max*y_max] = _idx3; */
						/* u_0_1[p_idx_x + p_idx_y*x_max + p_idx_z*x_max] = p_idx_x + p_idx_y*x_max + p_idx_z*x_max; */
						/* end debugging */
						
						/* printf("u_0_0[_idx0]=%f,u_0_0[_idx2]=%f,u_0_0[_idx1]=%f\n",u_0_0[_idx0],u_0_0[_idx2],u_0_0[_idx1]); */
						
					}
				}
			}
		}
	}
	*u_0_1_out = u_0_1;
}
