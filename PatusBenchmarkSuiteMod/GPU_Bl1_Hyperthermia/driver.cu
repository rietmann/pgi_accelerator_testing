#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include <cutil.h>

#include "patusrt.h"

// forward_decls -->
__global__ void initialize(float *  T_0_0, float *  T_0_1, float *  c_1_0, float *  c_2_0, float *  c_3_0, float *  c_4_0, float *  c_5_0, float *  c_6_0, float *  c_7_0, float *  c_8_0, float *  c_9_0, int x_max, int y_max, int z_max, int cbx);
__global__ void hyperthermia(float *  *  T_0_1_out, float *  T_0_0, float *  T_0_1, float *  c_1_0, float *  c_2_0, float *  c_3_0, float *  c_4_0, float *  c_5_0, float *  c_6_0, float *  c_7_0, float *  c_8_0, float *  c_9_0, int x_max, int y_max, int z_max, int cbx);

// <--


int main (int argc, char** argv)
{
	int i;
	cudaError_t res;
	
	// prepare grids
	// declare_grids -->
	float *  T_0_1_out;
	float *  T_0_0;
	float *  T_0_1;
	float *  c_1_0;
	float *  c_2_0;
	float *  c_3_0;
	float *  c_4_0;
	float *  c_5_0;
	float *  c_6_0;
	float *  c_7_0;
	float *  c_8_0;
	float *  c_9_0;
	if ((argc!=5))
	{
		printf("Wrong number of parameters. Syntax:\n%s <x_max> <y_max> <z_max> <cbx>\n", argv[0]);
		exit(-1);
	}
	int x_max = atoi(argv[1]);
	int y_max = atoi(argv[2]);
	int z_max = atoi(argv[3]);
	int cbx = atoi(argv[4]);
	// <--
	
	// allocate_grids -->
	T_0_0=((float * )malloc(((((x_max+4)*(y_max+4))*(z_max+4))*sizeof (float))));
	T_0_1=((float * )malloc(((((x_max+4)*(y_max+4))*(z_max+4))*sizeof (float))));
	c_1_0=((float * )malloc((((x_max*y_max)*z_max)*sizeof (float))));
	c_2_0=((float * )malloc((((x_max*y_max)*z_max)*sizeof (float))));
	c_3_0=((float * )malloc((((x_max*y_max)*z_max)*sizeof (float))));
	c_4_0=((float * )malloc((((x_max*y_max)*z_max)*sizeof (float))));
	c_5_0=((float * )malloc((((x_max*y_max)*z_max)*sizeof (float))));
	c_6_0=((float * )malloc((((x_max*y_max)*z_max)*sizeof (float))));
	c_7_0=((float * )malloc((((x_max*y_max)*z_max)*sizeof (float))));
	c_8_0=((float * )malloc((((x_max*y_max)*z_max)*sizeof (float))));
	c_9_0=((float * )malloc((((x_max*y_max)*z_max)*sizeof (float))));
	// <--
	
	
	// declare_GPU_grids -->
	float *  T_0_1_out_gpu;
	float *  T_0_0_gpu;
	float *  T_0_1_gpu;
	float *  c_1_0_gpu;
	float *  c_2_0_gpu;
	float *  c_3_0_gpu;
	float *  c_4_0_gpu;
	float *  c_5_0_gpu;
	float *  c_6_0_gpu;
	float *  c_7_0_gpu;
	float *  c_8_0_gpu;
	float *  c_9_0_gpu;
	dim3 thds(cbx, 1, 1);
	dim3 blks((x_max/cbx), (y_max*z_max), 1);
	// <--
	
	// allocate_GPU_grids -->
	cudaMalloc(((void *  * )( & c_9_0_gpu)), (((x_max*y_max)*z_max)*sizeof (float)));
	cudaMalloc(((void *  * )( & c_5_0_gpu)), (((x_max*y_max)*z_max)*sizeof (float)));
	cudaMalloc(((void *  * )( & c_8_0_gpu)), (((x_max*y_max)*z_max)*sizeof (float)));
	cudaMalloc(((void *  * )( & c_6_0_gpu)), (((x_max*y_max)*z_max)*sizeof (float)));
	cudaMalloc(((void *  * )( & c_1_0_gpu)), (((x_max*y_max)*z_max)*sizeof (float)));
	cudaMalloc(((void *  * )( & c_4_0_gpu)), (((x_max*y_max)*z_max)*sizeof (float)));
	cudaMalloc(((void *  * )( & c_3_0_gpu)), (((x_max*y_max)*z_max)*sizeof (float)));
	cudaMalloc(((void *  * )( & T_0_0_gpu)), ((((x_max+4)*(y_max+4))*(z_max+4))*sizeof (float)));
	cudaMalloc(((void *  * )( & T_0_1_out_gpu)), ((((x_max+4)*(y_max+4))*(z_max+4))*sizeof (float * )));
	cudaMalloc(((void *  * )( & c_2_0_gpu)), (((x_max*y_max)*z_max)*sizeof (float)));
	cudaMalloc(((void *  * )( & c_7_0_gpu)), (((x_max*y_max)*z_max)*sizeof (float)));
	cudaMalloc(((void *  * )( & T_0_1_gpu)), ((((x_max+4)*(y_max+4))*(z_max+4))*sizeof (float)));
	// <--
	
	// copy_grids_to_GPU -->
	cudaMemcpy(((void * )c_9_0_gpu), ((void * )c_9_0), (((x_max*y_max)*z_max)*sizeof (float)), cudaMemcpyHostToDevice);
	cudaMemcpy(((void * )c_5_0_gpu), ((void * )c_5_0), (((x_max*y_max)*z_max)*sizeof (float)), cudaMemcpyHostToDevice);
	cudaMemcpy(((void * )c_8_0_gpu), ((void * )c_8_0), (((x_max*y_max)*z_max)*sizeof (float)), cudaMemcpyHostToDevice);
	cudaMemcpy(((void * )c_6_0_gpu), ((void * )c_6_0), (((x_max*y_max)*z_max)*sizeof (float)), cudaMemcpyHostToDevice);
	cudaMemcpy(((void * )c_1_0_gpu), ((void * )c_1_0), (((x_max*y_max)*z_max)*sizeof (float)), cudaMemcpyHostToDevice);
	cudaMemcpy(((void * )c_4_0_gpu), ((void * )c_4_0), (((x_max*y_max)*z_max)*sizeof (float)), cudaMemcpyHostToDevice);
	cudaMemcpy(((void * )c_3_0_gpu), ((void * )c_3_0), (((x_max*y_max)*z_max)*sizeof (float)), cudaMemcpyHostToDevice);
	cudaMemcpy(((void * )T_0_0_gpu), ((void * )T_0_0), ((((x_max+4)*(y_max+4))*(z_max+4))*sizeof (float)), cudaMemcpyHostToDevice);
	cudaMemcpy(((void * )c_2_0_gpu), ((void * )c_2_0), (((x_max*y_max)*z_max)*sizeof (float)), cudaMemcpyHostToDevice);
	cudaMemcpy(((void * )c_7_0_gpu), ((void * )c_7_0), (((x_max*y_max)*z_max)*sizeof (float)), cudaMemcpyHostToDevice);
	cudaMemcpy(((void * )T_0_1_gpu), ((void * )T_0_1), ((((x_max+4)*(y_max+4))*(z_max+4))*sizeof (float)), cudaMemcpyHostToDevice);
	// <--
	
	
	// initialize_grids -->
	initialize<<<blks, thds>>>(T_0_0_gpu, T_0_1_gpu, c_1_0_gpu, c_2_0_gpu, c_3_0_gpu, c_4_0_gpu, c_5_0_gpu, c_6_0_gpu, c_7_0_gpu, c_8_0_gpu, c_9_0_gpu, x_max, y_max, z_max, cbx);
	// <--
	
	cudaThreadSynchronize ();
	res = cudaGetLastError ();
	if (res != cudaSuccess)
	{
		printf ("CUDA Error [Initialization]: %s.\n", cudaGetErrorString (res));
	}
	
	long nFlopsPerStencil = 16;
	long nGridPointsCount = 5 * ((x_max*y_max)*z_max);
	long nBytesTransferred = 5 * (((((((((((((x_max*y_max)*z_max)*sizeof (float))+(((x_max*y_max)*z_max)*sizeof (float)))+(((x_max*y_max)*z_max)*sizeof (float)))+(((x_max*y_max)*z_max)*sizeof (float)))+(((x_max*y_max)*z_max)*sizeof (float)))+(((x_max*y_max)*z_max)*sizeof (float)))+(((x_max*y_max)*z_max)*sizeof (float)))+(((x_max*y_max)*z_max)*sizeof (float)))+(((x_max*y_max)*z_max)*sizeof (float)))+((((x_max+4)*(y_max+4))*(z_max+4))*sizeof (float)))+(((x_max*y_max)*z_max)*sizeof (float)));
	
	// warm up
	// compute_stencil -->
	hyperthermia<<<blks, thds>>>(( & T_0_1_out_gpu), T_0_0_gpu, T_0_1_gpu, c_1_0_gpu, c_2_0_gpu, c_3_0_gpu, c_4_0_gpu, c_5_0_gpu, c_6_0_gpu, c_7_0_gpu, c_8_0_gpu, c_9_0_gpu, x_max, y_max, z_max, cbx);
	// <--
	
	cudaThreadSynchronize ();
	res = cudaGetLastError ();
	if (res != cudaSuccess)
	{
		printf ("CUDA Error [Stencil]: %s.\n", cudaGetErrorString (res));
	}
	
	// run the benchmark
	tic ();
	for (i = 0; i < 5; i++)
	{
		// compute_stencil -->
		hyperthermia<<<blks, thds>>>(( & T_0_1_out_gpu), T_0_0_gpu, T_0_1_gpu, c_1_0_gpu, c_2_0_gpu, c_3_0_gpu, c_4_0_gpu, c_5_0_gpu, c_6_0_gpu, c_7_0_gpu, c_8_0_gpu, c_9_0_gpu, x_max, y_max, z_max, cbx);
		// <--
		
		cudaThreadSynchronize ();
	}
	toc (nFlopsPerStencil, nGridPointsCount, nBytesTransferred);
	
	// free memory
	// deallocate_grids -->
	cudaFree(((void * )c_9_0_gpu));
	cudaFree(((void * )c_5_0_gpu));
	cudaFree(((void * )c_8_0_gpu));
	cudaFree(((void * )c_6_0_gpu));
	cudaFree(((void * )c_1_0_gpu));
	cudaFree(((void * )c_4_0_gpu));
	cudaFree(((void * )c_3_0_gpu));
	cudaFree(((void * )T_0_0_gpu));
	cudaFree(((void * )T_0_1_out_gpu));
	cudaFree(((void * )c_2_0_gpu));
	cudaFree(((void * )c_7_0_gpu));
	cudaFree(((void * )T_0_1_gpu));
	free(T_0_0);
	free(T_0_1);
	free(c_1_0);
	free(c_2_0);
	free(c_3_0);
	free(c_4_0);
	free(c_5_0);
	free(c_6_0);
	free(c_7_0);
	free(c_8_0);
	free(c_9_0);
	// <--
	
	
	cudaThreadExit ();
	return EXIT_SUCCESS;
}
