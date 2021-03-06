#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include <cutil.h>

#include "patusrt.h"

// forward_decls -->
__global__ void initialize(float *  u_0_0, float *  ux_1_0, float *  uy_2_0, float *  uz_3_0, float alpha, float beta, float gamma, int x_max, int y_max, int z_max);
__global__ void divergence(float *  *  u_0_0_out, float *  u_0_0, float *  ux_1_0, float *  uy_2_0, float *  uz_3_0, float alpha, float beta, float gamma, int x_max, int y_max, int z_max);

// <--


int main (int argc, char** argv)
{
	int i;
	cudaError_t res;
	
	// prepare grids
	// declare_grids -->
	float *  u_0_0_out;
	float *  u_0_0;
	float *  ux_1_0;
	float *  uy_2_0;
	float *  uz_3_0;
	if ((argc!=4))
	{
		printf("Wrong number of parameters. Syntax:\n%s <x_max> <y_max> <z_max>\n", argv[0]);
		exit(-1);
	}
	int x_max = atoi(argv[1]);
	int y_max = atoi(argv[2]);
	int z_max = atoi(argv[3]);
	// <--
	
	// allocate_grids -->
	u_0_0=((float * )malloc((((x_max*y_max)*z_max)*sizeof (float))));
	ux_1_0=((float * )malloc(((((x_max+4)*y_max)*z_max)*sizeof (float))));
	uy_2_0=((float * )malloc((((x_max*(y_max+4))*z_max)*sizeof (float))));
	uz_3_0=((float * )malloc((((x_max*y_max)*(z_max+4))*sizeof (float))));
	// <--
	
	
	// declare_GPU_grids -->
	float *  u_0_0_out_gpu;
	float *  u_0_0_gpu;
	float *  ux_1_0_gpu;
	float *  uy_2_0_gpu;
	float *  uz_3_0_gpu;
	dim3 thds(1, 1, 1);
	dim3 blks(x_max, (y_max*z_max), 1);
	// <--
	
	// allocate_GPU_grids -->
	cudaMalloc(((void *  * )( & u_0_0_gpu)), (((x_max*y_max)*z_max)*sizeof (float)));
	cudaMalloc(((void *  * )( & u_0_0_out_gpu)), (((x_max*y_max)*z_max)*sizeof (float * )));
	cudaMalloc(((void *  * )( & ux_1_0_gpu)), ((((x_max+4)*y_max)*z_max)*sizeof (float)));
	cudaMalloc(((void *  * )( & uy_2_0_gpu)), (((x_max*(y_max+4))*z_max)*sizeof (float)));
	cudaMalloc(((void *  * )( & uz_3_0_gpu)), (((x_max*y_max)*(z_max+4))*sizeof (float)));
	// <--
	
	// copy_grids_to_GPU -->
	cudaMemcpy(((void * )u_0_0_gpu), ((void * )u_0_0), (((x_max*y_max)*z_max)*sizeof (float)), cudaMemcpyHostToDevice);
	cudaMemcpy(((void * )ux_1_0_gpu), ((void * )ux_1_0), ((((x_max+4)*y_max)*z_max)*sizeof (float)), cudaMemcpyHostToDevice);
	cudaMemcpy(((void * )uy_2_0_gpu), ((void * )uy_2_0), (((x_max*(y_max+4))*z_max)*sizeof (float)), cudaMemcpyHostToDevice);
	cudaMemcpy(((void * )uz_3_0_gpu), ((void * )uz_3_0), (((x_max*y_max)*(z_max+4))*sizeof (float)), cudaMemcpyHostToDevice);
	// <--
	
	
	// initialize_grids -->
	initialize<<<blks, thds>>>(u_0_0_gpu, ux_1_0_gpu, uy_2_0_gpu, uz_3_0_gpu, 0.1, 0.2, 0.30000000000000004, x_max, y_max, z_max);
	// <--
	
	cudaThreadSynchronize ();
	res = cudaGetLastError ();
	if (res != cudaSuccess)
	{
		printf ("CUDA Error [Initialization]: %s.\n", cudaGetErrorString (res));
	}
	
	long nFlopsPerStencil = 8;
	long nGridPointsCount = 5 * ((x_max*y_max)*z_max);
	long nBytesTransferred = 5 * ((((((x_max*(y_max+4))*z_max)*sizeof (float))+((((x_max+4)*y_max)*z_max)*sizeof (float)))+(((x_max*y_max)*(z_max+4))*sizeof (float)))+(((x_max*y_max)*z_max)*sizeof (float)));
	
	// warm up
	// compute_stencil -->
	divergence<<<blks, thds>>>(( & u_0_0_out_gpu), u_0_0_gpu, ux_1_0_gpu, uy_2_0_gpu, uz_3_0_gpu, 0.4, 0.5, 0.6, x_max, y_max, z_max);
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
		divergence<<<blks, thds>>>(( & u_0_0_out_gpu), u_0_0_gpu, ux_1_0_gpu, uy_2_0_gpu, uz_3_0_gpu, 0.7, 0.7999999999999999, 0.8999999999999999, x_max, y_max, z_max);
		// <--
		
		cudaThreadSynchronize ();
	}
	toc (nFlopsPerStencil, nGridPointsCount, nBytesTransferred);
	
	// free memory
	// deallocate_grids -->
	cudaFree(((void * )u_0_0_gpu));
	cudaFree(((void * )u_0_0_out_gpu));
	cudaFree(((void * )ux_1_0_gpu));
	cudaFree(((void * )uy_2_0_gpu));
	cudaFree(((void * )uz_3_0_gpu));
	free(u_0_0);
	free(ux_1_0);
	free(uy_2_0);
	free(uz_3_0);
	// <--
	
	
	cudaThreadExit ();
	return EXIT_SUCCESS;
}
