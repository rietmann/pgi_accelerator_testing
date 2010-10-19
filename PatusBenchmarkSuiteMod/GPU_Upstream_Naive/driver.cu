#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include <cutil.h>

#include "patusrt.h"

// forward_decls -->
__global__ void initialize(double *  u_0_0, double *  u_0_1, double a, int x_max, int y_max, int z_max);
__global__ void upstream_5_3d(double *  *  u_0_1_out, double *  u_0_0, double *  u_0_1, double a, int x_max, int y_max, int z_max);

// <--


int main (int argc, char** argv)
{
	int i;
	cudaError_t res;
	
	// prepare grids
	// declare_grids -->
	double *  u_0_1_out;
	double *  u_0_0;
	double *  u_0_1;
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
	u_0_0=((double * )malloc(((((x_max+10)*(y_max+10))*(z_max+10))*sizeof (double))));
	u_0_1=((double * )malloc(((((x_max+10)*(y_max+10))*(z_max+10))*sizeof (double))));
	// <--
	
	
	// declare_GPU_grids -->
	double *  u_0_1_out_gpu;
	double *  u_0_0_gpu;
	double *  u_0_1_gpu;
	dim3 thds(1, 1, 1);
	dim3 blks(x_max, (y_max*z_max), 1);
	// <--
	
	// allocate_GPU_grids -->
	cudaMalloc(((void *  * )( & u_0_1_gpu)), ((((x_max+10)*(y_max+10))*(z_max+10))*sizeof (double)));
	cudaMalloc(((void *  * )( & u_0_0_gpu)), ((((x_max+10)*(y_max+10))*(z_max+10))*sizeof (double)));
	cudaMalloc(((void *  * )( & u_0_1_out_gpu)), ((((x_max+10)*(y_max+10))*(z_max+10))*sizeof (double * )));
	// <--
	
	// copy_grids_to_GPU -->
	cudaMemcpy(((void * )u_0_1_gpu), ((void * )u_0_1), ((((x_max+10)*(y_max+10))*(z_max+10))*sizeof (double)), cudaMemcpyHostToDevice);
	cudaMemcpy(((void * )u_0_0_gpu), ((void * )u_0_0), ((((x_max+10)*(y_max+10))*(z_max+10))*sizeof (double)), cudaMemcpyHostToDevice);
	// <--
	
	
	// initialize_grids -->
	initialize<<<blks, thds>>>(u_0_0_gpu, u_0_1_gpu, 0.1, x_max, y_max, z_max);
	// <--
	
	cudaThreadSynchronize ();
	res = cudaGetLastError ();
	if (res != cudaSuccess)
	{
		printf ("CUDA Error [Initialization]: %s.\n", cudaGetErrorString (res));
	}
	
	long nFlopsPerStencil = 22;
	long nGridPointsCount = 5 * ((x_max*y_max)*z_max);
	long nBytesTransferred = 5 * (((((x_max+10)*(y_max+10))*(z_max+10))*sizeof (double))+(((x_max*y_max)*z_max)*sizeof (double)));
	
	// warm up
	// compute_stencil -->
	upstream_5_3d<<<blks, thds>>>(( & u_0_1_out_gpu), u_0_0_gpu, u_0_1_gpu, 0.2, x_max, y_max, z_max);
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
		upstream_5_3d<<<blks, thds>>>(( & u_0_1_out_gpu), u_0_0_gpu, u_0_1_gpu, 0.30000000000000004, x_max, y_max, z_max);
		// <--
		
		cudaThreadSynchronize ();
	}
	toc (nFlopsPerStencil, nGridPointsCount, nBytesTransferred);
	
	// free memory
	// deallocate_grids -->
	cudaFree(((void * )u_0_1_gpu));
	cudaFree(((void * )u_0_0_gpu));
	cudaFree(((void * )u_0_1_out_gpu));
	free(u_0_0);
	free(u_0_1);
	// <--
	
	
	cudaThreadExit ();
	return EXIT_SUCCESS;
}
