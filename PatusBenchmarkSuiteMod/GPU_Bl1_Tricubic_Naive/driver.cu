#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include <cutil.h>

#include "patusrt.h"

// forward_decls -->
__global__ void initialize(double *  u_0_0, double *  u_0_1, double *  a_1_0, double *  b_2_0, double *  c_3_0, int x_max, int y_max, int z_max, int cbx);
__global__ void tricubic_interpolation(double *  *  u_0_1_out, double *  u_0_0, double *  u_0_1, double *  a_1_0, double *  b_2_0, double *  c_3_0, int x_max, int y_max, int z_max, int cbx);

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
	double *  a_1_0;
	double *  b_2_0;
	double *  c_3_0;
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
	u_0_0=((double * )malloc(((((x_max+6)*(y_max+6))*(z_max+6))*sizeof (double))));
	u_0_1=((double * )malloc(((((x_max+6)*(y_max+6))*(z_max+6))*sizeof (double))));
	a_1_0=((double * )malloc((((x_max*y_max)*z_max)*sizeof (double))));
	b_2_0=((double * )malloc((((x_max*y_max)*z_max)*sizeof (double))));
	c_3_0=((double * )malloc((((x_max*y_max)*z_max)*sizeof (double))));
	// <--
	
	
	// declare_GPU_grids -->
	double *  u_0_1_out_gpu;
	double *  u_0_0_gpu;
	double *  u_0_1_gpu;
	double *  a_1_0_gpu;
	double *  b_2_0_gpu;
	double *  c_3_0_gpu;
	dim3 thds(cbx, 1, 1);
	dim3 blks((x_max/cbx), (y_max*z_max), 1);
	// <--
	
	// allocate_GPU_grids -->
	cudaMalloc(((void *  * )( & a_1_0_gpu)), (((x_max*y_max)*z_max)*sizeof (double)));
	cudaMalloc(((void *  * )( & c_3_0_gpu)), (((x_max*y_max)*z_max)*sizeof (double)));
	cudaMalloc(((void *  * )( & u_0_0_gpu)), ((((x_max+6)*(y_max+6))*(z_max+6))*sizeof (double)));
	cudaMalloc(((void *  * )( & b_2_0_gpu)), (((x_max*y_max)*z_max)*sizeof (double)));
	cudaMalloc(((void *  * )( & u_0_1_out_gpu)), ((((x_max+6)*(y_max+6))*(z_max+6))*sizeof (double * )));
	cudaMalloc(((void *  * )( & u_0_1_gpu)), ((((x_max+6)*(y_max+6))*(z_max+6))*sizeof (double)));
	// <--
	
	// copy_grids_to_GPU -->
	cudaMemcpy(((void * )a_1_0_gpu), ((void * )a_1_0), (((x_max*y_max)*z_max)*sizeof (double)), cudaMemcpyHostToDevice);
	cudaMemcpy(((void * )c_3_0_gpu), ((void * )c_3_0), (((x_max*y_max)*z_max)*sizeof (double)), cudaMemcpyHostToDevice);
	cudaMemcpy(((void * )u_0_0_gpu), ((void * )u_0_0), ((((x_max+6)*(y_max+6))*(z_max+6))*sizeof (double)), cudaMemcpyHostToDevice);
	cudaMemcpy(((void * )b_2_0_gpu), ((void * )b_2_0), (((x_max*y_max)*z_max)*sizeof (double)), cudaMemcpyHostToDevice);
	cudaMemcpy(((void * )u_0_1_gpu), ((void * )u_0_1), ((((x_max+6)*(y_max+6))*(z_max+6))*sizeof (double)), cudaMemcpyHostToDevice);
	// <--
	
	
	// initialize_grids -->
	initialize<<<blks, thds>>>(u_0_0_gpu, u_0_1_gpu, a_1_0_gpu, b_2_0_gpu, c_3_0_gpu, x_max, y_max, z_max, cbx);
	// <--
	
	cudaThreadSynchronize ();
	res = cudaGetLastError ();
	if (res != cudaSuccess)
	{
		printf ("CUDA Error [Initialization]: %s.\n", cudaGetErrorString (res));
	}
	
	long nFlopsPerStencil = 318;
	long nGridPointsCount = 5 * ((x_max*y_max)*z_max);
	long nBytesTransferred = 5 * (((((((x_max*y_max)*z_max)*sizeof (double))+(((x_max*y_max)*z_max)*sizeof (double)))+((((x_max+6)*(y_max+6))*(z_max+6))*sizeof (double)))+(((x_max*y_max)*z_max)*sizeof (double)))+(((x_max*y_max)*z_max)*sizeof (double)));
	
	// warm up
	// compute_stencil -->
	tricubic_interpolation<<<blks, thds>>>(( & u_0_1_out_gpu), u_0_0_gpu, u_0_1_gpu, a_1_0_gpu, b_2_0_gpu, c_3_0_gpu, x_max, y_max, z_max, cbx);
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
		tricubic_interpolation<<<blks, thds>>>(( & u_0_1_out_gpu), u_0_0_gpu, u_0_1_gpu, a_1_0_gpu, b_2_0_gpu, c_3_0_gpu, x_max, y_max, z_max, cbx);
		// <--
		
		cudaThreadSynchronize ();
	}
	toc (nFlopsPerStencil, nGridPointsCount, nBytesTransferred);
	
	// free memory
	// deallocate_grids -->
	cudaFree(((void * )a_1_0_gpu));
	cudaFree(((void * )c_3_0_gpu));
	cudaFree(((void * )u_0_0_gpu));
	cudaFree(((void * )b_2_0_gpu));
	cudaFree(((void * )u_0_1_out_gpu));
	cudaFree(((void * )u_0_1_gpu));
	free(u_0_0);
	free(u_0_1);
	free(a_1_0);
	free(b_2_0);
	free(c_3_0);
	// <--
	
	
	cudaThreadExit ();
	return EXIT_SUCCESS;
}
