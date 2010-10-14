#include <stdio.h>
#include <stdlib.h>

#include <omp.h>
#include "patusrt.h"

// forward_decls -->
void initialize(float *  u_0_0, float *  u_0_1, int x_max, int y_max, int z_max);
void laplacian(float *  *  u_0_1_out, float *  u_0_0, float *  u_0_1, int x_max, int y_max, int z_max);

// <--


int main (int argc, char** argv)
{
	int i;
	
	// prepare grids
	// declare_grids -->
	float *  u_0_1_out;
	float *  u_0_0;
	float *  u_0_1;
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
	u_0_0=((float * )malloc(((((x_max+4)*(y_max+4))*(z_max+4))*sizeof (float))));
	u_0_1=((float * )malloc(((((x_max+4)*(y_max+4))*(z_max+4))*sizeof (float))));
	// <--
	
	
	// initialize
#pragma omp parallel private(i)
	{
		// initialize_grids -->
		initialize(u_0_0, u_0_1, x_max, y_max, z_max);
		// <--
		
	
	long nFlopsPerStencil = 7;
	long nGridPointsCount = 5 * ((x_max*y_max)*z_max);
	long nBytesTransferred = 5 * (((((x_max+4)*(y_max+4))*(z_max+4))*sizeof (float))+(((x_max*y_max)*z_max)*sizeof (float)));
	
	// warm up
#pragma omp barrier
	{
		// compute_stencil -->
		laplacian(( & u_0_1_out), u_0_0, u_0_1, x_max, y_max, z_max);
		// <--
		
	}
	
#pragma omp barrier

	// run the benchmark
#pragma omp master
	tic ();
	for (i = 0; i < 5; i++)
	{
		// compute_stencil -->
		laplacian(( & u_0_1_out), u_0_0, u_0_1, x_max, y_max, z_max);
		// <--
		
#pragma omp barrier
	}
#pragma omp master
	toc (nFlopsPerStencil, nGridPointsCount, nBytesTransferred);
	}

	// free memory
	// deallocate_grids -->
	free(u_0_0);
	free(u_0_1);
	// <--
	
	
	return EXIT_SUCCESS;
}
