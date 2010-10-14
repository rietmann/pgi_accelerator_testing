#include <stdio.h>
#include <stdlib.h>

#include <omp.h>
#include "patusrt.h"

// forward_decls -->
void initialize(float *  u_0_0, float *  ux_1_0, float *  uy_2_0, float *  uz_3_0, float alpha, float beta, float gamma, int x_max, int y_max, int z_max);
void divergence(float *  *  u_0_0_out, float *  u_0_0, float *  ux_1_0, float *  uy_2_0, float *  uz_3_0, float alpha, float beta, float gamma, int x_max, int y_max, int z_max);

// <--


int main (int argc, char** argv)
{
	int i;
	
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
	
	
	// initialize
#pragma omp parallel private(i)
	{
		// initialize_grids -->
		initialize(u_0_0, ux_1_0, uy_2_0, uz_3_0, 0.1, 0.2, 0.30000000000000004, x_max, y_max, z_max);
		// <--
		
	
	long nFlopsPerStencil = 8;
	long nGridPointsCount = 5 * ((x_max*y_max)*z_max);
	long nBytesTransferred = 5 * (((((((x_max+4)*y_max)*z_max)*sizeof (float))+(((x_max*(y_max+4))*z_max)*sizeof (float)))+(((x_max*y_max)*(z_max+4))*sizeof (float)))+(((x_max*y_max)*z_max)*sizeof (float)));
	
	// warm up
#pragma omp barrier
	{
		// compute_stencil -->
		divergence(( & u_0_0_out), u_0_0, ux_1_0, uy_2_0, uz_3_0, 0.4, 0.5, 0.6, x_max, y_max, z_max);
		// <--
		
	}
	
#pragma omp barrier

	// run the benchmark
#pragma omp master
	tic ();
	for (i = 0; i < 5; i++)
	{
		// compute_stencil -->
		divergence(( & u_0_0_out), u_0_0, ux_1_0, uy_2_0, uz_3_0, 0.7, 0.7999999999999999, 0.8999999999999999, x_max, y_max, z_max);
		// <--
		
#pragma omp barrier
	}
#pragma omp master
	toc (nFlopsPerStencil, nGridPointsCount, nBytesTransferred);
	}

	// free memory
	// deallocate_grids -->
	free(u_0_0);
	free(ux_1_0);
	free(uy_2_0);
	free(uz_3_0);
	// <--
	
	
	return EXIT_SUCCESS;
}
