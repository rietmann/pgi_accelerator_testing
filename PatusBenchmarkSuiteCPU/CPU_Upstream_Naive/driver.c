#include <stdio.h>
#include <stdlib.h>

#include <omp.h>
#include "patusrt.h"

// forward_decls -->
void initialize(double *  u_0_0, double *  u_0_1, double a, int x_max, int y_max, int z_max);
void upstream_5_3d(double *  *  u_0_1_out, double *  u_0_0, double *  u_0_1, double a, int x_max, int y_max, int z_max);

// <--


int main (int argc, char** argv)
{
	int i;
	
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
	
	
	// initialize
#pragma omp parallel private(i)
	{
		// initialize_grids -->
		initialize(u_0_0, u_0_1, 0.1, x_max, y_max, z_max);
		// <--

	long nFlopsPerStencil = 22;
	long nGridPointsCount = 5 * ((x_max*y_max)*z_max);
	long nBytesTransferred = 5 * (((((x_max+10)*(y_max+10))*(z_max+10))*sizeof (double))+(((x_max*y_max)*z_max)*sizeof (double)));

	// warm up
		// compute_stencil -->
		upstream_5_3d(( & u_0_1_out), u_0_0, u_0_1, 0.2, x_max, y_max, z_max);
		// <--
		
	
	// run the benchmark
#pragma omp master
	tic ();

	for (i = 0; i < 5; i++)
	{
		// compute_stencil -->
		upstream_5_3d(( & u_0_1_out), u_0_0, u_0_1, 0.30000000000000004, x_max, y_max, z_max);
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
