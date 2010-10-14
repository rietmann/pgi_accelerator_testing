#include <stdio.h>
#include <stdlib.h>

#include <omp.h>
#include "patusrt.h"

// forward_decls -->
void initialize(double *  u_0_0, double *  u_0_1, double *  a_1_0, double *  b_2_0, double *  c_3_0, int x_max, int y_max, int z_max);
void tricubic_interpolation(double *  *  u_0_1_out, double *  u_0_0, double *  u_0_1, double *  a_1_0, double *  b_2_0, double *  c_3_0, int x_max, int y_max, int z_max);

// <--


int main (int argc, char** argv)
{
	int i;
	
	// prepare grids
	// declare_grids -->
	double *  u_0_1_out;
	double *  u_0_0;
	double *  u_0_1;
	double *  a_1_0;
	double *  b_2_0;
	double *  c_3_0;
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
	u_0_0=((double * )malloc(((((x_max+6)*(y_max+6))*(z_max+6))*sizeof (double))));
	u_0_1=((double * )malloc(((((x_max+6)*(y_max+6))*(z_max+6))*sizeof (double))));
	a_1_0=((double * )malloc((((x_max*y_max)*z_max)*sizeof (double))));
	b_2_0=((double * )malloc((((x_max*y_max)*z_max)*sizeof (double))));
	c_3_0=((double * )malloc((((x_max*y_max)*z_max)*sizeof (double))));
	// <--
	
	
	// initialize
#pragma omp parallel private(i)
	{
		// initialize_grids -->
		initialize(u_0_0, u_0_1, a_1_0, b_2_0, c_3_0, x_max, y_max, z_max);
		// <--
#pragma omp barrier
	
	long nFlopsPerStencil = 318;
	long nGridPointsCount = 5 * ((x_max*y_max)*z_max);
	long nBytesTransferred = 5 * (((((((x_max*y_max)*z_max)*sizeof (double))+((((x_max+6)*(y_max+6))*(z_max+6))*sizeof (double)))+(((x_max*y_max)*z_max)*sizeof (double)))+(((x_max*y_max)*z_max)*sizeof (double)))+(((x_max*y_max)*z_max)*sizeof (double)));
	
	// warm up
#pragma omp parallel
	{
		// compute_stencil -->
		tricubic_interpolation(( & u_0_1_out), u_0_0, u_0_1, a_1_0, b_2_0, c_3_0, x_max, y_max, z_max);
		// <--
		
	}
	
#pragma omp barrier

	// run the benchmark
#pragma omp master
	tic ();

	for (i = 0; i < 5; i++)
	{
		// compute_stencil -->
		tricubic_interpolation(( & u_0_1_out), u_0_0, u_0_1, a_1_0, b_2_0, c_3_0, x_max, y_max, z_max);
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
	free(a_1_0);
	free(b_2_0);
	free(c_3_0);
	// <--
	
	
	return EXIT_SUCCESS;
}
