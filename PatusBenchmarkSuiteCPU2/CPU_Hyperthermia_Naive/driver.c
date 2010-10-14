#include <stdio.h>
#include <stdlib.h>

#include <omp.h>
#include "patusrt.h"

// forward_decls -->
void initialize(float *  T_0_0, float *  T_0_1, float *  c_1_0, float *  c_2_0, float *  c_3_0, float *  c_4_0, float *  c_5_0, float *  c_6_0, float *  c_7_0, float *  c_8_0, float *  c_9_0, int x_max, int y_max, int z_max);
void hyperthermia(float *  *  T_0_1_out, float *  T_0_0, float *  T_0_1, float *  c_1_0, float *  c_2_0, float *  c_3_0, float *  c_4_0, float *  c_5_0, float *  c_6_0, float *  c_7_0, float *  c_8_0, float *  c_9_0, int x_max, int y_max, int z_max);

// <--


int main (int argc, char** argv)
{
	int i;
	
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
	T_0_0=((float * )malloc(((((x_max+2)*(y_max+2))*(z_max+2))*sizeof (float))));
	T_0_1=((float * )malloc(((((x_max+2)*(y_max+2))*(z_max+2))*sizeof (float))));
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
	
	
	// initialize
#pragma omp parallel
	{
		// initialize_grids -->
		initialize(T_0_0, T_0_1, c_1_0, c_2_0, c_3_0, c_4_0, c_5_0, c_6_0, c_7_0, c_8_0, c_9_0, x_max, y_max, z_max);
		// <--
		
	}
	
	long nFlopsPerStencil = 16;
	long nGridPointsCount = 5 * ((x_max*y_max)*z_max);
	long nBytesTransferred = 5 * (((((((((((((x_max*y_max)*z_max)*sizeof (float))+((((x_max+2)*(y_max+2))*(z_max+2))*sizeof (float)))+(((x_max*y_max)*z_max)*sizeof (float)))+(((x_max*y_max)*z_max)*sizeof (float)))+(((x_max*y_max)*z_max)*sizeof (float)))+(((x_max*y_max)*z_max)*sizeof (float)))+(((x_max*y_max)*z_max)*sizeof (float)))+(((x_max*y_max)*z_max)*sizeof (float)))+(((x_max*y_max)*z_max)*sizeof (float)))+(((x_max*y_max)*z_max)*sizeof (float)))+(((x_max*y_max)*z_max)*sizeof (float)));
	
	// warm up
#pragma omp parallel
	{
		// compute_stencil -->
		hyperthermia(( & T_0_1_out), T_0_0, T_0_1, c_1_0, c_2_0, c_3_0, c_4_0, c_5_0, c_6_0, c_7_0, c_8_0, c_9_0, x_max, y_max, z_max);
		// <--
		
	}
	
	// run the benchmark
	tic ();
#pragma omp parallel private(i)
	for (i = 0; i < 5; i++)
	{
		// compute_stencil -->
		hyperthermia(( & T_0_1_out), T_0_0, T_0_1, c_1_0, c_2_0, c_3_0, c_4_0, c_5_0, c_6_0, c_7_0, c_8_0, c_9_0, x_max, y_max, z_max);
		// <--
		
#pragma omp barrier
	}
	toc (nFlopsPerStencil, nGridPointsCount, nBytesTransferred);
	
	// free memory
	// deallocate_grids -->
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
	
	
	return EXIT_SUCCESS;
}
