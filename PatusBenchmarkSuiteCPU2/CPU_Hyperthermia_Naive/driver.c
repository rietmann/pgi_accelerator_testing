#include <stdio.h>
#include <stdlib.h>

#include <omp.h>
#include "patusrt.h"

// forward_decls -->
void initialize(float *  T_0_0, float *  T_0_1, float *  c_1_0, float *  c_2_0, float *  c_3_0, float *  c_4_0, float *  c_5_0, float *  c_6_0, float *  c_7_0, float *  c_8_0, float *  c_9_0, int x_max, int y_max, int z_max);
void hyperthermia(float *  *  T_0_1_out, float *  T_0_0, float *  T_0_1, float *  c_1_0, float *  c_2_0, float *  c_3_0, float *  c_4_0, float *  c_5_0, float *  c_6_0, float *  c_7_0, float *  c_8_0, float *  c_9_0, int x_max, int y_max, int z_max,int t_max);
void hyperthermia_cpu(float *  *  T_0_1_out, float *  T_0_0, float *  T_0_1, float *  c_1_0, float *  c_2_0, float *  c_3_0, float *  c_4_0, float *  c_5_0, float *  c_6_0, float *  c_7_0, float *  c_8_0, float *  c_9_0, int x_max, int y_max, int z_max,int t_max);

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

	float *  T_0_1_out_cpu;
	float *  T_0_0_cpu;
	float *  T_0_1_cpu;
	float *  c_1_0_cpu;
	float *  c_2_0_cpu;
	float *  c_3_0_cpu;
	float *  c_4_0_cpu;
	float *  c_5_0_cpu;
	float *  c_6_0_cpu;
	float *  c_7_0_cpu;
	float *  c_8_0_cpu;
	float *  c_9_0_cpu;

	int iterations = 100;
	
	if ((argc<4))
	{
		printf("Wrong number of parameters. Syntax:\n%s <x_max> <y_max> <z_max>\n", argv[0]);
		exit(-1);
	}
	
	int x_max = atoi(argv[1]);
	int y_max = atoi(argv[2]);
	int z_max = atoi(argv[3]);
	if(argc==5)
	  iterations = atoi(argv[4]);
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

	T_0_0_cpu=((float * )malloc(((((x_max+2)*(y_max+2))*(z_max+2))*sizeof (float))));
	T_0_1_cpu=((float * )malloc(((((x_max+2)*(y_max+2))*(z_max+2))*sizeof (float))));
	c_1_0_cpu=((float * )malloc((((x_max*y_max)*z_max)*sizeof (float))));
	c_2_0_cpu=((float * )malloc((((x_max*y_max)*z_max)*sizeof (float))));
	c_3_0_cpu=((float * )malloc((((x_max*y_max)*z_max)*sizeof (float))));
	c_4_0_cpu=((float * )malloc((((x_max*y_max)*z_max)*sizeof (float))));
	c_5_0_cpu=((float * )malloc((((x_max*y_max)*z_max)*sizeof (float))));
	c_6_0_cpu=((float * )malloc((((x_max*y_max)*z_max)*sizeof (float))));
	c_7_0_cpu=((float * )malloc((((x_max*y_max)*z_max)*sizeof (float))));
	c_8_0_cpu=((float * )malloc((((x_max*y_max)*z_max)*sizeof (float))));
	c_9_0_cpu=((float * )malloc((((x_max*y_max)*z_max)*sizeof (float))));
	// <--
	
	
	// initialize
#pragma omp parallel
	{
		// initialize_grids -->
		initialize(T_0_0, T_0_1, c_1_0, c_2_0, c_3_0, c_4_0, c_5_0, c_6_0, c_7_0, c_8_0, c_9_0, x_max, y_max, z_max);
		initialize(T_0_0_cpu, T_0_1_cpu, c_1_0_cpu, c_2_0_cpu, c_3_0_cpu, c_4_0_cpu, c_5_0_cpu, c_6_0_cpu, c_7_0_cpu, c_8_0_cpu, c_9_0_cpu, x_max, y_max, z_max);
		// <--
		
	}
	
	long nFlopsPerStencil = 16;
	long nGridPointsCount = iterations * ((x_max*y_max)*z_max);
	long nBytesTransferred = iterations * (((((((((((((x_max*y_max)*z_max)*sizeof (float))+((((x_max+2)*(y_max+2))*(z_max+2))*sizeof (float)))+(((x_max*y_max)*z_max)*sizeof (float)))+(((x_max*y_max)*z_max)*sizeof (float)))+(((x_max*y_max)*z_max)*sizeof (float)))+(((x_max*y_max)*z_max)*sizeof (float)))+(((x_max*y_max)*z_max)*sizeof (float)))+(((x_max*y_max)*z_max)*sizeof (float)))+(((x_max*y_max)*z_max)*sizeof (float)))+(((x_max*y_max)*z_max)*sizeof (float)))+(((x_max*y_max)*z_max)*sizeof (float)));

	/* **************************************** benchmark pgi gpu accelerator **************************************** */
	
	// warm up
	
	// compute_stencil -->
	hyperthermia(( & T_0_1_out), T_0_0, T_0_1, c_1_0, c_2_0, c_3_0, c_4_0, c_5_0, c_6_0, c_7_0, c_8_0, c_9_0, x_max, y_max, z_max,iterations);
	// <--
	// run the benchmark
	tic ();
	// compute_stencil -->
	hyperthermia(( & T_0_1_out), T_0_0, T_0_1, c_1_0, c_2_0, c_3_0, c_4_0, c_5_0, c_6_0, c_7_0, c_8_0, c_9_0, x_max, y_max, z_max,iterations);
	// <--
	toc (nFlopsPerStencil, nGridPointsCount, nBytesTransferred);

	/* **************************************** benchmark naive cpu **************************************** */
	
	// warm up 

	// compute_stencil -->
	hyperthermia_cpu(( & T_0_1_out_cpu), T_0_0_cpu, T_0_1_cpu, c_1_0_cpu, c_2_0_cpu, c_3_0_cpu, c_4_0_cpu, c_5_0_cpu, c_6_0_cpu, c_7_0_cpu, c_8_0_cpu, c_9_0_cpu, x_max, y_max, z_max,iterations);
	// <--
	// run the benchmark
	tic ();
	// compute_stencil -->
	hyperthermia_cpu(( & T_0_1_out), T_0_0_cpu, T_0_1_cpu, c_1_0_cpu, c_2_0_cpu, c_3_0_cpu, c_4_0_cpu, c_5_0_cpu, c_6_0_cpu, c_7_0_cpu, c_8_0_cpu, c_9_0_cpu, x_max, y_max, z_max,iterations);
	// <--
	toc (nFlopsPerStencil, nGridPointsCount, nBytesTransferred);


	int error_count=0;
	for(i=0;i<(x_max)*(y_max)*(z_max);i++) {
	  if(abs(T_0_1_out[i] - T_0_1_out_cpu[i])>0.001) {
	    error_count++;
	    printf("%dth error encountered at u[%d]: |%f-%f|=%5.16f\n",error_count,i,T_0_1_out[i],T_0_1_out_cpu[i],abs(T_0_1_out[i] - T_0_1_out_cpu[i]));
	    if(error_count>30) {
	      printf("too many errors\n"); exit(1);
	    }

	  }
	}
	
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
	
	free(T_0_0_cpu);
	free(T_0_1_cpu);
	free(c_1_0_cpu);
	free(c_2_0_cpu);
	free(c_3_0_cpu);
	free(c_4_0_cpu);
	free(c_5_0_cpu);
	free(c_6_0_cpu);
	free(c_7_0_cpu);
	free(c_8_0_cpu);
	free(c_9_0_cpu);
	
	// <--
	
	
	return EXIT_SUCCESS;
}
