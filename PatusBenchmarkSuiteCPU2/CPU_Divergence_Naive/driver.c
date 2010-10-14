#include <stdio.h>
#include <stdlib.h>

#include <omp.h>
#include "patusrt.h"

// forward_decls -->
void initialize(float *  u_0_0, float *  ux_1_0, float *  uy_2_0, float *  uz_3_0, float alpha, float beta, float gamma, int x_max, int y_max, int z_max);
void divergence(float *  *  u_0_0_out, float *  u_0_0, float *  ux_1_0, float *  uy_2_0, float *  uz_3_0, float alpha, float beta, float gamma, int x_max, int y_max, int z_max,int t_max);
void divergence_cpu(float *  *  u_0_0_out, float *  u_0_0, float *  ux_1_0, float *  uy_2_0, float *  uz_3_0, float alpha, float beta, float gamma, int x_max, int y_max, int z_max,int t_max);
// <--


int main (int argc, char** argv)
{
	int i;
	int iterations = 100;
	// prepare grids
	// declare_grids -->
	float *  u_0_0_out;
	float *  u_0_0;
	float *  ux_1_0;
	float *  uy_2_0;
	float *  uz_3_0;
	float *  u_0_0_out_cpu;
	float *  u_0_0_cpu;
	float *  ux_1_0_cpu;
	float *  uy_2_0_cpu;
	float *  uz_3_0_cpu;
	if ((argc<4))
	{
		printf("Wrong number of parameters. Syntax:\n%s <x_max> <y_max> <z_max> <# of iterations>\n", argv[0]);
		exit(-1);
	}
	int x_max = atoi(argv[1]);
	int y_max = atoi(argv[2]);
	int z_max = atoi(argv[3]);
	if(argc==5)
	  iterations = atoi(argv[4]);
	// <--
	
	// allocate_grids -->
	u_0_0=((float * )malloc((((x_max*y_max)*z_max)*sizeof (float))));
	ux_1_0=((float * )malloc(((((x_max+2)*y_max)*z_max)*sizeof (float))));
	uy_2_0=((float * )malloc((((x_max*(y_max+2))*z_max)*sizeof (float))));
	uz_3_0=((float * )malloc((((x_max*y_max)*(z_max+2))*sizeof (float))));
	u_0_0_cpu=((float * )malloc((((x_max*y_max)*z_max)*sizeof (float))));
	ux_1_0_cpu=((float * )malloc(((((x_max+2)*y_max)*z_max)*sizeof (float))));
	uy_2_0_cpu=((float * )malloc((((x_max*(y_max+2))*z_max)*sizeof (float))));
	uz_3_0_cpu=((float * )malloc((((x_max*y_max)*(z_max+2))*sizeof (float))));
	// <--
	
	
	// initialize
	// initialize_grids -->
	initialize(u_0_0, ux_1_0, uy_2_0, uz_3_0, 0.1, 0.2, 0.30000000000000004, x_max, y_max, z_max);
	initialize(u_0_0_cpu, ux_1_0_cpu, uy_2_0_cpu, uz_3_0_cpu, 0.1, 0.2, 0.30000000000000004, x_max, y_max, z_max);
	// <--
	
	long nFlopsPerStencil = 8;
	long nGridPointsCount = iterations * ((x_max*y_max)*z_max);
	long nBytesTransferred = iterations * (((((((x_max+2)*y_max)*z_max)*sizeof (float))+(((x_max*(y_max+2))*z_max)*sizeof (float)))+(((x_max*y_max)*(z_max+2))*sizeof (float)))+(((x_max*y_max)*z_max)*sizeof (float)));
	

	/* *************************** PGI GPU-acc benchmark ********************* */

	
	// warm up
	
	{
		// compute_stencil -->
	  divergence(( & u_0_0_out), u_0_0, ux_1_0, uy_2_0, uz_3_0, 0.4, 0.5, 0.6, x_max, y_max, z_max,iterations);
		// <--
	}
	
	// run the benchmark
	tic ();
	{
		// compute_stencil -->
	  divergence(( & u_0_0_out), u_0_0, ux_1_0, uy_2_0, uz_3_0, 0.7, 0.7999999999999999, 0.8999999999999999, x_max, y_max, z_max,iterations);
		// <--	
	}
	toc (nFlopsPerStencil, nGridPointsCount, nBytesTransferred);

	
	/* *************************** ******************** ********************* */	

	/* *************************** Naive CPU Comparison ********************* */
	
	// warm up cpu comparison

	{
		// compute_stencil -->
	  divergence(( & u_0_0_out_cpu), u_0_0_cpu, ux_1_0_cpu, uy_2_0_cpu, uz_3_0_cpu, 0.4, 0.5, 0.6, x_max, y_max, z_max,iterations);
		// <--
	}
	
	// run the benchmark
	tic ();
	{
		// compute_stencil -->
	  divergence_cpu(( & u_0_0_out_cpu), u_0_0_cpu, ux_1_0_cpu, uy_2_0_cpu, uz_3_0_cpu, 0.7, 0.7999999999999999, 0.8999999999999999, x_max, y_max, z_max,iterations);
		// <--	
	}
	toc (nFlopsPerStencil, nGridPointsCount, nBytesTransferred);


	// checking "correctness" (assuming cpu version is correct)
	int error_count=0;
	for(i=0;i<(x_max)*(y_max)*(z_max);i++) {
	  if(abs(u_0_0_out[i] - u_0_0_out_cpu[i])>0.001) {
	    error_count++;
	    printf("%dth error encountered at u[%d]: |%f-%f|=%5.16f\n",error_count,i,u_0_0_out[i],u_0_0_out_cpu[i],abs(u_0_0_out[i] - u_0_0_out_cpu[i]));
	    if(error_count>30) {
	      printf("too many errors\n"); exit(1);
	    }

	  }
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
