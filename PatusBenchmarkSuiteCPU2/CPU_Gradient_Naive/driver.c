#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <omp.h>
#include "patusrt.h"

// forward_decls -->
void initialize(float *  u_0_0, float *  ux_1_0, float *  uy_2_0, float *  uz_3_0, float alpha, float beta, float gamma, int x_max, int y_max, int z_max);
void gradient(float *  *  ux_1_0_out, float *  *  uy_2_0_out, float *  *  uz_3_0_out, float *  u_0_0, float *  ux_1_0, float *  uy_2_0, float *  uz_3_0, float alpha, float beta, float gamma, int x_max, int y_max, int z_max,int t_max);
void gradient_cpu(float *  *  ux_1_0_out, float *  *  uy_2_0_out, float *  *  uz_3_0_out, float *  u_0_0, float *  ux_1_0, float *  uy_2_0, float *  uz_3_0, float alpha, float beta, float gamma, int x_max, int y_max, int z_max,int t_max);

// <--


int main (int argc, char** argv)
{
	int i;
	
	// prepare grids
	// declare_grids -->
	int iterations = 100;
	float *  ux_1_0_out;
	float *  uy_2_0_out;
	float *  uz_3_0_out;
	float *  ux_1_0_out_cpu;
	float *  uy_2_0_out_cpu;
	float *  uz_3_0_out_cpu;

	float *  u_0_0;
	float *  ux_1_0;
	float *  uy_2_0;
	float *  uz_3_0;
	float *  u_0_0_cpu;
	float *  ux_1_0_cpu;
	float *  uy_2_0_cpu;
	float *  uz_3_0_cpu;
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
	u_0_0=((float * )malloc(((((x_max+2)*(y_max+2))*(z_max+2))*sizeof (float))));
	ux_1_0=((float * )malloc((((x_max*y_max)*z_max)*sizeof (float))));
	uy_2_0=((float * )malloc((((x_max*y_max)*z_max)*sizeof (float))));
	uz_3_0=((float * )malloc((((x_max*y_max)*z_max)*sizeof (float))));
	u_0_0_cpu=((float * )malloc(((((x_max+2)*(y_max+2))*(z_max+2))*sizeof (float))));
	ux_1_0_cpu=((float * )malloc((((x_max*y_max)*z_max)*sizeof (float))));
	uy_2_0_cpu=((float * )malloc((((x_max*y_max)*z_max)*sizeof (float))));
	uz_3_0_cpu=((float * )malloc((((x_max*y_max)*z_max)*sizeof (float))));
	// <--
	
	
	// initialize
	printf("initializing\n");
	// initialize_grids -->
	initialize(u_0_0, ux_1_0, uy_2_0, uz_3_0, 0.1, 0.2, 0.30000000000000004, x_max, y_max, z_max);
	initialize(u_0_0_cpu, ux_1_0_cpu, uy_2_0_cpu, uz_3_0_cpu, 0.1, 0.2, 0.30000000000000004, x_max, y_max, z_max);
	// <--		
	printf("done initializing\n");
	
	long nFlopsPerStencil = 6;
	long nGridPointsCount = iterations * ((x_max*y_max)*z_max);
	long nBytesTransferred = iterations * (((((x_max+2)*(y_max+2))*(z_max+2))*sizeof (float))+(((((x_max*y_max)*z_max)*sizeof (float))+(((x_max*y_max)*z_max)*sizeof (float)))+(((x_max*y_max)*z_max)*sizeof (float))));
	

	/* ****************************** PGI GPU benchmark ****************************** */
	// warm up
	// compute_stencil -->
	gradient(( & ux_1_0_out), ( & uy_2_0_out), ( & uz_3_0_out), u_0_0, ux_1_0, uy_2_0, uz_3_0, 0.4, 0.5, 0.6, x_max, y_max, z_max,iterations);
	// <--
	printf("starting benchmark loop\n");
	// run the benchmark
	tic ();
	// compute_stencil -->
	gradient(( & ux_1_0_out), ( & uy_2_0_out), ( & uz_3_0_out), u_0_0, ux_1_0, uy_2_0, uz_3_0, 0.7, 0.7999999999999999, 0.8999999999999999, x_max, y_max, z_max,iterations);
		
	toc (nFlopsPerStencil, nGridPointsCount, nBytesTransferred);


	/* ****************************** naive single thread cpu benchmark ****************************** */
	// warm up

	// compute_stencil -->
	printf("starting naive cpu benchmark loop\n");
	gradient_cpu(( & ux_1_0_out_cpu), ( & uy_2_0_out_cpu), ( & uz_3_0_out_cpu), u_0_0_cpu, ux_1_0_cpu, uy_2_0_cpu, uz_3_0_cpu, 0.4, 0.5, 0.6, x_max, y_max, z_max,iterations);
	// <--
	// run the benchmark
	tic ();
	// compute_stencil -->
	gradient_cpu(( & ux_1_0_out_cpu), ( & uy_2_0_out_cpu), ( & uz_3_0_out_cpu), u_0_0_cpu, ux_1_0_cpu, uy_2_0_cpu, uz_3_0_cpu, 0.7, 0.7999999999999999, 0.8999999999999999, x_max, y_max, z_max,iterations);
		
	toc (nFlopsPerStencil, nGridPointsCount, nBytesTransferred);

	// checking "correctness" (assuming cpu version is correct)
	int error_count=0;
	for(i=0;i<(x_max)*(y_max)*(z_max);i++) {
	  if(fabs(ux_1_0_out[i] - ux_1_0_out_cpu[i])>0.001
	     && fabs(uy_2_0_out[i] - uy_2_0_out_cpu[i])>0.001
	     && fabs(uz_3_0_out[i] - uz_3_0_out_cpu[i])>0.001) {
	    error_count++;
	    printf("%dth error encountered at u[%d]: |%f-%f|=%5.16f\n",error_count,i,ux_1_0_out[i],ux_1_0_out_cpu[i],fabs(ux_1_0_out[i] - ux_1_0_out_cpu[i]));
	    if(error_count>30) {
	      printf("too many errors\n"); exit(1);
	    }

	  }
	}
	printf("No errors found in u_i\n");
	// free memory
	// deallocate_grids -->
	free(u_0_0);
	free(ux_1_0);
	free(uy_2_0);
	free(uz_3_0);
	free(u_0_0_cpu);
	free(ux_1_0_cpu);
	free(uy_2_0_cpu);
	free(uz_3_0_cpu);
	// <--
	
	
	return EXIT_SUCCESS;
}
