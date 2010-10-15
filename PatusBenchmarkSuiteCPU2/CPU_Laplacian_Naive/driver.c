#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <omp.h>
#include "patusrt.h"

// forward_decls -->
void initialize(float *  u_0_0, float *  u_0_1, int x_max, int y_max, int z_max);
void laplacian(float *  *  u_0_1_out, float *  u_0_0, float *  u_0_1, int x_max, int y_max, int z_max, int t_max);
void laplacian_cpu(float *  *  u_0_1_out, float *  u_0_0, float *  u_0_1, int x_max, int y_max, int z_max,int t_max);
// <--


int main (int argc, char** argv)
{
	int i;
	
	// prepare grids
	// declare_grids -->
	float *  u_0_1_out;
	float *  u_0_0;
	float *  u_0_1;
	float *  u_0_1_out_cpu;
	float *  u_0_0_cpu;
	float *  u_0_1_cpu;
	int iterations = 100;
	if ((argc<4))
	{
		printf("Wrong number of parameters. Syntax:\n%s <x_max> <y_max> <z_max>\n", argv[0]);
		exit(-1);
	}
	int x_max = atoi(argv[1]);
	int y_max = atoi(argv[2]);
	int z_max = atoi(argv[3]);
	if(argc==5) iterations = atoi(argv[4]);
	// <--
	
	// allocate_grids -->
	u_0_0=((float * )malloc(((((x_max+2)*(y_max+2))*(z_max+2))*sizeof (float))));
	u_0_1=((float * )malloc(((((x_max+2)*(y_max+2))*(z_max+2))*sizeof (float))));
	u_0_0_cpu=((float * )malloc(((((x_max+2)*(y_max+2))*(z_max+2))*sizeof (float))));
	u_0_1_cpu=((float * )malloc(((((x_max+2)*(y_max+2))*(z_max+2))*sizeof (float))));
	// <--
	
	
	// initialize

	// initialize_grids -->
	initialize(u_0_0, u_0_1, x_max, y_max, z_max);
	initialize(u_0_0_cpu, u_0_1_cpu, x_max, y_max, z_max);
	// <--
	
	long nFlopsPerStencil = 7;
	long nGridPointsCount = iterations * ((x_max*y_max)*z_max);
	long nBytesTransferred = iterations * (((((x_max+2)*(y_max+2))*(z_max+2))*sizeof (float))+(((x_max*y_max)*z_max)*sizeof (float)));
	

	// compute_stencil -->
	laplacian(( & u_0_1_out), u_0_0, u_0_1, x_max, y_max, z_max,iterations);
	// <--			
	// run the benchmark
	tic ();
	// compute_stencil -->
	laplacian(( & u_0_1_out), u_0_0, u_0_1, x_max, y_max, z_max,iterations);
	// <--
	toc (nFlopsPerStencil, nGridPointsCount, nBytesTransferred);

	/* **************************************** naive cpu single core **************************************** */
	// compute_stencil -->
	laplacian_cpu(( & u_0_1_out_cpu), u_0_0_cpu, u_0_1_cpu, x_max, y_max, z_max,iterations);
	// <--			
	// run the benchmark
	tic ();
	// compute_stencil -->
	laplacian_cpu(( & u_0_1_out_cpu), u_0_0_cpu, u_0_1_cpu, x_max, y_max, z_max,iterations);
	// <--
	toc (nFlopsPerStencil, nGridPointsCount, nBytesTransferred);

	// checking "correctness" (assuming cpu version is correct)
	int error_count=0;
	for(i=0;i<(x_max)*(y_max)*(z_max);i++) {
	  if(fabs(u_0_1_out[i] - u_0_1_out_cpu[i])>0.001) {
	    error_count++;
	    printf("%dth error encountered at u[%d]: |%f-%f|=%5.16f\n",error_count,i,u_0_1_out[i],u_0_1_out_cpu[i],fabs(u_0_1_out[i] - u_0_1_out_cpu[i]));
	    if(error_count>30) {
	      printf("too many errors\n"); exit(1);
	    }

	  }
	}
	// free memory
	// deallocate_grids -->
	free(u_0_0);
	free(u_0_1);
	free(u_0_0_cpu);
	free(u_0_1_cpu);
	// <--
	
	
	return EXIT_SUCCESS;
}
