#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "patusrt.h"

// forward_decls -->
void initialize(float *  u_0_0, float *  u_0_1, int x_max, int y_max, int z_max);
void laplacian(float *  *  u_0_1_out, float *  u_0_0, float *  u_0_1, int x_max, int y_max, int z_max);
void laplacian_cpu(float *  *  u_0_1_out, float *  u_0_0, float *  u_0_1, int x_max, int y_max, int z_max);

// <--


int main (int argc, char** argv)
{
	int i;
	int error_count;
	// prepare grids
	// declare_grids -->
	float *  u_0_1_out;
	float *  u_0_0;
	float *  u_0_1;
	float *  u_0_1_out_test;
	float *  u_0_0_test;
	float *  u_0_1_test;
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
	u_0_0_test=((float * )malloc(((((x_max+4)*(y_max+4))*(z_max+4))*sizeof (float))));
	u_0_1_test=((float * )malloc(((((x_max+4)*(y_max+4))*(z_max+4))*sizeof (float))));
	// <--
	
	
	// initialize
	initialize(u_0_0, u_0_1, x_max, y_max, z_max);
	initialize(u_0_0_test, u_0_1_test, x_max, y_max, z_max);
	
	long nFlopsPerStencil = 7;
	long nGridPointsCount = 100*((x_max*y_max)*z_max);
	long nBytesTransferred = 100*(((((x_max+4)*(y_max+4))*(z_max+4))*sizeof (float))+(((x_max*y_max)*z_max)*sizeof (float)));
	
	// warm up
	// compute_stencil -->
	laplacian(( & u_0_1_out), u_0_0, u_0_1, x_max, y_max, z_max);
	// <--
		
	// run the benchmark
	tic ();
	for (i = 0; i < 1; i++)
	{
		// compute_stencil -->
		laplacian(( & u_0_1_out), u_0_0, u_0_1, x_max, y_max, z_max);
		// <--
	}
	toc (nFlopsPerStencil, nGridPointsCount, nBytesTransferred);
	
	
	// warm up test
	// compute_stencil test -->
	laplacian_cpu(( & u_0_1_out_test), u_0_0_test, u_0_1_test, x_max, y_max, z_max);
	// <--
		
	// run the purely CPU benchmark test
	tic ();
	for (i = 0; i < 1; i++)
	{
		// compute_stencil -->
		laplacian_cpu(( & u_0_1_out_test), u_0_0_test, u_0_1_test, x_max, y_max, z_max);
		// <--
	}
	toc (nFlopsPerStencil, nGridPointsCount, nBytesTransferred);

	// checking "correctness" (assuming cpu version is correct)
	error_count=0;
	for(i=0;i<(x_max+4)*(y_max+4)*(z_max+4);i++) {
	  if(fabs(u_0_1_out[i] - u_0_1_out_test[i])>0.001) {
	    error_count++;
	    printf("%dth error encountered at u[%d]: |%f-%f|=%5.16f\n",error_count,i,u_0_1_out[i],u_0_1_out_test[i],fabs(u_0_1_out[i] - u_0_1_out_test[i]));
	    if(error_count>30) {
	      printf("too many errors\n"); exit(1);
	    }

	  }
	}
	
	// free memory
	// deallocate_grids -->
	free(u_0_0);
	free(u_0_1);
	free(u_0_0_test);
	free(u_0_1_test);
	// <--
	
	
	return EXIT_SUCCESS;
}
