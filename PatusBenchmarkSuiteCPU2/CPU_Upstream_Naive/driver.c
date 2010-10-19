#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <omp.h>
#include "patusrt.h"

// forward_decls -->
void initialize(double *  u_0_0, double *  u_0_1, double a, int x_max, int y_max, int z_max);
void upstream_5_3d(double *  *  u_0_1_out, double *  u_0_0, double *  u_0_1, double a, int x_max, int y_max, int z_max,int t_max);
void upstream_5_3d_cpu(double *  *  u_0_1_out, double *  u_0_0, double *  u_0_1, double a, int x_max, int y_max, int z_max, int t_max);

// <--


int main (int argc, char** argv)
{
	int i;
	
	// prepare grids
	// declare_grids -->
	double *  u_0_1_out;
	double *  u_0_0;
	double *  u_0_1;

	double *  u_0_1_out_cpu;
	double *  u_0_0_cpu;
	double *  u_0_1_cpu;

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
	u_0_0=((double * )malloc(((((x_max+5)*(y_max+5))*(z_max+5))*sizeof (double))));
	u_0_1=((double * )malloc(((((x_max+5)*(y_max+5))*(z_max+5))*sizeof (double))));

	u_0_0_cpu=((double * )malloc(((((x_max+5)*(y_max+5))*(z_max+5))*sizeof (double))));
	u_0_1_cpu=((double * )malloc(((((x_max+5)*(y_max+5))*(z_max+5))*sizeof (double))));
	// <--
	
	
	// initialize
	// initialize_grids -->
	initialize(u_0_0, u_0_1, 0.1, x_max, y_max, z_max);
	initialize(u_0_0_cpu, u_0_1_cpu, 0.1, x_max, y_max, z_max);
	// <--
		
	
	
	long nFlopsPerStencil = 22;
	long nGridPointsCount = iterations * ((x_max*y_max)*z_max);
	long nBytesTransferred = iterations * (((((x_max+5)*(y_max+5))*(z_max+5))*sizeof (double))+(((x_max*y_max)*z_max)*sizeof (double)));
	
	// warm up
	// compute_stencil -->
	upstream_5_3d(( & u_0_1_out), u_0_0, u_0_1, 0.2, x_max, y_max, z_max,iterations);
	// <--			
	
	// run the benchmark
	tic ();
	// compute_stencil -->
	upstream_5_3d(( & u_0_1_out), u_0_0, u_0_1, 0.30000000000000004, x_max, y_max, z_max,iterations);
	// <--
	toc (nFlopsPerStencil, nGridPointsCount, nBytesTransferred);

	/* **************************************** naive cpu version **************************************** */
	// warm up
	// compute_stencil -->
	upstream_5_3d_cpu(( & u_0_1_out_cpu), u_0_0_cpu, u_0_1_cpu, 0.2, x_max, y_max, z_max,iterations);
	// <--			
	
	// run the benchmark
	tic ();
	// compute_stencil -->
	upstream_5_3d_cpu(( & u_0_1_out_cpu), u_0_0_cpu, u_0_1_cpu, 0.30000000000000004, x_max, y_max, z_max,iterations);
	// <--
	toc (nFlopsPerStencil, nGridPointsCount, nBytesTransferred);

	int error_count=0;
	for(i=0;i<(x_max)*(y_max)*(z_max);i++) {
	  double error = fabs(u_0_1_out[i] - u_0_1_out_cpu[i]);
	  if(error>0.0001) {
	    error_count++;
	    printf("%dth error encountered at u[%d]: |%f-%f|=%f\n",error_count,i,u_0_1_out[i],u_0_1_out_cpu[i],error);
	    if(error_count>30) {
	      printf("too many errors\n"); exit(1);
	    }

	  }
	}
	printf("Error check successful\n");
	
	// free memory
	// deallocate_grids -->
	free(u_0_0);
	free(u_0_1);
	free(u_0_0_cpu);
	free(u_0_1_cpu);
	// <--
	
	
	return EXIT_SUCCESS;
}
