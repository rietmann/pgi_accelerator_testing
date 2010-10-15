#include <stdio.h>
#include <stdlib.h>

#include <omp.h>
#include "patusrt.h"

// forward_decls -->
void initialize(double *  u_0_0, double *  u_0_1, double *  a_1_0, double *  b_2_0, double *  c_3_0, int x_max, int y_max, int z_max);
void tricubic_interpolation(double *  *  u_0_1_out, double *  u_0_0, double *  u_0_1, double *  a_1_0, double *  b_2_0, double *  c_3_0, int x_max, int y_max, int z_max,int t_max);
void tricubic_interpolation_cpu(double *  *  u_0_1_out, double *  u_0_0, double *  u_0_1, double *  a_1_0, double *  b_2_0, double *  c_3_0, int x_max, int y_max, int z_max, int t_max);

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
	
	double *  u_0_1_out_cpu;
	double *  u_0_0_cpu;
	double *  u_0_1_cpu;
	double *  a_1_0_cpu;
	double *  b_2_0_cpu;
	double *  c_3_0_cpu;
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
	u_0_0=((double * )malloc(((((x_max+3)*(y_max+3))*(z_max+3))*sizeof (double))));
	u_0_1=((double * )malloc(((((x_max+3)*(y_max+3))*(z_max+3))*sizeof (double))));
	a_1_0=((double * )malloc((((x_max*y_max)*z_max)*sizeof (double))));
	b_2_0=((double * )malloc((((x_max*y_max)*z_max)*sizeof (double))));
	c_3_0=((double * )malloc((((x_max*y_max)*z_max)*sizeof (double))));
	
	u_0_0_cpu=((double * )malloc(((((x_max+3)*(y_max+3))*(z_max+3))*sizeof (double))));
	u_0_1_cpu=((double * )malloc(((((x_max+3)*(y_max+3))*(z_max+3))*sizeof (double))));
	a_1_0_cpu=((double * )malloc((((x_max*y_max)*z_max)*sizeof (double))));
	b_2_0_cpu=((double * )malloc((((x_max*y_max)*z_max)*sizeof (double))));
	c_3_0_cpu=((double * )malloc((((x_max*y_max)*z_max)*sizeof (double))));
	// <--
	
	
	// initialize
	// initialize_grids -->
	initialize(u_0_0, u_0_1, a_1_0, b_2_0, c_3_0, x_max, y_max, z_max);
	initialize(u_0_0_cpu, u_0_1_cpu, a_1_0_cpu, b_2_0_cpu, c_3_0_cpu, x_max, y_max, z_max);
	// <--
	
	long nFlopsPerStencil = 318;
	long nGridPointsCount = iterations * ((x_max*y_max)*z_max);
	long nBytesTransferred = iterations * ((((((((x_max+3)*(y_max+3))*(z_max+3))*sizeof (double))+(((x_max*y_max)*z_max)*sizeof (double)))+(((x_max*y_max)*z_max)*sizeof (double)))+(((x_max*y_max)*z_max)*sizeof (double)))+(((x_max*y_max)*z_max)*sizeof (double)));
	
	// warm up
	// compute_stencil -->
	tricubic_interpolation(( & u_0_1_out), u_0_0, u_0_1, a_1_0, b_2_0, c_3_0, x_max, y_max, z_max,iterations);
	// <--			
	
	// run the benchmark
	tic ();
	// compute_stencil -->
	tricubic_interpolation(( & u_0_1_out), u_0_0, u_0_1, a_1_0, b_2_0, c_3_0, x_max, y_max, z_max,iterations);
	// <--		
	toc (nFlopsPerStencil, nGridPointsCount, nBytesTransferred);


	/* **************************************** naive single core cpu **************************************** */
	// warm up
	// compute_stencil -->
	tricubic_interpolation_cpu(( & u_0_1_out_cpu), u_0_0_cpu, u_0_1_cpu, a_1_0_cpu, b_2_0_cpu, c_3_0_cpu, x_max, y_max, z_max,iterations);
	// <--			
	
	// run the benchmark
	tic ();
	// compute_stencil -->
	tricubic_interpolation_cpu(( & u_0_1_out_cpu), u_0_0_cpu, u_0_1_cpu, a_1_0_cpu, b_2_0_cpu, c_3_0_cpu, x_max, y_max, z_max,iterations);
	// <--		
	toc (nFlopsPerStencil, nGridPointsCount, nBytesTransferred);

	// checking "correctness" (assuming cpu version is correct)
	int error_count=0;
	for(i=0;i<(x_max)*(y_max)*(z_max);i++) {
	  if(abs(u_0_1_out[i] - u_0_1_out_cpu[i])>0.001) {
	    error_count++;
	    printf("%dth error encountered at u[%d]: |%f-%f|=%5.16f\n",error_count,i,u_0_1_out[i],u_0_1_out_cpu[i],abs(u_0_1_out[i] - u_0_1_out_cpu[i]));
	    if(error_count>30) {
	      printf("too many errors\n"); exit(1);
	    }

	  }
	}
	
	// free memory
	// deallocate_grids -->
	free(u_0_0);
	free(u_0_1);
	free(a_1_0);
	free(b_2_0);
	free(c_3_0);

	free(u_0_0_cpu);
	free(u_0_1_cpu);
	free(a_1_0_cpu);
	free(b_2_0_cpu);
	free(c_3_0_cpu);
	// <--
	
	
	return EXIT_SUCCESS;
}
