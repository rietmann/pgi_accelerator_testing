#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* compiled with */
/* pgcc -lm -Minfo -Mbounds -ta=nvidia -o test_pgi simple_kernel.c */

void test_doubly_nested_loop_pgi_acc(float *restrict u_out, float *restrict u, int x_max, int y_max) {

  int x,y;
#pragma acc region copyin(u[0:x_max*y_max]) copy(u_out[0:x_max*y_max])
  {
    /* because we use dynamically allocated memory, we had to
       explictly tell compiler that loops are indeed parallel. */
    /* this is an example stencil calculation, so we don't calculate
       for values on the edge of the 2D square, which is why the loop
       indeces go from 1:N-1 */
#pragma acc for independent
    for (y=1; y<y_max-1; y++) {
#pragma acc for independent
      for(x=1; x<x_max-1; x++) {
	
	/* Apparent computer BUG: this line does not give correct
	   results. It seems to delete x_max*y and just use x */

	int id = x + x_max*y;
	
	/* workaround: */
	/* If we split the previous line into two arithmetic operations, it works correctly */
	/* int id = x_max*y; */
	/* id = x + id;		 */
	
	/* for testing, we just store the calculated id in its self-calculated position. fails if "id" was generated in a multiply add all in one line */

	u_out[id] = id;

	/*	
	result should be:
	u[0] = -42; <-- first edge of square
	u[1] = -42; ...
	u[2] = -42; ...
	...
	u[9] =  9 <-- first non-edge point
	u[10]= 10 
	u[11]= 11
	...	
	u[x_max*y_max-1] = -42;
	*/	

	/* failure output looks like: (GTX280 @ CSCS.ch) */
	/* pgaccelinfo provided at bottom of code */
	
	/* u_out_pgi[0]=-42.00,
	/* u_out_pgi[1]=1.00, <-- should be -42, because id(x=1,y=1) = 1+8*1 = 9
	/* u_out_pgi[2]=2.00, 
	/* u_out_pgi[3]=3.00, 
	/* u_out_pgi[4]=4.00, 
	/* u_out_pgi[5]=5.00, 
	/* u_out_pgi[6]=6.00, <-- seems like id = x; (without the y*x_max addition)
	/* u_out_pgi[7]=-42.00,
	/* u_out_pgi[8]=-42.00,
	/* u_out_pgi[9]=-42.00,  
	
	   
	/* ideally we would like to compute a stencil type code like the 2D laplacian below */
	/* int idx_m1 = id-1; */
	/* int idx_p1 = id+1; */
	/* int idy_m1 = id-x_max; */
	/* int idy_p1 = id+x_max; */
	/* u_out[id] = -u[id] + 0.25*(u[idx_m1] + u[idx_p1] + u[idy_m1] + u[idy_p1]); */

      }
    
    }
  } 
}

void test_doubly_nested_loop_cpu(float *u_out, float *u, int x_max, int y_max) {

  int x,y;
  
  for (y=1; y<y_max-1; y++) {

    for(x=1;x<x_max-1; x++) {

      int id = x + x_max*y;
      int idx_m1 = id-1;
      int idx_p1 = id+1;
      int idy_m1 = id-x_max;
      int idy_p1 = id+x_max;	 

      /* test output corresponding to pgi-acc method above*/
      u_out[id] = id;

      /* a test laplacian stencil */
      /* u_out[id] = -u[id] + 0.25*(u[idx_m1] + u[idx_p1] + u[idy_m1] + u[idy_p1]); */

      
    }
    
  }
  
}

int main(int argc, char** argv) { 

  int x,y;
  
  int x_max = 8;
  int y_max = 8;

  float* u;
  float* u_out;
  float* u_cpu;
  float* u_out_cpu;
  
  /* allocate memory for pgi-acc function */
  u = (float*)malloc(x_max*(y_max+4)*sizeof(float));
  u_out = (float*)malloc(x_max*(y_max+4)*sizeof(float));

  /* allocate memory for simple cpu version */
  u_cpu = (float*)malloc(x_max*(y_max+4)*sizeof(float));
  u_out_cpu = (float*)malloc(x_max*(y_max+4)*sizeof(float));  


  /* initialize memory to some simple defaults */
  for(x=0;x<x_max*y_max;x++) {
    u[x] = 1.1;
    u_out[x] = -42;
    u_cpu[x] = 1.1;
    u_out_cpu[x] = -42;
  }

  printf("launching pgi kernel\n");  
  test_doubly_nested_loop_pgi_acc(u_out, u, x_max, y_max);
  
  printf("launching simple cpu kernel\n");
  test_doubly_nested_loop_cpu(u_out_cpu,u_cpu, x_max, y_max);
  
  printf("the pgi-gpu and cpu results (0-63)\n");
  for(x=0;x<64;x++) {
    printf("u_out_pgi[%d]=%2.2f, u_out_cpu[%d]=%2.2f\n",x,u_out[x],x,u_out_cpu[x]);
  }

  /* check for errors */
  /* double error=0; */
  /* int error_count=0; */
  /* for(x=0;x<x_max*y_max;x++) { */
  /*   error = fabs(u_out[x] - u_out_cpu[x]); */
  /*   if(error > 0.001) {       */
  /*     error_count++; */
  /*     printf("Error %d@u_out[%d]: |%f-%f|=%f\n",error_count, x, u_out[x],u_out_cpu[x],error); */
  /*     if(error_count > 30) { break;} */
  /*   } */
  /* } */
  
}

/*
$> pgaccelinfo

CUDA Driver Version:           3000

Device Number:                 0
Device Name:                   GeForce GTX 285
Device Revision Number:        1.3
Global Memory Size:            2147287040
Number of Multiprocessors:     30
Number of Cores:               240
Concurrent Copy and Execution: Yes
Total Constant Memory:         65536
Total Shared Memory per Block: 16384
Registers per Block:           16384
Warp Size:                     32
Maximum Threads per Block:     512
Maximum Block Dimensions:      512, 512, 64
Maximum Grid Dimensions:       65535 x 65535 x 1
Maximum Memory Pitch:          2147483647B
Texture Alignment:             256B
Clock Rate:                    1476 MHz
Initialization time:           8348 microseconds
Current free memory:           2092617728
Upload time (4MB):             1169 microseconds ( 992 ms pinned)
Download time:                 1878 microseconds ( 746 ms pinned)
Upload bandwidth:              3587 MB/sec (4228 MB/sec pinned)
Download bandwidth:            2233 MB/sec (5622 MB/sec pinned)
*/
