#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void laplacian_pgi(float *u_out, float *u, int x_max, int y_max) {

  int x,y;
#pragma acc region copyin(u[0:x_max*y_max-1]) copyout(u_out[0:x_max*y_max-1])
  {
#pragma acc for independent
    for (y=1; y<y_max-1; y++) {
#pragma acc for independent
      for(x=1;x<x_max-1; x++) {

	int id = x + x_max*y;
	int idx_m1 = id-1;
	int idx_p1 = id+1;
	int idy_m1 = id-x_max;
	int idy_p1 = id+x_max;

	u_out[id] = -u[id] + 0.25*(u[idx_m1] + u[idx_p1] + u[idy_m1] + u[idy_p1]);
      
      }
    
    }
  } 
}

void laplacian_cpu(float *u_out, float *u, int x_max, int y_max) {

  int x,y;
  
  for (y=1; y<y_max-1; y++) {

    for(x=1;x<x_max-1; x++) {

      int id = x + x_max*y;
      int idx_m1 = id-1;
      int idx_p1 = id+1;
      int idy_m1 = id-x_max;
      int idy_p1 = id+x_max;

      u_out[id] = -u[id] + 0.25*(u[idx_m1] + u[idx_p1] + u[idy_m1] + u[idy_p1]);
      
    }
    
  }
  
}

int main() {

  int x,y;
  
  int x_max = 8;
  int y_max = 8;

  float* u;
  float* u_out;
  float* u_cpu;
  float* u_out_cpu;

  u = (float*)malloc(x_max*y_max*sizeof(float));
  u_out = (float*)malloc(x_max*y_max*sizeof(float));
  u_cpu = (float*)malloc(x_max*y_max*sizeof(float));
  u_out_cpu = (float*)malloc(x_max*y_max*sizeof(float));

  for(x=0;x<x_max*y_max;x++) {
    u[x] = 1.1;
    u_out[x] = -42;
    u_cpu[x] = 1.1;
    u_out_cpu[x] = -42;
  }

  printf("launching pgi kernel\n");
  laplacian_pgi(u_out, u, x_max, y_max);
  printf("launching simple cpu kernel\n");
  laplacian_cpu(u_out_cpu,u_cpu, x_max, y_max);
  
  double error=0;
  int error_count=0;
  for(x=0;x<x_max*y_max;x++) {
    error = fabs(u_out[x] - u_out_cpu[x]);
    if(error > 0.001) {      
      error_count++;
      printf("Error %d@u_out[%d]: |%f-%f|=%f\n",error_count, x, u_out[x],u_out_cpu[x],error);
      if(error_count > 30) { exit(1); }
    }
  }
  
}
