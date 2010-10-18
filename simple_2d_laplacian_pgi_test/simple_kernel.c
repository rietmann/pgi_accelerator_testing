#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void laplacian_cpu(float *u_out, float *u, int x_max, int y_max);
void laplacian_pgi(float *u_out, float *u, int x_max, int y_max);

void laplacian_pgi(float *restrict u_out, float *restrict u, int x_max, int y_max) {

  int x,y;
#pragma acc region copyin(u[0:x_max*y_max]) copy(u_out[0:x_max*y_max])
  {
#pragma acc for independent
    for (y=1; y<y_max-1; y++) {
#pragma acc for independent
      for(x=1; x<x_max-1; x++) {
	
	/* int id = x + x_max*y; */
	int id = x_max*y;
	id = x + id;
	
	int idx_m1 = id-1;
	int idx_p1 = id+1;
	int idy_m1 = id-x_max;
	int idy_p1 = id+x_max;
	/* u_out[id] = -u[id] + 0.25*(u[idx_m1] + u[idx_p1] + u[idy_m1] + u[idy_p1]); */
	/* u_out[x + x_max*y] =  x + x_max*y; */
	u_out[id] = id;
	
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

      /* u_out[id] = -u[id] + 0.25*(u[idx_m1] + u[idx_p1] + u[idy_m1] + u[idy_p1]); */
      u_out[id] = id;
      
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
  float** u2;
  float** u_out2;
  
  
  u = (float*)malloc(x_max*(y_max+4)*sizeof(float));
  u_out = (float*)malloc(x_max*(y_max+4)*sizeof(float));
  u_cpu = (float*)malloc(x_max*(y_max+4)*sizeof(float));
  u_out_cpu = (float*)malloc(x_max*(y_max+4)*sizeof(float));  
  
  for(x=0;x<x_max*y_max;x++) {
    u[x] = 1.1;
    u_out[x] = -42;
    u_cpu[x] = 1.1;
    u_out_cpu[x] = -42;
  }

  /* printf("trying pgi_v2 w/ multi-dimensional array\n"); */
  /* laplacian_pgi_v2(u2,u_out2,x_max,y_max); */
  printf("launching pgi kernel\n");  
  laplacian_pgi(u_out, u, x_max, y_max);
  printf("launching simple cpu kernel\n");
  laplacian_cpu(u_out_cpu,u_cpu, x_max, y_max);


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

  printf("the first few cpu & gpu results\n");
  for(x=0;x<64;x++) {
    int x1 = floor(x/x_max);
    int x2 = x%x_max;
    printf("u_out_pgi[%d]=%2.2f, u_out_cpu[%d]=%2.2f\n",x,u_out[x],x,u_out_cpu[x]);
  }
  
}
