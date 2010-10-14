#include <stdlib.h>

#define DIM 4

int main(int argc, char* argv[])
{
   float* samples;
   int i, threads;

   threads = 1024;  /* number of threads */

   samples = (float*) malloc(sizeof(float) * DIM * threads);

   #pragma acc region copyin(samples[0:(DIM * threads)-1])
   {
      #pragma acc for
      for (i=0; i<threads; i++) {
         int j;
         
         for (j=0; j<3; j++) {
            float sum;
            int j;

            sum = samples[i * DIM + 0] + samples[i * DIM + 1] + samples[i * DIM + 2];
         }
      }
   }
}
