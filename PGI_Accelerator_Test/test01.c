#include <stdlib.h>

int main (int argc, char** argv)
{
	int i;
	float *a, *b;
	float c[1000];

	a = (float*) malloc (1000 * sizeof (float));
	b = (float*) malloc (1000 * sizeof (float));

	for (i = 0; i < 1000; i++)
		b[i] = 123;

	#pragma acc region copyin(b[0:999]) copyout(a[0:999])
	{
		for (i = 0; i < 1000; i++)
			a[i] = b[i] + 1;
	}

	#pragma acc region
	{
		for (i = 0; i < 1000; i++)
			c[i] = a[i] + 1;
	}

	return EXIT_SUCCESS;
}
