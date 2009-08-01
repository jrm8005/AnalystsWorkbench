#include <stdio.h>
#include <stdlib.h>

#define chatty 0

int main(int argc, char** argv)
{
	float expect, val;
	if (argc < 2) {
		perror("Usage: mse expected-val < values");
		exit(1);
	}
	
	expect = atoff(argv[1]);
	if (chatty)
		printf("%f\n", expect);
	
	while (scanf("%f ", &val) != EOF) {
		printf("%f ", (expect-val)*(expect-val));
	}
	printf("\n");
}
