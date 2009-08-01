#include <stdio.h>
#include <stdlib.h>

int main()
{
	int n=0;
	float sum=0, val=0;
	float mean=0, diff=0, sd=0;
	
	while(scanf("%f ", &val)!=EOF) {
		n+=1;
		sum += val;
	}
	
	mean = sum/n;
	
	fseek(stdin, 0, SEEK_SET);
	
	while(scanf("%f ", &val)!=EOF) {
		diff += (mean-val)*(mean-val);
	}
	
	diff /= n;
	sd = sqrt(diff);
		
	printf("%f %f %f\n", sum, mean, sd);
}
