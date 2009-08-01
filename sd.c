#include <stdio.h>
#include <stdlib.h>
#include <float.h>

int main(int argc, char** argv)
{
	int n=0;
	float sum=0, val=0, max=FLT_MIN, min=FLT_MAX;
	float mean=0, diff=0, sd=0;
	
	while(scanf("%f ", &val)!=EOF) {
		n+=1;
		sum += val;
		if(val > max)
			max = val;
		if(val < min)
			min = val;
	}
	
	mean = sum/n;
	
	fseek(stdin, 0, SEEK_SET);
	
	while(scanf("%f ", &val)!=EOF) {
		diff += (mean-val)*(mean-val);
	}
	
	diff /= n;
	sd = sqrt(diff);
		
	printf("%f %f %f %f %f\n", sum, max, min, mean, sd);
}