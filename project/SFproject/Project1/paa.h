#include <stdlib.h>
#include "config.h"

double * paa(double * original, int length)
{
	double * now = (double *)malloc(sizeof(double) * DIM);  //开辟数组空间，长度为DIM
	
	int i;
	for(i = 1;i <= DIM;i++)   
	{
		double num = 0;
		int j = length * (i - 1) / DIM + 1;
		int max = length * i / DIM;
		for(;j <= max;j++)
		{
			num += original[j-1];
		}
		now[i-1] = DIM * num / length;
		//cout << now[i-1] << endl;
	}
	return now;
}
