#include <stdio.h>
#include <stdlib.h>
#include <time.h>

double gaussRand (void)
{
   int i;
   double sum;

   sum = 0;
   for (i=1;i<=12;i++)
       sum+=((double)rand())/((double)RAND_MAX);        

   return (sum-6.0);
}

void main (int argc, char * argv[])
{
    int n=atoi(argv[1]);
    int i=0;
    
    srand((unsigned)clock()+(unsigned)time(NULL));
    
    for (i=0;i<n;i++)
	printf("%.5lf\n", gaussRand());
    
}
