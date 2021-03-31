/* test of central limit theorem */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>


int main (int argc, char * argv[])
{
    int n=0, N=0;
    double r=0.0;
    double hist[100];
    double b=0.01;
    int i, j;
    
    for (i=0;i<100;i++) hist[i]=0.0;
    
    n=atoi(argv[1]);
    N=atoi(argv[2]);
    
    printf("# test of central limit theorem; n=%i, N=%i\n", n, N);
    
    srand((unsigned)clock()+(unsigned)time(NULL));

    for (i=0;i<N;i++) /* number of sums */
    {
	r=0.0;
	for (j=0;j<n;j++) /* number of independent random numbers per sum */
	    r+=((double)rand())/((double)RAND_MAX);
	r/=n/2;  /* makes sure the final sum value is [0,2] */
	r-=1.0; /* makes sure the final sum value is [-1,1] */
/*	hist[(int)(r/b)]++; */
        printf("%.5lf\n", r);
    }
    
/*     r=0.0;
    for (i=0;i<100;i++)
	r+=hist[i];
    for (i=0;i<100;i++)
	printf("%.5lf %.5lf\n", i*b-0.5, hist[i]/r);
*/
}
