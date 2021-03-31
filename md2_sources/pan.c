#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
typedef struct VECTOR3
{
    double x;
    double y;
    double z;
} vector_type;
double f (double x)
{
    /* map function */
    return 0.5-0.5*cos(M_PI*x);
}
void main (int argc, char * argv[])
{
    int i = 0, j=0;
    int n, N;
    double F, f0=0, f1=1, x0=0, x1=1, x;
    vector_type c1={0, 0, 0}, c2={0, 0, 0}, d={0, 0, 0};
    short oneD = 0;
    
    for (i=1;i<argc;i++)
    {
	if (!strcmp(argv[i], "-n")) n = atoi(argv[++i]);
	else if (!strcmp(argv[i], "-N")) N = atoi(argv[++i]);
	else if (!strcmp(argv[i], "-c1"))
	{
	    c1.x = atof(argv[++i]);
	    c1.y = atof(argv[++i]);
	    c1.z = atof(argv[++i]);
	}
	else if (!strcmp(argv[i], "-c2"))
	{
	    c2.x = atof(argv[++i]);
	    c2.y = atof(argv[++i]);
	    c2.z = atof(argv[++i]);
	}
	else if (!strcmp(argv[i], "-1d"))
	    oneD = 1;
	else if (!strcmp(argv[i], "-x0"))
	    x0 = atof(argv[++i]);
	else if (!strcmp(argv[i], "-x1"))
	    x1 = atof(argv[++i]);
	else {printf("%s?\n", argv[i]);exit(0);}
    }
    if (x0 < 0.0 || x1 > 1.0) {printf("bad bounds\n");exit(0);}
    f0=f(x0);
    f1=f(x1);
/*    printf("n/N=%.5lf\n", (double)n/(double)N);*/
    x=((double)n/N*(x1-x0)+x0);
    F=f(x);
    F=(F-f0)/(f1-f0);
    d.x=c1.x+(F)*(c2.x-c1.x);
    d.y=c1.y+(F)*(c2.y-c1.y);
    d.z=c1.z+(F)*(c2.z-c1.z);
    
    printf("%.10lf", d.x);
    if (!oneD) printf(" %.10lf %.10lf", d.y, d.z);
    printf("\n");
    
}
