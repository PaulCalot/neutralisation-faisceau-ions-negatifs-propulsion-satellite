#include <stdio.h>
#include <math.h>

int main (int argc, char * argv[])
{
    double r=0.0;
    double dr=2*M_PI/100;
    int i;
    
    for (i=0;i<100;i++)
    {
	r=dr*i;
	printf("%.5lf %.5lf\n", cos(r), sin(r));
    }
}
