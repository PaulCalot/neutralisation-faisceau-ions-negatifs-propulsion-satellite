#include <stdio.h>
#include <math.h>
/*
 * intg.c -- uses the trapezoidal rule to integrate a function
 * supplied in the discrete form x_i, f(x_i).  Computes the maximum
 * value of f, the integral of f, and 'thickness' of f (integral/max).
 * 
 * two columns of input data can be read from stdin (no arguments)
 * or from an input file (single argument).  Output is to stdout.
 * MAXPTS can be redefined (below) if more than 2000 points are 
 * being considered.
 * 
 * (c) 1999 cam abrams
 */
typedef struct POINT
{
    double x, y;
} pt;

#define MAXPTS 2000
pt data[MAXPTS];
int nData;
char ln[255];
int main (int argc, char * argv[])
{
    FILE *fp=NULL;
    char *fn=NULL;
    int i;
    double s1, max=1.e-9;
    
    for (i=1;i<argc;i++)
    {
	if (argv[i][0]!='-') fn=argv[i];
    }

    if (!fn) fp=stdin;
    else fp=fopen(fn, "r");
    
    i=0;
    while (fgets(ln, 255, fp))
    {
	if (ln[0]!='#')
	{
	    sscanf(ln, "%lf %lf", &(data[i].x), &(data[i].y));
	    max=(data[i].y>max?data[i].y:max);
	    i++;
	}
    }
    nData=i;
    
    s1=0.0;
    for (i=0;i<nData-1;i++) s1+=(data[i+1].x-data[i].x)*(data[i].y+data[i+1].y);
    s1*=0.5;
    
    printf("# nData = %i; I = %.5lf; I/max = %.5lf\n", nData, s1, s1/max);
    
    if (fp!=stdin) fclose(fp);

}
