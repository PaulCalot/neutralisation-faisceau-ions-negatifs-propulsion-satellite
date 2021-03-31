#include <stdio.h>
#include <math.h>
/*
 * intg.c -- uses 2nd order finite difference to differentiate a function
 * supplied in the discrete form x_i, f(x_i).
 * 
 * Two columns of input data can be read from stdin (no arguments)
 * or from an input file (single argument).  Output is 2 columns, 
 * [x_i, f'(x_i)], to stdout.
 * MAXPTS can be redefined (below) if more than 5000 points are 
 * being considered.
 * 
 * (c) 1999 cam abrams
 */
typedef struct POINT
{
    double x, y;
} pt;

#define MAXPTS 5000
pt data[MAXPTS];
int nData;
char ln[255];
int main (int argc, char * argv[])
{
    FILE *fp=NULL;
    char *fn=NULL;
    int i;
    double diff, max=1.e-9;
    
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
    
    diff=0.0;
    for (i=0;i<nData;i++)
    {
	if (i==0) diff=(data[i].y-data[i+1].y)/(data[i].x-data[i+1].x);
	else if (i==(nData-1))
		  diff=(data[i-1].y-data[i].y)/(data[i-1].x-data[i].x);
	else diff=(data[i-1].y-data[i+1].y)/(data[i-1].x-data[i+1].x);
	printf("%.5lf %.5lf\n", data[i].x, diff);
    }
    
    printf("# nData = %i\n", nData);
    
    if (fp!=stdin) fclose(fp);

}
