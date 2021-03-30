/*
 * ydist.c    by cam abrams
 * (c) 1999 cfa
 * 
 * ydist.c is a simple code to compute a histogram distribution from a column 
 * of numbers.  It can optionally display (1) the histogram data and (2) the
 * stddev in that precedence. 
 * 
 * usage:
 * 
 * ydist [-bs #] [-ymin #] [-ymax #] [-f #] [-n #] [-| |file]
 * 
 * -ymin    allows user to specify minimum bound on histogram;
 *	    if not supplied by user, it is set to the minimum value of
 *	    of the data read in.
 * -ymin    allows user to specify maximum bound on histogram;
 *	    if not supplied by user, it is set to the maximum value of
 *	    of the data read in.
 * -bs	    allows user to specify bin size of histogram; default is 1.0.
 * -n	    allows user to specify the normalization constant, if it is
 *	    different from the count of data points.
 * -f	    bin growth factor, for binning logarithmic data; default is 1.0
 *	    (no growth)
 * - (or nothing)
            tells ydist.c that the input is a column of numbers from stdin, 
 *          terminated by ^D or EOF.
 * (filename) 
 *          tells ydist.c to open and read a single column from the file
 *          named 'filename'.
 * 
 * compile this as 'cc -o ydist ydist.c -lm'
 * 
 * 15Apr1999
 *
 */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#define MAXLINE 255
char scr_line_[MAXLINE];
#define MAXBINS 5000
double hist[5000], secav[5000];
#define MAXDATA 1000000
double data[MAXDATA], secdata[MAXDATA];

void usage (void)
{
    printf("usage:\n");
    printf("ydist {- | <fileName>} (options) \n");
    printf("Options:\n");
    printf("\t-ymin # -ymax #\n");
    printf("\t-n <norm const(calculated)>\n");
    printf("\t-bs <bin size (1.0)>\n");
    printf("\t-f <bin size growth factor (1.0=const. bin size)>\n");
    exit(-1);
}

int bin (double x, double bs, double Ymin, double f)
{
    int i=0;
    double sum=Ymin;
    
    while (i<MAXBINS&&x>sum) sum+=bs*pow(f,i++);
    if (i==MAXBINS) {printf("error bigbin\n");exit(0);}
    
    return i;
}

double binv (int b, double bs, double Ymin, double f)
{
    int i=0;
    double sum=0.0;
    
    for (i=0;i<b;i++) sum+=bs*pow(f,i);
    sum-=0.5*bs*pow(f,--i);
    return sum+Ymin;
}

char y1str[50];
void main (int argc, char * argv[])
{
    int i=0, b=0;
    char * fn=NULL;
    char * p=NULL;
    double sum_a=0.0, sumxx_a=0.0;
    double y1=0.0, y2=0.0;
    double totaly=0.0;
    double ymin=0.0, ymax=0.0, bval=0.0;
    double Ymin=-1.e9, Ymax=-1.e9;
    int nData=0;
    double binsize=1.0;
    double f=1.0;
    FILE * fp=stdin;    
    
    for (i=0;i<MAXBINS;i++) hist[i]=secav[i]=0.0;
    for (i=0;i<MAXDATA;i++) data[i]=secdata[i]=0.0;

    for (i=1;i<argc;i++)
    {
	if (argv[i][0] != '-') fn=argv[i];
	else if (!strcmp(argv[i], "-")) fp=stdin;
	else if (!strcmp(argv[i], "-bs")) binsize=atof(argv[++i]);
	else if (!strcmp(argv[i], "-f")) f=atof(argv[++i]);
	else if (!strcmp(argv[i], "-ymin")) Ymin=atof(argv[++i]);
	else if (!strcmp(argv[i], "-ymax")) Ymax=atof(argv[++i]);
	else if (!strcmp(argv[i], "-n")) totaly=atof(argv[++i]);
        else usage();
    }
    if (fn)
    {
        fp=fopen(fn, "r");
        if (!fp) {printf("Error. Could not open %s.\n", fn);exit(-1);}
    }
    
    i=0;
    while (fgets(scr_line_, MAXLINE, fp))
    {
	p=scr_line_;
	/* read the first number on the line */
	sscanf(p, "%s", &y1str);
	y1=atof(y1str);
	
	data[i]=y1;
	sum_a+=y1;
	sumxx_a+=y1*y1;
	ymax=(y1>ymax?y1:ymax);
	ymin=(y1<ymin?y1:ymin);
	/* advance the pointer until we hit a second number or EOR */
	p+=strlen(y1str);
	while(isspace(*(++p)));
	if (p) 
	{
	    sscanf(p, "%lf", &y2);
	    secdata[i]=y2;
	}
	i++;
    }
    if (fp!=stdin) fclose(fp);
    nData=i;
    if (!totaly) totaly=nData;
    if (Ymax==-1.e9) Ymax=ymax;
    if (Ymin==-1.e9) Ymin=ymin;
    printf("# mean: %.5lf",sum_a/i);
    printf(" sd: %.5lf", sqrt((sumxx_a-sum_a*sum_a/i)/i));
    printf(" min: %.5lf max: %.5lf", Ymin, Ymax);
    printf("\n");

    if (!binsize) {printf("error zerobinsize\n");exit(0);}
    
    printf("# histogram: nData=%i binsize(0)=%.4lf growth=%.3lf\n", 
	nData, binsize, f);
    if (f!=1.0) printf("# Warning: distribution function not normalized.\n");
    for (i=0;i<nData;i++) 
    {
	b=bin(data[i],binsize,Ymin,f);
	hist[b]++;
	secav[b]+=secdata[i];
    }
    for (i=0;i<MAXBINS;i++) 
    {
	secav[i]=hist[i]?secav[i]/hist[i]:0.0;
	hist[i]/=totaly;
    }
    for (i=0;i<MAXBINS&&(bval=binv(i, binsize, Ymin, f))<=Ymax;i++) 
    {
	if (bval>=Ymin&&bval<=Ymax)
	    fprintf(stdout, "%i %.5le %.5le %.5le %.5le\n", 
		i, bval, hist[i], (f==1.0?hist[i]/binsize:0.0), secav[i]);
    }
    
}
