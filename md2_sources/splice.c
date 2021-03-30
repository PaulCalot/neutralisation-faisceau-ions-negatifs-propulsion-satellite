#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>

/*
 * splice.c -- splices an arbitrary number (MAX 25) of parallel data files 
 * to produce a BIG data file.  All files must have the same
 * number of rows.  First column is assumed to be x-data,  and
 * should be the same in each file.
 * All fields are assumed to be floating-point values.
 * 
 * (c) 1999 cameron abrams
 * 
 */
 
char ln[5000], hln[5000];
#define MAXFILES 25
#define MAXFIELDS 50
FILE * fp[MAXFILES];
char * fn[MAXFILES];

int ln2row (char * p, double row[])
{
    char substr[255];
    int i=0;
    if (!p) return 0;
    while (*p&&sscanf(p,"%s",substr))
    {
	p+=strlen(substr);
	while(isspace(*(++p)));
	if (substr[strlen(substr)-1]==',') substr[strlen(substr)-1]='\0';
	if (strcmp(substr, "+/-")) row[i++]=atof(substr);
    }
    return i;
}

#ifndef SKIP_HEADER
#define SKIP_HEADER 1
#endif
short skipheader=SKIP_HEADER;
int main (int argc, char * argv[])
{
    int i=0, j=0;
    double buf[MAXFIELDS];
    double rowout[MAXFILES][MAXFIELDS];
    int nr[MAXFILES];
    int nFiles=argc-1;
    char * reading,  * p;
    char x[25], xstr[25];
    int dum;
    short header, okheader;
    for (i=1;i<argc;i++) fn[i-1]=argv[i];
    fprintf(stderr, "# splice: %i files\n", nFiles);
    
    for (i=0;i<nFiles;i++) fp[i]=fopen(fn[i], "r");
    reading=fn[0];
    header=okheader=0;
    while (reading)
    {
	for (i=0;i<MAXFILES;i++) for (j=0;j<MAXFIELDS;j++) rowout[i][j]=0.0;
	for (i=0;i<MAXFILES;i++) nr[i]=0;
	for (i=0;i<nFiles;i++)
	{
	    reading=fgets(ln, 5000, fp[i]);
	    if (reading&&ln[0]!='#'&&ln[0]!='%') 
	    {
		if (!okheader) okheader=1;
/*		printf("m:%s", ln);fflush(stdout); */
                p=ln;
                sscanf(p,"%s",x);
                p+=strlen(x);
                while(isspace(*(++p)));
		nr[i]=ln2row(p, buf);
		for (j=0;j<nr[i];j++) rowout[i][j]=buf[j];
	    }
           if (i==0&&ln[0]=='#') {sprintf(hln,"%s",ln);}
 	}
	if (!header&&okheader&&!skipheader)
	{
	    header=1;
	    printf("#\t");
	    for (i=0;i<nFiles;i++)
	    {
		printf("+----\t");
		for (j=0;j<(nr[i]/2-2);j++) printf("-----\t");
		printf("file%i", i+1);
		dum=1;
		if (nr[i]%2) dum=2;
		for (j=0;j<(nr[i]/2-dum);j++) printf("\t-----");
		printf("\t----+\t");
	    }
	    printf("\n");
	    hln[strlen(hln)-1]='\0';
	    sscanf(hln, "%s", xstr);
	    p=&(hln[0])+strlen(xstr);
	    while(isspace(*(++p)));
	    printf("%s\t", xstr);
	    for (i=0;i<nFiles;i++)
	    {
		printf("%s\t", p);
	    }
	    printf("\n");
	}
	if (reading&&ln[0]!='#'&&ln[0]!='%')
	{
	    printf("%s", x);
	    for (i=0;i<nFiles;i++)
		for (j=0;j<nr[i];j++) printf("\t%.3lf", rowout[i][j]);
	    printf("\n");
	}
    }
    for (i=0;i<nFiles;i++) fclose(fp[i]);
    
}
