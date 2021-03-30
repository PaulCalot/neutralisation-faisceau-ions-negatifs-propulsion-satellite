/* 
 * md series 2 
 * 
 * (c) 1999 cameron abrams
 * 
 * 
 *
 * abc.c: operates on output of rtd script to compute absolute indices
 * of all atoms from incident ions,  ultimately computing residence time
 * distributions for these atoms.
 *
 */

#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <math.h>

#define MXAPERION 10
#define MXY 100

typedef struct DATA
{
   int n, nci, ci[MXAPERION], nfi, fi[MXAPERION], ns, s[MXY];
} data_type;

void dt_out (FILE * fp, data_type * dt)
{
    int i=0;
    if (!fp || !dt) return;
    fprintf(fp, "n%i nci%i ", dt->n, dt->nci);
    for (i=0;i<dt->nci;i++) fprintf(fp, "c%i ", dt->ci[i]);
    fprintf(fp, "nfi%i ", dt->nfi);
    for (i=0;i<dt->nfi;i++) fprintf(fp, "f%i ", dt->fi[i]);
    fprintf(fp, "ns%i ", dt->ns);
    for (i=0;i<dt->ns;i++) fprintf(fp, "s%i ", dt->s[i]);
    fprintf(fp, "\n");
}

#define MAXBINS 2000
#define MAXLINES 4000
data_type data[MAXLINES];
int crtd[MAXBINS], frtd[MAXBINS];
char ln[1000];
#define NEXT_FIELD {j=strlen(substr);while(*p&&j--)p++;while(*p&&isspace(*(++p)));}
int ln2dt (char ln[], data_type * dt)
{
    char * p = ln;
    char substr[255];
    int i=0, j=0;
    if (!p) return 0;
    
    /* First field is the impact number */
    sscanf(p,"%s",substr);NEXT_FIELD;
    dt->n=atoi(substr);
    /* Second field is the number of carbons in ion */
    sscanf(p,"%s",substr);NEXT_FIELD;
    dt->nci=atoi(substr);
    for (i=0;i<dt->nci;i++)
    {
	sscanf(p,"%s",substr);NEXT_FIELD;
	dt->ci[i]=atoi(substr);
    }
    /* Next field is number of f's in ion */
    sscanf(p,"%s",substr);NEXT_FIELD;
    dt->nfi=atoi(substr);
    for (i=0;i<dt->nfi;i++)
    {
	sscanf(p,"%s",substr);NEXT_FIELD;
	dt->fi[i]=atoi(substr);
    }
    /* Next field is number of sputtered atoms */
    if (*p) 
    {
	sscanf(p,"%s",substr);NEXT_FIELD;
	dt->ns=atoi(substr);
	for (i=0;i<dt->ns;i++)
	{
	    sscanf(p,"%s",substr);NEXT_FIELD;
	    dt->s[i]=atoi(substr);
	}
    }
}

int main (int argc, char * argv[])
{
    int nData=0,i=0, j=0, k=0, l=0, dum1, dum2, id, cif,fif, ctot, ftot, 
        tid, ni, nk;
    data_type * dti,  * dtk;
    int cg=1;
    double csum, ccum, fsum, fcum;
    char * p;
    short dbg=0;
    for (i=0;i<MAXBINS;i++) crtd[i]=frtd[i]=0;
    
    for (i=1;i<argc;i++)
    {
	if (!strcmp(argv[i], "-cg")) cg=atoi(argv[++i]);
        else if (!strcmp(argv[i], "-D")) dbg=1;
    }
    
    while (fgets(ln,1000,stdin))
    {
	if (ln[0]!='#')
	{
	    ln2dt(ln, &(data[nData++]));
	    if (dbg) dt_out(stdout, &(data[nData-1]));
	}
    }    
    ni=0;
    
    for (i=0;i<nData;i++)
    {
/*        printf("#pre %i %i ", data[i].n, data[i].ci); 
	for (j=0;j<data[i].nc;j++)
	{
	    printf("%i ", data[i].cs[j]);
	}
	printf("\n"); */
	ni+=data[i].ns;
    }
    printf("# read data from %i lines\n", nData);
    printf("# total number of sputtered nonspecific atoms is %i\n", ni);
    
    
    for (i=0;i<nData;i++)
    {
	dti=&(data[i]);
	ni=dti->n;
        /* trace each ion carbon until a match is found in the
	 * sputtered ids of this trajectory or a later one */
	for (j=0;j<dti->nci;j++)
	{
	    id=dti->ci[j];
	    tid=id;
	    cif=0;
	    for (k=i;k<nData&&!cif;k++)
	    {
		dtk=&(data[k]);
		nk=dtk->n;
		for (l=0;l<dtk->ns&&dtk->s[l]!=id;l++) 
		    if (dtk->s[l]<id) tid--;
		if (dtk->s[l]==id) cif=1;
		else id=tid;
	    }
	    if (cif) 
	    {
		crtd[nk-ni+1]++;
	    }
	    else 
	    {
		crtd[0]++;
	    }
	}
       /* trace each ion fluorine until a match is found in the
	 * sputtered ids of this trajectory or a later one */
	for (j=0;j<dti->nfi;j++)
	{
	    id=dti->fi[j];
	    tid=id;
	    fif=0;
	    for (k=i;k<nData&&!fif;k++)
	    {
		dtk=&(data[k]);
		nk=dtk->n;
		for (l=0;l<dtk->ns&&dtk->s[l]!=id;l++) 
		    if (dtk->s[l]<id) tid--;
		if (dtk->s[l]==id) fif=1;
		else id=tid;
	    }
	    if (fif) 
	    {
		frtd[nk-ni+1]++;
	    }
	    else 
	    {
		frtd[0]++;
	    }
	}
    }
    /* output the rtd, course grained according to value of cg */
    ctot=0;
    for (i=1;i<MAXBINS;i++) ctot+=crtd[i];
    ftot=0;
    for (i=1;i<MAXBINS;i++) ftot+=frtd[i];
    printf("# (c,f) total **sputtered** %i %i, adsorbed %i %i, total %i %i\n", 
	ctot, ftot, crtd[0], frtd[0], ctot+crtd[0], ftot+frtd[0]);
    ccum=fcum=0.0;
    for (i=1;i<MAXBINS;i+=cg)
    {
	csum=0.0;
	for (j=0;j<cg;j++) csum+=((i+j)<MAXBINS)?(double)crtd[i+j]:0;
        ccum+=csum/cg/ctot; 
	fsum=0.0;
	for (j=0;j<cg;j++) fsum+=((i+j)<MAXBINS)?(double)frtd[i+j]:0;
        fcum+=fsum/cg/ftot; 
	if (csum||fsum) printf("%i %.4lf %.4lf %.6lf %.6lf\n", 
	    i, csum/cg, fsum/cg, ccum, fcum);
    }
}

