/* cfgrdf:  a post-processing program that computes rdfs
 * based on a series of static config files.  This is the series
 * 2 implementation.
 * 
 * (c) 1999 Cameron Abrams
 * University of California, Berkeley
 * Department of Chemical Engineering
 * 
 */

#include "atom.h"
#include "chem.h"
#include "cfg.h"
#include "genforce.h"
#include "tbt_sicf.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

double version_=2.00;
int build_=1;

extern atomPtr L_;
extern int nAtom_;
extern element per_table[];
extern int nElem_[];
extern pt Lr_, half_Lr_;
extern unit_type U_;
extern double T_, KE_;
extern ptPtr Tau_;
extern char cfg_infile_[];
void usage(void)
{
    printf("Usage:\n");
    printf("%%cfgrdf [list of cfg files] [keywords/values]\n", VERSION_NUMBER);
    printf("\tKeywords:\n");
    printf("\t-zlo\tspatial domain lower boundary (0.0)\n");
    printf("\t-zhi\tspatial domain upper boundary (100.0)\n");
    printf("\t-bs\thistogram binsize (0.1)\n");
    printf("\t-hhi\thistogram domain maximum (10.0)\n");
    exit(0);
}
#define MBINS 1000
#define MAXCFGS 500
void main (int argc, char * argv[])
{
    atomPtr a = NULL, b = NULL;
    pt TPT={0, 0, 0}, *tpt=&TPT;
    int i = 0, j = 0, k=0, n=0;
    FILE * fp = NULL;
    double hist[MBINS];
    char * cfgs[MAXCFGS];
    double zlo=0.0, zhi=50.0, bs=0.1, hhi=10.0, zmax, zmin, r, ru;
    int hn = 0, bin=0;
    double dens = 0.0, cnst=0.0, nideal=0.0, rijmax=0, vf=0.0;
    int nCfgs = 0;
    
    for(j=0;j<MBINS;j++) hist[j]=0.0;
    if (argc < 2) usage();
    for (i = 1; i < argc; i++)
    {
	if (!strcmp(argv[i], "-zlo")) zlo=atof(argv[++i]);
	else if (!strcmp(argv[i], "-zhi")) zhi=atof(argv[++i]);
	else if (!strcmp(argv[i], "-bs")) bs=atof(argv[++i]);
	else if (!strcmp(argv[i], "-vf")) vf=atof(argv[++i]);
	else if (!strcmp(argv[i], "-hhi")) hhi=atof(argv[++i]);
	else if (isalnum(argv[i][0])) cfgs[nCfgs++]=argv[i];
	else usage();
    }
    hn=(int)(hhi/bs);
    if (hn>=MBINS) {printf("hn too big\n");exit(0);}
    
    Chem_InitializePeriodicTable();
    pef_initialize();
    cfg_initialize();
    
    for (i=0;i<nCfgs;i++)
    {
	for (j=0;j<MBINS;j++) hist[j]=0;
	fp = fopen(cfgs[i], "r");
	if (!fp) {sleep(2);fopen(cfgs[i], "r");}
	if (!fp) {printf("error -- cannot find %s\n", cfgs[i]);exit(0);}
	strcpy(cfg_infile_, cfgs[i]);
	cfg_establish();
    
	/* shift so that lowest atoms are at z=0 */
	for (a=L_;a;a=a->next) a->pos->z+=half_Lr_.z;
	printf("#%s: Total %i T=%.2lf", cfgs[i], nAtom_, T_);
	fflush(stdout);
	for (k = 0; k < MAXNUMELEMENTS; k++)
	{
	    if (nElem_[k])
	    {
		printf(" %s %i", per_table[k].sym, nElem_[k]);
	    }
	}
	printf("\n");
	n=0;
	zmax=-1.e9;
	zmin=1.e9;
	rijmax=-1.e9;
	for (a=L_;a;a=a->next)
	{
	    zmax=(a->pos->z>zmax?a->pos->z:zmax);
	    zmin=(a->pos->z<zmin?a->pos->z:zmin);
	    if (a->pos->z>=zlo&&a->pos->z<=zhi)
	    {
		n++;
		for (b=a->next;b;b=b->next)
		{
		    if (b->pos->z>=zlo&&b->pos->z<=zhi)
		    {
			ptPtr_minimg(tpt, ptPtr_subtract(tpt, a->pos, b->pos));
			r=sqrt(ptPtr_sqabs(tpt));
			rijmax=(r>rijmax?r:rijmax);
			bin=(int)(r/bs);
			hist[bin]+=2;
		    }
		}
	    }
	}
	zmin=(zlo>zmin?zlo:zmin);
	dens=n/(Lr_.x*Lr_.y*(zmax-zmin+vf));
	cnst=4.0*M_PI*dens/3.0;
	printf("# n=%i, V=%.5lfx%.5lfx%.5lf=%.5lf, dens=%.5lf\n", 
	    n, Lr_.x, Lr_.y, (zmax-zmin), Lr_.x*Lr_.y*(zmax-zmin), dens);
	printf("# rijmax=%.5lf\n", rijmax);
	for (k=0;k<hn;k++) 
	{
	    r=k*bs;
	    ru=r+bs;
	    nideal=cnst*(pow(ru,3)-pow(r,3));
	    printf("%.5lf %.5lf %.5lf %.5lf\n", 
		(r+ru)/2.0, nideal, hist[k], hist[k]/(nideal*n));
	}
    }

}


