#include <stdio.h>

/*
 * This program simply reports the PE of the cluster
 * Si---F- - -F
 * (1) (2)   (3)
 * as a function of r23, where r12 is fixed at 1.6 A.
 * 
 */
 
#include "tbt_sicf.h"
#include "atom.h"
#include "cfg.h"
#include "chem.h"

double version_;
int build_;

extern char cfg_infile_[];
extern short pe_set_;
extern double PE_;

int main (int argc,  char * argv[])
{
    double r2x, r12;
    int i=0;
    atomPtr si=NULL, f1=NULL, f2=NULL;
    double del=0.05, r12i=4.0, r12f=1.0;
    Chem_InitializePeriodicTable();
    pef_initialize();
    cfg_initialize();

    strcpy(cfg_infile_, argv[1]);
    
    for (i=2;i<argc;i++)
    {
	if (!strcmp(argv[i], "-del")) del=atof(argv[++i]);
	else if (!strcmp(argv[i], "-r12i")) r12i=atof(argv[++i]);
	else if (!strcmp(argv[i], "-r12f")) r12f=atof(argv[++i]);
    }
    
    printf("%s\n", cfg_infile_);
    
    cfg_establish();
    
    si=cfg_addr(0);
    f1=cfg_addr(1);
    f2=cfg_addr(2);
    f2->pos->x=f1->pos->x+r12i;
    
    for (r12=r12i;r12>r12f;f2->pos->x-=del,r12=f2->pos->x-f1->pos->x)
    {
	pe_set_=0;
	cfg_setforce();
	printf("%.5lf %.5lf\n", r12, PE_);
	fflush(stdout);
    }
    
    cfg_clear();
    
}
