/* nbrof:  given an input cfg file and a number corresponding to 
 * the "atom of interest" in that cfg, outputs that atom and it's
 * neighbors in config format to a specified output file.
 *
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
    printf("%%nbrof (input cfg) (#of atom of interest) (output cfg)\n");
    exit(0);
}

void main (int argc, char * argv[])
{
    atomPtr a = NULL;
    nNodePtr np = NULL;
    int id;
    FILE * fp = NULL;
    char * inf=0L, * of=0L;
    
    if (argc < 2) usage();
    inf=argv[1];
    id=atoi(argv[2]);
    if (argc>3) of=argv[3];

    Chem_InitializePeriodicTable();
    pef_initialize();
    cfg_initialize();
    strcpy(cfg_infile_, inf);
    cfg_establish();
    
    printf("#%s: Total %i T=%.2lf\n", inf, nAtom_, T_);
    cfg_setforce();

    ptPtr_scalmult(Tau_, Tau_, 26.70167/cfg_v()); /* GPa */
    fflush(stdout);

    for (a=L_;a&&a->id!=id;a=a->next);
    if (a)
    {
	if (of) fp=fopen(of, "w");
	else fp=stdout;
	fprintf(fp,"#BoxSize.xyz\t%.14lf %.14lf %.14lf\n", 
	    Lr_.x, Lr_.y, Lr_.z);   
	cfg_atomout(a, fp);
	for (np=a->nList;np;np=np->next) 
	    if (atom_getfij(a, np->addr)==1.0&&
		atom_peij(a, np->addr) < 0.0)
		    cfg_atomout(np->addr, fp);
	if (fp!=stdout) fclose(fp);
    }
    else
    {
	printf("Error: could not find atom %i in %s.\n", id, inf);
	exit(0);
    }

}


