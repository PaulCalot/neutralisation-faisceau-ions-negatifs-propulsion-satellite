/*
 * md series 2 -- Molecular Dynamics of Plasma-Surface Chemistry
 * (c) 1999 Cam Abrams, Dept. of Chemical Engineering
 * Regents of the University of California, Berkeley
 * 
 * Module Name:  cfgfix
 * 
 * Module Description:
 * 	    usage:  cfgfix {cfg-file-name} {factor} {actioncode}
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "chem.h"
#include "cfg.h"
#include "genforce.h"
#include "tbt_sicf.h"
#include "point.h"

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
static FILE * fp;
int main (int argc, char * argv[])
{
    atomPtr a = NULL;
    double fracToFix = 0.0, factor = 0.0, topz=0.0, botz=0.0;
    int nAtom = 0, nAlreadyFixed = 0, nToFix = 0, i = 0;
    int actionCode = 0;
    
    if (argc < 2)
    {
	printf("Usage: ");
	printf("cfgFix <cfgFileName> <factor> <actionCode>\n");
	printf("\tcfgFileName\t\tname of file containing cfg data\n");
	printf("\tactionCode\t\t1 == Fix N=factor*Natoms atoms at bottom (default)\n");
	printf("\t          \t\t2 == Remove factor highest atoms\n");
	printf("\t          \t\t3 == Remove factor lowest atoms\n");
	printf("\t          \t\t4 == Remove atoms in factor highest Angstroms\n");
	printf("\t          \t\t5 == Remove atoms in factor lowest Angstroms\n");
	exit(-1);
    }
    strcpy(cfg_infile_, argv[1]);
    Chem_InitializePeriodicTable();
    pef_initialize();
    cfg_initialize();
    factor = atof(argv[2]);
    if (argv[3]) actionCode = atoi(argv[3]);
    else actionCode = 1;

    cfg_establish();

    atomList_zpossort(&L_);
    if (actionCode == 1)
    {
	fracToFix=factor;
	nToFix = (int)(fracToFix*nAtom_);
	if (nToFix)
	{
	    i=0;
	    for (a=L_;i<nToFix&&a;a=a->next)
	    {
		a->state=IS_FIXED;
		i++;
	    }
	    fprintf(stderr, "%% Fixed %i atoms\n", nToFix);
	}
    }
    else if (actionCode == 2)
    {
        atomList_reverse(&L_);
	fprintf(stderr, "%i hi atoms to be removed.\n", (int)factor);
	fflush(stderr);
	for (i = 0; i < (int)factor; i++)
	{
	    L_=L_->next;
	}
    }
    else if (actionCode == 3)
    {
	fprintf(stderr, "%i lo atoms to be removed.\n", (int)factor);
	for (i = 0; i < (int)factor; i++)
	{
	    L_=L_->next;
	}
    }
    else if (actionCode == 4)
    {
	fprintf(stderr, "Atoms in top %lf Angstroms to be removed.\n", factor);
        atomList_reverse(&L_);
	topz=L_->pos->z;
	while (topz-L_->pos->z>=factor) cfg_removeatom(&L_, NULL);
    }
    else if (actionCode == 5)
    {
	fprintf(stderr, "Atoms in bottom %lf Angstroms to be removed.\n", factor);
	fprintf(stderr, "Not implemented.\n");
    }
    else
    {
	printf("cfgFix:  actionCode [%i] is not recognized.\n", actionCode);
	exit(-1);
    }

    cfg_out(stdout);
    return 0;
       
}
