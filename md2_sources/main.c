/* 
 * Molecular dynamics simulations of plasma-surface chemistry
 *
 * (c) 1999 Cameron F. Abrams, 
 *	    Department of Chemical Engineering, 
 *	    University of California, Berkeley
 * 
 * and	    The Regents of the University of California, Berkeley
 *
 *
 * md series 2 (c) 1999 cfa
 *
 * main.c:  multiple trajectory driver code. 
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "args.h"
#include "atom.h"
#include "cfg.h"
#include "chem.h"
#include "domd.h"

#ifdef VERSION_NUMBER
double version_=VERSION_NUMBER;
#else
double version_=1.0;
#endif
#ifdef BUILD_NUMBER
int build_=BUILD_NUMBER;
#else
int build_=1;
#endif

extern short quiet_, cfg_newcrystal_;
extern int impact_num_, impact0_num_, impact1_num_, nAtom_, nElem_[], Ndt_;
extern double rt_, tlim1_, tlim0_;
extern element per_table[];

static int *i=&(impact_num_), 
	   *i0=&(impact0_num_),
	   *i1=&(impact1_num_);
static time_t t0, t1;
static int sdt, dt;
static double adt;
static int j, n;
void time_info (void)
{
    t1=time(NULL);
    sdt+=(dt=difftime(t1, t0));
    adt=(double)sdt/((*i)-(*i0)+1);
    if (!quiet_) 
    {
	printf("# Tr %i, N %i, tm %i s = %.3lf m, <tm/tr> %.3lf s = %.3lf m, ", 
		*i, nAtom_, dt, dt/60., adt, adt/60.); fflush(stdout);
	printf("<tm/st/a> %.3le s, ", adt/(Ndt_+1)/nAtom_); fflush(stdout);
	printf("Final El.#: ");
	n=0;
	for(j=0;j<MAXNUMELEMENTS;j++) 
	{
	    if (nElem_[j]) 
	    {
		printf("%s %i ", per_table[j].sym, nElem_[j]);
		n+=nElem_[j];
	    }
	}
	printf("Total %i\n", n);	
	fflush(stdout);
    }
}

int main (int argc, char * argv[])
{
    /* Handle all command line arguments */
    arg_handler(argc, argv, TRAP_BADWORDS);
    
    /* Say hello */
    if (!quiet_) 
	printf("# Welcome to md series 2 v %.2lf b %i (c) 1999 cfa\n",
	    version_, build_);

     /* Perform initializations */
    Chem_InitializePeriodicTable();
    pef_initialize();
    cfg_initialize();
    filenames_initialize();
    if (!quiet_) {cfg_showsizes();filenames_echo();}
    sdt=adt=0;
    
    /* Begin the trajectory loop. */
    for ((*i)=(*i0);(*i)<=(*i1);(*i)++)
    {
	if (!quiet_) {printf("#\n# Trajectory #%i (%i out of %i) begins here\n", 
	    *i, (*i)-(*i0)+1, (*i1)-(*i0)+1);fflush(stdout);}
	t0=time(NULL);

	/* Perform the i'th trajectory */
	do_md(*i);
	
	/* Output some timing/results info */
	time_info();
	
	/* Update the integration limits */
    	tlim1_=rt_+(tlim1_-tlim0_);
	tlim0_=rt_;
    }

    /* Say goodbye */
    if (!quiet_) printf("# Thank you for using md series 2.\n");
    return 0;
}
