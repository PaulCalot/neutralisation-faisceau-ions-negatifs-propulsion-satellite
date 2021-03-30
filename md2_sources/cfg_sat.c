/*
 * cfg_sat:  locates and saturates dangling bonds in a cfg.
 * 
 * (c) 1999 Cameron F. Abrams
 * and
 * Regents of the University of California, Berkeley
 */

#include "cfg_sat.h"

/* External globals */
extern atomPtr L_;  /* the configuration */
extern int nAtom_;  /* number of atoms in the configuration */
extern double PE_;  /* total potential energy */
extern double Eb_;  /* binding energy at specified setpoint T */
extern element per_table[]; /* periodic table */
extern short quiet_;

static atomPtr ip, jp;
static pt Rtmp, *r=&Rtmp;
static int iter;
static double Etrial, dE, Einit;
void cfg_sat (int nTrials, int * nSucc)
{
    short locate_tpos (ptPtr r, atomPtr a);
    int i, j, nT;
    
    nT=0;
    *nSucc=0;
    cfg_setforce();		/* initial call to set neighbor lists */
    printf("Pretrials:  PE = %.5le,  will perform %i trials\n", Einit=PE_,  nTrials);
    /* An nTrials of 0 means try every atom.  A finite nTrials
     * means only make that many trials. */
    for (ip=L_;ip&&(nTrials?nT<nTrials:1);ip=ip->next)
    {
	/* Locate an unoccupied tpos on atom ip */
	if ((ip->sym==C||ip->sym==Si)&&ip->state!=IS_FIXED)
	{
	    while (locate_tpos(r,ip))
	    {
		printf("Trial %i: Located tpos on atom %s_%i at %.5le %.5le %.5le\n", 
		    nT, per_table[ip->sym].sym, ip->id, r->x, r->y, r->z);
		jp=cfg_newatom();	/* create a new atom */
		jp->sym=F;		/* make it fluorine */
		ptPtr_copy(jp->pos, r); /* put it at position r */
		Einit=PE_;
		cfg_setforce();		/* compute total PE and all forces */
		/* Perform the conjugate gradient minimization
		 * allowing only the position of the new atom
		 * to vary. */
		cfg_cgmin(&iter,&Etrial,jp); 
		dE=Etrial-Einit;	/* compute change in potential energy */
		if (dE<-Eb_) 
		{
		    /* successful binding if change
		     * is lower than thermal binding energy */
		    (*nSucc)++;
		    printf("Trial %i successful, dE %.5le nAtom_ %i\n", 
			nT, dE, nAtom_);		
		    Einit+=dE;
		}
		else
		{
		    /* delete the new atom */
		    cfg_returnnewatom();
		    printf("\tTrial %i unsuccessful, dE %.5le nAtom_ %i\n", 
			nT, dE, nAtom_);		
		}
		nT++;
	    }
	}
    }
}

static nNodePtr np;
static int nn, db;
short locate_tpos (ptPtr r, atomPtr a)
{
    short locate_tpos3 (ptPtr r, atomPtr a);
    short locate_tpos2 (ptPtr r, atomPtr a);
    short locate_tpos1 (ptPtr r, atomPtr a);
    if (!a) return 0;
    /* count number of neighbors */
    nn=0;
    for (np=a->nList;np;np=np->next) nn++;
    db=4-nn;
    switch (db)
    {
	case 4:	printf("Error: atom %s_%i is solitary!\n", 
			per_table[a->sym].sym, a->id);
		return 0;
	case 3: return locate_tpos3(r, a);
	case 2: return locate_tpos2(r, a);
	case 1: return locate_tpos1(r, a);
	default: return 0;
    }
}

static pt TPT, *t=&TPT;
static double tl;
short locate_tpos1 (ptPtr r, atomPtr a)
{
    /* compute the centroid of the triangle formed by the
     * positions of the three neighboring atoms */
    ptPtr_clear(t);
    for (np=a->nList;np;np=np->next)
	ptPtr_add(t, t, np->addr->pos);
    ptPtr_scalmult(t, t, ONE_THIRD);
    /* compute vector from centroid to a->pos */
    ptPtr_subtract(t, a->pos, t);
    /* scale t to equilibrium bond length */
    tl=ptPtr_abs(t);
    ptPtr_scalmult(t, t, atom_reij(a->sym,F)/tl);
    /* r = a->pos + t sets final location of the tetrahedral position */
    ptPtr_add(r, t, a->pos);
    ptPtr_minimg(r, r);
    return 1;
}

short locate_tpos2 (ptPtr r, atomPtr a) {return 0;}
short locate_tpos3 (ptPtr r, atomPtr a) {return 0;}
