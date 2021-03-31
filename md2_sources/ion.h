/*
 * MD series 2 -- Molecular Dynamics of Plasma-Surface Chemistry
 * (c) 1999 
 * Cameron Abrams
 * and 
 * The Regents of the University of California, Berkeley
 * 
 * Department of Chemical Engineering
 * University of California, Berkeley	1995-1999
 *
 * Module Name:		ion
 * 
 */

#ifndef ION_H
#define ION_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "point.h"
#include "dblmat.h"
#include "atom.h"
#include "sicf_params.h"
#include "cfg.h"
#include "genforce.h"

#define MAXATOMSPERION	5
#ifndef MAXCHAR
#define MAXCHAR 255
#endif

typedef struct IONPOOL
{
    int spec;		/* species index */
    double wt;		/* weight */
} ionPool_t;

void ion_scanpool (void);

typedef struct IONINFO
{
    /* Parameters set by/computed from users specifications */
    int spec;		/* species index */
    short pool;		/* flag used to specify whether ion is
			 * selected randomly from pool (1), or
			 * is set explicitly by the user (0) */
    double E_i;		/* incident energy */
    double Th_i;	/* incident polar angle */
    double P_i;		/* incident azimuthal angle */
    short pRand;	/* flag: randomize P_i */
    nNodePtr mem;	/* member list built from bCards 
			 * donated by atom members */
    double mass;	/* ion mass */
    pt r0;		/* computed COM (center of mass) position */
    pt r0set;           /* desired COM (center of mass) position */
    short r0Set;        /* flag: user spec's COM position */
    short rExplic;	/* flag: ion atom positions input explicitly by user */
    pt v0;		/* initial COM (center of mass) velocity */
    short vExplic;	/* flag: ion velocity is input explicitly by user */
    double a1, a2, a3;	/* Euler orientation angles (for polyatomic ions) */
    short aRand;	/* flag: randomize angles */
    
    double imp_ef;	/* prescribed kinetic energy fraction at "impact";
			 * default is 0.75.  This is used for
			 * monoatomic ions and polyatomic ions which aren't
			 * energetic enough to dissociate.
			 * "Impact" is also defined as the point at which a
			 * polyatomic ion loses its first internal bond. */
    
    /* Dynamic variables */
    double k;	   	/* Current kinetic energy of all atoms in ion */
    double t;		/* Current (internal) potential energy */
    double tex;		/* Current potential energy via interaction with
			 * atoms external to the ion */
    pt rx;		/* COM position at impact time */
    pt vx;		/* COM velocity at impact time */
    pt r;		/* Current COM position */
    pt v;		/* Current COM velocity */
    pt Dr;		/* Total vector displacement from initial position */
    double impt;	/* impact time */
    int nBonds0;	/* initial number of intact bonds */
    int nBonds;		/* current number of intact bonds */
    short imp_;		/* impact flag: 
			 * imp_ == 0: ion has not yet impacted
			 * imp_ == 1: ion has impacted */
    /* Explicit atomic positions and velocities that may be user-specified */
    pt r0ex[MAXATOMSPERION];
    pt v0ex[MAXATOMSPERION];
    
    /* output info */
    int uint;		/* time step interval of updates */
    int oint;		/* time step interval of outputs */
    short ofex;		/* flag: output file existence */
    char outfile[MAXCHAR];
    char outfile_0[MAXCHAR];
} ion_type;

void ion_newion (void);	    /* introduces a new ion as specified by the 
			     * ionInfo structure */
void ion_lowerion (void);   /* recomputes the z-position of a newly introduced
			     * ion based on a new highest atom in the 
			     * active atomList L_ */

void ion_update (void);	    /* updates all dynamical data for the ion */
void ion_output (void);	    /* outputs all dynamical data for the ion */

int ion_a2i (char * a);
char * ion_i2a (int i);


#endif
