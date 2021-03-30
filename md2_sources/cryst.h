/* md series 2 (c) 1999 cam abrams
 * 
 * cryst.h -- header file for cryst.c
 */

#ifndef CRYST_H
#define CRYST_H

#include "point.h"

#define MAX_APERUC 16	    /* Maximum number of atoms per unit cell */

/* Number the unit cells: */
enum {DC_SILICON, /* Diamond Cubic Silicon, (100) face */ 
      OH_SILICON, /* Orthohexagonal Silicion, (111) face */
      DC_CARBON,  /* Diamond Cubic Carbon, (100) face */ 
      OH_CARBON,  /* Orthohexagonal Carbon, (111) face */
      NULL_UC};

typedef struct UC_DESC
{
    double a, b, c;	    /* lattice parameters */
    int ek[MAX_APERUC];	    /* element key */
    pt os[MAX_APERUC];	    /* (xyz)-offsets from center of unit cell
			     * as fractions of the lattice parameters */
} ucDesc_t;

void uc_a2i (char * ucStr);  /* Assigns uc_ based on the contents of
			      * string ucStr */
			     
void uc_initialize (void);  /* Computes lattice parameters from any
			     * user-supplied data */

int cryst_cnt (void);	    /* Based on the uc_ and nCell_ variables,
			     * returns the number of atoms in a simulation
			     * config of unit cell type uc with 
			     * nCell_.i * nCell_.j * nCell_.k unit
			     * cells */
void cryst_pos (void);	    /* Based on the nCell_ variable, assigns
			     * positions to a collection of atoms */
void cryst_vel (void);	    /* Based on the T_ variable, assigns
			     * thermal velocities to all atoms */

#endif
