/*
 * MD series -- Molecular Dynamics of Plasma-Surface Chemistry
 * (c) 1999 Cam Abrams, Dept. of Chemical Engineering
 * Regents of the University of California, Berkeley
 * 
 * Module Name:		atom 
 * Module Class:	base
 * 
 */

#ifndef ATOM_H
#define ATOM_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>
#include <string.h>
#include "point.h"
#include "chem.h"

/* on the SGI, FALSE and TRUE are undefined: no harm
 * in defining them here! */
#ifndef	TRUE
#define TRUE 1
#endif
#ifndef FALSE
#define FALSE 0
#endif

typedef struct theAtom *atomPtr;
typedef struct neighbor_list_node *nNodePtr;
 
/* nNode:  a Neighbor List entry */
typedef struct neighbor_list_node
{
    int		index;  /* index of this neighbor node */
    atomPtr	addr;   /* address of atom to which this node corresponds */
    nNodePtr	next;	/* pointer to next nNode in this nList */
    #ifndef TERSOFF_NEIGHBORS
    struct
    {
	unsigned int is_bound	:   1;
			/* this can be a bound neighbor (Tersoff-like) or
			 * either bound or unbound (Verlet-like) */
    } flags;
    nNodePtr	joint;	/* pointer to next nNode in a joint pair of nLists */
    #endif
} nNode;

enum {IS_UNDETERMINED, IS_FIXED, IS_SURF, IS_SPUTTERED, IS_UNCLEARED, 
      IS_PUSHTHRU, IS_BULK, IS_UNBOUND, IS_SCATTERED, IS_PROJECTILE, 
      ATOM_NSTATES};
#define XPLUS_ON    1
#define XMINUS_ON   2
#define YPLUS_ON    4
#define YMINUS_ON   8

typedef struct theAtom
{
    int		id;	/* unique index # of atom		*/
    sym_type	sym;	/* atomic symbol (a spec_table index)	*/
    ptPtr	pos;	/* pointer to pos_[id]			*/
    ptPtr	vel;	/* pointer to vel_[id]			*/
    ptPtr	hac;	/* pointer to hac_[id]			*/
    ptPtr	frc;	/* pointer to frc_[id]			*/
    double	mv2;	/* mass*velocity^2			*/
    double	q;	/* charge on/electron density at atom   */
    int		nNbrs;	/* number of neighbors on nList		*/
    nNodePtr	nList;	/* pointer to list of neighbors		*/
    nNodePtr	bCards; /* supply of business cards to "hand out"
			 * to neighbors	for their nLists	*/
    atomPtr	next;	/* ptr to the next atom			*/
    int		state;	/* see above				*/
    void	*c;	/* cluster designation	(a cNodePtr)	*/
    struct		/* bit-flags				*/
    {
	unsigned int is_RDF_a	:   1;	/* atom is RDF-active	*/
	unsigned int is_BAD_a	:   1;	/* atom is BAD-active	*/
	unsigned int is_PEF_a	:   1;	/* atom is PEF-dimer-active	*/
    } flags;
    short pbc;		/* 4-bit periodic boundary code.  This field
			 * is only used for sputtered/scattered atoms. */
} atom;

void atom_clear (atomPtr a);
void atom_copy (atomPtr a, atomPtr b);
void atom_out (FILE* fp, atomPtr a);
void atom_display (FILE* fp, atomPtr a);
void atom_becomeNeighbors (atomPtr a, atomPtr b, int * na, int * nb);
void atom_ungreet (atomPtr a);

void nNode_push (nNodePtr * L, nNodePtr a);
void atom_nreset (atomPtr a);
int nNode_n (nNodePtr L);
nNodePtr nNode_pop (nNodePtr * L);
nNodePtr nNode_remove (nNodePtr * L, nNodePtr a);
nNodePtr atom_bCardrecall (atomPtr v, atomPtr a);

void atomList_push (atomPtr a, atomPtr * L);
void atomList_remove (atomPtr a, atomPtr * L);
void atomList_indexsort (atomPtr * L);
void atomList_zpossort (atomPtr * L);
void atomList_reverse (atomPtr * L);
int atomList_sizeof (atomPtr L);
atomPtr atomList_pop (atomPtr * L);
atomPtr atomList_tailof (atomPtr L);

void atom_setState (atomPtr a, int State);
int atom_getState (atomPtr a);
char * atom_stateStr (atomPtr a);
char * atom_stateSymStr (atomPtr a);

#endif
