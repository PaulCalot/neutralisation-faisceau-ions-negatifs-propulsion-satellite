/*
 * MD series2 -- Molecular Dynamics of Plasma-Surface Chemistry
 * (c) 1999 
 * Cameron Abrams
 * and 
 * The Regents of the University of California, Berkeley
 * 
 * Department of Chemical Engineering
 * University of California, Berkeley	1995-1999
 * 
 * Module Name:	clust
 * 
 */

#ifndef CLUST_H
#define CLUST_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "atom.h"
#include "point.h"
#include "genforce.h"
#include "ion.h"

#define MAXNUMCLUSTERS  2000

enum {CS_BASE, CS_ION, CS_OK, CS_TRASH, CS_UNUSED, CS_NULLSTATE};
typedef struct cluster_list_node *cNodePtr;
typedef struct cluster_list_node
{
    int id;
    short state;
    nNodePtr mem;
    cNodePtr next;
} cNode;


void clust_clear (void);	/* resets all members of the static array */

void clust_clean (void);	/* cleans out isolated clusters according
				 * to an "evaporation" rule. */

void clust_cfgreport (void);	/* reports all clusters in the cluster list
				 * in cfg format which are cleaned from
				 * the simulation cfg */

/* Functions used for force-routine-implemented cluster identification */
void clust_enclust (atomPtr, atomPtr);
				/* Given two atoms that are bound, ensure
				 * that they belong to the same cluster. */
void clust_sweepsingles (void);	/* Searches the simulations cfg for any atoms
				 * that do not belong to a cluster. These
				 * are deemed "monomers" and each is put
				 * into its own cluster */
void clust_fixlist (void);	/* "Build" the linked list of clusters
				 * from viable members of the static array */

char  * clust_name (char * space, cNodePtr c);
ptPtr clust_com (cNodePtr c);
ptPtr clust_comv (cNodePtr c);
double clust_k (cNodePtr c); /* returns KE of cluster */
double clust_b (cNodePtr c); /* returns external binding energy of cluster */

#endif

