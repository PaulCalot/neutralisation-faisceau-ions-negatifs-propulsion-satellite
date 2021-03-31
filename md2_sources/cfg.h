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
 * Module Name:		cfg
 * 
 */

#ifndef CFG_H
#define CFG_H
#include <time.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include "point.h"
#include "chem.h"
#include "atom.h"
#include "args.h"
#include "units.h"
#include "clust.h"
#include "ion.h"

#define MAXNUMATOMS	MXA
#define MAXNUMNEIGHBORS 100
#define MAXCHAR		255
#ifndef MAXBIN
#define MAXBIN		500
#endif

void cfg_initialize (void);	    /* Initializes the entire cfg data segment */
void cfg_establish (void);	    /* 1. Reads in cfg from file.
				     * 2. Identifies clusters.
				     * 3. Introduces any new atoms (like an ion).
				     * 4. Cleans up the cfg (removes sputtered
				     *    clusters, etc., if requested). */
void cfg_report (void);		    /* Outputs the current cfg */
void cfg_snapout (int idnum);	    /* Writes a cfg snapshot to a specially-
				     * named file where the numerals of idnum
				     * replace ?'s in the given file name */

void cfg_atomout (atomPtr a, FILE *fp);
				    /* Outputs the atom a in cfg format */

void cfg_clear (void);		    /* zeros all elements in the cfg */
void cfg_nreset (void);		    /* resets all nbr lists */
void cfg_creset (void);		    /* resets all cluster ids */

double cfg_v (void);		    /* returns volume of cfg */

atomPtr cfg_addr (int i);	    /* returns the address of the i'th
				     * atom in the static atom array. */
				     
void cfg_removeatom (atomPtr * a, int * i);   
				    /* Removes an atom from the atomList
				     * interface to the static atom array
				     * and places it on the "trashHeap"
				     * interface. The data in the static
				     * arrays is not moved, and no atom
				     * addresses are changed. (*a) must be
				     * a legitimate atomPtr.  If the contents
				     * of (*a) point to an existing member
				     * of the static atom array, then this is
				     * the atom that is transferred. 
				     * If "i" is not null and "*a" is, then 
				     * the atom with array index "*i" in the 
				     * static array of atoms is transferred, 
				     * and the address of this atom is returned
				     * in (*a). */

void cfg_deleteatom (atomPtr * a);  /* deletes atom pointed to by *a. */  

atomPtr cfg_newatom (void);	    /* returns a pointer to the next available
				     * atom data structure */

void cfg_returnnewatom (void);   /* undoes latest invocation of cfg_newatom */


void cfg_showsizes (void);	    /* Shows sizes of all data segments */

void RDF_normalize_output_hist();   /* normalizes and outputs the RDF */
/* void cfg_generate (FILE *fp); */

#endif

