/*
 * MD series 2 -- Molecular Dynamics of Plasma-Surface Chemistry
 * (c) 1999 Cam Abrams, Dept. of Chemical Engineering
 * Regents of the University of California, Berkeley
 * 
 * Module Name:  genforce
 * 
 * Module Description:  generic force routine header
 */

#ifndef GENFORCE_H
#define GENFORCE_H

#include <stdio.h>
#include "atom.h"
#include "cfg.h"
#include "point.h"

void pef_initialize (void);
void pef_paramout (FILE* ofp);
void cfg_setforce (void);

double atom_peij (const atomPtr ip, const atomPtr jp);
double array_peij (int i, int j);
double atom_rcuij (sym_type is, sym_type js);
double atom_reij (sym_type is, sym_type js);
short atom_bound (const atomPtr ip, const atomPtr jp);

#endif
