#ifndef CFG_CGMIN_H
#define CFG_CGMIN_H

#include "atom.h"
#include "point.h"
#include "cfg.h"
#include "genforce.h"

/* Performs conjugate gradient energy minimization on a configuration
 * of atoms.  Adapted from 

 Numerical Recipes in C, Second Edition, 
 by William H. Press,  Saul A. Teukolsky,  William T. Vetterling, 
    and Brian P. Flannery, 
 Cambridge University Press, 
 The Pitt Building,  Trumpington St., Cambridge CB2 1RPm
 40th West 20th St. New York, NY,  10011-4211
 and 10 Stamford Road, Oakleigh, Victoria 3166, Australia
 (c) 1992
 
 Chapter 10. pp394-455.
 
 The algorithm used here is Polak-Ribiere minimization method.
 
 *iter returns number of iterations;
 *fret returns the function value at the minimum;
 A is a control argument;  if NULL, all atoms in the cfg
 are allowed to move to find the minimum energy; if it
 points to an atom, then ONLY that atom is allowed to move. 
 */
void cfg_cgmin (int * iter, double * fret, atomPtr A);

#endif
