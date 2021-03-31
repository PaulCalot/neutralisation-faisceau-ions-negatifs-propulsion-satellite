/*
 * MD series 2 -- Molecular Dynamics of Plasma-Surface Chemistry
 * (c) 1999 Cam Abrams, Dept. of Chemical Engineering
 * Regents of the University of California, Berkeley
 * 
 * Module Name:  tbt_sicf
 * 
 * Module Description:  implements the force routine and parameter assignment
 * and handling routines to utilize the Tersoff/Brenner Silicon-Carbon PEF, 
 * with the F modifications of Tanaka and myself, 
 * including the high-energy, near-separation cubic splines.
 * 
 * Adapted from "Empirical potentials for C-Si-H systems with application
 * to C60 interactions with Si crystal surfaces", Keith Beardmore and
 * Roger Smith,  Phil. Mag. A, 1996, Vol. 74, No. 6, 1439-1466.
 * 
 */
#ifndef TBT_SICF_H
#define TBT_SICF_H
#include "atom.h"

double atom_getrij (const atomPtr ip, const atomPtr jp);
double atom_getfij (const atomPtr ip, const atomPtr jp);
void cfg_frc2tau (void);


#endif

