/*
 * md series 2 -- Molecular Dynamics of Plasma-Surface Chemistry
 * (c) 1999 Cam Abrams, Dept. of Chemical Engineering
 * Regents of the University of California, Berkeley
 * 
 * Module Name:  units
 * 
 * Module Description:  defines several macros for implementation of
 * the unit systems.
 * 
 */

#ifndef UNITS_H
#define UNITS_H

#define J_PER_EV		1.6022e-19
#define J_PER_HARTREE		4.3598e-19
#define KB			1.38066e-23	/* J/K */
#define KB_EVPERK		8.61733003e-5	/* eV/K */
#define KG_PER_AMU		1.66054e-27
#define GPCC_PER_AMUPCUA	1.66		/* (g/cc)/(amu/cu-Ang) */
#define E2_4PIE0		14.40034276	/* eV-Ang */
#define	APVKMASS_PER_AMU	1.036477518e-4	/* (eV.ps2/Ang2 per amu) */
#define GPA_PER_EVPERCUANG	160.2099

#define	ONE_THIRD		0.33333333333333333333333333333333333
#define	ONE_SIXTH		0.16666666666666666666666666666666667
#define	ONE_SEVENTH		0.14285714285714285714285714285714286
#define	ONE_NINTH		0.11111111111111111111111111111111111

/* Stillinger-Weber Silicon-based units: */
#define ANGSTROMS_PER_SIGMA     2.0951
#define EV_PER_EPSILON          2.1678
#define PICOSECS_PER_TAU        0.076634   
#define G_PER_M                 4.66362659e-23 
#define AMU_PER_M               28.085 
#define GperCC_PER_MperCUSIG    5.071185526
#define KELVIN_PER_SWT          25156.73798
#define CUANG_per_CUSIGMA       9.1963241

#define E2_EP0                  3.170524711 /*e^2/4(pi)(ep0) in sig-ep*/
/* E2_EP0 is the unit charge squared divided by (4*PI*epsilon0),
 * or when multiplied by formal charges of two bodies, it becomes
 * the constant of proportionality between inverse distance and
 * potential energy.  
 * The value given is in units of SW-epsilon(energy).SW-sigma(distance). */

#define SW2SI_T(A)              (A)*KELVIN_PER_SWT
#define SW2SI_DENS(A)           (A)*GperCC_PER_MperCUSIG
#define SW2SI_E_eV(A)           (A)*EV_PER_EPSILON
#define SW2SI_TIME_PS(A)        (A)*PICOSECS_PER_TAU
#define SW2SI_MASS_G(A)         (A)*G_PER_M
#define SW2SI_LENGTH_ANG(A)     (A)*ANGSTROMS_PER_SIGMA
#define SW2SI_MASS_AMU(A)       (A)*AMU_PER_M


#endif



