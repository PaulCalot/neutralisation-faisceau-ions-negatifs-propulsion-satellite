/*
 * MD series -- Molecular Dynamics of Plasma-Surface Chemistry
 * (c) 1996-1998 Cam Abrams, Dept. of Chemical Engineering
 * Regents of the University of California, Berkeley
 * 
 * Module Name:         tricubic
 * Module Class:        BrennerPEF
 * Module Description:  routines for tricubic interpolation of the Fcc
 *                      function in the Brenner/Tanaka PEF
 * 
 * Note:  Junichi Tanaka did all of the hard work for this module.
 * I'm just borrowing it from him.  I have added only comments, and
 * changed one function name so that it is compatible with my version
 * of the Brenner/Tersoff/Tanaka PEF force routine.
 * 
 */
#ifndef TRICUBIC_H
#define TRICUBIC_H

#include <stdio.h>
#include <stdlib.h>

/* If the macro "BUILD_NUMBER" is not defined, then this compilation is
 * independent of the MD series, and the declarations and macros in
 * "bicubic.h" are not available here. */
#ifndef BUILD_NUMBER
#define	X1_NGRIDPOINTS	11
#define X2_NGRIDPOINTS	11
#define X1_NGRIDSQUARES (X1_NGRIDPOINTS-1)
#define X2_NGRIDSQUARES (X2_NGRIDPOINTS-1)
void * xmalloc(size_t n);
#else
/* bicubic.h contains
 *  (1) the declaration of xmalloc();
 *  (2) the macro definitions for X1 and X2 grid sizes
 */
#include "bicubic.h"
#endif

#define	X3_NGRIDPOINTS	11
#define X3_NGRIDSQUARES (X3_NGRIDPOINTS-1)

/* tcucof:  (adapted from "bcucof", Numerical Recipes, 2nd Ed., p 126)
 * Generates the 64-member coefficient matrix c_(ijk), i=1->4, j=1->4, k=1->4.
 * c_(ijk) is used in the "tcuint" function to interpolate a function
 * value at a specified off-lattice point.
 */
void tcucof(double y[],	/* function values at the 8 corners of the grid cube */
	    double y1[], double y2[], double y3[], /* 1st derivatives */
	    double y12[], double y23[], double y31[], /* cross derivatives */
	    double y123[], /* tertiary cross derivatives */
            double d1, double d2, double d3, /* grid cube dimensions */
	    double c[4][4][4]);

/* tcuint:  (adapted from "bcuint", Numerical Recipes, 2nd Ed., p 127)
 * A grid cube is defined by the six coordinates x1u, x1l, x2u, x2l, x3u, x3l.
 * The 4x4x4 matrix of coefficients c_(ijk) is passed in as c[][][].  The
 * point at which an interpolated function value is wished is
 * (x1, x2, x3).  The interpolated function value is returned in *ansy.
 * Its first derivative with respect to x1 is returned in *ansy1.
 * Its first derivative with respect to x2 is returned in *ansy2.
 * Its first derivative with respect to x3 is returned in *ansy3.
 */
void tcuint(double x1l, double x1u, double x2l, double x2u, 
     double x3l, double x3u, 
     double x1, double x2, double x3, 
     double *ansy, double *ansy1, double *ansy2, double *ansy3, 
     double c[4][4][4]);

/* tricubic_genCoef:  Given the values of the function and its
 * first derivatives and cross derivatives at each grid point, 
 * this function computes the set of all 4x4x4 coefficient matrices
 * c_(ijk) for each grid cube. */
void tricubic_genCoef (double y[X1_NGRIDPOINTS][X2_NGRIDPOINTS][X3_NGRIDPOINTS],
		       double y1[X1_NGRIDPOINTS][X2_NGRIDPOINTS][X3_NGRIDPOINTS], 
		       double y2[X1_NGRIDPOINTS][X2_NGRIDPOINTS][X3_NGRIDPOINTS],
		       double y3[X1_NGRIDPOINTS][X2_NGRIDPOINTS][X3_NGRIDPOINTS],
		       double y12[X1_NGRIDPOINTS][X2_NGRIDPOINTS][X3_NGRIDPOINTS], 
		       double y23[X1_NGRIDPOINTS][X2_NGRIDPOINTS][X3_NGRIDPOINTS], 
		       double y31[X1_NGRIDPOINTS][X2_NGRIDPOINTS][X3_NGRIDPOINTS], 
		       double y123[X1_NGRIDPOINTS][X2_NGRIDPOINTS][X3_NGRIDPOINTS], 
		       double (*****c)[4][4][4]);

/* Fcc_genCoef:  Initialize the set of 4x4x4 coefficient matrices for
 * the Fcc function.  This is the place where the spline knots are
 * specified. */
void Fcc_genCoef (void);

/* Fcc_tricubicint:  Tricubic interpolation of F_(cc) at point (x1,x2,x3).
 * The value of Fcc at the specified point is returned in *y,
 * d(Fcc)/d(x1) is returned in *y1, d(Fcc)/d(x2) is returned in *y2,
 * d(Fcc)/d(x3) is returned in *y3. */
void Fcc_tricubicint (double x1, double x2, double x3,
		      double *y, double *y1, double *y2, double *y3);

#endif
