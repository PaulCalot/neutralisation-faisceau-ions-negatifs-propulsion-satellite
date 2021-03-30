#ifndef POINT_H
#define POINT_H
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct POINT * ptPtr;
typedef struct POINT
{
    double x;
    double y;
    double z;
} pt;
typedef struct IPOINT * iptPtr;
typedef struct IPOINT
{
    int i;
    int j;
    int k;
} ipt;


double ptPtr_abs (ptPtr p);
double ptPtr_sqabs (ptPtr p);
double ptPtr_dotprod (ptPtr a, ptPtr b);

ptPtr ptPtr_create (double x, double y, double z);
ptPtr ptPtr_cross (ptPtr res, ptPtr a, ptPtr b);
ptPtr ptPtr_add (ptPtr res, ptPtr a, ptPtr b);
ptPtr ptPtr_subtract (ptPtr res, ptPtr a, ptPtr b);
ptPtr ptPtr_shift (ptPtr res, ptPtr a, double x, double y, double z);
ptPtr ptPtr_scalmult (ptPtr res, ptPtr a, double x);
ptPtr ptPtr_vectmult (ptPtr res, ptPtr a, ptPtr b);
ptPtr ptPtr_minimg (ptPtr res, ptPtr a);

void ptPtr_clear (ptPtr a);
void ptPtr_copy (ptPtr res, ptPtr a);
void ptPtr_swap (ptPtr a, ptPtr b);
void ptPtr_out (FILE * fp, ptPtr p);
void ptPtr_destroy (ptPtr p);

#endif

