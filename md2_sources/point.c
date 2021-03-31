#include "point.h"

double ptPtr_abs (ptPtr p)
{
    if (!p) return 0.0;
    return sqrt(p->x*p->x + p->y*p->y + p->z*p->z);
}

double ptPtr_sqabs (ptPtr p)
{
    if (!p) return 0.0;
    return (p->x*p->x + p->y*p->y + p->z*p->z);
}

double ptPtr_dotprod (ptPtr a, ptPtr b)
{
    if (!a || !b) return 0.0;
    return (a->x*b->x + a->y*b->y + a->z*b->z);
}

ptPtr ptPtr_create (double x, double y, double z)
{
    ptPtr p = malloc(sizeof(pt));
    if (p)
    {
	p->x = x;
	p->y = y;
	p->z = z;
    }
    else {printf("memalloc error: pt\n");exit(0);}
    return p;
}

ptPtr ptPtr_cross (ptPtr res, ptPtr a, ptPtr b)
{
    if (!a || !b) return res;
    if (res==a||res==b) 
    {
	printf("Error -- illegal overwrite in ptPtr_cross\n");
	return res;
    }
    if (!res) res = ptPtr_create(0.0, 0.0, 0.0);
    res->x = a->y*b->z - b->y*a->z;
    res->y = a->z*b->x - b->z*a->x;
    res->z = a->x*b->y - b->x*a->y;
    return res;
}

ptPtr ptPtr_add (ptPtr res, ptPtr a, ptPtr b)
{
    if (!a || !b) return res;
    if (!res) res = ptPtr_create(0.0, 0.0, 0.0);
    res->x = a->x+b->x;
    res->y = a->y+b->y;
    res->z = a->z+b->z;
    return res;
}

ptPtr ptPtr_subtract (ptPtr res, ptPtr a, ptPtr b)
{
    if (!a || !b) return res;
    if (!res) res = ptPtr_create(0.0, 0.0, 0.0);
    res->x = a->x-b->x;
    res->y = a->y-b->y;
    res->z = a->z-b->z;
    return res;
}

ptPtr ptPtr_shift (ptPtr res, ptPtr a, double x, double y, double z)
{
    if (!a) return res;
    if (!res) res = ptPtr_create(0.0, 0.0, 0.0);
    res->x=a->x+x;
    res->y=a->y+y;
    res->z=a->z+z;
    return res;
}

ptPtr ptPtr_scalmult (ptPtr res, ptPtr a, double x)
{
    if (!a) return res;
    if (!res) res = ptPtr_create(0.0, 0.0, 0.0);
    res->x=a->x*x;
    res->y=a->y*x;
    res->z=a->z*x;
    return res;
}

ptPtr ptPtr_vectmult (ptPtr res, ptPtr a, ptPtr b)
{
    if (!a || !b) return res;
    if (!res) res = ptPtr_create(0.0, 0.0, 0.0);
    res->x = a->x*b->x;
    res->y = a->y*b->y;
    res->z = a->z*b->z;
    return res;
}

extern pt Lr_, half_Lr_;
extern ipt per_;
ptPtr ptPtr_minimg (ptPtr res, ptPtr a)
{
    int shift = 0;
    /* shift is one of -1, 0, 1 */
    /* the following three adjustments implement the minimum image
     * convention for a set cube size (given by the global Lr_) */

    if (res!=a) ptPtr_copy(res, a);
    shift = per_.i*res->x/half_Lr_.x;
    if (shift == 2) shift = 1;
    if (shift == -2) shift = -1;
    if (shift) res->x -= Lr_.x * shift;
    shift = per_.j*res->y/half_Lr_.y;
    if (shift == 2) shift = 1;
    if (shift == -2) shift = -1;
    if (shift) res->y -= Lr_.y * shift;
    shift = per_.k*res->z/half_Lr_.z;
    if (shift == 2) shift = 1;
    if (shift == -2) shift = -1;
    if (shift) res->z -= Lr_.z * shift;
    
    return res;
}

void ptPtr_clear (ptPtr a)
{
    if (!a) return;
    a->x=a->y=a->z=0.0000000000000000L;
}

void ptPtr_copy (ptPtr res, ptPtr a)
{
    if (!a) return;
    if (!res) res = ptPtr_create(0.0, 0.0, 0.0);
    res->x=a->x;
    res->y=a->y;
    res->z=a->z;
}

void ptPtr_swap (ptPtr a, ptPtr b)
{
    double x;
    if (!a || !b) return;
    x=a->x;
    a->x=b->x;
    b->x=x;
    x=a->y;
    a->y=b->y;
    b->y=x;
    x=a->z;
    a->z=b->z;
    b->z=x;
}

void ptPtr_out(FILE * fp, ptPtr p)
{
    if (!fp || !p) return;
    fprintf(fp, "%.5le, %.5le, %.5le\n", p->x, p->y, p->z);
}

void ptPtr_destroy (ptPtr p)
{
    if (!p) return;
    free((ptPtr)p);
}
