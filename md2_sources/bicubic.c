/*
 * md series 2 -- Molecular Dynamics of Plasma-Surface Chemistry
 * (c) 1999 Cam Abrams, Dept. of Chemical Engineering
 * Regents of the University of California, Berkeley
 * 
 * Module Name:		bicubic
 * Module Class:	BrennerPEF
 * Module Description:	routines for bicubic interpolation of the Hij
 *			function in the Brenner/Tanaka PEF
 * 
 * Note:  Junichi Tanaka did all of the hard work for this module.
 * I'm just borrowing it from him.  I have added only comments, and
 * changed function names so that it is compatible with my version
 * of the Brenner/Tersoff/Tanaka PEF force routine.
 * 
 */
#include "bicubic.h"

static const int mp = (X1_NGRIDPOINTS);
static const int np = (X2_NGRIDPOINTS);
static const int ms = (X1_NGRIDSQUARES);
static const int ns = (X2_NGRIDSQUARES);

int BICUBIC_DIAG_=0;
short hcc_opt_=0;
short hcf_opt_=0;

static	double (***cCC)[4][4];
static	double (***cCH)[4][4];
static	double (***cCF)[4][4];
/* A ptr to...__|||    |____|
 *	         ||       |_______________________ 
 *               ||                               |
 *		 \|__ an "m"x"n" 2D array of...   |___ 4x4 2D arrays.
 *
 *  Here, we state that the 2D grid of points ("knots") which define the
 *  function we are interpolating is "m+1"x"n+1" and has, therefore, 
 *  m x n 4-membered grid squares.  Each one of these grid squares 
 *  has associated with it a 4x4 matrix of coefficients, c_(ij), 
 *  which are functions of the function values and its derivatives
 *  at the four gridpoints in that square.  This coefficient matrix
 *  c_(ij) is used in the bicubic interpolation routine.  We wish to
 *  initially compute *all* 4x4 matrices (c_(ij))_(mn) initially, and
 *  store them for use by the interpolation routine.  In order to store
 *  the m x n 4x4 arrays, we declare "c" as a "(pointer to a) 
 *  matrix of matrices",  or (*)(**c)[4][4], and dynamically allocate
 *  the "m" rows of "n" columns each using the xmalloc function
 *  given in Numerical Recipes.
 *
 */

void * xmalloc(size_t n)
{
    void *p;

    if ((p=malloc(n)) == NULL)
    {
	fprintf(stderr, "xmalloc:memory allocation\n");
	exit(1);
    }
    return p;
}

/* bcucof:  (Numerical Recipes, 2nd Ed., p 126)  Generates the
 * 16-member coefficient matrix c_(ij), i=1->4, j=1->4.  c_(ij)
 * is used in the "bcuint" function to interpolate a function
 * value at a specified off-lattice point.
 */
void bcucof(double y[], double y1[], double y2[], double y12[],
            double d1, double d2, double c[4][4])
{
    static int wt[16][16] = /* Weighting factors; pretabulated long long ago */
    {
	1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
	-3, 0, 0, 3, 0, 0, 0, 0,-2, 0, 0,-1, 0, 0, 0, 0,
	2, 0, 0,-2, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0,
	0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
	0, 0, 0, 0,-3, 0, 0, 3, 0, 0, 0, 0,-2, 0, 0,-1,
	0, 0, 0, 0, 2, 0, 0,-2, 0, 0, 0, 0, 1, 0, 0, 1,
	-3, 3, 0, 0,-2,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0,-3, 3, 0, 0,-2,-1, 0, 0,
	9,-9, 9,-9, 6, 3,-3,-6, 6,-6,-3, 3, 4, 2, 1, 2,
	-6 ,6,-6, 6,-4,-2, 2, 4,-3, 3, 3,-3,-2,-1,-1,-2,
	2,-2, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 2,-2, 0, 0, 0, 0, 0, 0,
	-6, 6,-6, 6,-3,-3, 3, 3,-4, 4, 2,-2,-2,-2,-1,-1,
	4,-4, 4,-4, 2, 2,-2,-2, 2,-2,-2, 2, 1, 1, 1, 1
    };
    int i, j, k, l;
    double xx, d1d2, cl[16], x[16];
    
    #if DIAG
    if (BICUBIC_DIAG_)
    {
	printf("bcucof begins\n");
	fflush(stdout);
    }
    #endif
    
    
    d1d2 = d1*d2;
    for (i = 0; i < 4; i++) /* Build a temporary vector */
    {
	x[i] = y[i];
	x[i+4] = y1[i]*d1;
	x[i+8] = y2[i]*d2;
	x[i+12] = y12[i]*d1d2;
    }
    for (i = 0; i < 16; i++) /* Multiply the matrix wt[][] by vector x[] and 
				store result in vector cl[] */
    {
	xx = 0.0;
	for (k = 0; k < 16; k++) xx += wt[i][k]*x[k];
	cl[i] = xx;
    }
    l = 0;
    for (i = 0; i < 4; i++) /* Place members of cl[] into their
			       proper locations in c[][]. */
	for (j = 0; j < 4; j++) c[i][j] = cl[l++];
 
    #if DIAG
    if (BICUBIC_DIAG_)
    {
	printf("bcucof ends\n");
	fflush(stdout);
    }
    #endif
 }

void bcuint(double x1l, double x1u, double x2l, double x2u, 
     double x1, double x2, double *ansy, double *ansy1, double *ansy2,
     double c[4][4])
{
    int i;
    double t, u, d1, d2;
    
    d1 = x1u - x1l;
    d2 = x2u - x2l;
    if (!d1 || !d2)
    {
	fprintf(stderr, "2D Cubic Interpolation: Bad gridpoint coords:\n");
	fprintf(stderr, "\tx1u(%.5e) x1l(%.5e) x2u(%.5e) x2l(%.5e)\n", 
	    x1u, x1l, x2u, x2l);
	fprintf(stderr, "Program exits\n");
	exit(0);
    }
    
    t = (x1-x1l)/d1;
    u = (x2-x2l)/d2;
    *ansy = (*ansy2) = (*ansy1) = 0.0;
    for (i = 3; i >= 0; i--)
    {
	*ansy = t*(*ansy) + ((c[i][3]*u + c[i][2])*u + c[i][1])*u + c[i][0];
	*ansy2 = t*(*ansy2) + (3.0*c[i][3]*u + 2.0*c[i][2])*u + c[i][1];
	*ansy1 = u*(*ansy1) + (3.0*c[3][i]*t + 2.0*c[2][i])*t + c[1][i];
    }
    *ansy1 /= d1;
    *ansy2 /= d2;
}

void bicubic_genCoef (double y[X1_NGRIDPOINTS][X2_NGRIDPOINTS],
		      double y1[X1_NGRIDPOINTS][X2_NGRIDPOINTS], 
		      double y2[X1_NGRIDPOINTS][X2_NGRIDPOINTS],
		      double y12[X1_NGRIDPOINTS][X2_NGRIDPOINTS], 
		      double (****c)[4][4])
{
    int i, j;
    double z[4], z1[4], z2[4], z12[4];
    
    #if DIAG
    if (BICUBIC_DIAG_)
    {
	printf("bicubic_genCoef begins\n");
	fflush(stdout);
    }
    #endif
    
    /* Initialize *c as a matrix of matrices */
    (*c) = (double (***)[4][4])xmalloc(ms*sizeof(double (**)[4][4]));
    for (i = 0; i < ms; i++)
    {
	(*c)[i] = (double (**)[4][4])xmalloc(ns*sizeof(double (*)[4][4]));
	for (j = 0; j < ns; j++)
	    ((*c)[i])[j] = (double (*)[4][4])xmalloc(16*sizeof(double));
    }
    
    /* For each grid square (i,j) compute the 4x4 matrix c_(ij).
     * 
     *		(4)------(3)
     *		 |        |
     *           |        |
     *           |        |
     *          (1)------(2)
     */
    for (i = 0; i < ms; i++)
    {
	for (j = 0; j < ns; j++)
	{
	    /* The array z[] holds the four function values at the
	     * corners of this grid square. */
	    z[0] = y[i][j];
	    z[1] = y[i+1][j];
	    z[2] = y[i+1][j+1];
	    z[3] = y[i][j+1];
	    /* The array z1[] holds the first derivative of the function
	     * with respect to x1 at each of the corners of this grid square.
	     */
	    z1[0] = y1[i][j];
	    z1[1] = y1[i+1][j];
	    z1[2] = y1[i+1][j+1];
	    z1[3] = y1[i][j+1];
	    /* The array z2[] holds the first derivative of the function
	     * with respect to x2 at each of the corners of this grid square.
	     */
	    z2[0] = y2[i][j];
	    z2[1] = y2[i+1][j];
	    z2[2] = y2[i+1][j+1];
	    z2[3] = y2[i][j+1];
	    /* The array z12[] holds the cross derivative of the function
	     * at each of the corners of this grid square. */
	    z12[0] = y12[i][j];
	    z12[1] = y12[i+1][j];
	    z12[2] = y12[i+1][j+1];
	    z12[3] = y12[i][j+1];
	    
	    /* Use bcucof() to compute the 4x4 coefficient matrix c_(ij) */
	    bcucof(z, z1, z2, z12, 1.0, 1.0, *((*c)[i][j]));
	}
    }
    #if DIAG
    if (BICUBIC_DIAG_)
    {
	printf("bicubic_genCoef ends\n");
	fflush(stdout);
    }
    #endif
}

void Hcc_genCoef (void)
{
    int i, j;
    double h[X1_NGRIDPOINTS][X2_NGRIDPOINTS]; 
    double h1[X1_NGRIDPOINTS][X2_NGRIDPOINTS]; 
    double h2[X1_NGRIDPOINTS][X2_NGRIDPOINTS]; 
    double h12[X1_NGRIDPOINTS][X2_NGRIDPOINTS];
    
    for (i = 0; i < X1_NGRIDPOINTS; i++)
	for (j = 0; j < X2_NGRIDPOINTS; j++)
	    h[i][j] = h1[i][j] = h2[i][j] = h12[i][j] = 0.0;

    /* Set values of Hcc at the various integer grid points. */
    /* Using brenner CC parameters. */
    #ifdef USE_BRENNER_HCC
    h[1][1] = -0.0226;
    h[2][0] = -0.0061;
    h[2][1] = 0.0160;
    h[3][0] = 0.0173;
    h[1][2] = 0.0149;
    #endif
    
    /* We allow the user to choose Tanaka Set I or
     * Tanaka Set II. */
    if (hcc_opt_==1)
    {
	/* These are the values J. Tanaka emailed to me on 14 Apr 1999. */
	h[1][1] = 0.0090;   /* hexafluorobenzene */
	h[2][0] = 0.0454;   /* C2F4 (I think) */
	h[2][1] = 0.1060;   /* perfluorocyclohexane */
	h[3][0] = 0.0891;   /* perfluoroethane */
    }
    else if (hcc_opt_==0) /* make this default! */
    {
	/* These are the values as they appear in a draft of Tanaka99b. */
	h[1][1] =  0.002;
	h[2][1] = -0.028;
	h[3][0] = -0.03783;
    }

    /* Derivative values by centered finite difference */
    h1[1][1] = 0.5*(h[2][1] - h[0][1]);
    h1[2][0] = 0.5*(h[3][0] - h[1][0]);

    /* Create the set of 4x4 coefficient matrices */
    bicubic_genCoef(h, h1, h2, h12, &cCC);
}

void Hcf_genCoef (void)
{
    int i, j;
    double h[X1_NGRIDPOINTS][X2_NGRIDPOINTS]; 
    double h1[X1_NGRIDPOINTS][X2_NGRIDPOINTS]; 
    double h2[X1_NGRIDPOINTS][X2_NGRIDPOINTS]; 
    double h12[X1_NGRIDPOINTS][X2_NGRIDPOINTS];
    
    #if DIAG
    if (BICUBIC_DIAG_)
    {
	printf("Hcf_genCoef begins\n");
	fflush(stdout);
    }
    #endif
    
    for (i = 0; i < X1_NGRIDPOINTS; i++)
	for (j = 0; j < X2_NGRIDPOINTS; j++)
	    h[i][j] = h1[i][j] = h2[i][j] = h12[i][j] = 0.0;

    /* Set values of Hcf at the various integer grid points. */
    /* Values that Junichi Tanaka determined. */
    /* We allow the user to choose the parameter set at runtime. */
    if (hcf_opt_==1)
    {
	/* These are the values J. Tanaka emailed to me on 14 Apr 1999. */
	h[0][1] =  0.10643;	    /* FCCF */
	h[0][2] =  0.1389;	    /* hexafluorobenzene */
	h[0][3] =  0.0;	    /* isobutane */
	h[1][0] =  0.06353;	    /* CF2 */
	h[1][1] =  0.1892;	    /* C2F4 */
	h[1][2] = -0.0021;	    /* perfluorocyclohexane */
	h[2][0] =  0.21895;	    /* CF3 */
	h[2][1] =  0.0104;	    /* C2F6 */
	h[3][0] =  0.01806;	    /* CF4 */
    }
    else if (hcf_opt_==0)
    {
	/* These are the values as they appear in a draft of Tanaka99b. */
	h[1][0] = 0.01582;	    /* CF2 */
	h[2][0] = 0.01819;	    /* CF3 */
	h[3][0] = -0.1277;	    /* CF4 */
	h[0][1] = -0.0002;	    /* CF b.e. in FCCF */
	h[0][2] = 0.0100;	    /* CF b.e. in C6F6 */
	h[1][1] = -0.0454;	    /* CF b.e. in F2CCF2 */
	h[2][1] = -0.1052;	    /* CF b.e. in F3CCF3 */
	h[1][2] = -0.055;	    /* CF b.e. in c-C6F12 */
    } 
    /* Derivative values by centered finite difference */
    h1[1][1] = 0.5*(h[2][1] - h[0][1]);
    h1[1][0] = 0.5*(h[2][0] - h[0][0]);
    h1[2][0] = 0.5*(h[3][0] - h[1][0]);
    h2[0][1] = 0.5*(h[0][2] - h[0][0]);
    h2[0][2] = 0.5*(h[0][3] - h[0][1]);
    h2[1][1] = 0.5*(h[1][2] - h[1][0]);
    
    /* Create the set of 4x4 coefficient matrices */
    bicubic_genCoef(h, h1, h2, h12, &cCF);
    #if DIAG
    if (BICUBIC_DIAG_)
    {
	printf("Hcf_genCoef ends\n");
	fflush(stdout);
    }
    #endif
}

void Hch_genCoef (void)
{
    int i, j;
    double h[X1_NGRIDPOINTS][X2_NGRIDPOINTS]; 
    double h1[X1_NGRIDPOINTS][X2_NGRIDPOINTS]; 
    double h2[X1_NGRIDPOINTS][X2_NGRIDPOINTS]; 
    double h12[X1_NGRIDPOINTS][X2_NGRIDPOINTS];
    
    for (i = 0; i < X1_NGRIDPOINTS; i++)
	for (j = 0; j < X2_NGRIDPOINTS; j++)
	    h[i][j] = h1[i][j] = h2[i][j] = h12[i][j] = 0.0;

    /* Set values of Hch at the various integer grid points. */
    h[0][1] =  -0.2479;
    h[0][2] =  -0.3221;
    h[0][3] =  -0.4460;
    h[1][0] =  -0.0984;
    h[1][1] =  -0.3344;
    h[1][2] =  -0.4449;
    h[2][0] =  -0.2878;
    h[2][1] =  -0.4438;
    h[3][0] =  -0.4507;
         
    /* Derivative values by centered finite difference */
    h1[1][1] = 0.5*(h[2][1] - h[0][1]);
    h1[1][0] = 0.5*(h[2][0] - h[0][0]);
    h1[2][0] = 0.5*(h[3][0] - h[1][0]);
    h2[0][2] = 0.5*(h[0][3] - h[0][1]);
    h2[1][1] = 0.5*(h[1][2] - h[1][0]);
    
    /* Create the set of 4x4 coefficient matrices */
    bicubic_genCoef(h, h1, h2, h12, &cCH);
}

void Hcc_bicubicint (double x1, double x2,
		     double *y, double *y1, double *y2)
{
    int i, j;
    double x1max = ms-1.e-10, x2max = ms-1.e-10;
    x1 = (x1 > x1max ? x1max : x1);
    x2 = (x2 > x2max ? x2max : x2);
    
    /* Perform the interpolation on a unit square:  This makes
     * "t" and "u" 1.0 in the bcuint routine and thus saves time
     * by simplifying multiplication operations. */
    i = x1;
    j = x2;
    bcuint(0.0, 1.0, 0.0, 1.0, x1-i, x2-j, y, y1, y2, *(cCC[i][j]));
}

void Hcf_bicubicint (double x1, double x2,
		     double *y, double *y1, double *y2)
{
    int i, j;
    double x1max = ms-1.e-10, x2max = ms-1.e-10;
    x1 = (x1 > x1max ? x1max : x1);
    x2 = (x2 > x2max ? x2max : x2);
    
    /* Perform the interpolation on a unit square:  This makes
     * "t" and "u" 1.0 in the bcuint routine and thus saves time
     * by simplifying multiplication operations. */
    i = x1;
    j = x2;
    bcuint(0.0, 1.0, 0.0, 1.0, x1-i, x2-j, y, y1, y2, *(cCF[i][j]));
}

void Hch_bicubicint (double x1, double x2,
		     double *y, double *y1, double *y2)
{
    int i, j;
    double x1max = ms-1.e-10, x2max = ms-1.e-10;
    x1 = (x1 > x1max ? x1max : x1);
    x2 = (x2 > x2max ? x2max : x2);
    
    /* Perform the interpolation on a unit square:  This makes
     * "t" and "u" 1.0 in the bcuint routine and thus saves time
     * by simplifying multiplication operations. */
    i = x1;
    j = x2;
    bcuint(0.0, 1.0, 0.0, 1.0, x1-i, x2-j, y, y1, y2, *(cCH[i][j]));
}

