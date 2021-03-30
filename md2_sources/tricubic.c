/*
 * md series 2 -- Molecular Dynamics of Plasma-Surface Chemistry
 * (c) 1999 Cam Abrams, Dept. of Chemical Engineering
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
#include "tricubic.h"

static const int mp = (X1_NGRIDPOINTS);
static const int np = (X2_NGRIDPOINTS);
static const int ms = (X1_NGRIDSQUARES);
static const int ns = (X2_NGRIDSQUARES);
static const int lp = (X3_NGRIDPOINTS);
static const int ls = (X3_NGRIDSQUARES);

int TRICUBIC_DIAG_=0;
short fcc_opt_=0;

static  double (****cCC)[4][4][4];
/* A ptr to...__||||    |_______|
 *	         |||         |_________________________ 
 *               |||                                   |
 *		 \|/				       |
 *		  |____an "m"x"n"x"l" 3D array of...   |___ 4x4x4 3D arrays.
 *
 *  Here, we state that the 3D grid of points ("knots") which define the
 *  function we are interpolating is "m+1"x"n+1"x"l+1" and has, therefore, 
 *  m x n x l 8-membered grid cubes.  Each one of these grid cubes 
 *  has associated with it a 4x4x4 matrix of coefficients, c_(ijk), 
 *  which are functions of the function values and its derivatives
 *  at the eight gridpoints in that cube.  This coefficient matrix
 *  c_(ijk) is used in the tricubic interpolation routine.  We wish to
 *  initially compute *all* 4x4x4 matrices (c_(ijk))_(mnl) initially, and
 *  store them for use by the interpolation routine.  In order to store
 *  the m x n x l 4x4x4 arrays, we declare "c" as a "(pointer to a) 
 *  3D-matrix of 3D-matrices",  or (*)(***c)[4][4][4], and dynamically
 *  allocate the "m" superrows of "n" columns of "l" rows each using the
 *  xmalloc function adapted from the bicubic analog in Numerical Recipes.
 *
 */

void tcucof(double y[8], double y1[8], double y2[8], double y3[8],
	    double y12[8], double y23[8], double y31[8], double y123[8],
            double d1, double d2, double d3, double c[4][4][4])
{
    /* Declare the 4,096 tabulated weighting factors in a [64]x[64]
     * array. */
    #include "tricof.h"
    
    int l, k, j, i;
    double xx, cl[64], x[64], d1d2, d2d3, d3d1, d1d2d3;
    
    d1d2 = d1*d2;
    d2d3 = d2*d3;
    d3d1 = d3*d1;
    d1d2d3 = d1d2*d3;
    
    #if DIAG
    if (TRICUBIC_DIAG_)
    {
	printf("tcucof begins\n");
	fflush(stdout);
    }
    #endif
    
    /* Build the temporary vector */
    for (i = 0; i < 8; i++)
    {
	#if DIAG
	if (TRICUBIC_DIAG_)
	{
	    printf("%i y=%.5le ", i, y[i]);fflush(stdout);
	    printf("y1=%.5le ", y1[i]);fflush(stdout);
	    printf("y2=%.5le ", y2[i]);fflush(stdout);
	    printf("y3=%.5le ", y3[i]);fflush(stdout);
	    printf("y12=%.5le ", y12[i]);fflush(stdout);
	    printf("y23=%.5le ", y23[i]);fflush(stdout);
	    printf("y31=%.5le ", y31[i]);fflush(stdout);
	    printf("y123=%.5le\n", y123[i]);fflush(stdout);
	}
	#endif
	x[i]	=   y[i];
	x[i+8]	=   y1[i]*d1;
	x[i+16]	=   y2[i]*d2;
	x[i+24]	=   y3[i]*d3;
	x[i+32]	=   y12[i]*d1d2;
	x[i+40]	=   y23[i]*d2d3;
	x[i+48]	=   y31[i]*d3d1;
	x[i+56]	=   y123[i]*d1d2d3;
    }
    
    /* Matrix multiply the weighting factors with the temporary vector,
     * and store the result in a second vector. */
    for (i = 0; i < 64; i++)
    {
	xx = 0.0;
	for (k = 0; k < 64; k++) xx += twt[i][k]*x[k];
	cl[i] = xx;
    }
    
    /* Unpack the second temporary vector into the 4x4x4 coefficient matrix. */
    l = 0;
    for (i = 0; i < 4; i++)
	for (j = 0; j < 4; j++)
	    for (k = 0; k < 4; k++)
		c[i][j][k] = cl[l++];
    #if DIAG
    if (TRICUBIC_DIAG_)
    {
	printf("tcucof ends\n");
	fflush(stdout);
    }
    #endif
}

void tcuint(double x1l, double x1u, double x2l, double x2u, 
     double x3l, double x3u, 
     double x1, double x2, double x3, 
     double *ansy, double *ansy1, double *ansy2, double *ansy3, 
     double c[4][4][4])
{
    int i, j, k;
    double t, u, v, d1, d2, d3;
    
    d1 = x1u - x1l;
    d2 = x2u - x2l;
    d3 = x3u - x3l;
    if (!d1 || !d2 || !d3)
    {
	fprintf(stderr, "3D Cubic Interpolation: Bad gridpoint coords:\n");
	fprintf(stderr, "\tx1u(%.5e) x1l(%.5e) x2u(%.5e) x2l(%.5e) x3u(%.5e) x3l(%.5e)\n", 
	    x1u, x1l, x2u, x2l, x3u, x3l);
	fprintf(stderr, "Program exits\n");
	exit(0);
    }
    t = (x1 - x1l)/d1;
    u = (x2 - x2l)/d2;
    v = (x3 - x3l)/d3;
    *ansy = (*ansy1) = (*ansy2) = (*ansy3) = 0.0;
    
    for (i = 3; i >= 0; i--)
    {
        *ansy = t*(*ansy)
	    + (((((c[i][3][3]*v + c[i][3][2])*v + c[i][3][1])*v + c[i][3][0])*u
	    +   (((c[i][2][3]*v + c[i][2][2])*v + c[i][2][1])*v + c[i][2][0]))*u
	    +   (((c[i][1][3]*v + c[i][1][2])*v + c[i][1][1])*v + c[i][1][0]))*u
	    +   (((c[i][0][3]*v + c[i][0][2])*v + c[i][0][1])*v + c[i][0][0]);
	*ansy1 = u*(*ansy1)
	    + ((((c[3][i][3]*v + c[3][i][2])*v + c[3][i][1])*v + c[3][i][0])*3*t
	    +  (((c[2][i][3]*v + c[2][i][2])*v + c[2][i][1])*v + c[2][i][0])*2)*t
	    +  (((c[1][i][3]*v + c[1][i][2])*v + c[1][i][1])*v + c[1][i][0]);
	*ansy2 = t*(*ansy2)
	    + (3*(((c[i][3][3]*v + c[i][3][2])*v + c[i][3][1])*v + c[i][3][0])*u
	    +  2*(((c[i][2][3]*v + c[i][2][2])*v + c[i][2][1])*v + c[i][2][0]))*u
	    +    (((c[i][1][3]*v + c[i][1][2])*v + c[i][1][1])*v + c[i][1][0]);
	*ansy3 = t*(*ansy3)
	    + ((((c[i][3][3]*u + c[i][2][3])*u + c[i][1][3])*u + c[i][0][3])*3*v
	    +  (((c[i][3][2]*u + c[i][2][2])*u + c[i][1][2])*u + c[i][0][2])*2)*v
	    +  (((c[i][3][1]*u + c[i][2][1])*u + c[i][1][1])*u + c[i][0][1]);
    }
    *ansy1 /= d1;
    *ansy2 /= d2;
    *ansy3 /= d3;    
}

void tricubic_genCoef (double y[X1_NGRIDPOINTS][X2_NGRIDPOINTS][X3_NGRIDPOINTS],
		       double y1[X1_NGRIDPOINTS][X2_NGRIDPOINTS][X3_NGRIDPOINTS], 
		       double y2[X1_NGRIDPOINTS][X2_NGRIDPOINTS][X3_NGRIDPOINTS],
		       double y3[X1_NGRIDPOINTS][X2_NGRIDPOINTS][X3_NGRIDPOINTS],
		       double y12[X1_NGRIDPOINTS][X2_NGRIDPOINTS][X3_NGRIDPOINTS], 
		       double y23[X1_NGRIDPOINTS][X2_NGRIDPOINTS][X3_NGRIDPOINTS], 
		       double y31[X1_NGRIDPOINTS][X2_NGRIDPOINTS][X3_NGRIDPOINTS], 
		       double y123[X1_NGRIDPOINTS][X2_NGRIDPOINTS][X3_NGRIDPOINTS], 
		       double (*****c)[4][4][4])
{
    int i, j, k, s, ii, jj, kk;
    double z[8], z1[8], z2[8], z3[8], z12[8], z23[8], z31[8], z123[8];
    int ip[][3] = {{0,0,0},{1,0,0},{1,1,0},{0,1,0},{0,0,1},{1,0,1},
                   {1,1,1},{0,1,1}};

    /* The array ip[][3] determines the order of visitation of the 8 vertices
     * of the cube:
     * 
     * 
     *              (8)------(7)
     *              /|       /|
     *             / |      / |		 Directions:
     *		(4)------(3)  |
     *		 |  (5)---|--(6)               /\    __ 
     *           |  /     |  /                 ||    //|
     *           | /      | /                  ||   //
     *          (1)------(2)             1===>  2  3 
     * 
     */

    #if DIAG
    if (TRICUBIC_DIAG_)
    {
	printf("tricubic_genCoef begins\n");
	fflush(stdout);
    }
    #endif

    /* Initialize *c as a matrix of matrices */
    (*c) = (double (****)[4][4][4])xmalloc(ls*sizeof(double (***)[4][4][4]));
    for (i = 0; i < lp; i++)
    {
	(*c)[i] = (double (***)[4][4][4])xmalloc(ms*sizeof(double (**)[4][4][4]));
	for (j = 0; j < mp; j++)
	{
	    (*c)[i][j] = (double (**)[4][4][4])xmalloc(ns*sizeof(double (*)[4][4][4]));
	    for (k = 0; k < np; k++)
		(*c)[i][j][k] = (double (*)[4][4][4])xmalloc(64*sizeof(double));
	}
    }

    for (i = 0; i < ls; i++)
    {
	for (j = 0; j < ms; j++)
	{
	    for (k = 0; k < ns; k++)
	    {
		#if DIAG
		if (TRICUBIC_DIAG_>1)
		{
		    printf("%i %i %i\n", i, j, k);
		    fflush(stdout);
		}
		#endif
		for (s = 0; s < 8; s++)
		{
		    ii = i + ip[s][0];
		    jj = j + ip[s][1];
		    kk = k + ip[s][2];
		    if (ii>=X1_NGRIDPOINTS) {printf("error big_i\n");exit(0);}
		    if (jj>=X2_NGRIDPOINTS) {printf("error big_j\n");exit(0);}
		    if (kk>=X3_NGRIDPOINTS) {printf("error big_k\n");exit(0);}
		    if (ii<0) {printf("error little_i\n");exit(0);}
		    if (jj<0) {printf("error little_j\n");exit(0);}
		    if (kk<0) {printf("error little_k\n");exit(0);}
		    z   [s] = y   [ii][jj][kk];
		    z1  [s] = y1  [ii][jj][kk];
		    z2  [s] = y2  [ii][jj][kk];
		    z3  [s] = y3  [ii][jj][kk];
		    z12 [s] = y12 [ii][jj][kk];
		    z23 [s] = y23 [ii][jj][kk];
		    z31 [s] = y31 [ii][jj][kk];
		    z123[s] = y123[ii][jj][kk];
		}
	        /* Use tcucof() to compute the 4x4x4 coefficient matrix c_(ijk) */
		tcucof(z,z1,z2,z3,z12,z23,z31,z123,1.0,1.0,1.0,*((*c)[i][j][k]));
	    }
	}
    }
    #if DIAG
    if (TRICUBIC_DIAG_)
    {
	printf("tricubic_genCoef ends\n");
	fflush(stdout);
    }
    #endif
}

void Fcc_genCoef (void)
{
    int i, j, k;
    double fcc[X1_NGRIDPOINTS][X2_NGRIDPOINTS][X3_NGRIDPOINTS];
    double fcc1[X1_NGRIDPOINTS][X2_NGRIDPOINTS][X3_NGRIDPOINTS];
    double fcc2[X1_NGRIDPOINTS][X2_NGRIDPOINTS][X3_NGRIDPOINTS];
    double fcc3[X1_NGRIDPOINTS][X2_NGRIDPOINTS][X3_NGRIDPOINTS];
    double fcc12[X1_NGRIDPOINTS][X2_NGRIDPOINTS][X3_NGRIDPOINTS];
    double fcc23[X1_NGRIDPOINTS][X2_NGRIDPOINTS][X3_NGRIDPOINTS];
    double fcc31[X1_NGRIDPOINTS][X2_NGRIDPOINTS][X3_NGRIDPOINTS];
    double fcc123[X1_NGRIDPOINTS][X2_NGRIDPOINTS][X3_NGRIDPOINTS];

    #if DIAG
    if (TRICUBIC_DIAG_)
    {
	printf("Fcc_genCoef begins\n");
	fflush(stdout);
    }
    #endif

    for (i = 0; i < lp; i++)
	for (j = 0; j < mp; j++)
	    for (k = 0; k < np; k++)
		fcc[i][j][k] = 
		fcc1[i][j][k] = fcc2[i][j][k] = fcc3[i][j][k] =
		fcc12[i][j][k] = fcc23[i][j][k] = fcc31[i][j][k] =
		fcc123[i][j][k] = 0.0;

    /* Set values of the Fcc function at the integer locations:
     * These values are from Brenner90, and represent values
     * for the hydrocarbon potential.  
     * fcc[a][b][c]:  a = "coordination" of c_i excluding j, 
     *		      b = "coordination" of c_j excluding j, 
     *		      c = Nij_conj.
     */
    
    /* We allow the user to choose which paramter set at runtime */
    
    if (fcc_opt_==1)
    {
	/* These are the values J. Tanaka emailed to me
	 * on 14 Apr 1999. */
	fcc[1][1][1] = -0.02880;   /* FCCF */
	fcc[2][2][1] =  0.04150;    /* C2F4 */
	fcc[1][2][1] = -0.0900;       /* FC-CF2 */
	fcc[1][3][1] = fcc[1][3][1] =  0.0;
	fcc[0][3][1] = fcc[0][3][1] =  0.0;
	fcc[0][2][2] = 0.0;
	fcc[0][2][1] = 0.0;
	fcc[0][1][1] = -0.02882;   /* C2F */
	fcc[1][1][2] = 0.0;
	fcc[2][3][1] = fcc[2][3][2] = -0.0363;  /* Evac for diamond */
	fcc[1][2][2] = -0.0243;		    /* Evac for graphite */
    }
    else if (fcc_opt_==0)
    {
	/* These are the values as they appear in a draft
	 * form of Tanaka99b.  Quite different. */
	fcc[0][1][1] = -0.0412;		/* C2F */
	fcc[0][1][2] = 0.0484;		/* C3 */
	fcc[0][2][1] = -0.0683;		/* CC dbl b. in C=CF2 */
	fcc[1][1][1] = -0.0635;		/* CC tpl b. in FCCF */
	fcc[1][1][2] = -0.0127;		/* C4 */
	fcc[1][2][1] = -0.1690;		/* CC b.e. in C2F3 */
	fcc[2][2][1] = -0.1248;		/* CC b.e. in F2CCF2 */
	fcc[2][3][1] = fcc[2][3][2] = -0.0363;  /* Evac for diamond */
	fcc[1][2][2] = -0.0243;		    /* Evac for graphite */
    }
    /* Enforce the symmetry of Fcc */
    for (i=0;i<lp;i++)
    {
	for (j=0;j<mp;j++)
	{
	    for(k=0;k<np;k++)
	    {
		if (k>2) 
		{
		    fcc[i][j][k] =fcc[i][j][2];
		    fcc1[i][j][k]=fcc1[i][j][2];
		    fcc2[i][j][k]=fcc2[i][j][2];
		    fcc3[i][j][k]=0.0;
		}
		if (i>j) fcc[i][j][k]=fcc[j][i][k];
	    }
	}
    }

    /* Differentiate, 2nd order except 8 corners */

    for (k=0;k<np;k++)
    {
	fcc12[0][0][k] =(fcc[1][1][k]-fcc[1][0][k]-fcc[0][1][k]+fcc[0][0][k]);
	fcc12[ls][0][k]=(fcc[ls][1][k]-fcc[ls][0][k]-fcc[ls-1][1][k]+fcc[ls-1][0][k]);
	fcc12[0][ms][k]=(fcc[1][ms][k]-fcc[1][ms-1][k]-fcc[0][ms][k]+fcc[0][ms-1][k]);
	fcc12[ls][ms][k]=(fcc[ls][ms][k]-fcc[ls][ms-1][k]-fcc[ls-1][ms][k]+fcc[ls-1][ms-1][k]);
	for(i=1;i<ls;i++)
	{
	    fcc12[i][0][k] =
		(-3*(fcc[i+1][0][k]-fcc[i-1][0][k])+
		4*(fcc[i+1][1][k]-fcc[i-1][1][k])-(fcc[i+1][2][k]-fcc[i-1][2][k]))/4;
	    fcc12[i][ms][k] =
		-(-3*(fcc[i+1][ms][k]-fcc[i-1][ms][k])+
		4*(fcc[i+1][ms-1][k]-fcc[i-1][ms-1][k])-(fcc[i+1][ms-2][k]-fcc[i-1][ms-2][k]))/4;
	}
	for (j=1;j<ms;j++)
	{
	    fcc12[0][j][k] =
		(-3*(fcc[0][j+1][k]-fcc[0][j-1][k])+
		4*(fcc[1][j+1][k]-fcc[1][j-1][k])-(fcc[2][j+1][k]-fcc[2][j-1][k]))/4;
	    fcc12[ls][j][k] =
		-(-3*(fcc[ls][j+1][k]-fcc[ls][j-1][k])+
		4*(fcc[ls-1][j+1][k]-fcc[ls-1][j-1][k])-(fcc[ls-2][j+1][k]-fcc[ls-2][j-1][k]))/4;
	}
    }
    for (i=0;i<lp;i++)
    {
	fcc23[i][0][0] =(fcc[i][1][1]-fcc[i][1][0]-fcc[i][0][1]+fcc[i][0][0]);
	fcc23[i][ms][0]=(fcc[i][ms][1]-fcc[i][ms][0]-fcc[i][ms-1][1]+fcc[i][ms-1][0]);
	fcc23[i][0][ns]=(fcc[i][1][2]-fcc[i][1][1]-fcc[i][0][2]+fcc[i][0][1]);
	fcc23[i][ms][ns]=(fcc[i][ms][2]-fcc[i][ms][1]-fcc[i][ms-1][2]+fcc[i][ms-1][1]);
	for (k=1;k<ns;k++)
	{
	    fcc23[i][0][k] =
		(-3*(fcc[i][0][k+1]-fcc[i][0][k-1])+
		4*(fcc[i][1][k+1]-fcc[i][1][k-1])-(fcc[i][2][k+1]-fcc[i][2][k-1]))/2;
	    fcc23[i][ms][k] =
		-(-3*(fcc[i][ms][k+1]-fcc[i][ms][k-1])+
		4*(fcc[i][ms-1][k+1]-fcc[i][ms-1][k-1])-(fcc[i][ms-2][k+1]-fcc[i][ms-2][k-1]))/2;
	}
	for (j=1;j<ms;j++)
	{
	    fcc23[i][j][0] =
		(-3*(fcc[i][j+1][0]-fcc[i][j-1][0])+
		4*(fcc[i][j+1][1]-fcc[i][j-1][1])-(fcc[i][j+1][2]-fcc[i][j-1][2]));
	    fcc23[i][j][ns] =
		-(-3*(fcc[i][j+1][ns]-fcc[i][j-1][ns])+
		4*(fcc[i][j+1][ns-1]-fcc[i][j-1][ns-1])-(fcc[i][j+1][ns-2]-fcc[i][j-1][ns-2]));
	}
    }

    for (j=0;j<mp;j++)
    {
	fcc31[0][j][0] = (fcc[1][j][1]-fcc[1][j][0]-fcc[0][j][1]+fcc[0][j][0]);
	fcc31[ls][j][0] = (fcc[ls][j][1]-fcc[ls][j][0]-fcc[ls-1][j][1]+fcc[ls-1][j][0]);
	fcc31[0][j][ns] = (fcc[1][j][ns]-fcc[1][j][ns-1]-fcc[0][j][ns]+fcc[0][j][ns-1]);
	fcc31[ls][j][ns] = (fcc[ls][j][ns]-fcc[ls][j][ns-1]-fcc[ls-1][j][ns]+fcc[ls-1][j][ns-1]);
	for (k=1;k<ns;k++)
	{
	    fcc31[0][j][k] =
		(-3*(fcc[0][j][k+1]-fcc[0][j][k-1])+
		4*(fcc[1][j][k+1]-fcc[1][j][k-1])-(fcc[2][j][k+1]-fcc[2][j][k-1]))/2;
	    fcc31[ls][j][k] =
		-(-3*(fcc[ls][j][k+1]-fcc[ls][j][k-1])+
		4*(fcc[ls-1][j][k+1]-fcc[ls-1][j][k-1])-(fcc[ls-2][j][k+1]-fcc[ls-2][j][k-1]))/2;
	}
	for (i=1;i<ls;i++)
	{
	    fcc31[i][j][0] = 
		(-3*(fcc[i+1][j][0]-fcc[i-1][j][0])+
		4*(fcc[i+1][j][1]-fcc[i-1][j][1])-(fcc[i+1][j][2]-fcc[i-1][j][2]))/4;
	    fcc31[i][j][ns] =
		-(-3*(fcc[i+1][j][ns]-fcc[i-1][j][ns])+
		4*(fcc[i+1][j][ns-1]-fcc[i-1][j][ns-1])-(fcc[i+1][j][ns-2]-fcc[i-1][j][ns-2]))/4;
	}
    }

    for (k=0;k<np;k++)
    {
	for (j=0;j<mp;j++)
	{
	    fcc1[0][j][k] =(-3*fcc[0][j][k]+4*fcc[1][j][k]-fcc[2][j][k])/2;
	    fcc1[ls][j][k]=( 3*fcc[ls][j][k]-4*fcc[ls-1][j][k]+fcc[ls-2][j][k])/2; 
	    for (i=0;i<ls;i++)
	    {
		if (j==0) fcc2[i][j][k]=fcc[i][j+1][k]-fcc[i][j][k];
		if (j==ms) fcc2[i][j][k]=-(-3*fcc[i][j][k]+4*fcc[i][j-1][k]-fcc[i][j-2][k])/2;
		if (k==0) fcc3[i][j][k]=(-3*fcc[i][j][k]+4*fcc[i][j][k+1]-fcc[i][j][k+2])/2;
		if (k==ns) fcc3[i][j][k]=-(-3*fcc[i][j][k]+4*fcc[i][j][k-1]-fcc[i][j][k-2])/2; 
		if (i>0&&i<ls) fcc1[i][j][k]=(fcc[i+1][j][k]-fcc[i-1][j][k])/2;
		if (j>0&&j<ms) fcc2[i][j][k]=(fcc[i][j+1][k]-fcc[i][j-1][k])/2;
		if (k>0&&j<ns) fcc3[i][j][k]=(fcc[i][j][k+1]-fcc[i][j][k-1])/2;
		if (i>0&&i<ls&&j>0&&j<ms) fcc12[i][j][k] =
		    (fcc[i+1][j+1][k]-fcc[i-1][j+1][k]-fcc[i+1][j-1][k]+fcc[i-1][j-1][k]);
		if (j>0&&j<ms&&k>0&&k<ns) fcc23[i][j][k] =
		    (fcc[i][j+1][k+1]-fcc[i][j-1][k+1]-fcc[i][j+1][k-1]+fcc[i][j-1][k-1]);
		if (k>0&&k<ns&&i>0&&i<ls) fcc31[i][j][k] =
		    (fcc[i+1][j][k+1]-fcc[i-1][j][k+1]-fcc[i+1][j][k-1]+fcc[i-1][j][k-1]); 
	    }
	}
    }

    tricubic_genCoef(fcc, fcc1, fcc2, fcc3, fcc12, fcc23, fcc31, fcc123, &cCC);
    
    #if DIAG
    if(TRICUBIC_DIAG_)
    {
	printf("#\tHere are the 'uninterpolated' inputs:\n");
	printf("#\tfcc[1][1][1](0,1,2,3) = %.5le, %.5le, %.5le, %.5le\n", 
	    fcc[1][1][1], fcc1[1][1][1], fcc2[1][1][1],  fcc3[1][1][1]);
	printf("#\tfcc[2][2][1](0,1,2,3) = %.5le, %.5le, %.5le, %.5le\n", 
	    fcc[2][2][1], fcc1[2][2][1], fcc2[2][2][1],  fcc3[2][2][1]);
	printf("#\tfcc[1][2][1](0,1,2,3) = %.5le, %.5le, %.5le, %.5le\n", 
	    fcc[1][2][1], fcc1[1][2][1], fcc2[1][2][1],  fcc3[1][2][1]);
	printf("#\tfcc[1][3][1](0,1,2,3) = %.5le, %.5le, %.5le, %.5le\n", 
	    fcc[1][3][1], fcc1[1][3][1], fcc2[1][3][1],  fcc3[1][3][1]);
	printf("#\tfcc[1][3][2](0,1,2,3) = %.5le, %.5le, %.5le, %.5le\n", 
	    fcc[1][3][2], fcc1[1][3][2], fcc2[1][3][2],  fcc3[1][3][2]);
	printf("#\tjust to name a few.\n");
    }
    #endif
    /* Interpolate all integer-indexed values */
/*    for (i = 0; i < lp; i++)
    {
	for (j = 0; j < mp; j++)
	{
	    for(k = 0; k < np; k++)
	    {
		tcuint(0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 
		       0.0, 0.0, 0.0, 
		       &(fcc[i][j][k]), &(fcc1[i][j][k]), 
		       &(fcc2[i][j][k]), &(fcc3[i][j][k]), *(cCC[i][j][k]));
	    }
	}
    }
  */  #if DIAG
    if (TRICUBIC_DIAG_)
    {
	printf("#\tHere are the 'interpolated' inputs:\n");
	printf("#\tfcc[1][1][1](0,1,2,3) = %.5le, %.5le, %.5le, %.5le\n", 
	    fcc[1][1][1], fcc1[1][1][1], fcc2[1][1][1],  fcc3[1][1][1]);
	printf("#\tfcc[2][2][1](0,1,2,3) = %.5le, %.5le, %.5le, %.5le\n", 
	    fcc[2][2][1], fcc1[2][2][1], fcc2[2][2][1],  fcc3[2][2][1]);
	printf("#\tfcc[1][2][1](0,1,2,3) = %.5le, %.5le, %.5le, %.5le\n", 
	    fcc[1][2][1], fcc1[1][2][1], fcc2[1][2][1],  fcc3[1][2][1]);
	printf("#\tfcc[1][3][1](0,1,2,3) = %.5le, %.5le, %.5le, %.5le\n", 
	    fcc[1][3][1], fcc1[1][3][1], fcc2[1][3][1],  fcc3[1][3][1]);
	printf("#\tfcc[1][3][2](0,1,2,3) = %.5le, %.5le, %.5le, %.5le\n", 
	    fcc[1][3][2], fcc1[1][3][2], fcc2[1][3][2],  fcc3[1][3][2]);
	printf("#\tjust to name a few.\n");
    }
    #endif
    #if DIAG
    if (TRICUBIC_DIAG_)
    {
	printf("Fcc_genCoef ends\n");
	fflush(stdout);
    }
    #endif
}

void Fcc_tricubicint (double x1, double x2, double x3,
		      double *y, double *y1, double *y2, double *y3)
{
    int i, j, k;
    double x1max = ls-1e-10, x2max = ms-1e-10, x3max = ns-1e-10;
    
    
    x1 = (x1 > x1max ? x1max : x1);
    x2 = (x2 > x2max ? x2max : x2);
    x3 = (x3 > x3max ? x3max : x3);
    i = x1;
    j = x2;
    k = x3;
    
    /* if this is a spline knot, don't do the interpolation */
/*    if (i==x1 && j==x2 && k==x3)
    {
	*y = fcc[i][j][k];
	*y1 = fcc1[i][j][k];
	*y2 = fcc2[i][j][k];
	*y3 = fcc3[i][j][k];
    }
    else
*/    /* Perform the interpolation on a unit cube:  This makes
     * "t", "u", and "v" 1.0 in the tcuint routine and thus saves time
     * by simplifying multiplication operations. */
	tcuint(0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 
	   x1-i, x2-j, x3-k, 
	   y, y1, y2, y3, *(cCC[i][j][k]));
}

