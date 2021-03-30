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
#define PAIR_POTENTIAL 0    /* Tersoff/Brenner/Tanaka is NOT a pair potential */
#include <math.h>
#include "genforce.h"	/* generic force routine header */
#include "tbt_sicf.h"	/* declares specific functions available with this PEF */
#include "sicf_params.h"/* PEF parameters set as macros */
#include "bicubic.h"	/* header for bicubic interpolation routines */
#include "tricubic.h"	/* header for tricubic interpolation routines */
#include "clust.h"

double PE_=0.0;			/* total potential energy */
double xPE_=0.0;		/* old PE calculation from 2+3 body PE */
static pt TAU__={0, 0, 0};	/* diagonal elements of the stress tensor */
ptPtr Tau_=&TAU__;		/* homogeneous stress tensor components */

short pe_set_=0;		/* lock */
int bond_ctrl_=2;		/* bond definition control parameter */

short paramset_a_=0;		/* default it to use parameter set B (newer) */
extern short hcc_opt_;
extern short hcf_opt_;
extern short fcc_opt_;


extern	atomPtr L_;		/* the configuration */
extern	element	per_table[];	/* the periodic table */

int	SIC_DIAG2_;
int	SIC_DIAG3_;
#if DIAG
extern	int	DIAG_;
static  int	D2_local, D3_local;
extern	int	D3_local;
extern	int	IMG_TRACE_;
extern	int	N_DIAG_;
#endif

#define MA	MAXNUMATOMS
#define MN	MAXNUMNEIGHBORS
/* Global storage arrays: */
static	pt	rij_hat[MA][MN];    /* pair separation unit vectors */
static	double	rij_[MA][MN];	    /* pair separation scalars */
static	double	rij2_[MA][MN];	    /* pair square separation scalars */
static	double	fij_[MA][MN];	    /* fij (smooth cutoffs) */
static	double  fij1_[MA][MN];	    /* First derivative of fij */
static	double	Va_ij_[MA][MN];	    /* Attractive Morse */
static	double	Vr_ij_[MA][MN];	    /* Repuslive Morse */
static	double	dVa_drij_[MA][MN];  /* First derivative of Va */
static	double	dVr_drij_[MA][MN];  /* First derivative of Vr */
static	double	Phi_ij_[MA][MN];    /* Potential energy of pair (i,j) */
static	double	dtPhi_ij_[MA][MN];  /* Time rate of change of potential 
				     * energy of pair (i,j) */
	double  phi_rate_tol_;
static	double	Fik_[MA][MN];	    /* Fik (conj. smooth cutoffs) */
static	double  Fik1_[MA][MN];	    /* F'ik (conj. smooth cutoff gradients) */
static	pt	frc_ij_[MA][MN];    /* force on i due to j */
static	double  NiF_[MA];	    /* NiF (coord. of i by F atoms) */
static	double  NiC_[MA];	    /* NiC (coord. of i by C atoms) */
static	double  NiSi_[MA];	    /* NiSi (coord. of i by Si atoms) */

/* Arrays reused in each atom pair PE/Force computation: */
static	double	sijk_[MN];	    /* cos(theta_ijk) */
static	double	sjik_[MN];	    /* cos(theta_jik) */
static	double	gijk_[MN];	    /* gijk */
static	double	gjik_[MN];	    /* gjik */
static	double	gijk1_[MN];	    /* g'ijk */
static	double	gjik1_[MN];	    /* g'jik */
static	double	omijk_[MN];	    /* omega_ijk */
static	double	omjik_[MN];	    /* omega_jik */
static	double	epijk_[MN];	    /* epsilon_ijk */
static	double	epjik_[MN];	    /* epsilon_jik */

static	pt	RijT_ = {0.0, 0.0, 0.0};
static  ptPtr   pRijT_ = &RijT_;

/* Global 2body parameter variables: */
void TwoBody_SetParameters (sym_type is, sym_type js);
static	double Aij, lamij;		/* Repulsive Morse parameters */
static	double Bij, muij;		/* Attractive Morse parameters */
static	double Ms_ra, Ms_rb;		/* Moliere-spline r_a and r_b */
static	double Ms_a, Ms_b, Ms_c, Ms_s;	/* Moliere-spline a, b, c, d */
static	double M_a;			/* Moliere a */
static	double M_C;			/* Moliere C */
static	double rijCu1_;			/* R(1) in fij */
static	double rijCu2_;			/* R(2) in fij */
static	double rijCu2sq_;		/* square of R(2) in fij */
static	double PI_rijCu12_;		/* pi * [ R(2) - R(1) ] in fij */
static	short  moliere_only;		/* set to 1 if the pair interaction is
					 * described by the Moliere repuslive
					 * potential ONLY; set to 0 otherwise. */
/* Global 3body parameter variables for the b_ij functions: */
void ThreeBody_SetParameters (sym_type i, sym_type j, sym_type k);
static double g_a, g_c, g_d, g_h, g_c2, g_d2, g_c2d2, g_ac2, g_ac2d2, g_p1_, g_p2_;
static double rij_e, rik_e, rjk_e;
static double eta, delta, alpha, beta;

/* Fij_:  computes smooth cutoff function fij and its first
 * derivative f'ij at separation rij.  fij is retured in *fij, 
 * and f'ij in *f1ij. */
void Fij_ (double rij, double * fij, double * f1ij);

/* g_costh: computes bond angle function "g" and it's first
 * derivative "g'" at cos(theta) = "s" and returns them in
 * *g0 and *g1 respectively. */
void g_costh (double s, sym_type i, double * g0, double * g1);

/* omega_: computes bond length competition function "omega" 
 * and the factor "epsilon" in it's first derivative "omega'"
 * at bond lengths rij and rik.  (omega' = epsilon*omega)
 * omega is returned in *om_ and epsilon in *ep_. */
void omega_ (double rij, double rik, double * om_, double * ep_);

/* Fconj:  computes conjugate strength function F (the one that appears
 * in the Nij_conj summation) and its first derivative F' at N_(ik) = n.
 * Returns F in *F0 and F' in *F1. */
void Fconj (double n, double * F0, double * F1);

/* Moliere:  Computes Moliere-type interatomic energy Vr and its
 * first derivative Vr' at separation rij.  Vr is returned in *Vr0
 * and Vr' is returned in *Vr1. */
void Moliere (double rij, double * Vr0, double * Vr1, double * Vr2);

/* MolSpline:  Computes Spline-type interatomic energy Vr and its
 * first derivative Vr' at separation rij.  Vr is returned in *Vr0
 * and Vr' is returned in *Vr1. */
void MolSpline (double rij, double * Vr0, double * Vr1);

/* Morse functions compute Morse-type interatomic energy V and
 * its first derivative V' at separation rij.  fij is the smooth
 * cutoff function f(rij). V is returned in *V0, and V' in *V1. 
 * If Va2 (Vr2) is not NULL, the 2nd derivative of Va (Vr) is 
 * computed and returned in *Va2 (*Vr2) ONLY IF fij is 1.0. */
void MorseVa (double rij, double fij, double fij1, double * Va0,
	      double * Va1, double * Va2);
void MorseVr (double rij, double fij, double fij1, double * Vr0,
	      double * Vr1, double * Vr2);

/* Fcc_tester is for debugging purposes only.  It is a 3D gaussian 
 * which is well-behaved over all of 3space, as are its 3 first 
 * derivatives.  Fcc_tester is used to make sure that the forces are 
 * tracked correctly by the mainForce routine independent of any errors 
 * present in the tricubic interpolation of Fcc. */

double fcc_a_t = 0.05;
double fcc_b_t = 0.01;
void Fcc_tester (double x1, double x2, double x3, double * y, 
		 double * y1, double * y2, double * y3)
{
    *y = fcc_a_t*exp(-fcc_b_t*((x1-1)*(x1-1) + (x2-1)*(x2-1) + (x3-1)*(x3-1)));
    *y1 = (*y)*(-fcc_b_t*2*(x1-1));
    *y2 = (*y)*(-fcc_b_t*2*(x2-1));
    *y3 = (*y)*(-fcc_b_t*2*(x3-1));
}

double hcc_a_t = 0.05;
double hcc_b_t = 0.01;
void Hcc_tester (double x1, double x2, double * y, double * y1, double * y2)
{
    *y = hcc_a_t*exp(-hcc_b_t*((x1-1)*(x1-1) + (x2-1)*(x2-1)));
    *y1 = (*y)*(-hcc_b_t*2*(x1-1));
    *y2 = (*y)*(-hcc_b_t*2*(x2-1));
}

short pef_initialized_=0;
extern short quiet_;
void pef_initialize (void)
{

    if (!quiet_) printf("# md series 2 pef initializer (c) 1999 cfa\n");

    /* establish which parameter set is being used; default is B */
    if (paramset_a_)
    {
	/* setting these three flags to 1 tells the genCoef()'s to
	 * use parameter set A (the older one) rather than
	 * parameter set B (the newer one). */
	hcc_opt_=1;
	hcf_opt_=1;
	fcc_opt_=1;
	if (!quiet_) printf("# Using parameter set A\n");
    }

    /* Set up the bi- and tricubic interpolation functions */
    #if DIAG
    if (SIC_DIAG2_)
    {
	printf("Calling Hcc_genCoef()...\n");
	fflush(stdout);
    }
    #endif
    Hcc_genCoef();
    #if DIAG
    if (SIC_DIAG2_)
    {
	printf("Calling Hch_genCoef()...\n");
	fflush(stdout);
    }
    #endif
    Hch_genCoef();
    #if DIAG
    if (SIC_DIAG2_)
    {
	printf("Calling Hcf_genCoef()...\n");
	fflush(stdout);
    }
    #endif
    Hcf_genCoef();
    #if DIAG
    if (SIC_DIAG2_)
    {
	printf("Calling Fcc_genCoef()...\n");
	fflush(stdout);
    }
    #endif
    Fcc_genCoef();
    
    #if DIAG
    if (SIC_DIAG2_)
    {
	printf("Done setting up PEF\n");
	fflush(stdout);
    }
    #endif

    pef_initialized_=1;
}

void pef_paramout (FILE* ofp)
{
    fprintf(ofp, "Parameters of the TBT SiCF Potential:\n");
    fprintf(ofp, "(not yet implemented)\n");
}


/* cfg_setforce is the main driver function for the
 * force routine.  Its tasks are:
 * (1) clear all scratch arrays and initialize global energy variables;
 * (2) call the precomputations routine;
 * (3) call the main force routine.
 */
static double twoBodyEnergy = 0.0;
static double threeBodyEnergy = 0.0;
void cfg_setforce (void)
{
    char * d_h = "tbt_sicf::cfg_setforce";
    void Precomputations (void);   
    void mainForce (void);
    void cfg_frc2hac (void);
    int i, j;

    if (!pef_initialized_) {printf("error: pef not initialized.\n");exit(0);}
    if (pe_set_) 
    {printf("# Warning: set_force() already called for current cfg\n");return;}
    
    /* Clear all the scratch arrays. */
    for (i = 0; i < MA; i++)
    {
	NiF_[i] = NiC_[i] = NiSi_[i] = 0.0;
	for (j = 0; j < MN; j++)
	{
	    ptPtr_clear(&(rij_hat[i][j]));
	    ptPtr_clear(&(frc_ij_[i][j]));
	    rij_[i][j] = rij2_[i][j] = fij_[i][j] = fij1_[i][j] =
		Va_ij_[i][j] = Vr_ij_[i][j] = 
		dVr_drij_[i][j] = dVa_drij_[i][j] = 
		Fik_[i][j] = Fik1_[i][j] = dtPhi_ij_[i][j] = 0.0;
	}
    }

    /* initialize GLOBAL total potential energy */
    PE_ = 0.0;
    twoBodyEnergy = 0.0;
    threeBodyEnergy = 0.0;
    clust_clear();

    #if DIAG
    D2_local = SIC_DIAG2_;
    D3_local = SIC_DIAG3_;
    
    if (D2_local || D3_local)
	printf("%s Calling Precomputations routine...\n", d_h);
    #endif

    /* Call the precomputations function. */
    Precomputations();

    #if DIAG
    if (D2_local || D3_local)
	printf("%s Calling Main Force routine...\n", d_h);
    #endif

    /* Call the main force routine */
    mainForce();

    /* compute the total PE_ from the two- and three-body PE */
    xPE_ = twoBodyEnergy + threeBodyEnergy;
    
    /* printf("%.5le\t%.5le\t%.5le\n", twoBodyEnergy, threeBodyEnergy, PE_); */
    
    /* convert forces into half-accelerations */
    cfg_frc2hac();

    /* compute the force-term of the stress tensor */
    cfg_frc2tau();

    #if DIAG
    if (D2_local || D3_local)
    {
	printf("%s Two body energy = [%.5le] or [%.2lf%%] of Total.\n", 
	    d_h, twoBodyEnergy, twoBodyEnergy/PE_*100.0);
	printf("%s Three body energy = [%.5le] or [%.2lf%%] of Total.\n", 
	    d_h, threeBodyEnergy, threeBodyEnergy/PE_*100.0);
	printf("%s Total energy = [%.5le].\n", d_h, PE_);
	printf("%s ends.\n", d_h);
    }
    #endif
    
    clust_fixlist();

/*    {
	atomPtr ip;
	nNodePtr jp;
	for (ip=L_;ip;ip=ip->next)
	{
	    i=ip->id;
	    for (jp=ip->nList;jp;jp=jp->next)
	    {
		j=jp->addr->id;
		if (fabs(dtPhi_ij_[i][j])>phi_rate_tol_)
		{
		    printf("pair %i %i dphi=%.5le\n", i, j, dtPhi_ij_[i][j]);
		}
	    }
	}
    }
  */ 
    pe_set_=1;
    
}

/* Scratch variables for the main force routine.  Some of these
 * are "fetched" by the force routine from the precomputed arrays. */
static	double Vr, Va;		    /* Morse repulsive and attractive PEFs */
static  double dVa_drij;	    /* First derivative of Va */
static  double dVr_drij;	    /* First derivative of Vr */
static  double Va_2;		    /* Va/2.0 */
static  double Ni, Nj, Nk;	    /* smooth coordination of i, j, k */
static  double NiF, NiC, NiSi;	    /* F, C, and Si coordination of i */
static  double NjF, NjC, NjSi;	    /* F, C, and Si coordination of j */
static  double NkF, NkC, NkSi;	    /* F, C, and Si coordination of k */
static  double Nij_F, Nji_F;	    /* F!=j coor. of i, and F!=i coor. of j */
static  double Nij_C, Nji_C;	    /* C!=j coor. of i, and C!=i coor. of j */
static  double Nij_Si, Nji_Si;	    /* Si!=j coor. of i, and Si!=i coor. of j */
static  double Nij, Nji;	    /* coor. of i excl. j, coor. of j excl. i */
static  double Nik, Nki;	    /* coor. of i excl. k, coor. of k excl. i */
static  double Njk, Nkj;	    /* coor. of j excl. k, coor. of k excl. j */
static  double Fik, Fki;	    /* F(N_(ik)) and F(N_(ki)) respectively */
static  double Fik1, Fki1;	    /* F'(N_(ik)) and F'(N_(ki)) respectively */
static  double Fjk, Fkj;	    /* F(N_(jk)) and F(N_(kj)) respectively */
static  double Fjk1, Fkj1;	    /* F'(N_(jk)) and F'(N_(kj)) respectively */
static  double Fil, Fli;	    /* F(N_(il)) and F(N_(li)) respectively */
static  double Fil1, Fli1;	    /* F'(N_(il)) and F'(N_(li)) respectively */
static  double Fjl, Flj;	    /* F(N_(jl)) and F(N_(lj)) respectively */
static  double Fjl1, Flj1;	    /* F'(N_(jl)) and F'(N_(lj)) respectively */
static  double Nij_conj;	    /* Conjugated coor. of C-C bond (i,j) */
static  double rij, rik, rjk;	    /* scalar separation of (i,j), (i,k), & (j,k) */
static  double rij2, rik2, rjk2;    /* scalar squared separation of (i,j), (i,k), & (j,k) */
static  double fij, fij1;	    /* (i,j) smooth cutoff and 1st derivative */
static  double fik, fik1;	    /* (i,k) smooth cutoff and 1st derivative */
static  double fil, fil1;	    /* (i,l) smooth cutoff and 1st derivative */
static  double fjk, fjk1;	    /* (j,k) smooth cutoff and 1st derivative */
static  double fjl, fjl1;	    /* (j,l) smooth cutoff and 1st derivative */
static  double fkl, fkl1;	    /* (k,l) smooth cutoff and 1st derivative */
static  double sijk, sjik;	    /* cos(theta_ijk) and cos(theta_jik) */
static  double gijk, gijk1;	    /* bond angle function for (i,j,k) and 1st der. */
static  double gjik, gjik1;	    /* bond angle function for (j,i,k) and 1st der. */
static  double omijk, epijk;	    /* omega_ijk and epsilon_ijk (the exp. func.) */
static  double omjik, epjik;	    /* omega_jik and epsilon_jik (the exp. func.) */
static  double xi_ij, xi_ji;	    /* Bond competition function for i->j, i<-j 
				     * This is "zeta" in Brenner's definition, but
				     * I can't write a "zeta" on paper so I called
				     * it "xi". */
static  double Hij, Hij1, Hij2;	    /* BO correction function & 2 1st der. for i->j */
static  double Hji, Hji1, Hji2;	    /* BO correction function & 2 1st der. for i<-j */
static  double Fcc, Fcc1, Fcc2, Fcc3;
				    /* CC correction function & 3 1st der. */
				    
short  USE_FCC = 1;	    /* flag that can be controlled by the debugging
			     * routines. */
short  USE_HIJ = 1;	    /* flag that can be controlled by the debugging
			     * routines. */
short  HIJ_TESTER = 0;	    /* flag that can be controlled by the debugging
			     * routines. */
short  FCC_TESTER = 0;	    /* flag that can be controlled by the debugging
			     * routines. */
short	USE_HIE_CORE=1;	    /* flag for controlling whether to use the
			     * Beardmore/Smith-type high energy repulsive
			     * core functions -- default is on */


static  double b_ij, b_ji;	    /* bond order function for i->j and i<-j */
static  double b_ij_bar;	    /* Avg BO correction function for i<->j */
static  double beta_ij, beta_ji;    /* beta-coefficients for i->j and i<-j */

static  double fknl__[8];	    /* Three body force kernel magnitudes */
static  pt F2_ij, F2_ji;	    /* Two-body force j->i and i->j */
static  pt F3_ij, F3_ji;	    /* Three-body force j->i and i->j */
static  pt F3k_ij, F3k_ji;	    /* Force on (i,j) neighbor k */

static  const ptPtr pF2_ij = &F2_ij;/* pF2_ij always pts to F2_ij */
static  const ptPtr pF2_ji = &F2_ji;

static  ptPtr Rij_hatp, Rik_hatp, Rjk_hatp, Rkl_hatp;
			    /* Pointers to unit vector separations for
			     * (i, j), (i, k), (j, k), and (k, l). */
			     
static pt TMP_PT;		    /* temporary 3-vector */

void Precomputations (void)
{
    atomPtr ip = NULL, jp = NULL;
    sym_type is, js;
    char * d_h = "tbt_sicf::Precomputations";
    int i = 0, j = 0, fatal = 0;
    int nij = 0, nji = 0;
 
    #if DIAG
    D2_local = SIC_DIAG2_;
    if (D2_local) 
    {
	printf("\n");
	printf("%s diagnostics\n", d_h);
    }
    #endif
    
    
    /* Clear all force vectors. */
    for (ip = L_; ip; ip = ip->next) ptPtr_clear(ip->frc);
    
    /* Clear all scratch variables */
    Va = dVa_drij = Vr = dVr_drij = fij = fij1 = 0.0;
    
    /* Reset the neighbor lists */
    cfg_nreset();

    /* Examine each unique pair of atoms in the system. */
    
    for (ip=L_;ip->next;ip=ip->next)
    {
	i = ip->id;
	is = ip->sym;
	#if DIAG
	if (D2_local)
	{
	    printf("2b: atom %s_%i\n", per_table[is].sym, i);
	    fflush(stdout);
	}
	#endif
	for (jp = ip->next; jp; jp = jp->next)
	{
	    j = jp->id;
	    js = jp->sym;
	    
	    /* Assign all global two-body parameter values. */
	    TwoBody_SetParameters(is, js);

	    #if DIAG
	    if (D2_local > 1)
	    {
		printf("(%s(%i)?%s(%i))", 
			per_table[ip->sym].sym, i, 
			per_table[jp->sym].sym, j);
		fflush(stdout);
	    }
	    #endif

	    /* Compute the separation vector Rij and its squared 
	     * magnitude.  (We'll take the square root only if it's
	     * necessary.) */
	    
	    ptPtr_minimg(pRijT_, ptPtr_subtract(pRijT_, ip->pos, jp->pos));
	    rij2 = ptPtr_sqabs(pRijT_);
	    
	    #if DIAG
	    if (D2_local > 1)
	    {
		printf("\t%s(%i)???%s(%i), rij=%.5le\n", 
		    per_table[is].sym, i, per_table[js].sym, j, 
		    sqrt(rij2));
	    }
	    #endif
	    #ifdef RDF
	    /* Update the RDF histogram if requested */
	    if (RDF_
		&& ip->flags.is_RDF_a && jp->flags.is_RDF_a)
	    {
		rdf_thisBin = (int)(sqrt(rij2)/RDF_dR_) + 1;
		if (rdf_thisBin <= MAXBIN) RDF_Hist_[rdf_thisBin] += 2;
	    }
	    #endif
	    #if DIAG
	    if (D2_local > 2)
	    {
		printf("%s %s(%i)(i)->pos = %.5le, %.5le, %.5le\n", 
		    d_h, per_table[is].sym, i, ip->pos->x, ip->pos->y, 
			ip->pos->z);
		printf("%s %s(%i)(j)->pos = %.5le, %.5le, %.5le\n", 
		    d_h, per_table[js].sym, j, jp->pos->x, jp->pos->y, 
			jp->pos->z);

	    }
	    #endif

	    /* Compare the squared separation between (i, j) with the
	     * appropriate cutoff value.  If it is within cutoff, 
	     * perform operations to neighboring pair (i, j). */
	     
	    if (rij2 < rijCu2sq_)
	    {
		nij = nji = 0;
		#if DIAG
		if (D2_local > 1)
		{
		    printf("\t%s(%i)<->%s(%i):  exchanging bCards.\n", 
			per_table[is].sym, i, per_table[js].sym, j);
		}
		#endif
		
		atom_becomeNeighbors(ip, jp, &nij, &nji);
		fatal = (nij>=MN||nji>=MN);
		if (fatal)
		{
		    fprintf(stderr, "FATAL ERROR:  too many neighbors\n");
		    fprintf(stderr, "%s_%i\t%s_%i\n", 
			per_table[is].sym, i, per_table[js].sym, j);
		    fprintf(stderr, "%i\t%i\n", nij, nji);
		    fprintf(stderr, "MN is %i\n", 
			MN);
		    exit(0);
		}
		
		#if DIAG
		if (D2_local > 1)
		{
		    printf("\t%s(%i)<->%s(%i):  exchanged bCards.\n", 
			per_table[is].sym, i, per_table[js].sym, j);
		    fflush(stdout);
		}
		#endif
		
		/* Perform all pair-wise computations and store the results
		 * in the global arrays for use in the main computation 
		 * loop. */
		
		rij = rij_[i][nij] = rij_[j][nji] = sqrt(rij2);
		rij2_[i][nij] = rij2_[j][nji] = rij2;
		ptPtr_scalmult(pRijT_, pRijT_, 1.0/rij);
		ptPtr_copy(&(rij_hat[i][nij]), pRijT_);
		ptPtr_copy(&(rij_hat[j][nji]), 
			   ptPtr_scalmult(&TMP_PT, pRijT_, -1.0));
		
		#if DIAG
		if (D2_local)
		{
		    printf("\t%s(%i)<->%s(%i), nij=%i, nji=%i, rij=%.5le\n", 
			per_table[is].sym, i, per_table[js].sym, j, nij, nji, rij);
		    fflush(stdout);
		}
		#endif
		
		/* If this is a bonding (non-moliere) interaction:
		 *  (1) compute the smooth cutoff and its first derivative;
		 *  (2) compute the attractive Morse potential and its
		 *      first derivative.
		 */
		if (!moliere_only)
		{
		    Fij_ (rij, &fij, &fij1);
			if (fij<0.0){printf("err fij in pc f(%s_%i,%s_%i)=%.5le\n", 
			per_table[is].sym, i, per_table[js].sym, j, fij);exit(0);}
		    MorseVa(rij, fij, fij1, &Va, &dVa_drij, NULL);
		}
		else
		{
		    fij = fij1 = Va = dVa_drij = 0.0;
		}
		
		
		/* Vr is given by one of three functions, depending on
		 * rij: 
		 *  if (rij < Ms_ra)	    Vr = Moliere(rij)
		 *  else if (rij < Ms_rb)   Vr = Spline(rij)
		 *  else		    Vr = fij*Aij*exp(-lamij*rij)
		 */
		if (moliere_only || (rij < Ms_ra && USE_HIE_CORE))
		    Moliere(rij, &Vr, &dVr_drij, NULL);
		else if (rij < Ms_rb && USE_HIE_CORE)
		    MolSpline(rij, &Vr, &dVr_drij);
		else
		    MorseVr(rij, fij, fij1, &Vr, &dVr_drij, NULL);

		/* Store the values of fij, dfij/drij, 
		 * Va, Vr, dVa/drij, and dVr/drij in
		 * the proper global arrays. */
		 
		fij_[i][nij] = fij_[j][nji] = fij;
		fij1_[i][nij] = fij1_[j][nji] = fij1;
		Va_ij_[i][nij] = Va_ij_[j][nji] = Va;
		dVa_drij_[i][nij] = dVa_drij_[j][nji] = dVa_drij;
		Vr_ij_[i][nij] = Vr_ij_[j][nji] = Vr;
		dVr_drij_[i][nij] = dVr_drij_[j][nji] = dVr_drij;
		
		/* Now that fij is established, update the elemental
		 * coordinations of atoms 'i' and 'j'. */
		 
		if (jp->sym == F)			    NiF_[i] += fij;
		else if (jp->sym == C)			    NiC_[i] += fij;
		else if (jp->sym == Si || jp->sym == Si_0)  NiSi_[i] += fij;
		else if (jp->sym == Ar) {} /* do nothing */
		else{fprintf(stderr,"%s:%s?\n",d_h,per_table[jp->sym].sym);}

		if (ip->sym == F)			    NiF_[j] += fij;
		else if (ip->sym == C)			    NiC_[j] += fij;
		else if (ip->sym == Si || ip->sym == Si_0)  NiSi_[j] += fij;
		else if (ip->sym == Ar) {} /* do nothing */
		else{fprintf(stderr,"%s:s?\n",d_h,per_table[ip->sym].sym);}

		#if DIAG
		if (D2_local > 1)
		{
		    printf("\t%s(%i)---%s(%i) --> DATA STORAGE\n", 
			per_table[ip->sym].sym, i, 
			per_table[jp->sym].sym, j);
		    printf("\tR-hat(%i,%i) = %.5le, %.5le, %.5le\n", 
			i, j, rij_hat[i][nij].x, rij_hat[i][nij].y, rij_hat[i][nij].z);
		    printf("\tR-hat(%i,%i) = %.5le, %.5le, %.5le\n", 
			j, i, rij_hat[j][nji].x, rij_hat[j][nji].y, rij_hat[j][nji].z);
		    printf("\tr(%i,%i) = %.10f\n", i, j, rij_[i][nij]);
		    printf("\tr2(%i,%i) = %.10f\n", i, j, rij2_[i][nij]);
		    printf("\tf(%i,%i) = %.10f\n", i, j, fij_[i][nij]);
		    printf("\tf'(%i,%i) = %.10f\n", i, j, fij1_[i][nij]);
		    printf("\tVa(%i,%i) = %.10f, dVa_drij = %.10f\n",
			i, j, Va, dVa_drij);
		    printf("\tVr(%i,%i) = %.10f, dVr_drij = %.10f\n",
			i, j, Vr, dVr_drij);
		}
		#endif
	    }
	}
    } /* end for each atom (*ip) in AtomList L */

    #if DIAG
    if (D2_local)
    {
	printf("%s ====== End precomputations ======\n\n\n", d_h); 
    }
    #endif
} /* end sic::Precomputations */

static void atom_enforce (atomPtr ip, atomPtr jp, ptPtr frc, ptPtr f_ij, ptPtr f_ji)
{		    
    ptPtr_add(ip->frc, ip->frc, frc);
    ptPtr_add(f_ij, f_ij, frc);
    ptPtr_scalmult(frc, frc, -1.0);
    ptPtr_add(jp->frc, jp->frc, frc);
    ptPtr_add(f_ji, f_ji, frc);
}		    

static double a__, af1__, af0__, g0om__, g1om__, af0g1om__, 
	      acc1__, acc2__, acc3__, acc3ik__, xiH_ij__, xiH_ji__, dxi_ijk__, thisPhij;
static ptPtr f_ij__, f_ji__, f_ik__, f_ki__, f_jk__, f_kj__, f_kl__, f_lk__;
char * mainForce_d_h_ = "tbt_sicf::mainForce";
void mainForce (void)
{
    atomPtr ip = NULL, jp = NULL, kp = NULL, lp = NULL;
    nNodePtr nip = NULL, njp = NULL, nkp = NULL, nlp = NULL, rare_nkp = NULL;
    nNodePtr nikp=NULL, nklp=NULL, njkp=NULL;
    char * d_h = mainForce_d_h_;
    int i = 0, j = 0, k = 0, p = 0;
    int nij = 0, nji = 0;
    int nik = 0, njk = 0;
    int nil = 0, njl = 0;
    int nkl = 0, nki = 0, nkj = 0, rare_nik = 0, rare_njk = 0, nlk=0;
    sym_type is, js, ks;
    short cc_ij = 0;  /* set to one if (i,j) is a c-c bond */
    short allrep = 0; /* set to one if (i,j) is a fully repulsive interaction */
    
    #if DIAG
    D3_local = SIC_DIAG3_;
    if (D3_local)
    {
	printf("%s begins.\n", d_h);
	fflush(stdout);
    }
    #endif

    Ni = Nj = Nk = 
    NiF = NiC = NiSi =
    NjF = NjC = NjSi =
    NkF = NkC = NkSi =
    Nij_F = Nji_F = Nij_C = Nji_C = Nij_Si = Nji_Si =
    Nij = Nji = Nik = Nki = Njk = Nkj = Nij_conj =
    rij = fij = fij1 = rik = fik = fik1 = rjk = fjk = fjk1 = fkl1 =
    Fik = Fik1 = Fjk = Fjk1 = Fil = Fil1 = Fjl = Fjl1 = 
    xi_ij = xi_ji = b_ij = b_ji = 
    beta_ij = beta_ji = 
    Hij = Hji = Hij1 = Hij2 = Hji1 = Hji2 =
    sijk = sjik = omijk = omjik = epijk = epjik = 
    gijk = gijk1 = gjik = gjik1 = 
    Fcc = Fcc1 = Fcc2 = Fcc3 = 0.0;
    Rij_hatp = Rik_hatp = Rjk_hatp = Rkl_hatp = NULL;
    f_ij__ = f_ji__ = f_ik__ = f_ki__ = f_jk__ = f_kj__ = f_kl__ = f_lk__ = NULL;
 
    for (ip = L_; ip->next; ip = ip->next)
    {
	i = ip->id;
	is = ip->sym;

	Nij_F = NiF = NiF_[i];
	Nij_C = NiC = NiC_[i];
	Nij_Si = NiSi = NiSi_[i];
	
	Ni = NiF + NiC + NiSi;

	#if DIAG
	if (D2_local)
	{
	    printf("%s Atom %s(%i)\n", d_h, per_table[is].sym, i);
	    for (njp=ip->nList;njp;njp=njp->next)
		printf("%i:%s_%i ", njp->index, 
		    per_table[njp->addr->sym].sym, njp->addr->id);
	    printf("\n");
	    fflush(stdout);
	}
	#endif

	/* Identify the unique pair of atoms (i,j):  The nList was
	 * constructed so that the atoms are represented in *descending*
	 * order on the list.  Therefore, the list is "divided" into two
	 * segments:  (1) the first segment contains references to atoms
	 * neighboring the list owner ("i") whose id's are *greater* than
	 * the list owner's,  and (2) the second segment contains references
	 * to atoms neighboring the list owner whose id's are *less* than
	 * the list owner's.  To ensure that only *unique* pairs of atoms
	 * are identified (i.e., we don't want to visit pair (j, i) after
	 * we've already visited pair (i, j) since they're the same), we
	 * identify a pair as (1) the list owner and (2) a list member whose
	 * id is greater than that of the list owner.  So, we visit the list
	 * members sequentially, each one defining a unique pair with the list
	 * owner, until we visit a list member whose id is less than the list
	 * owner's.  In this fashion we visit each unique pair exactly once. */
	for (njp = ip->nList; njp && njp->addr->id > i; njp = njp->next)
	{
	    jp = njp->addr;
	    j = jp->id;
	    ptPtr_clear(pF2_ij);
	    ptPtr_clear(pF2_ji);
	    for (k = 0; k < MN; k++) 
	    {	
		sijk_[k] = gijk_[k] = gijk1_[k] = omijk_[k] = epijk_[k] = 0.0;
		sjik_[k] = gjik_[k] = gjik1_[k] = omjik_[k] = epjik_[k] = 0.0;
	    }
	
	    js = jp->sym;
	    TwoBody_SetParameters(is, js);
	    if (is == C && js == C) cc_ij = 1;
	    else cc_ij = 0;
	    if (chem_isInert(ip->sym) || chem_isInert(jp->sym)) allrep = 1;
	    else allrep = 0;
	    nij = njp->index;
	    rij = rij_[i][nij];
	    rij2 = rij2_[i][nij];
	    Rij_hatp = &(rij_hat[i][nij]);
	    f_ij__ = &(frc_ij_[i][nij]);
	
	    #if DIAG
	    if (D3_local)
	    {
		printf("%i:%s_%i\n", njp->index, per_table[jp->sym].sym, jp->id);
		fflush(stdout);
	    }
	    #endif
	
	    fij = fij_[i][nij];
	    fij1 = fij1_[i][nij];

	    Nji_F = NjF = NiF_[j];
	    Nji_C = NjC = NiC_[j];
	    Nji_Si = NjSi = NiSi_[j];
	    
	    Nj = NjF + NjC + NjSi;

	    Va = Va_ij_[i][nij];
	    Vr = Vr_ij_[i][nij];
	    dVa_drij = dVa_drij_[i][nij];
	    dVr_drij = dVr_drij_[i][nij];
	    Va_2 = Va/2.0;
	
	    /* Compute Nij, Nji, Nij(F), Nij(Si), Nij(C),
	     * Nji(F), Nji(Si), Nji(C) */
	    Nij = Ni - fij;
	    Nji = Nj - fij;
	    /* Set Nij_x to their initial values */
	    Nij_F = NiF;
	    Nij_C = NiC;
	    Nij_Si = NiSi;
	    if (js == F) Nij_F -= fij;
	    else if (js == Si) Nij_Si -= fij;
	    else if (js == C) Nij_C -= fij;
	    if (is == F) Nji_F -= fij;
	    else if (is == Si) Nji_Si -= fij;
	    else if (is == C) Nji_C -= fij;
	
	    #if DIAG
	    if (D3_local)
	    {
		printf("Searching for i(%s_%i) in j's(%s_%i) nbrs...\n", 
		    per_table[is].sym, i, per_table[jp->sym].sym, jp->id);
		printf("ip %x jp %x ip->nList %x jp->nList %x\n", 
		    ip, jp, ip->nList, jp->nList);
	    
		for (nip = jp->nList; nip && nip->addr != ip; nip = nip->next)
		{
		    printf("%i:%s_%i ", nip->index, 
			    per_table[nip->addr->sym].sym, nip->addr->id);
		}
		printf("\n");
		fflush(stdout);
	    }
	    #endif
	    
	    /* Assign nip so that it points to the nNode on j's neighbor list
	     * that corresponds to atom i. */
	    for (nip = jp->nList; nip && nip->addr != ip; nip = nip->next);

	    if (!nip)
	    {
		printf("%s --> strangeness in triplet loop:\n", d_h);
		printf("\tI(%s_%i)->J(%s_%i) but J!->I.\n", 
		    per_table[is].sym, i, per_table[js].sym, j);
		fflush(stdout);
		exit(0);
	    }
	    nji = nip->index;
	    nip = NULL;
	    f_ji__ = &(frc_ij_[j][nji]);

	    #if DIAG
	    if (D3_local > 1)
	    {
		printf("\n\n\tBonded Pair %s(%i)--%s(%i).\n", 
		    per_table[is].sym, i, per_table[js].sym, j);
		printf("\t(j)%s(%i) is Neighbor #%i to (i)%s(%i).\n", 
		    per_table[js].sym, j, nij, per_table[is].sym, i);
		printf("\t(i)%s(%i) is Neighbor #%i to (j)%s(%i).\n", 
		    per_table[is].sym, i, nji, per_table[js].sym, j);
		printf("\tDATA RETRIEVAL:\n");
		printf("\tR-hat(%i,%i) = %.5le, %.5le, %.5le\n", 
		    i, j, Rij_hatp->x, Rij_hatp->y, Rij_hatp->z);
		printf("\tr(%i,%i) = %.10f\n", i, j, rij_[i][nij]);
		printf("\tf(%i,%i) = %.10f\n", i, j, fij_[i][nij]);
		printf("\tf'(%i,%i) = %.10f\n", i, j, fij1_[i][nij]);
		printf("\tVa(%i,%i) = %.10f, dVa_drij = %.10f\n",
		    i, j, Va, dVa_drij);
		printf("\tVr(%i,%i) = %.10f, dVr_drij = %.10f\n",
		    i, j, Vr, dVr_drij);
		printf("\tN_F(%i) = %.5le, N_F(%i) = %.5le\n", 
		    i, NiF, j, NjF);
		printf("\tN_F(%i,%i) = %.5le, N_F(%i,%i) = %.5le\n", 
		    i, j, Nij_F, j, i, Nji_F);
		printf("\tN_C(%i) = %.5le, N_C(%i) = %.5le\n", 
		    i, NiC, j, NjC);
		printf("\tN_C(%i,%i) = %.5le, N_C(%i,%i) = %.5le\n", 
		    i, j, Nij_C, j, i, Nji_C);
		printf("\tN_Si(%i) = %.5le, N_Si(%i) = %.5le\n", 
		    i, NiSi, j, NjSi);
		printf("\tN_Si(%i,%i) = %.5le, N_Si(%i,%i) = %.5le\n", 
		    i, j, Nij_Si, j, i, Nji_Si);
	    }
	    #endif
		
	    /* Compute Hij and derivatives if i is C and j is C or F */
	    Hij = Hij1 = Hij2 = Hji = Hji1 = Hji2 = 0.0;
	    if (cc_ij && USE_HIJ)
	    {
		if (HIJ_TESTER) 
		{
		    Hcc_tester(Nij_F, Nij_C+Nij_Si, &Hij, &Hij1, &Hij2);
		    Hcc_tester(Nji_F, Nji_C+Nji_Si, &Hji, &Hji1, &Hji2);
		}
		else
		{
		    Hcc_bicubicint(Nij_F, Nij_C+Nij_Si, &Hij, &Hij1, &Hij2);
		    Hcc_bicubicint(Nji_F, Nji_C+Nji_Si, &Hji, &Hji1, &Hji2);
		}
		#if DIAG
		if (D3_local > 1)
		{
		    printf("\tDATA COMPUTED AND STORED:\n");
		    printf("\tHcc(%.5le,%.5le)_(%i->%i) = %.5le\n", 
			Nij_F, Nij_C+Nij_Si, i, j, Hij);
		    printf("\tHcc(%.5le,%.5le)_(%i->%i) = %.5le\n", 
		    Nji_F, Nji_C+Nji_Si, j, i, Hji);
		}
		#endif
	    }
	    else if (is == C && js == F && USE_HIJ)
	    {
		if (HIJ_TESTER) Hcc_tester(Nij_F, Nij_C+Nij_Si, &Hij, &Hij1, &Hij2);
		else Hcf_bicubicint(Nij_F, Nij_C+Nij_Si, &Hij, &Hij1, &Hij2);
		#if DIAG
		if (D3_local > 1)
		{
		    printf("\tDATA COMPUTED AND STORED:\n");
		    printf("\tHcf(%.5le,%.5le)_(%i->%i) = %.5le, (1)%.5le, (2)%.5le\n", 
			Nij_F, Nij_C+Nij_Si, i, j, Hij, Hij1, Hij2);
		}
		#endif
	    }
	    else if (is == F && js == C && USE_HIJ)
	    {
		if (HIJ_TESTER) Hcc_tester(Nji_F, Nji_C+Nji_Si, &Hji, &Hji1, &Hji2);
		else Hcf_bicubicint(Nji_F, Nji_C+Nji_Si, &Hji, &Hji1, &Hji2);
		#if DIAG
		if (D3_local > 1)
		{
		    printf("\tDATA COMPUTED AND STORED:\n");
		    printf("\tHcf(%.5le,%.5le)_(%i->%i) = %.5le, (1)%.5le, (2)%.5le\n", 
		    Nji_F, Nji_C+Nji_Si, j, i, Hji, Hji1, Hji2);
		}
		#endif
	    }
	    if (Hij!=Hij||Hij==HUGE_VAL) {printf("error hij in mf\n");exit(0);}
	    if (Hij1!=Hij1||Hij1==HUGE_VAL) {printf("error hij1 in mf\n");exit(0);}
	    if (Hij2!=Hij2||Hij2==HUGE_VAL) {printf("error hij1 in mf\n");exit(0);}
	    if (Hji!=Hji||Hji==HUGE_VAL) {printf("error hji in mf\n");exit(0);}
	    if (Hji1!=Hji1||Hji1==HUGE_VAL) {printf("error hji1 in mf\n");exit(0);}
	    if (Hji2!=Hji2||Hji2==HUGE_VAL) {printf("error hji1 in mf\n");exit(0);}

	    /* Three Body precomputations and storage */
	
	    xi_ij = 0.0;
	    Nij_conj = 1.0;
	    b_ij = 1.0;
	    beta_ij = 0.0;
	    #if THREE_BODY
	    /* "allrep" is 1 if either "i" or "j" is an inert atom (e.g. argon) */
	    if (!allrep) 
	    {
		#if DIAG
		if (D3_local > 1)
		{
		    printf("\tLoop over neighbors of 'i' %s_%i\n", 
			per_table[is].sym, i);
		    fflush(stdout);
		}
		#endif
		for (nkp = ip->nList; nkp; nkp = nkp->next)
		{
		    kp = nkp->addr;
		    if (kp != jp && !chem_isInert(kp->sym))
		    {
			k = kp->id;
			ks = kp->sym;
			ThreeBody_SetParameters(is, js, ks);
			#if DIAG
			if (D3_local > 1)
			{
			    printf("\tThreebody parameters, (%s-%s-%s), (%i-%i-%i)\n", 
				per_table[is].sym, per_table[js].sym, per_table[ks].sym, 
				i, j, k);
			}
			#endif
			nik = nkp->index;
			rik = rij_[i][nik];
			rik2 = rij2_[i][nik];
			fik = fij_[i][nik];
			if (fik<0.0){printf("err fik(dxi) in mf f(%s_%i,%s_%i)=%.5le\n", 
			per_table[is].sym, i, per_table[ks].sym, k, fik);exit(0);}
			Nik = Ni - fik;
			Nk = NiF_[k] + NiC_[k] + NiSi_[k];
			Nki = Nk - fik;
			
			/* Compute rjk2 for sijk computation*/
			Rjk_hatp = &TMP_PT;
			ptPtr_minimg(Rjk_hatp, 
			    ptPtr_subtract(Rjk_hatp, jp->pos, kp->pos));
			rjk2 = ptPtr_sqabs(Rjk_hatp);
		    
			/* Compute sijk = cos(theta_ijk) */
			sijk = (rij2 + rik2 - rjk2)/(2*rij*rik);
		    
			/* Compute gijk and its first derivative */
			g_costh(sijk, is, &gijk, &gijk1);
			if(gijk<0.0){printf("error gijk in mf g(%s_%i,%s_%i, %s_%i)=%.5le\n", 
			per_table[is].sym, i, per_table[js].sym, j, per_table[ks].sym, k, gijk);exit(0);}
			
			/* Compute omega_ijk and epsilon_ijk */
			omega_(rij, rik, &omijk, &epijk);
			if(omijk<0.0){printf("error omijk in mf om(%s_%i,%s_%i, %s_%i)=%.5le\n", 
			per_table[is].sym, i, per_table[js].sym, j, per_table[ks].sym, k, omijk);exit(0);}
		    
			/* Store the computed data */
			sijk_[nik] = sijk;
			gijk_[nik] = gijk;
			gijk1_[nik] = gijk1;
			omijk_[nik] = omijk;
			epijk_[nik] = epijk;
		    
			/* Tally xi_ij */
			dxi_ijk__ = fik*gijk*omijk;
			if(dxi_ijk__<0.0){printf("error dxi_ijk in mf dxi_ijk(%s_%i,%s_%i, %s_%i)=%.5le\n", 
			per_table[is].sym, i, per_table[js].sym, j, per_table[ks].sym, k, dxi_ijk__);exit(0);}
			xi_ij += dxi_ijk__;
			
			/* If i--j is a carbon-carbon bond, update Nij_conj */
			if (cc_ij && USE_FCC)
			{
			    /* Compute F_ik and F'_ik. */
			    Fconj(Nki, &(Fik_[i][nik]), &(Fik1_[i][nik]));
			    if (ks == C) Nij_conj += fik*Fik_[i][nik];			
			}
		    
			#if DIAG
			if (D3_local > 1)
			{
			    printf("\tk--ij: Triplet %s(%i)(i->%i) - - %s(%i)--%s(%i)\n",
				per_table[ks].sym, k, nik,  
				per_table[is].sym, i,
				per_table[js].sym, j);
			    printf("\tDATA RETRIEVED(ik):\n");
			    printf("\tr(%i,%i)(%.5le) r2(%i,%i)(%.5le) f(%i,%i)(%.5le)\n", 
				i, k, rik, i, k, rik2, i, k, fik);
			    printf("\tDATA COMPUTED BUT NOT STORED:\n");
			    printf("\tR(%i,%i) = %.5le, %.5le, %.5le\n", 
				j, k, Rjk_hatp->x, Rjk_hatp->y, Rjk_hatp->z);
			    printf("\tr2(%i,%i)(%.5le), r(%i,%i)(%.5le)\n", 
				j, k, rjk2, j, k, sqrt(rjk2));
			    printf("\tDATA COMPUTED AND STORED(ijk):\n");
			    printf("\ts(%i,%i,%i)(%.5le) g(%i,%i,%i)(%.5le) g'(%i,%i,%i)(%.5le)\n", 
				i, j, k, sijk, i, j, k, gijk, i, j, k, gijk1);
			    printf("\tom(%i,%i,%i)(%.5le) epijk(%i,%i,%i)(%.5le)\n", 
				i, j, k, omijk, i, j, k, epijk);
			    printf("\txi(%i,%i,%i) = %.5le\n", i, j, k, fik*gijk*omijk);
			    if (cc_ij && USE_FCC)
			    {
				printf("\tN(%i,%i)(%.5le) F(%i,%i)(%.5le) F1(%i,%i)(%.5le)\n", 
				    i, k, Nik, i, k, Fik_[i][nik], i, k, Fik1_[i][nik]);
				printf("\tN(%i,%i)_conj -> %.5le\n", i, j, Nij_conj);
			    }
			    fflush(stdout);
			}
			#endif
		    }
		}
	   
		/* Bond order computation for i->j.  "eta" and "delta"
		 * were globally set by ThreeBody_SetParameters for
		 * the triplet subtended by atom "i". */
	    
		xiH_ij__ = xi_ij + Hij;
		if (eta<1.0&&xiH_ij__<0.0) {printf("error xiH_ij in mf xi(%s_%i-%s_%i)(%.6le) Hij(%.6le)\n", 
		per_table[is].sym, i, per_table[js].sym, j, xi_ij, Hij);exit(0);}
		if (xiH_ij__) b_ij += pow(xiH_ij__, eta);
		b_ij = pow(b_ij, -delta);
		
		if (xiH_ij__) beta_ij = -delta*eta*pow(b_ij, (delta+1)/delta)
		      *(eta == 1 ? 1.0 : pow(xiH_ij__, eta-1.0));
		      
		if (beta_ij!=beta_ij||beta_ij==HUGE_VAL) {printf("error beta_ij in mf\n");exit(0);}
	    }
	    #endif /* if THREE_BODY */
	    
	    #if DIAG
	    if (D3_local > 1)
	    {
		printf("\tBond order (%i->%i) is %.6le\n", i, j, b_ij);
		printf("\txi_(%i,%i)(%.5le), H(%i,%i)(%.5le), eta(%.5le), delta(%.5le)\n", 
		    i, j, xi_ij, i, j, Hij, eta, delta);
		printf("\tbeta(%i->%i) is %.6le\n", i, j, beta_ij);
	    }
	    #endif
	    
	    xi_ji = 0.0;
	    b_ji = 1.0;
	    beta_ji = 0.0;
	    #if THREE_BODY
	    if (!allrep)
	    {
		#if DIAG
		if (D3_local > 1)
		{
		    printf("\tLoop over neighbors of 'j' %s_%i\n", 
			per_table[js].sym, j);
		    fflush(stdout);
		}
		#endif	    
		for (nkp = jp->nList; nkp; nkp = nkp->next)
		{
		    kp = nkp->addr;
		    if (kp != ip && !chem_isInert(kp->sym))
		    {
			k = kp->id;
			ks = kp->sym;
			ThreeBody_SetParameters(js, is, ks);
			#if DIAG
			if (D3_local > 1)
			{
			    printf("\tThreebody parameters, (%s-%s-%s), (%i-%i-%i)\n", 
				per_table[js].sym, per_table[is].sym, per_table[ks].sym, 
				j, i, k);
			}
			#endif
			njk = nkp->index;
			rjk = rij_[j][njk];
			rjk2 = rij2_[j][njk];
			fjk = fij_[j][njk];
			Njk = Nj - fjk;
			Nk = NiF_[k] + NiC_[k] + NiSi_[k];
			Nkj = Nk - fjk;
		    
			/* Compute Rik, rik2 and rik */
			Rik_hatp = &TMP_PT;
			ptPtr_minimg(Rik_hatp, 
			    ptPtr_subtract(Rik_hatp, ip->pos, kp->pos));
			rik2 = ptPtr_sqabs(Rik_hatp);
		    
			/* Compute sjik = cos(theta_jik) */
			sjik = (rij2 + rjk2 - rik2)/(2*rij*rjk);
		    
			/* Compute gjik and its first derivative */
			g_costh(sjik, js, &gjik, &gjik1);
			    
			/* Compute omega_jik and epsilon_jik */
			omega_(rij, rjk, &omjik, &epjik);
		    
			/* Store the computed data */
			sjik_[njk] = sjik;
			gjik_[njk] = gjik;
			gjik1_[njk] = gjik1;
			omjik_[njk] = omjik;
			epjik_[njk] = epjik;

			/* Tally xi_ji */
			xi_ji += fjk*gjik*omjik;

			/* If j--i is a carbon-carbon bond, update Nij_conj */
			if (cc_ij && USE_FCC)
			{
			    /* Compute F_jk and F'_jk. */
			    Fconj(Nkj, &(Fik_[j][njk]), &(Fik1_[j][njk]));
			    if (ks == C) Nij_conj += fjk*Fik_[j][njk];			
			}
		    
			#if DIAG
			if (D3_local > 1)
			{
			    printf("\tk--ji: Triplet %s(%i)(j->%i) - - %s(%i)--%s(%i)\n",
				per_table[ks].sym, k, njk, 
				per_table[js].sym, j,
				per_table[is].sym, i);
			    printf("\tDATA RETRIEVED(jk):\n");
			    printf("\tr(%i,%i)(%.5le) r(%i,%i)2(%.5le) f(%i,%i)(%.5le)\n", 
				j, k, rjk, j, k, rjk2, j, k, fjk);
			    printf("\tDATA COMPUTED BUT NOT STORED:\n");
			    printf("\tR(%i,%i) = %.5le, %.5le, %.5le\n", 
				i, k, Rik_hatp->x, Rik_hatp->y, Rik_hatp->z);
			    printf("\tr2(%i,%i)(%.5le), r(%i,%i)(%.5le)\n",
				i, k, rik2, i, k, sqrt(rik2));
			    printf("\tDATA COMPUTED AND STORED(jik):\n");
			    printf("\ts(%i,%i,%i)(%.5le) g(%i,%i,%i)(%.5le) g'(%i,%i,%i)(%.5le)\n", 
				j, i, k, sjik, j, i, k, gjik, j, i, k, gjik1);
			    printf("\tom(%i,%i,%i)(%.5le) ep(%i,%i,%i)(%.5le)\n", 
				j, i, k, omjik, j, i, k, epjik);
			    printf("\txi(%i,%i,%i) = %.5le\n", j, i, k, fjk*gjik*omjik);
			    if (cc_ij && USE_FCC)
			    {
				printf("\tN(%i,%i)(%.5le) F(%i,%i)(%.5le) F1(%i,%i)(%.5le)\n", 
				    j, k, Njk, j, k, Fik_[j][njk], j, k, Fik1_[j][njk]);
				printf("\tN(%i,%i)_conj -> %.5le\n", i, j, Nij_conj);
			    }
			    fflush(stdout);
			}
			#endif
		    }
		}

		/* Bond order computation for j->i.  "eta" and "delta"
		 * were globally set by ThreeBody_SetParameters for
		 * the triplet subtended by atom "j". */

		xiH_ji__ = xi_ji + Hji;
		if (eta<1.0&&xiH_ji__<0.0) {printf("error xiH_ji in mf xi(%s_%i-%s_%i)(%.6le) Hji(%.6le)\n", 
		per_table[js].sym, j, per_table[is].sym, i, xi_ji, Hji);exit(0);}
		if (xiH_ji__) b_ji += pow(xiH_ji__, eta);
		b_ji = pow(b_ji, -delta);
	    
		if (xiH_ji__) beta_ji = -delta*eta*pow(b_ji, (delta+1)/delta)
		      *(eta == 1 ? 1.0 : pow(xiH_ji__, eta-1.0));
		      
		if (beta_ji!=beta_ji||beta_ji==HUGE_VAL) {printf("error beta_ji in mf\n");exit(0);}
	    }
	    #endif /* if THREEBODY */
	    
	    #if DIAG
	    if (D3_local > 1)
	    {
		printf("\tBond order (%i->%i) is %.6le\n", j, i, b_ji);
		printf("\txi_(%i,%i)(%.5le), H(%i,%i)(%.5le), eta(%.5le), delta(%.5le)\n", 
		    j, i, xi_ji, j, i, Hji, eta, delta);
		printf("\tbeta(%i->%i) is %.6le\n", j, i, beta_ji);
	    }
	    #endif
	    
	    /* If i--j is a carbon-carbon bond, compute Fcc and its
	     * three first derivatives. */
	    Fcc=Fcc1=Fcc2=Fcc3=0.0;
	    if (cc_ij && USE_FCC)
	    {
		if (!FCC_TESTER) Fcc_tricubicint(Nij, Nji, Nij_conj, &Fcc, &Fcc1, &Fcc2, &Fcc3);
		if (FCC_TESTER) Fcc_tester(Nij, Nji, Nij_conj, &Fcc, &Fcc1, &Fcc2, &Fcc3);
		if (Fcc!=Fcc||Fcc==HUGE_VAL) {printf("error fcc in mf\n");exit(0);}
		if (Fcc1!=Fcc1||Fcc1==HUGE_VAL) {printf("error fcc1 in mf\n");exit(0);}
		if (Fcc2!=Fcc2||Fcc2==HUGE_VAL) {printf("error fcc2 in mf\n");exit(0);}
		if (Fcc3!=Fcc3||Fcc3==HUGE_VAL) {printf("error fcc3 in mf\n");exit(0);}
		#if DIAG
		if (D3_local > 1)
		{
		    printf("\tFcc(%i,%i)[%.5le, %.5le, %.5le] = %.5le\n", 
			i, j, Nij, Nji, Nij_conj, Fcc);
		    printf("\tFcc1 = %.5le, Fcc2 = %.5le, Fcc3 = %.5le\n", 
			Fcc1, Fcc2, Fcc3);
		}
		#endif
		acc1__ = Va_2*Fcc1;
		acc2__ = Va_2*Fcc2;
		acc3__ = Va_2*Fcc3;
	    }

	    /* Average bond order computation */
	    b_ij_bar = 0.5*(b_ij + b_ji + Fcc);
	    
	    #if DIAG
	    if (D3_local > 1)
	    {
		printf("\n\tb_ij_bar\t = 0.5*(%.5e + %.5e + %.5e)\n", 
		    b_ij, b_ji, Fcc);
		printf("\t\t\t = %.7e\n", b_ij_bar);
	    }
	    #endif
	    
	    /* Compute and store Phi_ij_ */
	    thisPhij=Vr-b_ij_bar*Va;
	    PE_+=thisPhij;
	    dtPhi_ij_[i][nij] = thisPhij-Phi_ij_[i][nij]; 
	    dtPhi_ij_[j][nji] = thisPhij-Phi_ij_[j][nji]; 
	    Phi_ij_[i][nij] = Phi_ij_[j][nji] = thisPhij;
	    
	    /* Apply the requested bonding definition parameter 
	     * to decide whether these two atoms belong in the
	     * same cluster. */
	    if (bond_ctrl_==1) 
	    {
		if (fij_[i][nij]>0) clust_enclust(ip, jp);
	    }
	    else if (bond_ctrl_==2)
	    {
		if (fij_[i][nij]==1.0&&Phi_ij_[i][nij]<0.0)
		    clust_enclust(ip, jp);
	    }
	    /* Continue tally of twoBodyEnergy and threeBodyEnergy */
	    twoBodyEnergy += Vr - Va;
	    threeBodyEnergy += -(b_ij_bar - 1.0)*Va;
	    
	    #if DIAG
	    if (D3_local > 1)
	    {
		printf("\n\tb_(%i,%i)_bar = %.5le, phi_(%i,%i) = %.5le\n", 
		    i, j, b_ij_bar, i, j, Phi_ij_[i][nij]);
	    }
	    #endif

	    /* Two body forces:  vector force F2_ij is computed
	     * from dVr_drij, dVa_drij, b_ij_bar, and the unit
	     * vector Rij-hat.  Force on atom "i" is incremented by
	     * F2_ij and force on atom "j" is incremented by 
	     * -F2_ij. */
	    ptPtr_scalmult(pF2_ij, Rij_hatp, -dVr_drij+b_ij_bar*dVa_drij);
	    ptPtr_add(ip->frc, ip->frc, pF2_ij);
	    ptPtr_add(f_ij__, f_ij__, pF2_ij);
	    ptPtr_scalmult(pF2_ji, pF2_ij, -1.0);
	    ptPtr_add(jp->frc, jp->frc, pF2_ji);
	    ptPtr_add(f_ji__, f_ji__, pF2_ji);
	    
	    #if DIAG
	    if (D3_local > 1)
	    {
		printf("\t2Body forces:\n");
		printf("\tF2_(%i<-%i) = -F2_(%i<-%i) = %.5le, %.5le, %.5le\n", 
		    i, j, j, i, pF2_ij->x, pF2_ij->y, pF2_ij->z);
		printf("%s 2body forces (%i,%i) complete.  Entering 3body force loop.\n\n", d_h, 
		    i, j);
	    }
	    #endif
	    
	    /* Three Body force computations */
	    #if THREE_BODY
	    a__ = Va_2*beta_ij;
	    if (!allrep) for (nkp = ip->nList; nkp; nkp = nkp->next)
	    {
		kp = nkp->addr;
		if (kp != jp && !chem_isInert(kp->sym))
		{
		    /* Fetch stored precomputations. */
		    
		    k	    = kp->id;
		    ks	    = kp->sym;
		    nik	    = nkp->index;
		    rik	    = rij_[i][nik];
		    rik2    = rij2_[i][nik];
		    fik	    = fij_[i][nik];
		    fik1    = fij1_[i][nik];
		    Fik	    = Fik_[i][nik];
		    Fik1    = Fik1_[i][nik];
		    Nik	    = Ni - fik;
		    sijk    = sijk_[nik];	    
		    gijk    = gijk_[nik];	    
		    gijk1   = gijk1_[nik];	    
		    omijk   = omijk_[nik];	    
		    epijk   = epijk_[nik];
		    Rik_hatp = &(rij_hat[i][nik]);
		    f_ik__  = &(frc_ij_[i][nik]);
		    
		    for (nikp=kp->nList;nikp&&nikp->addr!=ip;nikp=nikp->next);
		    if (nikp) nki=nikp->index;
		    else {printf("error: no parity, i, k\n");exit(0);}
		    f_ki__  = &(frc_ij_[k][nki]);
		    
		    #if DIAG
		    if (D3_local > 1)
		    {
			printf("\n\tk--ij: Triplet %s(%i)(i->%i) - - %s(%i)--%s(%i)\n",
			    per_table[ks].sym, k, nik,  
			    per_table[is].sym, i,
			    per_table[js].sym, j);
			printf("\tDATA RETRIEVAL:\n");
			printf("\tr(%i,%i)(%.5le) r(%i,%i)2(%.5le) f(%i,%i)(%.5le) f'(%i,%i)(%.5le)\n", 
			    i, k, rik, i, k, rik2, i, k, fik, i, k, fik1);
			printf("\ts(%i,%i,%i)(%.5le) g(%i,%i,%i)(%.5le) g'(%i,%i,%i)(%.5le)\n", 
			    i, j, k, sijk, i, j, k, gijk, i, j, k, gijk1);
			printf("\tom(%i,%i,%i)(%.5le) ep(%i,%i,%i)(%.5le)\n", 
			    i, j, k, omijk, i, j, k, epijk);
			fflush(stdout);
		    }
		    #endif
		    
		    /* Compute force kernel factors. */
		    
		    af1__	= a__*fik1;
		    af0__	= a__*fik;
		    g0om__	= gijk*omijk;
		    g1om__	= gijk1*omijk;
		    af0g1om__	= af0__*g1om__;

		    #if DIAG
		    if (D3_local > 1)
		    {
			printf("\tForce kernel factors, (%i - %i<->%i):\n", k, i, j);
			printf("\ta(%.5le) af0(%.5le), af1(%.5le)\n", a__, af0__, af1__);
			printf("\tg0om(%.5le), g1om(%.5le), af0g1om(%.5le)\n",
			    g0om__, g1om__, af0g1om__);
		    }
		    #endif
		    
		    /* Compute force kernel magnitudes. */
		    
		    fknl__[0]=af1__*(g0om__+Hij1*(ks==F)+Hij2*(ks==C||ks==Si));
		    fknl__[1]=af0g1om__/rik;
		    fknl__[3]=af0g1om__/rij;
		    fknl__[5]=af0__*g0om__*epijk;
		    fknl__[2]=-fknl__[3]*sijk;
		    fknl__[4]=-fknl__[1]*sijk;
		    fknl__[6]=-fknl__[5];
		    
		    #if DIAG
		    if (D3_local > 1)
		    {
			printf("\tForce kernel magnitudes, (%i - %i<->%i):\n", k, i, j);
			for (p=0;p<8;p++) printf("\tfknl[%i] = %.5le\n", p, fknl__[p]);
		    }
		    #endif
		    
		    /* Use the TMP_PT to compute each force
		     * kernel and add it to the proper atom's force
		     * vector. */
		    		    
		    /* Kernels 0, 4, and 6 add and then multiply
		     * Rik-hat;  resulting force vector adds to i's
		     * force and subtracts from k's force. */
		    fknl__[0]+=fknl__[4]+fknl__[6];
		    ptPtr_scalmult(&TMP_PT, Rik_hatp, fknl__[0]);
		    atom_enforce (ip, kp, &TMP_PT, f_ik__, f_ki__);

		    /* Kernel 1 multiplies Rij-hat; it adds to i's force,
		     * and subtracts from k's force. */
		    ptPtr_scalmult(&TMP_PT, Rij_hatp, fknl__[1]);
		    atom_enforce (ip, kp, &TMP_PT, f_ik__, f_ki__);
		    
		    /* Kernels 2 and 5 add, then multiply 
		     * Rij-hat; resulting force vector adds to i's force,
		     * and subtracts from j's force. */
		    fknl__[2]+=fknl__[5];
		    ptPtr_scalmult(&TMP_PT, Rij_hatp, fknl__[2]);
		    atom_enforce (ip, jp, &TMP_PT, f_ij__, f_ji__);
		    
		    /* Kernel 3 multiplies Rik-hat; it adds to i's force,
		     * and subtracts from j's force. */
		    ptPtr_scalmult(&TMP_PT, Rik_hatp, fknl__[3]);
		    atom_enforce (ip, jp, &TMP_PT, f_ij__, f_ji__);
		    
		    /* Forces if (ij) is a C-C bond:
		     * Reusing force kernel magnitudes 0-3. */
		    if (cc_ij && USE_FCC)
		    {
			fknl__[0] = acc1__ * fik1;
			
			/* Kernel 0 multiplies Rik-hat; it adds to i's
			 * force and subtracts from k's force. */
			ptPtr_scalmult(&TMP_PT, Rik_hatp, fknl__[0]);
			atom_enforce (ip, kp, &TMP_PT, f_ik__, f_ki__);
			
			if (ks == C)  /* k-i-j potentially conjugated */
			{
			    
			    fknl__[1] = acc3__ * fik1 * Fik;
			    
			    /* Kernel 1 multiplies Rik-hat; it adds to i's
			     * force and subtracts from k's force. */
			    ptPtr_scalmult(&TMP_PT, Rik_hatp, fknl__[1]);
			    atom_enforce (ip, kp, &TMP_PT, f_ik__, f_ki__);

			    /* The next kernel should only be non-zero in the rare event
			     * when k is a neighbor of BOTH i and j. The values of fjk
			     * and Fjk1 must be fetched from storage, but the ordinate of
			     * k on j's neighbor list must be determined first. */
			    for (rare_nkp = jp->nList; 
				 rare_nkp && rare_nkp->addr != kp; rare_nkp = rare_nkp->next);
			    
			    if (rare_nkp)
			    {
				#if DIAG
				if (D3_local)
				fprintf(stdout, "Strangeness:  %s_%i is a neighbor of both %s_%i and %s_%i\n", 
				    per_table[ks].sym, kp->id, per_table[is].sym, ip->id, 
				    per_table[js].sym, jp->id);
				#endif
				rare_njk = rare_nkp->index;
				fjk = fij_[j][rare_njk];
				Fjk1 = Fik1_[j][rare_njk];

				fknl__[2] = acc3__ * fjk * Fjk1 * fik1;

				/* Kernel 2 multiplies Rik-hat; it adds to i's
				 * force and subtracts from k's force. */
				ptPtr_scalmult(&TMP_PT, Rik_hatp, fknl__[2]);
				atom_enforce (ip, kp, &TMP_PT, f_ik__, f_ki__);
			    }
			    
			    /* Second-nearest neighbor forces:  Consider all atoms
			     * bound to neighbor k which are not atom i. */
			    acc3ik__ = acc3__ * fik * Fik1;
			    for (nlp = kp->nList; nlp; nlp = nlp->next)
			    {
				lp = nlp->addr;
				if (lp != ip && lp != jp)
				{
				    nkl = nlp->index;
				    fkl1 = fij1_[k][nkl];
				    Rkl_hatp = &(rij_hat[k][nkl]);
				    f_kl__ = &(frc_ij_[k][nkl]);
				    for (nklp=lp->nList;nklp&&nklp->addr!=kp;nklp=nklp->next);
				    if (nklp) nlk=nklp->index;
				    else {printf("error: no parity, k, l (i)\n");exit(0);}
				    f_lk__ = &(frc_ij_[k][nlk]);
				    fknl__[3] = acc3ik__ * fkl1;
				    
				    ptPtr_scalmult(&TMP_PT, Rkl_hatp, fknl__[3]);
				    atom_enforce (kp, lp, &TMP_PT, f_kl__, f_lk__);
				}
			    }
			}
		    }
		}
	    }
		    
	    a__ = Va_2*beta_ji;
	    if (!allrep) for (nkp = jp->nList; nkp; nkp = nkp->next)
	    {
		kp = nkp->addr;
		if (kp != ip && !chem_isInert(kp->sym))
		{

		    /* Fetch stored precomputations. */

		    k	    = kp->id;
		    ks	    = kp->sym;
		    njk	    = nkp->index;
		    rjk	    = rij_[j][njk];
		    rjk2    = rij2_[j][njk];
		    fjk	    = fij_[j][njk];
		    fjk1    = fij1_[j][njk];
		    Fjk	    = Fik_[j][njk];
		    Fjk1    = Fik1_[j][njk];
		    Njk	    = Nj - fjk;
		    sjik    = sjik_[njk];
		    gjik    = gjik_[njk];
		    gjik1   = gjik1_[njk];
		    omjik   = omjik_[njk];
		    epjik   = epjik_[njk];
		    Rjk_hatp = &(rij_hat[j][njk]);
		    f_jk__  = &(frc_ij_[j][njk]);
		    
		    for (njkp=kp->nList;njkp&&njkp->addr!=jp;njkp=njkp->next);
		    if (njkp) nkj=njkp->index;
		    else {printf("error: no parity, j, k\n");exit(0);}
		    f_kj__  = &(frc_ij_[k][nkj]);
		    
		    #if DIAG
		    if (D3_local > 1)
		    {
			printf("\tk--ji: Triplet %s(%i)(j->%i) - - %s(%i)--%s(%i)\n",
			    per_table[ks].sym, k, njk, 
			    per_table[js].sym, j,
			    per_table[is].sym, i);
			printf("\tDATA RETRIEVAL:\n");
			printf("\tr(%i,%i)(%.5le) r(%i,%i)2(%.5le) f(%i,%i)(%.5le) f'(%i,%i)(%.5le)\n", 
			    j, k, rjk, j, k, rjk2, j, k, fjk, j, k, fjk1);
			printf("\ts(%i,%i,%i)(%.5le) g(%i,%i,%i)(%.5le) g'(%i,%i,%i)(%.5le)\n", 
			    j, i, k, sjik, j, i, k, gjik, j, i, k, gjik1);
			printf("\tom(%i,%i,%i)(%.5le) ep(%i,%i,%i)(%.5le)\n", 
			    j, i, k, omjik, j, i, k, epjik);
			fflush(stdout);
		    }
		    #endif
		    
		    /* Compute force kernel factors. */

		    af1__	= a__*fjk1;
		    af0__	= a__*fjk;
		    g0om__	= gjik*omjik;
		    g1om__   	= gjik1*omjik;
		    af0g1om__	= af0__*g1om__;
		    
		    #if DIAG
		    if (D3_local > 1)
		    {
			printf("\tForce kernel factors, (%i - %i<->%i):\n", k, j, i);
			printf("\ta(%.5le) af0(%.5le), af1(%.5le)\n", a__, af0__, af1__);
			printf("\tg0om(%.5le), g1om(%.5le), af0g1om(%.5le)\n",
			    g0om__, g1om__, af0g1om__);
		    }
		    #endif
		    
		    /* Compute force kernel magnitudes. */

		    fknl__[0]=af1__*(g0om__+Hji1*(ks==F)+Hji2*(ks==C||ks==Si));
		    fknl__[1]=-af0g1om__/rjk;
		    fknl__[3]=af0g1om__/rij;
		    fknl__[5]=-af0__*g0om__*epjik;
		    fknl__[2]=fknl__[3]*sjik;
		    fknl__[4]=fknl__[1]*sjik;
		    fknl__[6]=fknl__[5];

		    #if DIAG
		    if (D3_local > 1)
		    {
			printf("\tForce kernel magnitudes, (%i - %i<->%i):\n", k, j, i);
			for (p=0;p<7;p++) printf("\tfknl[%i] = %.5le\n", p, fknl__[p]);
		    }
		    #endif
		    
		    /* Use the TMP_PT to compute each force
		     * kernel and add it to the proper atom's force
		     * vector. */
		    
		    /* Kernels 0, 4, and 6 add and then multiply
		     * Rjk-hat;  resulting force vector adds to j's
		     * force and subtracts from k's force. */
		    fknl__[0]+=fknl__[4]+fknl__[6];
		    ptPtr_scalmult(&TMP_PT, Rjk_hatp, fknl__[0]);
		    atom_enforce (jp, kp, &TMP_PT, f_jk__, f_kj__);
		    
		    /* Kernel 1 multiplies Rij-hat; it adds to j's force,
		     * and subtracts from k's force. */
		    ptPtr_scalmult(&TMP_PT, Rij_hatp, fknl__[1]);
		    atom_enforce (jp, kp, &TMP_PT, f_jk__, f_kj__);
		    
		    /* Kernels 2 and 5 add, then multiply
		     * Rij-hat; reulting force vecotor
		     * adds to j's force, and subtracts from i's force. */
		    fknl__[2]+=fknl__[5];
		    ptPtr_scalmult(&TMP_PT, Rij_hatp, fknl__[2]);
		    atom_enforce (jp, ip, &TMP_PT, f_ji__, f_ij__);
		    
		    /* Kernel 3 multiplies Rjk-hat; it adds to j's force,
		     * and subtracts from i's force. */
		    ptPtr_scalmult(&TMP_PT, Rjk_hatp, fknl__[3]);
		    atom_enforce (jp, ip, &TMP_PT, f_ji__, f_ij__);
		    		    		    
		    /* Forces if (ij) is a C-C bond:
		     * Reusing force kernel magnitudes 0-3. */
		    if (cc_ij && USE_FCC)
		    {
			fknl__[0] = acc2__ * fjk1;
			
			/* Kernel 0 multiplies Rjk-hat; it adds to j's
			 * force and subtracts from k's force. */
			ptPtr_scalmult(&TMP_PT, Rjk_hatp, fknl__[0]);
			atom_enforce (jp, kp, &TMP_PT, f_jk__, f_kj__);
			
			if (ks == C)  /* k-j-i potentially conjugated */
			{
			    
			    fknl__[1] = acc3__ * fjk1 * Fjk;
			    
			    /* Kernel 1 multiplies Rjk-hat; it adds to j's
			     * force and subtracts from k's force. */
			    ptPtr_scalmult(&TMP_PT, Rjk_hatp, fknl__[1]);
			    atom_enforce (jp, kp, &TMP_PT, f_jk__, f_kj__);

			    /* The next kernel should only be non-zero in the rare event
			     * when k is a neighbor of BOTH j and i. The values of fik
			     * and Fik1 must be fetched from storage, but the ordinate of
			     * k on i's neighbor list must be determined first. */
			    for (rare_nkp = ip->nList; 
				 rare_nkp && rare_nkp->addr != kp; rare_nkp = rare_nkp->next);
			    
			    if (rare_nkp)
			    {
				#if DIAG
				if (D3_local)
				fprintf(stdout, "Strangeness:  %s_%i is a neighbor of both %s_%i and %s_%i\n", 
				    per_table[ks].sym, kp->id, per_table[js].sym, jp->id, 
				    per_table[is].sym, ip->id);
				#endif
				rare_njk = rare_nkp->index;
				fik = fij_[i][rare_njk];
				Fik1 = Fik1_[i][rare_njk];

				fknl__[2] = acc3__ * fik * Fik1 * fjk1;

				/* Kernel 2 multiplies Rjk-hat; it adds to j's
				 * force and subtracts from k's force. */
				ptPtr_scalmult(&TMP_PT, Rjk_hatp, fknl__[2]);
				atom_enforce (jp, kp, &TMP_PT, f_jk__, f_kj__);
			    }
			    
			    /* Second-nearest neighbor forces:  Consider all atoms
			     * bound to neighbor k which are not atom j. */
			    acc3ik__ = acc3__ * fjk * Fjk1;
			    for (nlp = kp->nList; nlp; nlp = nlp->next)
			    {
				lp = nlp->addr;
				if (lp != jp && lp != ip)
				{
				    nkl = nlp->index;
				    fkl1 = fij1_[k][nkl];
				    Rkl_hatp = &(rij_hat[k][nkl]);
				    f_kl__ = &(frc_ij_[k][nkl]);
				    for (nklp=lp->nList;nklp&&nklp->addr!=kp;nklp=nklp->next);
				    if (nklp) nlk=nklp->index;
				    else {printf("error: no parity, k, l (j)\n");exit(0);}
				    f_lk__ = &(frc_ij_[k][nlk]);
				    fknl__[3] = acc3ik__ * fkl1;
				    
				    ptPtr_scalmult(&TMP_PT, Rkl_hatp, fknl__[3]);
				    atom_enforce (kp, lp, &TMP_PT, f_kl__, f_lk__);
				}
			    }
			}
		    }	    /* end if cc_ij */
		}	    /* end if kp != ip and kp is inert */
	    }		    /* end for k neighbor of j */

	    #endif /* if THREEBODY */
	}		    /* end for pair i,j */
    }			    /* end for atom i */
    
    /* Since clusters so far have only been formed from pairs of bound
     * atoms, we need to make sure each isolated single atom gets
     * put in its own cluster. */
    clust_sweepsingles();

} /* end tbt_sicf::mainForce */

static double g_hs_, g_hs2_, g_dhs2_, g_dhs22_;
void g_costh (double s, sym_type i, double * g0, double * g1)
{
    *g0 = 1.0;
    *g1 = 0.0;
    
    #if DIAG
    if (SIC_DIAG2_)
    {
	printf("\tg(costh): i=%s, s=%.5le\n", per_table[i].sym, s);
	printf("\th(%.5le) a(%.5le) c2(%.5le) d2(%.5le) p1(%.5le) p2(%0.5le)\n", 
	    g_h, g_a, g_c2, g_d2, g_p1_, g_p2_);
    }
    #endif
    
    g_hs_ = (g_h-s);
    g_hs2_ = g_hs_*g_hs_;
    if (i == C)
    {
	g_dhs2_ = g_d2 + g_hs2_;
	g_dhs22_ = g_dhs2_*g_dhs2_;
	/* defined by 3body parameter selector:
	 * g_p1_ = g_a + g_p2/g_d2
	 * g_p2_ = g_a*g_c2 */
	*g0 = g_p1_ - g_p2_/(g_dhs2_);
	*g1 = -2*g_p2_*g_hs_/g_dhs22_;
    }
    else if (i == Si)
    {
	*g0 = g_c + g_d*g_hs2_;
	*g1 = -2*g_d*g_hs_;
    }
    else if (i == F)
    {
	*g0 = g_c;
	*g1 = 0.0;
    }
    #if DIAG
    if (SIC_DIAG2_)
    {
	printf("\tg0 = %.5le, g1 = %0.5le\n", *g0, *g1);
	printf("\texiting g(costh)\n");
    }
    #endif
}

static const double ns___ = .5625000000;   /* 9/16 */
static const double os___ = .0625000000;   /* 1/16 */
static double maR__ = 0.0;
void Fij_ (double rij, double * fij, double * f1ij)
{
    *fij = *f1ij = 0.0;
    if (rij <= rijCu1_) *fij = 1.0;
    else if (PI_rijCu12_ && rij < rijCu2_)
    {
	#if BRENNER_FIJ
	/* Brenner/Tersoff smooth cutoff: */
	*fij = 0.5*(1.0+cos(PI_rijCu12_*(rij-rijCu1_)));
	*f1ij = -0.5*PI_rijCu12_*sin(PI_rijCu12_*(rij-rijCu1_));
	#else
	/* Murty/Atwater smoother smooth cutoff: */
	maR__ = PI_rijCu12_*(rij - 0.5*(rijCu1_+rijCu2_));
	*fij = 0.5 - ns___*sin(maR__) - os___*sin(3*maR__);
	*f1ij = -ns___*PI_rijCu12_*cos(maR__)-os___*3*PI_rijCu12_*cos(3*maR__);
	#endif
	if (*fij<0.0) *fij = *f1ij = 0.0;
	if (*fij>1.0) 
	{
	    *fij = 1.0;
	    *f1ij = 0.0;
	}
    }
}

void Ep_ (double r, double ra, double rb, double * ep, double * ep1)
{
    maR__ = M_PI/(rb - ra)*(r - 0.5*(ra + rb));
    *ep = 0.5 - ns___*sin(maR__) - os___*sin(3*maR__);
    *ep1 = -ns___*M_PI/(rb - ra)*cos(maR__)-os___*3*M_PI/(rb - ra)*cos(3*maR__);
}

static double rdiff, rdiff1, rdiff2;
void omega_ (double rij, double rik, double * om_, double * ep_)
{
    *om_ = 1.0;
    *ep_ = 0.0;
    if (alpha && beta)
    {
	rdiff = (rij-rij_e) - (rik-rik_e);
	
	#if DIAG
	if (D3_local > 2)
	{
	    printf("\t\t\tomega: rdiff = %.5e\n", rdiff);
	}
	#endif
	
	if (beta == 1.0) 
	{
	    rdiff1 = rdiff;
	    *ep_ = alpha;
	}
	else
	{
	    rdiff1 = pow(rdiff, beta);
	    *ep_ = alpha*beta*pow(rdiff, beta-1.0);
	}
	if (rdiff1 > -50.0) *om_ = alpha ? exp(alpha*rdiff1) : 1.0;
	else *om_ = 0.0;
    }

}

/* Fconj:  conjugation smooth cutoff function F0 and its first derivative
 * F1 are computed given coordination n. */

void Fconj (double n, double * F0, double * F1)
{
    if (n <= 2.0)
    {
	*F0 = 1.0;
	*F1 = 0.0;
    }
    else if (n < 3.0)
    {
	*F0 = (1.0+cos(M_PI*(n-2.0)))/2.0;
	*F1 = -M_PI*sin(M_PI*(n-2.0))/2.0;
    }
    else
    {
	*F0 = *F1 = 0.0;
    }
}

short no_spline_offset = 0;
static const double m_c[3] = {0.35, 0.55, 0.1};
static const double m_d[3] = {0.3, 1.2, 6.2};
static const double m_cd[3] = {0.105, 0.66, 0.62};
static const double m_cdd[3] = {0.0315, 0.792, 3.844};
static double rij_am, m_efac, m_v1_s1, m_v1_s2, m_v1_s3, mc_efac;
void Moliere (double rij, double * Vr0, double * Vr1, double * Vr2)
{
    int i = 0;
    
    *Vr0 = *Vr1 = m_v1_s1 = m_v1_s2 = m_v1_s3 =  0.0;
    rij_am = rij/M_a;	/* dimensionless rij for the Moliere function */
    
    for (i = 0; i < 3; i++)
    {
	m_efac = exp(-m_d[i]*rij_am);
	mc_efac = m_c[i]*m_efac;
	*Vr0 += mc_efac;
	m_v1_s1 += mc_efac;
	m_v1_s2 += m_cd[i]*m_efac;
	if (Vr2) m_v1_s3 += m_cdd[i]*m_efac;
    }
    *Vr0 *= M_C/rij_am;
    if (!no_spline_offset) *Vr0 += Ms_s;
    m_v1_s1 /= rij;
    m_v1_s2 /= M_a;
    *Vr1 = m_v1_s1 + m_v1_s2;
    *Vr1 *= -M_C/rij_am;
    if (Vr2)
    {
	m_v1_s1 *= 2/rij;
	m_v1_s2 *= 2/rij;
	m_v1_s3 /= (M_a*M_a);
	*Vr2 = m_v1_s1 + m_v1_s2 + m_v1_s3;
	*Vr2 *= M_C/rij_am;
    }

    #if DIAG
    if (D2_local > 2)
    {
	printf("\t\tMOLIERE: Vr(%.5e) = %.5le, dVr_drij = %.5le\n", 
	    rij, *Vr0, *Vr1);
	fflush(stdout);
    }
    #endif  
}

void MolSpline (double rij, double * Vr0, double * Vr1)
{
    *Vr0 = *Vr1 = 0.0;
    mc_efac = exp(Ms_a*rij + Ms_b);
    *Vr0 = mc_efac + (no_spline_offset ? 0 : Ms_c);
    *Vr1 = Ms_a*mc_efac;
    #if DIAG
    if (D2_local > 2)
    {
	printf("\t\t\tVr: MOLSPLINE: Vr(%.5e) = %.5e, dVr = %.5e\n", 
	    rij, *Vr0, *Vr1);
	fflush(stdout);
    }
    #endif
    
}

void MorseVa (double rij, double fij, double fij1, double * Va0, 
	      double * Va1, double * Va2)
{
    *Va0 = *Va1 = 0.0;
    if (fij > 0.0)
    {
	*Va0 = fij*Bij*exp(-muij*rij);
	*Va1 = (fij ? (*Va0)*(fij1/fij - muij) : 0.0);
	if (Va2) *Va2 = *Va0*muij;  /* not using cutoff */
    }
}

void MorseVr (double rij, double fij, double fij1, double * Vr0,
	      double * Vr1, double * Vr2)
{
    *Vr0 = *Vr1 = 0.0;
    if (fij > 0.0)
    {
	*Vr0 = fij*Aij*exp(-lamij*rij);
	*Vr1 = (fij ? (*Vr0)*(fij1/fij - lamij) : 0.0);    
	if (Vr2) *Vr2 = *Vr0*lamij;  /* not using cutoff */
    }
}

short USE_RCU=1;
#define BIG_VALUE1  5
#define BIG_VALUE2  5.3
#define BIG_VALUE1_2 25
#define BIG_PI_12 10.47197551
void TwoBody_SetParameters (sym_type is, sym_type js)
{
    Aij = Bij = lamij = muij = 
    Ms_ra = Ms_rb = Ms_a = Ms_b = Ms_c = Ms_s =
    M_a = M_C =
    rijCu1_ = rijCu2_ = rijCu2sq_ = PI_rijCu12_ = 0.0;
    moliere_only = 0;
    if ((is == Si || is == Si_0) && (js == Si) || (js == Si_0))
    {
	Aij = A_SiSi;
	Bij = B_SiSi;
	lamij = LAMBDA_SiSi;
	muij = MU_SiSi;
	Ms_ra = MS_RASiSi;
	Ms_rb = MS_RBSiSi;
	Ms_a = MS_ASiSi;
	Ms_b = MS_BSiSi;
	Ms_c = MS_CSiSi;
	Ms_s = MS_SSiSi;
	M_a = M_ASiSi;
	M_C = M_CSiSi;
	rijCu1_ = USE_RCU?RSISI_CU1:BIG_VALUE1;
	rijCu2_ = USE_RCU?RSISI_CU2:BIG_VALUE2;
	rijCu2sq_ = USE_RCU?RSISI_CU2SQ:BIG_VALUE1_2;
	PI_rijCu12_ = USE_RCU?PI_RSISI_CU12:BIG_PI_12;
    }
    else if (is == C && js == C)
    {
	Aij = A_CC;
	Bij = B_CC;
	lamij = LAMBDA_CC;
	muij = MU_CC;
	Ms_ra = MS_RACC;
	Ms_rb = MS_RBCC;
	Ms_a = MS_ACC;
	Ms_b = MS_BCC;
	Ms_c = MS_CCC;
	Ms_s = MS_SCC;
	M_a = M_ACC;
	M_C = M_CCC;
	rijCu1_ = USE_RCU?RCC_CU1:BIG_VALUE1;
	rijCu2_ = USE_RCU?RCC_CU2:BIG_VALUE2;
	rijCu2sq_ = USE_RCU?RCC_CU2SQ:BIG_VALUE1_2;
	PI_rijCu12_ = USE_RCU?PI_RCC_CU12:BIG_PI_12;
    }
    else if (is == F && js == F)
    {
	Ms_ra = MS_RAFF;
	Ms_rb = MS_RBFF;
	M_a = M_AFF;
	M_C = M_CFF;
	if (!paramset_a_)
	{
	    Aij = A_FF;
	    Bij = B_FF;
	    lamij = LAMBDA_FF;
	    muij = MU_FF;
	    Ms_a = MS_AFF;
	    Ms_b = MS_BFF;
	    Ms_c = MS_CFF;
	    Ms_s = MS_SFF;
	}
	else
	{
	    Aij = A_FF__A;
	    Bij = B_FF__A;
	    lamij = LAMBDA_FF__A;
	    muij = MU_FF__A;
	    Ms_a = MS_AFF__A;
	    Ms_b = MS_BFF__A;
	    Ms_c = MS_CFF__A;
	    Ms_s = MS_SFF__A;
	}
	rijCu1_ = USE_RCU?RFF_CU1:BIG_VALUE1;
	rijCu2_ = USE_RCU?RFF_CU2:BIG_VALUE2;
	rijCu2sq_ = USE_RCU?RFF_CU2SQ:BIG_VALUE1_2;
	PI_rijCu12_ = USE_RCU?PI_RFF_CU12:BIG_PI_12;
    }
    else if (((is == Si || is == Si_0) && js == C) 
	    || (is == C && (js == Si || js == Si_0)))
    {
	Aij = A_SiC;
	Bij = B_SiC;
	lamij = LAMBDA_SiC;
	muij = MU_SiC;
	Ms_ra = MS_RASiC;
	Ms_rb = MS_RBSiC;
	Ms_a = MS_ASiC;
	Ms_b = MS_BSiC;
	Ms_c = MS_CSiC;
	Ms_s = MS_SSiC;
	M_a = M_ASiC;
	M_C = M_CSiC;
	rijCu1_ = USE_RCU?RSIC_CU1:BIG_VALUE1;
	rijCu2_ = USE_RCU?RSIC_CU2:BIG_VALUE2;
	rijCu2sq_ = USE_RCU?RSIC_CU2SQ:BIG_VALUE1_2;
	PI_rijCu12_ = USE_RCU?PI_RSIC_CU12:BIG_PI_12;
    }
    else if ((is == F && js == C) || (is == C && js == F))
    {
	Ms_ra = MS_RACF;
	Ms_rb = MS_RBCF;
	M_a = M_ACF;
	M_C = M_CCF;
	if (!paramset_a_)
	{
	    Aij = A_CF;
	    Bij = B_CF;
	    lamij = LAMBDA_CF;
	    muij = MU_CF;
	    Ms_a = MS_ACF;
	    Ms_b = MS_BCF;
	    Ms_c = MS_CCF;
	    Ms_s = MS_SCF;
	}
	else
	{
	    Aij = A_CF__A;
	    Bij = B_CF__A;
	    lamij = LAMBDA_CF__A;
	    muij = MU_CF__A;
	    Ms_a = MS_ACF__A;
	    Ms_b = MS_BCF__A;
	    Ms_c = MS_CCF__A;
	    Ms_s = MS_SCF__A;
	}
	rijCu1_ = USE_RCU?RCF_CU1:BIG_VALUE1;
	rijCu2_ = USE_RCU?RCF_CU2:BIG_VALUE2;
	rijCu2sq_ = USE_RCU?RCF_CU2SQ:BIG_VALUE1_2;
	PI_rijCu12_ = USE_RCU?PI_RCF_CU12:BIG_PI_12;
    }
    else if (((is == Si || is == Si_0) && js == F) 
	    || (is == F && (js == Si || js == Si_0)))
    {
	Aij = A_SiF;
	Bij = B_SiF;
	lamij = LAMBDA_SiF;
	muij = MU_SiF;
	Ms_ra = MS_RASiF;
	Ms_rb = MS_RBSiF;
	Ms_a = MS_ASiF;
	Ms_b = MS_BSiF;
	Ms_c = MS_CSiF;
	Ms_s = MS_SSiF;
	M_a = M_ASiF;
	M_C = M_CSiF;
	rijCu1_ = USE_RCU?RSIF_CU1:BIG_VALUE1;
	rijCu2_ = USE_RCU?RSIF_CU2:BIG_VALUE2;
	rijCu2sq_ = USE_RCU?RSIF_CU2SQ:BIG_VALUE1_2;
	PI_rijCu12_ = USE_RCU?PI_RSIF_CU12:BIG_PI_12;
    }
    else if (((is == Si || is == Si_0) && js == Ar) 
	    || (is == Ar && (js == Si || js == Si_0)))
    {
	moliere_only = 1;
	Ms_s = MS_SSiAr;
	M_a = M_ASiAr;
	M_C = M_CSiAr;
	Ms_ra = MS_RASiAr;
	Ms_rb = MS_RBSiAr;
	rijCu2sq_ = RSIAr_CU2SQ;
    }
    else if ((is == C && js == Ar) 
	    || (is == Ar && js == C))
    {
	moliere_only = 1;
	Ms_s = MS_SCAr;
	M_a = M_ACAr;
	M_C = M_CCAr;
	Ms_ra = MS_RACAr;
	Ms_rb = MS_RBCAr;
	rijCu2sq_ = RCAr_CU2SQ;
    }
    else if ((is == F && js == Ar) 
	    || (is == Ar && js == F))
    {
	moliere_only = 1;
	Ms_s = MS_SFAr;
	M_a = M_AFAr;
	M_C = M_CFAr;
	Ms_ra = MS_RAFAr;
	Ms_rb = MS_RBFAr;
	rijCu2sq_ = RFAr_CU2SQ;
    }
    else 
	Aij = Bij = lamij = muij = 
	Ms_ra = Ms_rb = Ms_a = Ms_b = Ms_c = Ms_s = M_a = M_C =
	rijCu1_ = rijCu2_ = rijCu2sq_ = PI_rijCu12_ = 0.0;
    
}

/* ThreeBody_SetParameters:  applies the appropriate selection criteria
 * to choose parameters for the b_ij functions based on the
 * elements in the triplet (i, j, k).  It assigns these values
 * to the global parameter variables declared above. 
 */
 
void ThreeBody_SetParameters (sym_type i, sym_type j, sym_type k)
{

    #if DIAG
    if (SIC_DIAG2_ > 1)
    {
	printf("\tThreeBody_SetParameters: (i%s, j%s, k%s)\n", 
	    per_table[i].sym, per_table[j].sym, per_table[k].sym);
	   
    }
    #endif
    
    g_a = g_c = g_c2 = g_ac2 = g_d = g_d2 = g_c2d2 = g_h
        = g_p1_ = g_p2_ = eta = delta = alpha = beta
	= rij_e = rik_e = rjk_e = 0.0;
    
    if (i == C)
    {
	if (j == C && (k == C || k == Si || k == F))
	{
	    g_a	    = a_Cbrenner;
	    g_c	    = c_Cbrenner;
	    g_c2    = c2_Cbrenner;
	    g_ac2   = ac2_Cbrenner;
	    g_d	    = d_Cbrenner;
	    g_d2    = d2_Cbrenner;
	    g_c2d2  = c2d2_Cbrenner;
	    g_ac2d2 = ac2d2_Cbrenner;
	    g_h	    = h_Cbrenner;
	    g_p1_   = p1_Cbrenner;
	    g_p2_   = p2_Cbrenner;
	    eta	    = eta_Cbrenner;
	    delta   = delta_Cbrenner;
	    alpha   = 0.0;
	    beta    = 0.0;
	    rij_e   = RCC_E;
	    if (k == C) rjk_e = rik_e = RCC_E;
	    else if (k == Si) rjk_e = rik_e = RSIC_E;
	    else if (k == F) rjk_e = rik_e = paramset_a_?RCF_E__A:RCF_E;
	    #if DIAG
	    if (SIC_DIAG2_ > 1)
	    {
		printf("\t(%s, %s, %s) set.\n", 
		    per_table[i].sym, per_table[j].sym, per_table[k].sym);
	    }
	    #endif
	}
	else if (j == Si && (k == C || k == Si || k == F))
	{
	    g_a	    = a_Ctersoff;
	    g_c	    = c_Ctersoff;
	    g_c2    = c2_Ctersoff;
	    g_ac2   = c2_Ctersoff;
	    g_d	    = d_Ctersoff;
	    g_d2    = d2_Ctersoff;
	    g_c2d2  = c2d2_Ctersoff;
	    g_ac2d2 = ac2d2_Ctersoff;
	    g_h	    = h_Ctersoff;
	    g_p1_   = p1_Ctersoff;
	    g_p2_   = p2_Ctersoff;
	    eta	    = eta_Ctersoff;
	    delta   = delta_Ctersoff;
	    alpha   = 0.0;
	    beta    = 0.0;
	    rij_e   = RSIC_E;
	    if (k == C) 
	    {
		rik_e = RCC_E;
		rjk_e = RSIC_E;
	    }
	    else if (k == Si)
	    {
		rik_e = RSIC_E;
		rjk_e = RSISI_E;
	    }
	    else if (k == F)
	    {
		rik_e = paramset_a_?RCF_E__A:RCF_E;
		rjk_e = RSIF_E;
	    }
	    #if DIAG
	    if (SIC_DIAG2_ > 1)
	    {
		printf("\t(%s, %s, %s) set.\n", 
		    per_table[i].sym, per_table[j].sym, per_table[k].sym);
	    }
	    #endif
	}
	else if (j == F && (k == C || k == Si))
	{
	    g_a	    = a_Cbrenner;
	    g_c	    = c_Cbrenner;
	    g_c2    = c2_Cbrenner;
	    g_ac2   = ac2_Cbrenner;
	    g_d	    = d_Cbrenner;
	    g_d2    = d2_Cbrenner;
	    g_c2d2  = c2d2_Cbrenner;
	    g_ac2d2 = ac2d2_Cbrenner;
	    g_h	    = h_Cbrenner;
	    g_p1_   = p1_Cbrenner;
	    g_p2_   = p2_Cbrenner;
	    eta	    = eta_Cbrenner;
	    delta   = delta_Cbrenner;
	    alpha   = 0.0;
	    beta    = 0.0;
	    rij_e   = paramset_a_?RCF_E__A:RCF_E;
	    if (k == C) {rik_e = RCC_E; rjk_e = paramset_a_?RCF_E__A:RCF_E;}
	    else if (k == Si)
	    {
		rik_e = RSIC_E;
		rjk_e = RSIF_E;
	    }
	    #if DIAG
	    if (SIC_DIAG2_ > 1)
	    {
		printf("\t(%s, %s, %s) set.\n", 
		    per_table[i].sym, per_table[j].sym, per_table[k].sym);
	    }
	    #endif
	}
	else if (j == F && k == F)
	{
	    g_a	    = a_Cbrenner;
	    g_c	    = c_Cbrenner;
	    g_c2    = c2_Cbrenner;
	    g_ac2   = ac2_Cbrenner;
	    g_d	    = d_Cbrenner;
	    g_d2    = d2_Cbrenner;
	    g_c2d2  = c2d2_Cbrenner;
	    g_ac2d2 = ac2d2_Cbrenner;
	    g_h	    = h_Cbrenner;
	    g_p1_   = p1_Cbrenner;
	    g_p2_   = p2_Cbrenner;
	    eta	    = eta_Cbrenner;
	    delta   = delta_Cbrenner;
	    alpha   = alpha_Cbrenner;
	    if (paramset_a_) alpha = alpha_Cbrenner__A;
	    beta    = beta_Cbrenner;
	    rij_e   = rik_e = paramset_a_?RCF_E__A:RCF_E;
	    rjk_e   = paramset_a_?RFF_E__A:RFF_E;
	    #if DIAG
	    if (SIC_DIAG2_ > 1)
	    {
		printf("\t(%s, %s, %s) set.\n", 
		    per_table[i].sym, per_table[j].sym, per_table[k].sym);
	    }
	    #endif
	}
	else
	{
	    fprintf(stdout, "FATAL ERROR:  Triplet (%s %s %s) not recognized by 3body setter.\n", 
		per_table[i].sym, per_table[j].sym, per_table[k].sym);
	    fprintf(stderr, "FATAL ERROR:  Triplet (%s %s %s) not recognized by 3body setter.\n", 
		per_table[i].sym, per_table[j].sym, per_table[k].sym);
	    exit(0);
	}
    }
    else if (i == Si)
    {
	if (j == C && (k == C || k == Si))
	{
	    g_c	    = c_Sitersoff;
	    g_d	    = d_Sitersoff;
	    g_h	    = h_Sitersoff;
	    eta	    = eta_Sitersoff;
	    delta   = delta_Sitersoff;
	    alpha   = 0.0;
	    beta    = 0.0;
	    rij_e   = RSIC_E;
	    if (k == C)
	    {
		rik_e = RSIC_E;
		rjk_e = RCC_E;
	    }
	    else if (k == Si)
	    {
		rik_e = RSISI_E;
		rjk_e = RSIC_E;
	    }
	    #if DIAG
	    if (SIC_DIAG2_ > 1)
	    {
		printf("\t(%s, %s, %s) set.\n", 
		    per_table[i].sym, per_table[j].sym, per_table[k].sym);
	    }
	    #endif
	}
	else if (j == C && k == F)
	{
	    eta	    = eta_Sitersoff;
	    delta   = delta_Sitersoff;
	    g_c	    = c_Simurty;
	    g_d	    = d_Simurty;
	    g_h	    = h_Simurty;
	    alpha   = alpha_Simurty;
	    beta    = beta_Simurty;
	    #ifdef TERSOFF_GCOS_SI
	    g_c	    = c_Sitersoff;
	    g_d	    = d_Sitersoff;
	    g_h	    = h_Sitersoff;
	    alpha   = alpha_Ftanaka;
	    if (paramset_a_) alpha = alpha_Ftanaka__A;
	    beta    = beta_Ftanaka;
	    #endif
	    
	    rij_e   = RSIC_E;
	    rik_e   = RSIF_E;
	    rjk_e   = paramset_a_?RCF_E__A:RCF_E;
	    #if DIAG
	    if (SIC_DIAG2_ > 1)
	    {
		printf("\t(%s, %s, %s) set.\n", 
		    per_table[i].sym, per_table[j].sym, per_table[k].sym);
	    }
	    #endif
	}
	else if (j == Si && k == C)
	{
	    g_c	    = c_Sitersoff;
	    g_d	    = d_Sitersoff;
	    g_h	    = h_Sitersoff;
	    eta	    = eta_Sitersoff;
	    delta   = delta_Sitersoff;
	    alpha   = 0.0;
	    beta    = 0.0;
	    rij_e   = RSISI_E;
	    rik_e   = RSIC_E;
	    rjk_e   = RSIC_E;	    
	    #if DIAG
	    if (SIC_DIAG2_ > 1)
	    {
		printf("\t(%s, %s, %s) set.\n", 
		    per_table[i].sym, per_table[j].sym, per_table[k].sym);
	    }
	    #endif
	}
	else if (j == Si && k == Si)
	{
	    g_c	    = c_Sitersoff;
	    g_d	    = d_Sitersoff;
	    g_h	    = h_Sitersoff;
	    eta	    = eta_Sitersoff;
	    delta   = delta_Sitersoff;
	    alpha   = alpha_Sitersoff;
	    beta    = beta_Sitersoff;
	    rij_e   = rik_e = rjk_e = RSISI_E;
	    #if DIAG
	    if (SIC_DIAG2_ > 1)
	    {
		printf("\t(%s, %s, %s) set.\n", 
		    per_table[i].sym, per_table[j].sym, per_table[k].sym);
	    }
	    #endif
	}
	else if (j == Si && k == F)
	{
	    g_c	    = c_Simurty;
	    g_d	    = d_Simurty;
	    g_h	    = h_Simurty;
	    alpha   = alpha_Simurty;
	    beta    = beta_Simurty;	    
	    #ifdef TERSOFF_GCOS_SI
	    g_c	    = c_Sitersoff;
	    g_d	    = d_Sitersoff;
	    g_h	    = h_Sitersoff;
	    alpha   = alpha_Ftanaka;
	    if (paramset_a_) alpha = alpha_Ftanaka__A;
	    beta    = beta_Ftanaka;
	    #endif

	    eta	    = eta_Sitersoff;
	    delta   = delta_Sitersoff;
	    rij_e   = RSISI_E;
	    rik_e   = RSIF_E;
	    rjk_e   = RSIF_E;	    
	    #if DIAG
	    if (SIC_DIAG2_ > 1)
	    {
		printf("\t(%s, %s, %s) set.\n", 
		    per_table[i].sym, per_table[j].sym, per_table[k].sym);
	    }
	    #endif
	}
	else if (j == F && (k == C || k == Si || k == F))
	{
	    g_c	    = c_Simurty;
	    g_d	    = d_Simurty;
	    g_h	    = h_Simurty;
	    eta	    = eta_Simurty;
	    delta   = delta_Simurty;
	    alpha   = alpha_Simurty;
	    beta    = beta_Simurty;
	    #ifdef TERSOFF_GCOS_SI
	    g_c	    = c_Sitersoff;
	    g_d	    = d_Sitersoff;
	    g_h	    = h_Sitersoff;
	    eta	    = eta_Ftanaka;
	    delta   = delta_Ftanaka;
	    alpha   = alpha_Ftanaka;
	    #ifdef alpha_Ftanaka__A
	    if (paramset_a_) alpha = alpha_Ftanaka__A;
	    #endif
	    beta    = beta_Ftanaka;
	    #endif
	    
	    rij_e   = RSIF_E;
	    rik_e   = (k == C ? RSIC_E :
		      (k == Si  ? RSISI_E : RSIF_E));
	    rjk_e   = (k == C ? paramset_a_?RCF_E__A:RCF_E :
		      (k == Si  ? RSIF_E : paramset_a_?RFF_E__A:RFF_E));
	    #if DIAG
	    if (SIC_DIAG2_ > 1)
	    {
		printf("\t(%s, %s, %s) set.\n", 
		    per_table[i].sym, per_table[j].sym, per_table[k].sym);
	    }
	    #endif
	}
	else
	{
	    fprintf(stdout, "FATAL ERROR:  Triplet (%s %s %s) not recognized by 3body setter.\n", 
		per_table[i].sym, per_table[j].sym, per_table[k].sym);
	    fprintf(stderr, "FATAL ERROR:  Triplet (%s %s %s) not recognized by 3body setter.\n", 
		per_table[i].sym, per_table[j].sym, per_table[k].sym);
	    exit(0);
	}
    }
    else if (i == F)
    {
	if (j == C && (k == C || k == Si))
	{
	    eta	    = eta_Ftanaka;
	    delta   = delta_Ftanaka;
	    g_c	    = c_Ftanaka;
	    if (paramset_a_) g_c = c_Ftanaka__A;
	    g_d	    = d_Ftanaka;
	    alpha   = alpha_Ftanaka;
	    if (paramset_a_) alpha = alpha_Ftanaka__A;
	    beta    = beta_Ftanaka;
	    rij_e   = paramset_a_?RCF_E__A:RCF_E;
	    rik_e   = (k == C ? paramset_a_?RCF_E__A:RCF_E : RSIF_E);
	    rjk_e   = (k == C ? RCC_E : RSIC_E);
	    #if DIAG
	    if (SIC_DIAG2_ > 1)
	    {
		printf("\t(%s, %s, %s) set.\n", 
		    per_table[i].sym, per_table[j].sym, per_table[k].sym);
	    }
	    #endif
	}
	else if (j == C && k == F)
	{
	    eta	    = eta_Ftanaka;
	    delta   = delta_Ftanaka;
	    g_c	    = c_Ftanaka;
	    if (paramset_a_) g_c = c_Ftanaka__A;
	    g_d	    = d_Ftanaka;
	    alpha   = alpha_Ftanaka;
	    if (paramset_a_) alpha = alpha_Ftanaka__A;
	    beta    = beta_Ftanaka;
	    rij_e   = paramset_a_?RCF_E__A:RCF_E;
	    rik_e   = paramset_a_?RFF_E__A:RFF_E;
	    rjk_e   = paramset_a_?RCF_E__A:RCF_E;
	    #if DIAG
	    if (SIC_DIAG2_ > 1)
	    {
		printf("\t(%s, %s, %s) set.\n", 
		    per_table[i].sym, per_table[j].sym, per_table[k].sym);
	    }
	    #endif
	}
	else if (j == Si && k == C)
	{
	    eta	    = eta_Ftanaka;
	    delta   = delta_Ftanaka;
	    g_c	    = c_Ftanaka;
	    if (paramset_a_) g_c = c_Ftanaka__A;
	    g_d	    = d_Ftanaka;
	    alpha   = alpha_Ftanaka;
	    #ifdef alpha_Ftanaka__A
	    if (paramset_a_) alpha = alpha_Ftanaka__A;
	    #endif
	    beta    = beta_Ftanaka;
	    rij_e   = RSIF_E;
	    rik_e   = paramset_a_?RCF_E__A:RCF_E;
	    rjk_e   = RSIC_E;
	    #if DIAG
	    if (SIC_DIAG2_ > 1)
	    {
		printf("\t(%s, %s, %s) set.\n", 
		    per_table[i].sym, per_table[j].sym, per_table[k].sym);
	    }
	    #endif
	}
	else if (j == Si && k == Si)
	{
	    eta	    = eta_Ftanaka;
	    delta   = delta_Ftanaka;
	    g_c	    = c_Ftanaka;
	    if (paramset_a_) g_c = c_Ftanaka__A;
	    g_d	    = d_Ftanaka;
	    alpha   = alpha_Ftanaka;
	    #ifdef alpha_Ftanaka__A
	    if (paramset_a_) alpha = alpha_Ftanaka__A;
	    #endif
	    beta    = beta_Ftanaka;
	    rij_e   = rik_e = RSIF_E;
	    rjk_e   = RSISI_E;
	    #if DIAG
	    if (SIC_DIAG2_ > 1)
	    {
		printf("\t(%s, %s, %s) set.\n", 
		    per_table[i].sym, per_table[j].sym, per_table[k].sym);
	    }
	    #endif
	}
	else if (j == Si && k == F)
	{
	    eta	    = eta_Ftanaka;
	    delta   = delta_Ftanaka;
	    g_c	    = c_Ftanaka;
	    if (paramset_a_) g_c = c_Ftanaka__A;
	    g_d	    = d_Ftanaka;
	    alpha   = alpha_Ftanaka;
	    if (paramset_a_) alpha = alpha_Ftanaka__A;
	    beta    = beta_Ftanaka;
	    rij_e   = RSIF_E;
	    rik_e   = paramset_a_?RFF_E__A:RFF_E;
	    rjk_e   = RSIF_E;
	    #if DIAG
	    if (SIC_DIAG2_ > 1)
	    {
		printf("\t(%s, %s, %s) set.\n", 
		    per_table[i].sym, per_table[j].sym, per_table[k].sym);
	    }
	    #endif
	}
	else if (j == F && (k == C || k == Si || k == F))
	{
	    eta	    = eta_Ftanaka;
	    delta   = delta_Ftanaka;
	    g_c	    = c_Ftanaka;
	    if (paramset_a_) g_c = c_Ftanaka__A;
	    g_d	    = d_Ftanaka;
	    alpha   = alpha_Ftanaka;
	    #ifdef alpha_Ftanaka__A
	    if (paramset_a_) alpha = alpha_Ftanaka__A;
	    #endif
	    beta    = beta_Ftanaka;
	    rij_e   = paramset_a_?RFF_E__A:RFF_E;
	    rik_e   = (k == C ? paramset_a_?RCF_E__A:RCF_E :
		      (k == Si  ? RSIF_E : paramset_a_?RFF_E__A:RFF_E));
	    rjk_e   = rik_e;
	    #if DIAG
	    if (SIC_DIAG2_ > 1)
	    {
		printf("\t(%s, %s, %s) set.\n", 
		    per_table[i].sym, per_table[j].sym, per_table[k].sym);
	    }
	    #endif
	}
	else
	{
	    fprintf(stdout, "FATAL ERROR:  Triplet (%s %s %s) not recognized by 3body setter.\n", 
		per_table[i].sym, per_table[j].sym, per_table[k].sym);
	    fprintf(stderr, "FATAL ERROR:  Triplet (%s %s %s) not recognized by 3body setter.\n", 
		per_table[i].sym, per_table[j].sym, per_table[k].sym);
	    exit(0);
	}
    }
    else
    {
	fprintf(stdout, "FATAL ERROR:  Triplet (%s %s %s) not recognized by 3body setter.\n", 
	    per_table[i].sym, per_table[j].sym, per_table[k].sym);
	fprintf(stderr, "FATAL ERROR:  Triplet (%s %s %s) not recognized by 3body setter.\n", 
	    per_table[i].sym, per_table[j].sym, per_table[k].sym);
	exit(0);
    }
}

double atom_peij (const atomPtr ip, const atomPtr jp)
{
    nNodePtr p = NULL;
    if (!ip || !jp) return 0.0;
    for (p = ip->nList; p && p->addr != jp; p = p->next);
    if (!p) return 0.0;
    else return Phi_ij_[ip->id][p->index];
}

double atom_rcuij (sym_type is, sym_type js)
{
    TwoBody_SetParameters(is,js);
    #if DIAG
    if (SIC_DIAG2_)
	printf("AtomPair_CutoffDistance: %s %s\n", 
	    per_table[is].sym, per_table[js].sym);
    #endif
    return (rijCu2_ ? rijCu2_ : sqrt(rijCu2sq_));
}

double atom_reij (sym_type is, sym_type js)
{
    ThreeBody_SetParameters(is,js,Si);
    return (rij_e);
}

short atom_bound (const atomPtr ip, const atomPtr jp)
{
    return (short)(atom_getfij(ip, jp) > 0.0 && atom_peij(ip, jp) < 0.0);
}

double atom_getrij (const atomPtr ip, const atomPtr jp)
{
    nNodePtr p = NULL;
    if (!ip || !jp) return 0.0;
    for (p = ip->nList; p && p->addr != jp; p = p->next);
    if (!p) return 0.0;
    else return rij_[ip->id][p->index];
}

double array_peij (int i, int j)
{
    if (i>MAXNUMATOMS) return 0;
    if (j> MAXNUMNEIGHBORS) return 0;
    return Phi_ij_[i][j];
}

double atom_getfij (const atomPtr ip, const atomPtr jp)
{
    nNodePtr p = NULL;
    if (!ip || !jp) return 0.0;
    for (p = ip->nList; p && p->addr != jp; p = p->next);
    if (!p) return 0.0;
    else return fij_[ip->id][p->index];
}

void cfg_frc2hac (void)
{
    atomPtr p=NULL;
    for (p=L_;p;p=p->next)
    {
	ptPtr_scalmult(p->hac, p->frc, 0.5/per_table[p->sym].mass);
    }
}

static pt tTau={0, 0, 0};
void cfg_frc2tau (void)
{
    atomPtr ip=NULL;
    nNodePtr jp=NULL;
    int ji=0, i=0;
    
    ptPtr_clear(Tau_);
    for (ip=L_;ip;ip=ip->next)
    {
        ptPtr_clear(&tTau);
	for (jp=ip->nList;jp;jp=jp->next)
	{
	    ji=jp->index;
	    ptPtr_scalmult(&tTau, &(rij_hat[i][ji]), rij_[i][ji]);
	    ptPtr_vectmult(&tTau, &(frc_ij_[i][ji]), &tTau);
	    ptPtr_add(Tau_, Tau_, &tTau);
	}
    }
    ptPtr_scalmult(Tau_, Tau_, 0.5);
}
