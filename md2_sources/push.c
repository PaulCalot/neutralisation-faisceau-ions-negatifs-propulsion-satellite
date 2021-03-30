#include "push.h"

/* The integrator */

/* Global variables needed by integrator, initialized to default values.
 * These can be altered by user input in the setup file, and/or command line
 * arguments.  See "args.c". */

double dt_=1.0e-3;	    /* time step */
short dtvar_=1;		    /* turns on/off time step adaptive variation */
double Lc_=0.01;	    /* max position del allowed per time step */
double maxdt_=5.0e-3;	    /* maximum allowed time step value */
double mindt_=1.0e-3;	    /* maximum allowed time step value */
double tlim0_=0.0;	    /* integration time lower limit */
double tlim1_=0.0;	    /* integration time upper limit */
int t0_=0;		    /* integration time step lower limit */
int Ndt_=0;		    /* explicit number of time steps */
int cfg_outint_=0;	    /* number of time steps between in-run cfg
			     * outputs */
int data_int_=1;	    /* number of time steps between integration
			     * data output */
			     
/* Global variables computed by this module:
 * (a) local variables */

static double dt2_=0.0;	    /* square of time step */

/* (b) variables used by other modules: */

double Eo_=0.0;		    /* initial total energy */
double E_=0.0;		    /* total energy */
int	t_=0;		    /* time step index */
double  rt_=0.0;	    /* run time value */
short quiet_=0;		    /* flag to supress message/header output */
short ke_set_=0;	    /* lock */

double bhb_Tset_=0.0;	    /* Berendsen-type heat bath temperature set point */
double bhb_tau_=0.0;	    /* Berendsen-type heat bath rise time, in same time
			     * units as the integration time step */
int bhb_start_=0;	    /* time step at which to begin heat bath */
double bmb_Pset_=0.0;	    /* Berendsen-type momentum (pressure)
			     * bath Pressure set point, in GPa */
double bmb_tau_=0.0;	    /* Berendsen-type momentum (pressure)
			     * bath rise time, in same time
			     * units as the integration time step */
int bmb_start_=0;	    /* time step at which to begin momentum bath */

/* externals: */

extern const atomPtr L_;    /* the configuration */
extern int nAtom_, nFixed_; /* number of atoms and number of fixed atoms */
extern ipt per_;	    /* Periodic boundary condition switches */
extern pt Lr_;		    /* box size */
extern pt half_Lr_;	    /* 1/2 box size */
extern double KE_;	    /* kinetic energy */
extern double sKE_;	    /* portion of kinetic energy that is not due
			     * to the unbound incident ion */
extern double PE_;	    /* potential energy */
extern double T_;	    /* temperature */
extern ptPtr Tau_;	    /* diagonals of stress tensor */
extern unit_type U_;	    /* unit system */
extern element per_table[]; /* periodic table */
extern short pe_set_;	    /* lock */

extern ion_type ion_;

/* functions used only in this module */
static void cfg_push1 (void);	/* Step 1 of Velocity Verlet */
static void cfg_push2 (void);	/* Step 2 of Velocity Verlet */

static void header_out (FILE *);
static void data_out (FILE *, short);
static void data_sum (void);
void cfg_push (void)
{
    FILE * efp=stdout;
    if (!dt2_) dt2_=dt_*dt_;
    
    if (!quiet_) 
    {
	printf("# md series 2 cfg pusher (c) 1999 cfa\n");
	fflush(stdout);
    }
    
    t_=t0_;
    rt_=tlim0_;
    if (!tlim1_) tlim1_=tlim0_+Ndt_;
    cfg_setforce();
    cfg_push2();
    E_=Eo_=KE_+PE_;
    header_out(efp);
    data_sum();
    data_out(efp,1);
    if (cfg_outint_) {cfg_snapout(t_);}
    if (ion_.spec!=-1 && ion_.uint &&!(t_%ion_.uint)) {ion_update();}
    if (ion_.spec!=-1 && ion_.oint &&!(t_%ion_.oint)) {ion_output();}

    while ((!Ndt_&&tlim1_&&rt_<=tlim1_)||(Ndt_&&(t_-t0_)<Ndt_))
    {
	t_++;
	rt_+=dt_;
	cfg_push1();	    /* Step 1 of Velocity Verlet */
	cfg_setforce();	    /* Force routine */
	cfg_push2();	    /* Step 2 of Velocity Verlet */
	E_=KE_+PE_;
	data_sum();
	if (!(t_%data_int_)) data_out(efp,0);
	if (cfg_outint_ &&!(t_%cfg_outint_))	{cfg_snapout(t_);}
	if (ion_.spec!=-1 && ion_.uint &&!(t_%ion_.uint)) {ion_update();}
	if (ion_.spec!=-1 && ion_.oint &&!(t_%ion_.oint)) {ion_output();}
    }
    
    if (cfg_outint_) {cfg_snapout(tlim1_);}
    if (efp!=stdout&&efp!=stderr) fclose(efp);
    if (!quiet_) {printf("# md series 2 cfg pusher completed\n");fflush(stdout);}
}

static void header_out (FILE * fp)
{
    if (!fp||quiet_) return;
    fprintf(fp, "#t\ttime\tKE\t\tPE\t\tE\t\tEdrift\t\t%%Edrift");
    fprintf(fp, "\t\tT\tnClst\tnAicl?\tPlat\n");
    fprintf(fp, "#[#]\t[ps]\t[eV]\t\t[eV]\t\t[eV]\t\t[eV]\t\t[%%]");
    fprintf(fp, "\t\t[K]\t[#]\t[0/1]\t[GPa]\n");
}

static double Ed_=0.0;
static double Edp_=0.0;
static double sumKE_=0.0;
static double sumPE_=0.0;
static double sumE_=0.0;
static double sumEd_=0.0;
static double sumEdp_=0.0;
static double sumT_=0.0;
static double sumTauLat_=0.0;
static void data_sum (void)
{
    Ed_=E_-Eo_;
    Edp_=fabs(Ed_/Eo_)*100.0;
    sumKE_+=KE_;    
    sumPE_+=PE_;    
    sumE_+=E_;    
    sumEd_+=Ed_;    
    sumEdp_+=Edp_;    
    sumT_+=T_;  
    sumTauLat_+=0.5*(per_.i*Tau_->x+per_.j*Tau_->y);
}
extern int nClust_;
extern int nAic_;
static int norm=1;
static void data_out (FILE * fp, short init)
{
    
    if (!fp) return;
    
    if (init) norm=1;
    else norm=data_int_;
    fprintf(fp, "%i\t%.4lf\t%.7le\t%.7le\t%.7le\t%.7le\t%.7le\t%.2lf\t%i\t%i\t%5lf\n", 
	         t_, rt_,   sumKE_/norm,   sumPE_/norm,
		    sumE_/norm,    sumEd_/norm,    
		    sumEdp_/norm,  sumT_/norm, 
	         nClust_, nAic_==nAtom_, sumTauLat_/norm);
    fflush(fp);
    sumKE_=sumPE_=sumE_=sumEd_=sumEdp_=sumT_=sumTauLat_=0.0;
}

static atomPtr ip = NULL;
static pt v = {0, 0, 0}, *vp=&v;
static pt h = {0, 0, 0}, *hp=&h;
static void cfg_push1 (void)
{
    for (ip=L_;ip;ip=ip->next)
    {
	if (ip->state != IS_FIXED)
	{
	    ptPtr_scalmult(vp, ip->vel, dt_);
	    ptPtr_scalmult(hp, ip->hac, dt2_);
	    ptPtr_add(ip->pos, ip->pos, vp);
	    ptPtr_add(ip->pos, ip->pos, hp);
	    ptPtr_minimg(ip->pos, ip->pos);
	    ptPtr_scalmult(hp, ip->hac, dt_);
	    ptPtr_add(ip->vel, ip->vel, hp);
	    
	}
    }
    /* Atomic positions have changed -> new cfg -> current pe info no longer
     * valid -> set pe_set_ lock variable to 0. */
    pe_set_=0;
    /* Atomic velocities have changed -> current ke no longer valid ->
     * set ke_set_ lock variable to 0. */
    ke_set_=0;
}

static double mass;
static double s;
static double smax=-1.e9;
static double vsk_fac, lsk_fac, Pcalc;
static int n;
static short sgn;
static void cfg_push2 (void)
{
    if (ke_set_) {printf("# KE already set for current cfg.\n");return;}
    KE_=sKE_=0.0;
    n=0;
    for (ip=L_;ip;ip=ip->next)
    {
	if (ip->state!=IS_FIXED)
	{
	    ptPtr_scalmult(hp, ip->hac, dt_);
	    ptPtr_add(ip->vel, ip->vel, hp);

	    s=ptPtr_sqabs(ip->vel);
            mass=per_table[ip->sym].mass;
	    ip->mv2=mass*s;
	    KE_+=ip->mv2;
	    if (ip->state!=IS_PROJECTILE&&ip->state!=IS_UNCLEARED)
	    {
		n++;
		sKE_+=ip->mv2;
	    }
	    smax=(s>smax?s:smax);
	    
	    /* Increment the normal stresses vector with the Ideal Gas Pressure */
	    if (ip->state!=IS_PROJECTILE&&ip->state!=IS_UNCLEARED)
	    {
		ptPtr_vectmult(vp, ip->vel, ip->vel);
		ptPtr_scalmult(vp, vp, 2*per_table[ip->sym].mass);
		ptPtr_add(Tau_, Tau_, vp);
	    }
	}
    }
    
    /* Convert the normal stresses vector from eV/cu.Ang to GPa */
    if (U_!=APVK) fprintf(stderr, "Warning: Tau in wrong units\n");
    ptPtr_scalmult(Tau_, Tau_, GPA_PER_EVPERCUANG/cfg_v());
    /* Compute the stress in the film.  per_ is the vector
     * of periodicities: per_.x is 1 if dir x is periodic.
     * An aperiodic (free) boundary cannot contribute to
     * the isotropic stress. */
    Pcalc=(per_.i||per_.j||per_.k)?
	  (per_.i*Tau_->x+per_.j*Tau_->y+per_.k*Tau_->z)/(per_.i+per_.j+per_.k):
	  0.0;

    KE_*=0.5;
    sKE_*=0.5;
    if (n)
    {
	if (U_==APVK) T_=sKE_/1.5/n/KB_EVPERK;
	else T_=sKE_/1.5/n/KB;
    }
    else
    {
	T_=KE_/1.5/nAtom_/KB_EVPERK;
    }
    if (dtvar_)
    {
	s=Lc_/sqrt(smax);
	dt_=(s<maxdt_)? s : maxdt_;
	dt_=(s>mindt_)? dt_ : mindt_;
	dt2_=dt_*dt_;
    }
    /* implement the berendsen heat bath */
    if (t_>=bhb_start_&&bhb_Tset_&&bhb_tau_)
    {
	vsk_fac = 1. + dt_/bhb_tau_*(bhb_Tset_/T_-1.0);
	vsk_fac = sqrt(vsk_fac);
	for (ip=L_;ip;ip=ip->next)
	    if (ip->state!=IS_FIXED&&ip->state!=IS_PROJECTILE
		&&ip->state!=IS_UNCLEARED) 
		ptPtr_scalmult(ip->vel, ip->vel, vsk_fac);
    }
    /* implement the momentum (pressure) bath */
    if (t_>=bmb_start_&&bmb_tau_)
    {
	lsk_fac = 1. - dt_/bmb_tau_*(bmb_Pset_-Pcalc);
	sgn=lsk_fac<0?1:0;
	lsk_fac = (sgn?-1:1)*pow(fabs(lsk_fac), 0.3333333333);
	vp->x=per_.i?lsk_fac:1.0;
	vp->y=per_.j?lsk_fac:1.0;
	vp->z=per_.k?lsk_fac:1.0;
	ptPtr_vectmult(&Lr_, &Lr_, vp);
	ptPtr_scalmult(&half_Lr_, &Lr_, 0.5);
	for (ip=L_;ip;ip=ip->next)
	    ptPtr_vectmult(ip->pos, ip->pos, vp);
    }
    
    ke_set_=1;
}
