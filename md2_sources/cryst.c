/* md series 2 (c) 1999 cam abrams
 * 
 * cryst.c -- defines unit cell parameters for several
 * simple crystal structures
 */
 
#include "cryst.h"
#include "atom.h"
#include "cfg.h"

extern ipt nCell_;
extern element per_table[];
int uc_=0;
char * ucNames_[NULL_UC]=
{
    "DC_Si","OH_Si","DC_C","OH_C"
};

int ucOcc_[]=
{
    8, 12, 8, 12
};
void uc_a2i (char * ucStr)
{
    uc_=0;
    if (!ucStr) return;
    for (uc_=0;uc_<NULL_UC&&strcmp(ucStr, ucNames_[uc_]);uc_++);
    if (uc_==NULL_UC) {printf("Error: %s is not a valid unit cell\n", 
	ucStr);exit(0);}
}
int cryst_cnt (void)
{
    if (uc_<0||uc_>=NULL_UC) 
	{printf("Error: uc # %i not valid\n", uc_);exit(0);}
    return nCell_.i*nCell_.j*nCell_.k*ucOcc_[uc_];
}

/* Some important constants */
#define SQRT_2		1.41421356237309504880    
#define SQRT_2_2	0.70710678118654752440    
#define SQRT_3		1.73205080756887729352
#define SQRT_6		2.44948974278317809819
#define SQRT_6_2	1.22474487139158904909
#define THD		0.33333333333333333333
#define TTHD		0.66666666666666666667
#define SXTH		0.16666666666666666667
#define TWTH		0.08333333333333333333
#define FVTWTH		0.41666666666666666667
#define SVTWTH		0.58333333333333333333

/* Converts diamond-cubic (DC) a0 into the orthohexagonal (OH) a,b,c */
#define a0_abc(A0,A,B,C) {(A)=(A0)*SQRT_2_2; (B)=(A0)*SQRT_6_2; (C)=(A0)*SQRT_3;}

/* Note: all lengths are in Angstroms */

ucDesc_t u_[NULL_UC]=
{
/* DC_Si */ {	5.434176017741, 5.434176017741, 5.434176017741,\
		/* The above values for a,b,c give 0.00000 GPa
		 * of stress in a 4x4x4 uc (512 atom) cfg with
		 * (xy) peridocity and 2 layers of fixed atoms */
		{Si,Si,Si,Si,Si,Si,Si,Si,XX,XX,XX,XX,XX,XX,XX,XX}, 
		{{-0.25,0.25,0.25},{0.25,-0.25,0.25},
		 {-0.5,0,0},{0,0.5,0},
		 {-0.25,-0.25,-0.25},{0.25,0.25,-0.25},
		 {-0.5,-0.5,-0.5},{0,0,-0.5},
		 {0,0,0},{0,0,0},{0,0,0},{0,0,0},
		 {0,0,0},{0,0,0},{0,0,0},{0,0,0}}
	    }, 
/* OH_Si */ {	3.842542712306, 6.655479207967, 9.412268960000,
		{Si,Si,Si,Si,Si,Si,Si,Si,Si,Si,Si,Si,XX,XX,XX,XX},
		{{-0.5,-0.5,-0.5},{0,0,-0.5},
		 {0,-THD,-SXTH},{-0.5,SXTH,-SXTH},
		 {-0.5,-SXTH,SXTH},{0,THD,SXTH},
		 {-0.5,-0.5,-0.25},{0,0,-0.25},
		 {0,-THD,TWTH},{-0.5,SXTH,TWTH},
		 {-0.5,-SXTH,FVTWTH},{0,THD,FVTWTH},
		 {0,0,0},{0,0,0},{0,0,0},{0,0,0}}
	    }, 
/* DC_C */ {	3.501352430020,3.501352430020,3.501352430020,
		/* The above values for a,b,c give 0.000000 GPa
		 * of stress in a 6x6x6 uc (1728 atom) cfg with
		 * (xy) peridocity and 360 fixed atoms */
		{C,C,C,C,C,C,C,C,XX,XX,XX,XX,XX,XX,XX,XX}, 
		{{-0.25,0.25,0.25},{0.25,-0.25,0.25},
		 {-0.5,0,0},{0,0.5,0},
		 {-0.25,-0.25,-0.25},{0.25,0.25,-0.25},
		 {-0.5,-0.5,-0.5},{0,0,-0.5},
		 {0,0,0},{0,0,0},{0,0,0},{0,0,0},
		 {0,0,0},{0,0,0},{0,0,0},{0,0,0}}
	    }, 
/* OH_C */ {	2.527336331266, 4.377474933568, 6.190684420000,
		/* The above values for a,b,c give 0.000000 GPa
		 * of stress in a 8x4x2 uc (768 atom) cfg with
		 * (xy) peridocity corresponding
		 * to a CC bond length of 1.5476711 Angstroms */
		{C,C,C,C,C,C,C,C,C,C,C,C,XX,XX,XX,XX}, 
		{{-0.5,-0.5,-0.5},{0,0,-0.5},
		 {0,-THD,-SXTH},{-0.5,SXTH,-SXTH},
		 {-0.5,-SXTH,SXTH},{0,THD,SXTH},
		 {-0.5,-0.5,-0.25},{0,0,-0.25},
		 {0,-THD,TWTH},{-0.5,SXTH,TWTH},
		 {-0.5,-SXTH,FVTWTH},{0,THD,FVTWTH},
		 {0,0,0},{0,0,0},{0,0,0},{0,0,0}}
	    } 
};

double xl_Re_=0.0;	/* atomic spacing -- if user-specified,
			 * crystal lattice parameters are 
			 * computed using it, otherwise they
			 * are taken from default values */
extern short quiet_;
void uc_initialize (void)
{
    double a0;
    if (xl_Re_)
    {
	/* 4*Re=sqrt(3)*a0 for the dc unit cell */
	a0=4*xl_Re_/SQRT_3;
	if (!quiet_) 
	    printf("# Crystal Re=%.10lf => DC a0 = %.10lf\n", xl_Re_, a0);
	if (uc_==OH_SILICON||uc_==OH_CARBON)
	{
	    a0_abc(a0,u_[uc_].a,u_[uc_].b,u_[uc_].c);
	}
	else u_[uc_].a=u_[uc_].b=u_[uc_].c=a0;
	if (!quiet_) 
	    printf("# Crystal %s Lattice Params: %.12lf %.12lf %.12lf\n", 
		ucNames_[uc_], u_[uc_].a, u_[uc_].b, u_[uc_].c);
	fflush(stdout);
    }
}

extern atomPtr L_;  /* The configuration */
extern int nAtom_, nFixed_;
extern pt Lr_, half_Lr_;
double ff_=0.2;
extern ipt per_;
void cryst_pos (void)
{
    int i=0, iz, iy, ix, ic;
    double z, zmin, zmax;
    pt c_l={0, 0, 0};
    pt TMP, *tpt=&TMP;
    c_l.x=u_[uc_].a;
    c_l.y=u_[uc_].b;
    c_l.z=u_[uc_].c;
    Lr_.x = nCell_.i*c_l.x;
    Lr_.y = nCell_.j*c_l.y;
    Lr_.z = nCell_.k*c_l.z;
    ptPtr_scalmult(&half_Lr_,&Lr_,0.50);
    
    printf("# a0 %.12lf, b0 %.12lf, c0 %.12lf\n", c_l.x, c_l.y, c_l.z);
    fflush(stdout);
    
    for (i=0;i<ucOcc_[uc_];i++)
    {
	L_[i].id=i;
	L_[i].next=((i+1)<nAtom_?&(L_[i+1]):NULL);
	L_[i].sym=u_[uc_].ek[i];
	u_[uc_].os[i].z+=0.5;
	if (u_[uc_].os[i].z>=0.5) u_[uc_].os[i].z-=1.0;
	ptPtr_copy(L_[i].pos, ptPtr_vectmult(tpt, &c_l, &(u_[uc_].os[i])));
    }
    for (iz=0;iz<nCell_.k;iz++)
    {
	for (iy=0;iy<nCell_.j;iy++)
	{
	    for (ix=0;ix<nCell_.i;ix++)
	    {
		if (ix||iy||iz) /* Cell 000 was assigned outside loop */
		{
		    for (ic=0;ic<ucOcc_[uc_];ic++)
		    {
			L_[i].id=i;
			L_[i].next=((i+1)<nAtom_?&(L_[i+1]):NULL);
			L_[i].sym=L_[ic].sym;
			L_[i].pos->x=L_[ic].pos->x+c_l.x*ix;
			L_[i].pos->y=L_[ic].pos->y+c_l.y*iy;
			L_[i].pos->z=L_[ic].pos->z+c_l.z*iz;
			i++;
		    }
		}
	    }
	}
    }
    /* shift center of box to origin  -- right now, origin is at center
     * of the first unit cell */
    zmin=1.e9;
    zmax=-1.e9;
    for (i=0;i<nAtom_;i++)
    {
	L_[i].pos->x-=half_Lr_.x-c_l.x/2;
	L_[i].pos->y-=half_Lr_.y-c_l.y/2;
	L_[i].pos->z-=half_Lr_.z-c_l.z/2;
	z=L_[i].pos->z;
	if (z<zmin) zmin=z;
	if (z>zmax) zmax=z;
    }
    nFixed_=0;
    if (!per_.k)
	for (i=0;i<nAtom_;i++)
	{
	    z=L_[i].pos->z;
	    if (z < zmin+ff_*(zmax-zmin)) 
	    {
		L_[i].state=IS_FIXED;
		nFixed_++;
	    }
	}

}

double gaussRand (void)
{
   int i;
   double sum;

   sum = 0;
   for (i=1;i<=12;i++)
       sum+=((double)rand())/((double)RAND_MAX);	

   return (sum-6.0);
}

#include <time.h>
extern double bhb_Tset_;
extern unit_type U_;
void cryst_vel (void)
{
    double summv2=0.0;	/* sum of squares of velocities */
    double setmv2=0.0;	/* calculated average v-squared */
    double sFac=0.0;	/* scaling factor */
    atomPtr p=NULL;
    int i=0;
    double Tcalc=0.0;
    pt V={0, 0, 0}, sV={0, 0, 0};
    double T=bhb_Tset_;

    srand((unsigned)clock()+(unsigned)time(NULL));
    setmv2=3.0*nAtom_*T*(1.0*(U_==STILLWEB)+(U_==APVK)*KB_EVPERK);
    printf("# Desired T %.3lf K\n", T);
    if (T)
    {
	for (p=L_;p;p=p->next)
	{
	    for (i=0;i<3;i++) *(&(p->vel->x)+i)=gaussRand();
	    ptPtr_add(&sV, &sV, p->vel);
	    /* update the running sum of square velocities */
	    summv2+=per_table[p->sym].mass*ptPtr_sqabs(p->vel);
	}
	ptPtr_scalmult(&V, &sV, 1.0/nAtom_);
	
	/* generate scaling factor */
	sFac=sqrt(setmv2/summv2);
	
	sV.x=sV.y=sV.z=summv2=0.0;
    
	/* scale and adjust the assigned velocity components appropriately */
	for (p=L_;p;p=p->next)
	{
	    ptPtr_subtract(p->vel, p->vel, &V);
	    ptPtr_scalmult(p->vel, p->vel, sFac);
	    ptPtr_add(&sV, &sV, p->vel);
	    summv2+=per_table[p->sym].mass*ptPtr_sqabs(p->vel);
	}
	 
	Tcalc=summv2/3.0/nAtom_/(1.0*(U_==STILLWEB)+(U_==APVK)*KB_EVPERK);
	printf("# Tcalc %.3lf K, %.1f%% of %.3lf.\n", Tcalc, Tcalc/T*100, T);
    }
    else
    {
	for (p=L_;p;p=p->next)
	{
	    for (i=0;i<3;i++) *(&(p->vel->x)+i)=0.0;
	}
    }
}

