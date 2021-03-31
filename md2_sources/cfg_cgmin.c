#include "cfg_cgmin.h"
#define ARMAX (3*MAXNUMATOMS)
#define ITMAX 200
#define EPS 1.e-10
/* External globals */
extern atomPtr L_;  /* the configuration */
extern int nAtom_;  /* number of atoms in the configuration */
extern double PE_;  /* total potential energy */
extern short quiet_;
/* Globals owned by this module */
extern double ftol_; /* convergence criterion */
int iter_;	    /* iteration counter for the outer convergence loop */
double * r_=NULL;
double * ng_=NULL;
double r0_[ARMAX];
double xi_[ARMAX];
extern short pe_set_;

static int nCalls=0;
static void header_out (FILE * fp)
{
    if (!fp||quiet_) return;
    fprintf(fp, "#iter\tPE\tTol\tnCalls\n");
    fprintf(fp, "#[#]\t[eV]\n");
}
static void data_out (FILE * fp, int i, double f, double t)
{
    if (!fp) return;
    
    fprintf(fp, "%i\t%.6le\t%.6le\t%i\n", i, f, t, nCalls);
    fflush(fp);
}

static int n=0;	/* n is assigned the number of degrees of freedom
		 * in cfg_cgmin */
double func (double lam)
{
    int i=0;
    if (!r_) {printf("error 1 in cgmin:g\n");exit(0);}
    for (i=0;i<n;i++) r_[i]=r0_[i]+lam*xi_[i];
    pe_set_=0;
    cfg_setforce();
    nCalls++;
    return PE_;
}

void cfg_cgmin (int * iter, double * fret, atomPtr A){
  int i=0, j=0, its=0;
  double gg, gam, fp, dgg, lam=0.0;
  double g[ARMAX], h[ARMAX];
  void linmin (double * lam, double * fret, double (*func)(double));
  nCalls=0;
  if (!A){ /* no single atom specified; minimize by allowing
	     * all atomic positions to vary */
    n=3*nAtom_;
    A=L_;
  }
  else n=3;	/* a single atom was specified; minimize by only
		 * allowing this single atom's position to vary */
  for (i=0;i<n;i++) xi_[i]=0.0;
  
  r_=(double*)&(A->pos->x);
  for (i=0;i<n;i++) r0_[i]=r_[i];
  ng_=(double*)&(A->frc->x);
  fp=func(0.0);
  for (j=0;j<n;j++) {
    xi_[j]=ng_[j];    /* xi is initially the negative gradient direction */
    g[j]=h[j]=xi_[j];
  }
  header_out(stdout);
  data_out(stdout, 0, fp, 0);
  for (its=0;its<ITMAX;its++) {
    *iter=its+1;
    linmin(&lam, fret, func);
    data_out(stdout, *iter, *fret, 2.0*fabs(*fret-fp));
    /* r_[] is minimum in dir. xi_[] */
    /* Convergence criterion: */
    if (2.0*fabs(*fret-fp) <= ftol_*(fabs(*fret)+fabs(fp)+EPS)) return;
    fp=*fret;
    /* Compute the conjugate gradient direction */
    dgg=gg=0.0;
    for (j=0;j<n;j++){
      gg+=g[j]*g[j];
      dgg+=(-ng_[j]+g[j])*-ng_[j];
    }
    if (gg == 0.0) return;
    gam=dgg/gg;
    for (j=0;j<n;j++){
      g[j]=ng_[j];
      xi_[j]=h[j]=g[j]+gam*h[j];
    }
  }
}

void linmin (double * lam, double * fret, double (*func)(double)){
  void mnbrak(double * ax, double * bx, double * cx, 
	      double * fa, double * fb, double * fc, double (*func)(double));
  double brent(double ax, double bx,  double cx,  double (*func)(double), 
	       double tol, double *xmin);
	       
  double xx, ax, bx, fx, fb, fa; int i;
  for (i=0;i<n;i++) r0_[i]=r_[i];
  ax=0.0;
  xx=1.0;
  /* All calls of function g(lam) in mnbrak and brent are done with respect
   * to r0 begin the "current" point in space. */
  mnbrak(&ax, &xx, &bx, &fa, &fx, &fb, func);
  *fret=brent(ax, xx, bx, func, ftol_, lam);
  /* At this point, the final call to g in brent 
   * has set r_[] as the new point in space that is a minimum.  
   * ng_[] are the negative gradients at this point.  If necessary, when
   * we reenter linmin, this point will become the new r0. 
   */
}

#define SHFT(a, b, c, d)    (a)=(b);(b)=(c);(c)=(d);
#define GOLD 1.618034
#define GLIMIT 100.0
#define TINY 1.e-20
static double dmaxarg1, dmaxarg2;
#define DMAX(a, b) (dmaxarg1=(a), dmaxarg2=(b), (dmaxarg1) > (dmaxarg2) ? (dmaxarg1) : (dmaxarg2))
#define SIGN(a, b) ((b) > 0.0 ? fabs(a) : -fabs(a)) 
void mnbrak(double * ax, double * bx, double * cx, 
	        double * fa, double * fb, double * fc, double (*func)(double)){
  double ulim, u, r, q, fu, dum;
  *fa=(*func)(*ax);
  *fb=(*func)(*bx);
  //printf("%f %f\n",*fa,*fb);
  if (*fb>*fa){
    SHFT(dum, *ax, *bx, dum);
    SHFT(dum, *fb, *fa, dum);
  }
  *cx=(*bx)+GOLD*(*bx-*ax);
  *fc=(*func)(*cx);
  while (*fb > *fc){
    r=(*bx-*ax)*(*fb-*fc);
    q=(*bx-*cx)*(*fb-*fa);
    u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/(2.0*SIGN(DMAX(fabs(q-r), TINY), q-r));
    ulim=(*bx)+GLIMIT*(*cx-*bx);
    if ((*bx-u)*(u-*cx) > 0.0){
      fu=(*func)(u);
      if (fu < *fc){
	*ax=(*bx);
	*bx=u;
	*fa=(*fb);
	*fb=fu;
	return;
      }
      else if (fu > *fb){
	*cx=u;
	*fc=fu;
	return;
      }
      u=(*cx)+GOLD*(*cx-*bx);
      fu=(*func)(u);
    }
    else if ((*cx-u)*(u-ulim) > 0.0){
      fu=(*func)(u);
      if (fu < *fc){
	SHFT(*bx, *cx, u, *cx+GOLD*(*cx-*bx));
	SHFT(*fb, *fc, fu, (*func)(u));
      }
    }
    else if ((u-ulim)*(ulim-*cx) >= 0.0){
      u=ulim;
      fu=(*func)(u);
    }
    else{
      u=(*cx)+GOLD*(*cx-*bx);
      fu=(*func)(fu);
    }
    SHFT(*ax, *bx, *cx, u);
    SHFT(*fa, *fb, *fc, fu);
  }
}

#define CGOLD 0.3819660
double brent(double ax, double bx,  double cx,  double (*func)(double), 
	       double tol, double *xmin)
{
    int iter;
    double a, b, d, etemp, fu, fv, fw, fx, p, q, r, tol1, tol2, u, v, w, x, xm;
    double e=0.0;
    
    a=(ax<cx?ax:cx);
    b=(ax>cx?ax:cx);

    x=w=v=bx;
    fw=fv=fx=(*func)(x);
    for (iter=0;iter<ITMAX;iter++)
    {
	xm=0.5*(a+b);
	tol2=2.0*(tol1=tol*fabs(x)+EPS);
	if (fabs(x-xm)<=(tol2-0.5*(b-a)))
	{
	    *xmin=x;
	    return fx;
	}
	if (fabs(e) > tol1)
	{
	    r=(x-w)*(fx-fv);
	    q=(x-v)*(fx-fw);
	    p=(x-v)*q-(x-w)*r;
	    q=2.0*(q-r);
	    if (q>0.0) p=-p;
	    q=fabs(q);
	    etemp=e;
	    e=d;
	    if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
		d=CGOLD*(e=(x >= xm ? a-x : b-x));
	    else
	    {
		d=p/q;
		u=x+d;
		if (u-a < tol2 || b-u < tol2)
		    d=SIGN(tol1, xm-x);
	    }
	}
	else
	{
	    d=CGOLD*(e=(x >= xm ? a-x : b-x));
	}
	u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1, d));
	fu=(*func)(u);
	if (fu<=fx)
	{
	    if (u >= x) a=x; else b=x;
	    SHFT(v, w, x, u);
	    SHFT(fv, fw, fx, fu);
	}
	else
	{
	    if (u < x) a=u; else b=u;
	    if (fu <= fw || w == x)
	    {
		v=w; w=u; fv=fw; fw=fu;
	    }
	    else if (fu <= fv || v == x || v == w)
	    {
		v=u; fv=fu;
	    }
	}
    }
    printf("Warning: too many iterations in brent\n");
    *xmin=x;
    return fx;
}
