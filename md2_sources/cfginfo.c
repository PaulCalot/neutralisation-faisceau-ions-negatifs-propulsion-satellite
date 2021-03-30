/* cfginfo:  a post-processing program that computes information
 * based on a config file.  This is the series 2 implementation.
 * 
 * (c) 1999 Cameron Abrams
 * University of California, Berkeley
 * Department of Chemical Engineering
 * 
 */

#include "atom.h"
#include "chem.h"
#include "cfg.h"
#include "genforce.h"
#include "tbt_sicf.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

double version_=2.00;
int build_=1;

extern atomPtr L_;
extern int nAtom_;
extern element per_table[];
extern int nElem_[];
extern pt Lr_, half_Lr_;
extern unit_type U_;
extern double T_, KE_;
extern ptPtr Tau_;
extern char cfg_infile_[];
void usage(void)
{
    printf("Usage:\n");
    printf("%% cfginfo [list of cfg files] [keywords/values]\n");
    printf("\tKeywords:\n");
    printf("\t-hlo\thistogram domain lower boundary (0.0)\n");
    printf("\t-hhi\thistogram domain upper boundary (50.0)\n");
    printf("\t-hdel\thistogram domain discretization (1.0)\n");
    exit(0);
}
#define MBINS 100

int main (int argc, char * argv[])
{
    atomPtr a = NULL, b = NULL;
    nNodePtr np = NULL, mp=NULL;
    int i = 0, j = 0, k=0, c=0, m=0;
    FILE * fp = NULL;
    double hist[MAXNUMELEMENTS][MBINS];
    double dens[MAXNUMELEMENTS][MBINS];
    short sptags[MAXNUMELEMENTS], shift=1, nominZ=0;
    double Z_z[MAXNUMELEMENTS][MBINS], thisZ=0.0, thisCZ=0.0;
    double dbdens[MAXNUMELEMENTS][MBINS];
    double sumhist[MAXNUMELEMENTS][MBINS];
    double sumdens[MAXNUMELEMENTS][MBINS];
    double sumZ_z[MAXNUMELEMENTS][MBINS];
    double sumdbdens[MAXNUMELEMENTS][MBINS];
    double bonddens[MAXNUMELEMENTS][MAXNUMELEMENTS][MBINS];
    double sumbonddens[MAXNUMELEMENTS][MAXNUMELEMENTS][MBINS];
    short bondtags[MAXNUMELEMENTS][MAXNUMELEMENTS];
    char * cfg=0L;
    int dangles[5];
    double hlo=0.0, hhi=50.0, hdel=1.0, binvol;
    int hn = 0, bin=0, nbin=0;
    pt TDV={0, 0, 0};
    double rij_l = 0.0, totN=0, totX=0;
    short ff = 0;
    double totdens = 0.0, thisV=0;
    short massdens=0, c_asc=0;
    sym_type nsym, asym;
    int nNbrs=0;
    double tBnd=0;
    int nn, nnf;
    double Qn[6]={0, 0, 0, 0, 0, 0}, Qsum, minZ, rcu_=1.97;
    /* count species: Si, SiF,SiF2,SiF3,SiF4,C,CF,CF2,CF3,CF4 */
    double sicfx[10]={0,  0,  0,   0,   0,   0,0, 0,  0,  0}, sicfxsum;
    
    printf("# cfginfo v 2.0 (c) 1999 Cameron Abrams\n");fflush(stdout);
    
    for(i=0;i<MAXNUMELEMENTS;i++) for(j=0;j<MAXNUMELEMENTS;j++) bondtags[i][j]=0;
    for(i=0;i<MAXNUMELEMENTS;i++) for(j=0;j<MAXNUMELEMENTS;j++) 
	for (k=0;k<MBINS;k++) bonddens[i][j][k]=sumbonddens[i][j][k]=0.0;
    for(i=0;i<MAXNUMELEMENTS;i++) for(j=0;j<MBINS;j++) 
	hist[i][j]=dbdens[i][j]=sumdbdens[i][j]=
	sumdens[i][j]=sumZ_z[i][j]=sumhist[i][j]=dens[i][j]=Z_z[i][j]=0.0;
    for(i=0;i<MAXNUMELEMENTS;i++) sptags[i]=0;
    for(i=0;i<5;i++) dangles[i]=0;
    for(i=0;i<10;i++) sicfx[i]=0.0;
    for(i=0;i<6;i++) Qn[i]=0.0;
    
    if (argc<2) usage();
    for (i=1;i<argc;i++)
    {
	if (!strcmp(argv[i],"-hlo")) hlo=atof(argv[++i]);
	else if (!strcmp(argv[i],"-hhi")) hhi=atof(argv[++i]);
	else if (!strcmp(argv[i],"-hdel")) hdel=atof(argv[++i]);
	else if (!strcmp(argv[i],"-noshift")) shift=0;
	else if (!strcmp(argv[i],"-massdens")) massdens=1;
	else if (!strcmp(argv[i],"-nominz")) nominZ=1;
	else if (!strcmp(argv[i],"-ff")) ff=1;
	else if (!strcmp(argv[i],"-")||isalnum(argv[i][0])
	        ||(argv[i][0]=='.'&&argv[i][1]=='.'))
	    cfg=argv[i];
	else usage();
    }
    hn = (int)((hhi-hlo)/hdel);
    if (hn >= MBINS) {printf("hn too big\n");exit(0);}
    
    Chem_InitializePeriodicTable();
    pef_initialize();
    cfg_initialize();
    
    /* Hard code the bonding tags so these bond densities are always output */
    bondtags[Si][Si]=1;
    bondtags[Si][C]=1;
    bondtags[Si][Csp1]=1;
    bondtags[Si][Csp2]=1;
    bondtags[Si][Csp3]=1;
    bondtags[Si][F]=1;
    bondtags[C][C]=1;
    bondtags[C][Csp1]=1;
    bondtags[C][Csp2]=1;
    bondtags[C][Csp3]=1;
    bondtags[F][C]=1;
    bondtags[Csp1][Csp1]=1;
    bondtags[Csp1][Csp2]=1;
    bondtags[Csp1][Csp3]=1;
    bondtags[F][Csp1]=1;
    bondtags[Csp2][Csp2]=1;
    bondtags[Csp2][Csp3]=1;
    bondtags[F][Csp2]=1;
    bondtags[Csp3][Csp3]=1;
    bondtags[F][Csp3]=1;
    sptags[C]=sptags[Si]=sptags[F]=1;
    
    for(i=0;i<MAXNUMELEMENTS;i++) 
	for(j=0;j<MAXNUMELEMENTS;j++) 
	    if (j<i) bondtags[i][j]=bondtags[j][i];

    for (k = 0;k<MAXNUMELEMENTS;k++) for (j=0;j<MBINS;j++) hist[k][j]=0;
    for (k = 0;k<MAXNUMELEMENTS;k++) for (j=0;j<MBINS;j++) Z_z[k][j]=0.0;
    for (j=0;j<MBINS;j++) for (k=0;k<MAXNUMELEMENTS;k++) 
	for (m=0;m<MAXNUMELEMENTS;m++) bonddens[k][m][j]=0.0;
    strcpy(cfg_infile_, cfg);
    cfg_establish();

    binvol=Lr_.x*Lr_.y*hdel;
    if (nominZ) {printf("# zmin is off\n"); minZ=-99.999;}
    else minZ=0.5*Lr_.z;
	
    /* shift so that lowest atoms are at z=0 */
    if (shift) for (a=L_;a;a=a->next) a->pos->z+=half_Lr_.z;
    printf("# %s: Total %i T=%.2lf", 
	strcmp(cfg,"-")?cfg:"stdin", nAtom_, T_);
    fflush(stdout);
    for (k = 0; k < MAXNUMELEMENTS; k++)
    {
	if (nElem_[k])
	{
	    printf(" %s %i", per_table[k].sym, nElem_[k]);
	}
    }
    printf("\n");

    cfg_setforce();

    ptPtr_scalmult(Tau_, Tau_, GPA_PER_EVPERCUANG/cfg_v()); /* GPa */
    printf("# bin volume = %.5lf cuA, = %.5le cc\n", binvol, binvol*1.e-24);
    printf("# %s: V=%.5le cuA, Pxy(0)=%.5le GPa\n", 
	 strcmp(cfg,"-")?cfg:"stdin", cfg_v(),
	 0.5*(Tau_->x+Tau_->y));
    fflush(stdout);
    for (nn=0;nn<6;nn++) Qn[nn]=0;
    
    for (a=L_;a;a=a->next)
    {
	bin=(int)((a->pos->z-hlo)/hdel);
	thisZ=0.0;
	nNbrs=0;
	for (np=a->nList;np;np=np->next) 
	{
	    thisZ+=(atom_peij(a,np->addr)<0.0?atom_getfij(a,np->addr):0);
	    nNbrs++;
	}
	if (asym!=F&&a->state!=IS_FIXED) dangles[4-(nNbrs<5?nNbrs:4)]++;
	asym=a->sym;
	if (asym==C)
	{
	    nn=0;
	    nnf=0;
	    for (np=a->nList;np;np=np->next)
	    {
		if (atom_getrij(a,np->addr)<=rcu_) 
		{
		    nn++;
		    if (np->addr->sym==F) nnf++;
		}
	    }
	    if (nn<6&&a->pos->z>minZ) Qn[nn]++;
	    if (nnf<5&&nnf>0&&a->pos->z>minZ) sicfx[nnf+5]++; 
	    if (thisZ > 2.5 && thisZ < 3.5) asym=Csp2;
	    else if (thisZ < 2.5) asym=Csp1;
	    else asym=Csp3;
	}
	if (asym==Si)
	{
	    nn=0;
	    nnf=0;
	    for (np=a->nList;np;np=np->next)
	    {
		if (atom_getrij(a,np->addr)<=rcu_) 
		{
		    nn++;
		    if (np->addr->sym==F) nnf++;
		}
	    }
	    if (nnf<5&&nnf>0&&a->pos->z>minZ) sicfx[nnf]++; 
	}
	nNbrs=0;
	for (np=a->nList;np;np=np->next)
	{
	    if (np->addr->id>a->id)
	    {
		nbin=(int)((np->addr->pos->z-hlo)/hdel);
		nsym=np->addr->sym;
		if (nsym==C)
		{
		    thisCZ=0.0;
		    for (mp=np->addr->nList;mp;mp=mp->next)
			thisCZ+=(atom_peij(a,np->addr)<0.0?
				 atom_getfij(a,np->addr):0);

		    if (thisCZ > 2.5 && thisCZ < 3.5) nsym=Csp2;
		    else if (thisCZ < 2.5) nsym=Csp1;
		    else nsym=Csp3;
		}
		tBnd=(atom_getfij(a,np->addr)==1.0&&atom_peij(a,np->addr)<0.0?0.5:0);
		bonddens[asym][nsym][bin]+=tBnd;
		if (nsym!=asym) bonddens[nsym][asym][bin]+=tBnd;
		bonddens[asym][nsym][nbin]+=tBnd;
		if (nsym!=asym) bonddens[nsym][asym][nbin]+=tBnd;
		c_asc=0;
		if (nsym==Csp1||nsym==Csp2||nsym==Csp3) {nsym=C;c_asc=1;}
		if (asym==Csp1||asym==Csp2||asym==Csp3) {asym=C;c_asc=1;}
		if (c_asc)
		{
		    bonddens[asym][nsym][bin]+=tBnd;
		    if (nsym!=asym) bonddens[nsym][asym][bin]+=tBnd;
		    bonddens[asym][nsym][nbin]+=tBnd;
		    if (nsym!=asym) bonddens[nsym][asym][nbin]+=tBnd;
		}
	    }
	    nNbrs++;
	}
	thisV=per_table[a->sym].v-nNbrs;
	hist[a->sym][bin]++;
	Z_z[a->sym][bin]+=thisZ;
	dbdens[a->sym][bin]+=thisV>0?thisV:0;
    }
    
    Qsum=0.0;
    for (nn=0;nn<6;nn++) Qsum+=Qn[nn];
    printf("# Q:\n#\t");
    for (nn=0;nn<6;nn++) printf("%i\t", nn);
    printf("\n# ");
    for (nn=0;nn<6;nn++) printf("%.5lf\t", Qsum?Qn[nn]/Qsum:0.0);
    printf("\n");

    sicfxsum=0.0;
    for (nn=0;nn<5;nn++) sicfxsum+=sicfx[nn];
    if (sicfxsum)
    {
	printf("# SiFx:\n#\t");
	for (nn=0;nn<5;nn++) printf("%i\t", nn);
	printf("\n# ");
	for (nn=0;nn<5;nn++) printf("%.5lf %.5lf\t", sicfx[nn], sicfxsum?sicfx[nn]/sicfxsum:0.0);
	printf("\n");
    }
    sicfxsum=0.0;
    for (nn=5;nn<10;nn++) sicfxsum+=sicfx[nn];
    if (sicfxsum)
    {
	printf("# CFx:\n#\t");
	for (nn=5;nn<10;nn++) printf("%i\t", nn-5);
	printf("\n# ");
	for (nn=5;nn<10;nn++) printf("%.5lf\t", sicfxsum?sicfx[nn]/sicfxsum:0.0);
	printf("\n");
    }

    for(j=0;j<MBINS;j++) 
    {
	for(k=0;k<MAXNUMELEMENTS;k++) 
	{
	    if (sptags[k])
	    {
		dens[k][j]=hist[k][j]/binvol;
		dbdens[k][j]/=binvol;
		for (m=0;m<MAXNUMELEMENTS;m++) bonddens[k][m][j]/=binvol;
		sumhist[k][j]+=hist[k][j];
		sumdens[k][j]+=dens[k][j];
		sumZ_z[k][j]+=(hist[k][j]?Z_z[k][j]/hist[k][j]:0.0);
		sumdbdens[k][j]+=dbdens[k][j];
		for (m=0;m<MAXNUMELEMENTS;m++) 
		    sumbonddens[k][m][j]+=(k==m?1:2)*bonddens[k][m][j];
	    }
	}
    }
    cfg_clear();
	    

    c=1;
    printf("#%i-depth/A\t", c++);
    for (i=0;i<MAXNUMELEMENTS;i++)
	if (sptags[i]) printf("%i-#_%s\t%i-!7q!6_%s\t%i-x_%s\t%i-Z_%s\t%i-dbd_%s\t",
	    c++, per_table[i].sym, c++, per_table[i].sym, c++, 
	    per_table[i].sym,  c++, per_table[i].sym,  c++, per_table[i].sym);
    printf("%i-!7q!6!itot!n\t%i-sumX\t", c++, c++);
    for (i=0;i<MAXNUMELEMENTS;i++) for(j=i;j<MAXNUMELEMENTS;j++)
	if (bondtags[i][j]) printf("%i-%s-%s\t", 
	    c++, per_table[i].sym, per_table[j].sym);
    printf("\n");
    for (j=0;j<hn;j++)
    {
	totdens = 0.0;
	printf("%.3lf\t", (j+0.5)*hdel);
	totN = 0.0;
	for (i=0;i<MAXNUMELEMENTS;i++) if (sptags[i]) totN += sumhist[i][j];
	if (!totN) totN = 1.0;
	totX = 0.0;
	for (i=0;i<MAXNUMELEMENTS;i++) if (sptags[i]) totX += sumhist[i][j]/totN;
/*	if (totX != 1.0) {printf("error x\n");exit(0);}*/
	for (i=0;i<MAXNUMELEMENTS;i++)
	{
	    if (sptags[i])
	    {	
		if (massdens) sumdens[i][j] *= 1.66*per_table[i].mass/APVKMASS_PER_AMU;
		printf("%.4lf\t%.5lf\t", sumhist[i][j], sumdens[i][j]);
		thisV=per_table[i].v-sumZ_z[i][j];
		printf("%.5lf\t%.5lf\t%.5lf\t", 
		    sumhist[i][j]/totN,sumZ_z[i][j],sumdbdens[i][j] );
		totdens += sumdens[i][j];
	    }
	}
	printf("%.5lf\t%.5lf\t", totdens, totX);
	for (i=0;i<MAXNUMELEMENTS;i++) for(k=i;k<MAXNUMELEMENTS;k++)
	    if (bondtags[i][k]) printf("%.5lf\t", sumbonddens[i][k][j]);
	printf("\n");
    }
    printf("# 1db %i, 2db %i, 3db %i, 4db %i, total %i\n", 
	dangles[1], dangles[2], dangles[3], dangles[4], 
	dangles[1]+dangles[2]+dangles[3]+dangles[4]);
    return(0);
    
}


