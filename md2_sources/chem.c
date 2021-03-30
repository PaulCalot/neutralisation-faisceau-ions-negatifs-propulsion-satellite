/*
 * MD series -- Molecular Dynamics of Plasma-Surface Chemistry
 * (c) 1996-1998 Cam Abrams, Dept. of Chemical Engineering
 * Regents of the University of California, Berkeley
 * 
 * Module Name:		chem
 * Module Class:	Base
 * Module Description:  chem defines and initializes the periodic
 * table data structure used by other modules in the series.
 * 
 */

#include <string.h>
#include "chem.h"
#include "units.h"

element per_table[MAXNUMELEMENTS];
short per_table_initialized_=0;
unit_type U_=APVK;
extern short quiet_;

char * UnitTypeLabels[NULL_UT]=
{
    "APVK", "Stillinger-Weber", "Lennard-Jones"
};

unit_type chem_a2ut (char * us)
{
    unit_type u=APVK;
    if (!strcmp(us, "SW") ||
	!strcmp(us, "StillWeb") ||
	!strcmp(us, "Stillinger-Weber")) u=STILLWEB;
    else if (!strcmp(us, "LJ") ||
	!strcmp(us, "LenJones") ||
	!strcmp(us, "Lennard-Jones")) u=LENJONES;
    return u;
}
char * chem_ut2a (unit_type u)
{
    return UnitTypeLabels[u];
}
#if DIAG
int CHEM_DIAG_;
#endif
static char * cheminp="chem.inp";
static FILE * fp=NULL;
#ifndef MAXCHAR
#define MAXCHAR 255
#endif
static char ln[MAXCHAR];
static char dumsym[3];
static double dumR, dumr, dumg, dumb;

void Chem_InitializePeriodicTable (void)
{
    int i;
    #if DIAG
    char * d_h = "chem::Chem_InitializePeriodicTable";
    #endif
    
    if (!quiet_) printf("# md series 2 periodic table initializer (c) 1999 cfa\n");
    
    for (i = 0; i < MAXNUMELEMENTS; i++)
    {
	per_table[i].sym[0] = '\0';
	per_table[i].mass = 0.00;
	per_table[i].name[0] = '\0';
	per_table[i].color.r = 0.0;
	per_table[i].color.g = 0.0;
	per_table[i].color.b = 0.0;
	per_table[i].radius = 0.0;
	per_table[i].z = 0;
	per_table[i].flags.is_recyclable = 0;
    }

    per_table[Si].flags.is_recyclable = 1;
    per_table[O].flags.is_recyclable  = 1;
    per_table[Cu].flags.is_recyclable  = 1;
    
    per_table[Si].mass	    =	MASS_SILICON;
    per_table[Si_0].mass    =	MASS_SILICON;
    per_table[O].mass	    =	MASS_OXYGEN;
    per_table[He].mass	    =	MASS_HELIUM;
    per_table[Ne].mass	    =	MASS_NEON;
    per_table[Ar].mass	    =	MASS_ARGON;
    per_table[Kr].mass	    =	MASS_KRYPTON;
    per_table[Xe].mass	    =	MASS_XENON;
    per_table[F].mass	    =	MASS_FLUORINE;
    per_table[Cl].mass	    =	MASS_CHLORINE;
    per_table[Cl_S].mass    =	MASS_CHLORINE;
    per_table[Br].mass	    =	MASS_BROMINE;
    per_table[Cu].mass	    =	MASS_COPPER;
    per_table[H].mass	    =	MASS_HYDROGEN;
    per_table[C].mass	    =	MASS_CARBON;

    per_table[Si].z	    =	Z_Si;
    per_table[Si_0].z	    =	Z_Si;
    per_table[O].z	    =	Z_O;
    per_table[He].z	    =	Z_He;
    per_table[Ne].z	    =	Z_Ne;
    per_table[Ar].z	    =	Z_Ar;
    per_table[Kr].z	    =	Z_Kr;
    per_table[Xe].z	    =	Z_Xe;
    per_table[F].z	    =	Z_F;
    per_table[Cl].z	    =	Z_Cl;
    per_table[Cl_S].z	    =	Z_Cl;
    per_table[Br].z	    =	Z_Br;
    per_table[Cu].z	    =	Z_Cu;
    per_table[H].z	    =	Z_H;
    per_table[C].z	    =	Z_C;
        
    per_table[Si].v         =   V_Si;
    per_table[Si_0].v       =   V_Si;
    per_table[O].v          =   V_O;
    per_table[He].v         =   V_He;
    per_table[Ne].v         =   V_Ne;
    per_table[Ar].v         =   V_Ar;
    per_table[Kr].v         =   V_Kr;
    per_table[Xe].v         =   V_Xe;
    per_table[F].v          =   V_F;
    per_table[Cl].v         =   V_Cl;
    per_table[Cl_S].v       =   V_Cl;
    per_table[Br].v         =   V_Br;
    per_table[Cu].v         =   V_Cu;
    per_table[H].v          =   V_H;
    per_table[C].v          =   V_C;
        
    strcpy(per_table[Si].name,	    "Silicon");
    strcpy(per_table[Si_0].name,    "Silicon");
    strcpy(per_table[O].name,	    "Oxygen");
    strcpy(per_table[Ar].name,	    "Argon");
    strcpy(per_table[He].name,	    "Helium");
    strcpy(per_table[Xe].name,	    "Xenon");
    strcpy(per_table[Ne].name,	    "Neon");
    strcpy(per_table[Kr].name,	    "Krypton");
    strcpy(per_table[F].name,	    "Fluorine");
    strcpy(per_table[Cl].name,	    "Chlorine");
    strcpy(per_table[Cl_S].name,    "Chlorine");
    strcpy(per_table[Br].name,	    "Bromine");
    strcpy(per_table[Cu].name,	    "Copper");
    strcpy(per_table[H].name,	    "Hydrogen");
    strcpy(per_table[C].name,	    "Carbon");
        
    strcpy(per_table[XX].sym,	    "XX");
    strcpy(per_table[Si_0].sym,     "S0");
    strcpy(per_table[Si].sym,	    "Si");
    strcpy(per_table[O].sym,	    "O");
    strcpy(per_table[Ar].sym,	    "Ar");
    strcpy(per_table[He].sym,	    "He");
    strcpy(per_table[Xe].sym,	    "Xe");
    strcpy(per_table[Ne].sym,	    "Ne");
    strcpy(per_table[Kr].sym,	    "Kr");
    strcpy(per_table[F].sym,	    "F");
    strcpy(per_table[Cl].sym,	    "Cl");
    strcpy(per_table[Cl_S].sym,     "Cl");
    strcpy(per_table[Br].sym,	    "Br");
    strcpy(per_table[Cu].sym,	    "Cu");
    strcpy(per_table[H].sym,	    "H");
    strcpy(per_table[C].sym,	    "C");
    strcpy(per_table[Csp1].sym,	    "C1");
    strcpy(per_table[Csp2].sym,	    "C2");
    strcpy(per_table[Csp3].sym,	    "C3");

    /* Radii below are in SW sigmas (don't ask me why).
     * They are converted to angstroms if required below. 
     * These are the default values.  If the file "chem.inp"
     * exists in the directory in which this is run, it 
     * is read immediately after these defaults are set
     * and any values it contains supercede these defaults. */

    per_table[Si].color.r = 0.0;
    per_table[Si].color.g = 0.15;
    per_table[Si].color.b = 0.8;
    per_table[Si].radius = 0.6;

    per_table[H].color.r = 0.9;
    per_table[H].color.g = 0.9;
    per_table[H].color.b = 0.9;
    per_table[H].radius = 0.2;

    per_table[C].color.r = 0.0;
    per_table[C].color.g = 1.0;
    per_table[C].color.b = 0.1;
    per_table[C].radius = 0.35;

    per_table[Si_0].color.r = 0.0;
    per_table[Si_0].color.g = 0.3;
    per_table[Si_0].color.b = 0.7;
    per_table[Si_0].radius = 0.5;

    per_table[Ar].color.r = 0.9;
    per_table[Ar].color.g = 0.9;
    per_table[Ar].color.b = 0.0;
    per_table[Ar].radius = 0.6;

    per_table[Ne].color.r = 1.00;
    per_table[Ne].color.g = .431;
    per_table[Ne].color.b = .706;
    per_table[Ne].radius = 0.3;

    per_table[O].color.r = 0.8;
    per_table[O].color.g = 0.7;
    per_table[O].color.b = 0.6;
    per_table[O].radius = 0.55;

    per_table[F].color.r = 1.0;
    per_table[F].color.g = 0.01;
    per_table[F].color.b = 0.01;
    per_table[F].radius = 0.4;

    per_table[Cl].color.r = 0.1;
    per_table[Cl].color.g = 0.95;
    per_table[Cl].color.b = 0.1;
    per_table[Cl].radius = 0.75;

    per_table[Cl_S].color.r = 0.1;
    per_table[Cl_S].color.g = 0.95;
    per_table[Cl_S].color.b = 0.1;
    per_table[Cl_S].radius = 0.75;

    per_table[Br].color.r = 0.1;
    per_table[Br].color.g = 0.95;
    per_table[Br].color.b = 0.1;
    per_table[Br].radius = 0.75;

    per_table[Cu].color.r = 0.5;
    per_table[Cu].color.g = 0.2;
    per_table[Cu].color.b = 0.2;
    per_table[Cu].radius = 0.61;

    fp=fopen(cheminp, "r");
    if (fp)
    {
	while (fgets(ln, MAXCHAR, fp))
	{
	    if (ln[0]!='#') 
	    {
		sscanf(ln, "%s %lf %lf %lf %lf\n", 
		    dumsym, &dumR, &dumr, &dumg, &dumb);
		for (i=0;strcmp(dumsym,per_table[i].sym)&&i<MAXNUMELEMENTS;i++);
		if (i<MAXNUMELEMENTS)
		{
		    per_table[i].color.r=dumr;
		    per_table[i].color.g=dumg;
		    per_table[i].color.b=dumb;
		    per_table[i].radius=dumR;
		}
	    }
	}
	fclose(fp);
    }
    
    /* Implement the unit system conventions */
    switch(U_)
    {
	case APVK:	for (i = 0; i < MAXNUMELEMENTS; i++)
			{
			    per_table[i].radius *= ANGSTROMS_PER_SIGMA;
			    per_table[i].mass *= APVKMASS_PER_AMU;
			}
			break;
	case STILLWEB:	for (i = 0; i < MAXNUMELEMENTS; i++)
			    per_table[i].mass /= MASS_SILICON;
			per_table[Si].mass = 1.0;
			per_table[Si_0].mass = 1.0;
			break;
	case LENJONES:	for (i = 0; i < MAXNUMELEMENTS; i++)
			    per_table[i].mass /= MASS_ARGON;
			per_table[Ar].mass = 1.0;
			break;
	default:	for (i = 0; i < MAXNUMELEMENTS; i++)
			{
			    per_table[i].radius *= ANGSTROMS_PER_SIGMA;
			    per_table[i].mass *= APVKMASS_PER_AMU;
			}
			break;
    }
    per_table_initialized_=1;
}

double Chem_PeriodicTable_MassOfSym (char sym[])
{
    int i = 0;
    double mass = 0.0;
    if (!per_table_initialized_) {printf("MOS: PT not initialized.\n");exit(0);}
    for (i = 0; !mass && i < MAXNUMELEMENTS; i++)
    {
	if (!strcmp(sym, per_table[i].sym)) mass = per_table[i].mass;
    }
    
    return mass;
}

sym_type Chem_PeriodicTable_SymOfSym (char sym[])
{
    int i=0;
    #if DIAG
    char * d_h="chem::Chem_PeriodicTable_SymOfSym";
    #endif
    sym_type s=0;
    if (!per_table_initialized_) {printf("SOS: PT not initialized.\n");exit(0);}
    #if DIAG
    if (CHEM_DIAG_)
	printf("%s searching for [%s] in per_table...", d_h, sym);
    #endif
    for (i=0;!s&&i<MAXNUMELEMENTS;i++)
	if (!strcmp(sym,per_table[i].sym)) s=i;
    #if DIAG
    if (CHEM_DIAG_)
	printf("returning %i[%s]\n", s, per_table[s].sym);
    #endif
    return s;
}

sym_type Chem_PeriodicTable_SymOfZ (int z)
{
    int i = 0;
    #if DIAG
    char * d_h = "chem::Chem_PeriodicTable_SymOfSym";
    #endif
    sym_type s = 0;
    #if DIAG
    if (CHEM_DIAG_)
	printf("%s searching for z=[%i] in per_table...", d_h, z);
    #endif
    if (!per_table_initialized_) {printf("SOZ: PT not initialized.\n");exit(0);}
    for (i = 0; !s && i < MAXNUMELEMENTS; i++)
    {
	if (z == per_table[i].z) s = i;
    }
    
    return s;
}

double Chem_MassOfCompound (int elemCnt[])
{
    int i = 0;
    double mass = 0.0;
    if (!per_table_initialized_) {printf("MOC: PT not initialized.\n");exit(0);}
    for (i = 0; i < MAXNUMELEMENTS; i++)
    {
	if (elemCnt[i]) mass += per_table[i].mass;
    }
    return mass;
}

char * Chem_NameCompound (char name[], int elemCnt[])
{
    int i = 0;
    char tmpStr[12];
    char tmpStr2[12];
    if (!per_table_initialized_) {printf("NC: PT not initialized.\n");exit(0);}
    name[0] = '\0';
    for (i = 0; i < MAXNUMELEMENTS; i++)
    {
	if (elemCnt[i])
	{
	    sprintf(tmpStr, "%s", per_table[i].sym);
	    if (elemCnt[i] > 1) 
	    {
		sprintf(tmpStr2, "!I%i!N", elemCnt[i]);
		strcat(tmpStr, tmpStr2);
	    }
	    strcat(name, tmpStr);
	}
    }
    return name;
}

void Pair_SwapMembers (element_pair p)
{
    sym_type tmp;
    
    tmp = p.a;
    p.a = p.b;
    p.b = tmp;
}

int chem_isInert (sym_type s)
{
    return (s==He||s==Ne||s==Ar||s==Kr||s==Xe);
}

int chem_valence (sym_type s)
{
    if (s>=MAXNUMELEMENTS||s<0) return -1;
    return (per_table[s].v);
}
