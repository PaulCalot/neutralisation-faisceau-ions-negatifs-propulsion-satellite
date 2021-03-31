/*
 * md  series 2 -- Molecular Dynamics of Plasma-Surface Chemistry
 * (c) 1999 Cam Abrams, Dept. of Chemical Engineering
 * Regents of the University of California, Berkeley
 * 
 * Module Name:  chem 
 * 
 * Module Description:  chem defines and initializes the periodic
 * table data structure used by other modules in the series.
 * 
 */

#ifndef CHEM_H
#define CHEM_H

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define MAXNUMELEMENTS      26

#define MASS_HELIUM	4.002
#define MASS_NEON	20.179
#define MASS_ARGON	39.948
#define MASS_KRYPTON	83.80
#define MASS_XENON      131.30
#define MASS_OXYGEN     15.999
#define MASS_SILICON    28.085
#define MASS_FLUORINE   18.998
#define MASS_CHLORINE   35.453
#define MASS_BROMINE	79.904
#define MASS_COPPER	63.546
#define MASS_HYDROGEN	1.008
#define MASS_CARBON	12.011

#define XX		5
#define Si_0    	0
#define Si      	1
#define F       	2
#define H		3
#define C		4
#define O       	7 
#define He      	10
#define Ne      	11
#define Kr      	12
#define Ar      	13
#define Xe      	14
#define Br      	15
#define Cu		16
#define Cl      	22
#define Cl_S      	25
#define Csp1          17
#define Csp2          18
#define Csp3          19

#define Z_Si      	14
#define Z_F       	9
#define Z_O       	8 
#define Z_Ar      	18
#define Z_Br      	35
#define Z_Cl      	17
#define Z_He      	2
#define Z_Ne      	10
#define Z_Kr      	36
#define Z_Xe      	54
#define Z_Cu		29
#define Z_H		1
#define Z_C		6

#define V_Si      	4
#define V_F       	1
#define V_O       	2 
#define V_Ar      	0
#define V_Br      	1
#define V_Cl      	1
#define V_He      	0
#define V_Ne      	0
#define V_Kr      	0
#define V_Xe      	0
#define V_Cu		6
#define V_H		1
#define V_C		4

typedef int sym_type;
typedef enum {APVK, STILLWEB, LENJONES, NULL_UT} unit_type;
unit_type chem_a2ut (char *);
char * chem_ut2a (unit_type);

typedef struct ELEMENT
{   /* record of data on an element */
    char    name[15];       /* name of element */
    double  mass;           /* atomic mass */
    int	    z;		    /* atomic number */
    int	    v;		    /* # valence electrons */
    char    sym[3];         /* atomic symbol */
    double  radius;	    /* radius (for rendering software) */
    struct
    { 
	double r, g, b; 
    } color;		    /* color (for rendering software) */
    struct
    {
	unsigned int is_recyclable  :	1;
    } flags;
} element;	

typedef struct ELEMENTPAIR
{
    sym_type a;
    sym_type b;
} element_pair;

/* Chem_InitializePeriodicTable:  initializes the global element array
 * per_table[]. */
void Chem_InitializePeriodicTable (void);

double Chem_PeriodicTable_MassOfSym (char sym[]);
double Chem_MassOfCompound (int elemCnt[]);

sym_type Chem_PeriodicTable_SymOfSym (char sym[]);
sym_type Chem_PeriodicTable_SymOfZ (int z);

char * Chem_NameCompound (char name[], int elemCnt[]);

void Pair_SwapMembers (element_pair p);

int chem_isInert (sym_type s);
int chem_valence (sym_type s);


#endif

