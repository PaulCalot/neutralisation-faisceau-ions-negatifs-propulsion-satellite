/* md series 2 */
#include "ion.h"
#define SEL_AB(A,B)   (paramset_a_?A:B)
#define s60    8.66025403784438600000e-01   /* sin(60) */
#define s19_45 3.32984122349684280000e-01   /* sin(19.45) */
#define s70_55 9.42932433561923020000e-01   /* sin(70.55) */

extern atomPtr L_;
extern double Zhi_;
extern sym_type symZhi_;
extern element per_table[];
extern unit_type U_;	    /* unit system denoter */
extern pt Lr_, half_Lr_;    /* boxsize,1/2 boxsize */
extern double KE_;
extern short quiet_;
extern short paramset_a_;
typedef struct sym_dvec
{
    sym_type sym;
    pt r;
} sym_pt_type;
typedef struct ion_int_coords
{
    sym_type base;
    int n;
    const sym_pt_type * x;
} iic_type;

/* Ion type ID's:				    |  ion classes: */
enum {I_He, I_C, I_F, I_Ne, I_Si, I_Ar, I_Kr, I_Xe, /* atoms */
      I_CC, I_FF, I_SiSi,			    /* homodimers AA */
      I_CF, I_SiF, I_SiC,			    /* heterodimers AX */
      I_CF2, I_SiF2, I_SiC2, I_Si2C,		    /* AX2 */
      I_CF3, I_SiF3, I_SiC3, I_Si3C,		    /* AX3 */
      I_CF4, I_SiF4, I_SiC4, I_Si4C,		    /* AX4 */
      NULL_ION};

char * IonIDstr[NULL_ION] =
{
    "He", "C", "F", "Ne", "Si", "Ar", "Kr", "Xe", 
    "CC", "FF", "SiSi", 
    "CF", "SiF", "SiC", 
    "CF2", "SiF2", "SiC2", "Si2C", 
    "CF3", "SiF3", "SiC3", "Si3C", 
    "CF4", "SiF4", "SiC4", "Si4C" 
};

char * IonIDstr_alt[NULL_ION] =
{
    "HE", "C", "F", "NE", "SI", "AR", "KR", "XE", 
    "CC", "FF", "SISI", 
    "CF", "SIF", "SIC", 
    "CF2", "SIF2", "SIC2", "SI2C", 
    "CF3", "SIF3", "SIC3", "SI3C", 
    "CF4", "SIF4", "SIC4", "SI4C" 
};

/* Internal atomic coordinates (ic's) for all ion types.  
 * Coordinates are in reference to the "base" atom 
 * in the molecule, which resides at (0, 0, 0).
 * Coordinates and atoms are represented by the
 * sym_pt_type data type, which is a structure
 * containing the symbol and a pt.
 * Molecules are represented by the iic_type, 
 * in which the base atom is denoted by the sym_type
 * 'base', and the number of *other* atoms in the
 * molecule is 'n', and the symbols and coordinates
 * of these atoms are in an array pointed to by 'x'.
 * These arrays are declared as 'const sym_pt_type XX_coord[]'
 * where XX is the molecule name. */

/* Dimers are vertically oriented by default */
const sym_pt_type CC_coord[] = {{C,0,0,RCC_E}};
const sym_pt_type FF_coord[] = {{F,0,0,RFF_E}};
const sym_pt_type SiSi_coord[] = {{Si,0,0,RSISI_E}};
const sym_pt_type SiF_coord[] = {{F,0,0,RSIF_E}};
const sym_pt_type CF_coord[] = {{F,0,0,RCF_E}};
const sym_pt_type SiC_coord[] = {{C,0,0,RSIC_E}};

/* AX2's are linear in the x direction */
const sym_pt_type CF2_coord[]={{F, -RCF_E, 0,0},{F, RCF_E, 0,0}}; /* CF2 */
const sym_pt_type SiF2_coord[]={{F, -RSIF_E,0,0},{F, RSIF_E,0,0}}; /* SiF2 */
const sym_pt_type SiC2_coord[]={{C, -RSIC_E,0,0},{C, RSIC_E,0,0}}; /* SiC2 */
const sym_pt_type Si2C_coord[]={{Si,-RSIC_E,0,0},{Si,RSIC_E,0,0}}; /* Si2C */

/* AX3's are planar in the xy plane with ax(1) bond in the +x direction */
const sym_pt_type CF3_coord[]={{F,RCF_E,0,0},
    {F,-RCF_E*0.5,RCF_E*s60,0}, 
    {F,-RCF_E*0.5,-RCF_E*s60,0}}; /* CF3 */
const sym_pt_type SiF3_coord[]={{F,RSIF_E,0,0},
    {F,-RSIF_E*0.5,RSIF_E*s60,0}, 
    {F,-RSIF_E*0.5,-RSIF_E*s60,0}}; /* SiF3 */
const sym_pt_type SiC3_coord[]={{C,RSIC_E,0,0},
    {C,-RSIC_E*0.5,RSIC_E*s60,0}, 
    {C,-RSIC_E*0.5,-RSIC_E*s60,0}}; /* SiC3 */
const sym_pt_type Si3C_coord[]={{Si,RSIC_E,0,0},
    {Si,-RSIC_E*0.5,RSIC_E*s60,0}, 
    {Si,-RSIC_E*0.5,-RSIC_E*s60,0}}; /* Si3C */

/* AX4's are tetragonal with ax(1) bond in +z dir and ax(2) bond in +x dir */
const sym_pt_type CF4_coord[]={{F,0,0,RCF_E},{F,s70_55*RCF_E,0,-s19_45*RCF_E}, 
    {F,-0.5*s70_55*RCF_E,s60*s70_55*RCF_E,-s19_45*RCF_E}, 
    {F,-0.5*s70_55*RCF_E,-s60*s70_55*RCF_E,-s19_45*RCF_E}}; /* CF4 */
const sym_pt_type SiF4_coord[]={{F,0,0,RSIF_E},
    {F,s70_55*RSIF_E,0,-s19_45*RSIF_E}, 
    {F,-0.5*s70_55*RSIF_E,s60*s70_55*RSIF_E,-s19_45*RSIF_E}, 
    {F,-0.5*s70_55*RSIF_E,-s60*s70_55*RSIF_E,-s19_45*RSIF_E}}; /* SiF4 */
const sym_pt_type SiC4_coord[]={{C,0,0,RSIC_E},
    {C,s70_55*RSIC_E,0,-s19_45*RSIC_E}, 
    {C,-0.5*s70_55*RSIC_E,s60*s70_55*RSIC_E,-s19_45*RSIC_E}, 
    {C,-0.5*s70_55*RSIC_E,-s60*s70_55*RSIC_E,-s19_45*RSIC_E}}; /* SiC4 */
const sym_pt_type Si4C_coord[]={{Si,0,0,RSIC_E},
    {Si,s70_55*RSIC_E,0,-s19_45*RSIC_E}, 
    {Si,-0.5*s70_55*RSIC_E,s60*s70_55*RSIC_E,-s19_45*RSIC_E}, 
    {Si,-0.5*s70_55*RSIC_E,-s60*s70_55*RSIC_E,-s19_45*RSIC_E}}; /* Si4C */

const iic_type IonInternalCoords[NULL_ION] =
{
    {He, 0, 0L},
    {C,  0, 0L},
    {F,  0, 0L},
    {Ne, 0, 0L},
    {Si, 0, 0L},
    {Ar, 0, 0L},
    {Kr, 0, 0L},
    {Xe, 0, 0L},
    {C,  1, CC_coord},      /* CC */
    {F,  1, FF_coord},      /* FF */
    {Si, 1, SiSi_coord},    /* SiSi */
    {C,  1, CF_coord},      /* CF */
    {Si, 1, SiF_coord},     /* SiF */
    {Si, 1, SiC_coord},     /* SiC */
    {C,  2, CF2_coord},     /* CF2 */
    {Si, 2, SiF2_coord},    /* SiF2 */
    {Si, 2, SiC2_coord},    /* SiC2 */
    {C,  2, Si2C_coord},    /* Si2C */
    {C,  3, CF3_coord},     /* CF3 */
    {Si, 3, SiF3_coord},    /* SiF3 */
    {Si, 3, SiC3_coord},    /* SiC3 */
    {C,  3, Si3C_coord},    /* Si3C */
    {C,  4, CF4_coord},     /* CF4 */
    {Si, 4, SiF4_coord},    /* SiF4 */
    {Si, 4, SiC4_coord},    /* SiC4 */
    {C,  4, Si4C_coord}     /* Si4C */
};

/* Parameter set A declarations for ion internal coordinates; parallel
 * to above code. */
const sym_pt_type CC_coord__A[] = {{C,0,0,RCC_E}};
const sym_pt_type FF_coord__A[] = {{F,0,0,RFF_E__A}};
const sym_pt_type SiSi_coord__A[] = {{Si,0,0,RSISI_E}};
const sym_pt_type SiF_coord__A[] = {{F,0,0,RSIF_E}};
const sym_pt_type CF_coord__A[] = {{F,0,0,RCF_E__A}};
const sym_pt_type SiC_coord__A[] = {{C,0,0,RSIC_E}};
const sym_pt_type CF2_coord__A[]={{F, -RCF_E__A, 0,0},{F, RCF_E__A, 0,0}};
const sym_pt_type SiF2_coord__A[]={{F, -RSIF_E,0,0},{F, RSIF_E,0,0}};
const sym_pt_type SiC2_coord__A[]={{C, -RSIC_E,0,0},{C, RSIC_E,0,0}};
const sym_pt_type Si2C_coord__A[]={{Si,-RSIC_E,0,0},{Si,RSIC_E,0,0}};
const sym_pt_type CF3_coord__A[]={{F,RCF_E__A,0,0},
    {F,-RCF_E__A*0.5,RCF_E__A*s60,0}, 
    {F,-RCF_E__A*0.5,-RCF_E__A*s60,0}}; /* CF3 */
const sym_pt_type SiF3_coord__A[]={{F,RSIF_E,0,0},
    {F,-RSIF_E*0.5,RSIF_E*s60,0}, 
    {F,-RSIF_E*0.5,-RSIF_E*s60,0}}; /* SiF3 */
const sym_pt_type SiC3_coord__A[]={{C,RSIC_E,0,0},
    {C,-RSIC_E*0.5,RSIC_E*s60,0}, 
    {C,-RSIC_E*0.5,-RSIC_E*s60,0}}; /* SiC3 */
const sym_pt_type Si3C_coord__A[]={{Si,RSIC_E,0,0},
    {Si,-RSIC_E*0.5,RSIC_E*s60,0}, 
    {Si,-RSIC_E*0.5,-RSIC_E*s60,0}}; /* Si3C */
const sym_pt_type CF4_coord__A[]={{F,0,0,RCF_E__A},{F,s70_55*RCF_E__A,0,-s19_45*RCF_E__A}, 
    {F,-0.5*s70_55*RCF_E__A,s60*s70_55*RCF_E__A,-s19_45*RCF_E__A}, 
    {F,-0.5*s70_55*RCF_E__A,-s60*s70_55*RCF_E__A,-s19_45*RCF_E__A}}; /* CF4 */
const sym_pt_type SiF4_coord__A[]={{F,0,0,RSIF_E},
    {F,s70_55*RSIF_E,0,-s19_45*RSIF_E}, 
    {F,-0.5*s70_55*RSIF_E,s60*s70_55*RSIF_E,-s19_45*RSIF_E}, 
    {F,-0.5*s70_55*RSIF_E,-s60*s70_55*RSIF_E,-s19_45*RSIF_E}}; /* SiF4 */
const sym_pt_type SiC4_coord__A[]={{C,0,0,RSIC_E},
    {C,s70_55*RSIC_E,0,-s19_45*RSIC_E}, 
    {C,-0.5*s70_55*RSIC_E,s60*s70_55*RSIC_E,-s19_45*RSIC_E}, 
    {C,-0.5*s70_55*RSIC_E,-s60*s70_55*RSIC_E,-s19_45*RSIC_E}}; /* SiC4 */
const sym_pt_type Si4C_coord__A[]={{Si,0,0,RSIC_E},
    {Si,s70_55*RSIC_E,0,-s19_45*RSIC_E}, 
    {Si,-0.5*s70_55*RSIC_E,s60*s70_55*RSIC_E,-s19_45*RSIC_E}, 
    {Si,-0.5*s70_55*RSIC_E,-s60*s70_55*RSIC_E,-s19_45*RSIC_E}}; /* Si4C */

const iic_type IonInternalCoords__A[NULL_ION] =
{
    {He, 0, 0L},
    {C,  0, 0L},
    {F,  0, 0L},
    {Ne, 0, 0L},
    {Si, 0, 0L},
    {Ar, 0, 0L},
    {Kr, 0, 0L},
    {Xe, 0, 0L},
    {C,  1, CC_coord__A},      /* CC */
    {F,  1, FF_coord__A},      /* FF */
    {Si, 1, SiSi_coord__A},    /* SiSi */
    {C,  1, CF_coord__A},      /* CF */
    {Si, 1, SiF_coord__A},     /* SiF */
    {Si, 1, SiC_coord__A},     /* SiC */
    {C,  2, CF2_coord__A},     /* CF2 */
    {Si, 2, SiF2_coord__A},    /* SiF2 */
    {Si, 2, SiC2_coord__A},    /* SiC2 */
    {C,  2, Si2C_coord__A},    /* Si2C */
    {C,  3, CF3_coord__A},     /* CF3 */
    {Si, 3, SiF3_coord__A},    /* SiF3 */
    {Si, 3, SiC3_coord__A},    /* SiC3 */
    {C,  3, Si3C_coord__A},    /* Si3C */
    {C,  4, CF4_coord__A},     /* CF4 */
    {Si, 4, SiF4_coord__A},    /* SiF4 */
    {Si, 4, SiC4_coord__A},    /* SiC4 */
    {C,  4, Si4C_coord__A}     /* Si4C */
};

/* Declare and initialize a pool of 10 ion types */
ionPool_t pool[] =
{
    {-1, 0.0}, {-1, 0.0}, {-1, 0.0}, {-1, 0.0}, {-1, 0.0}, 
    {-1, 0.0}, {-1, 0.0}, {-1, 0.0}, {-1, 0.0}, {-1, 0.0} 
};
char ionPool_filename_[20]="ionpool.dat";

/* Incident ion information data structure */
ion_type ion_ = 
/* type,pool,E_i,Th_i,P_i */
{-1, 0, 0.0, 0.0, 0.0, 
/* pRand, mem, mass, r0, r0set, r0Set, rExplic, v0, vExplic */
1, NULL, 0.0, {0, 0, 0}, {0, 0, 0}, 0, 0, {0, 0, 0}, 0, 
/* a1,a2,a3,aRand */
0.0, 0.0, 0.0, 1,
/* imp_ef, k,t,tex,rx,vx,r,v,Dr*/
0.75, 0.0, 0.0, 0.0, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, 
/* impt,nBonds0,nBonds,imp_*/
0.0, 0, 0, 0,
/* r0ex[] */
{{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0}}, 
/* v0ex[] */
{{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0}}, 
/* uint,oint,ofex */
1, 1, 0, "", "ion/####.ion"};

static char ln[255];
char dumion[10];
double dumwt=0.0;
int nIons=0;
extern short cfg_ionctrl_;
void ion_scanpool (void)
{
    FILE * fp=NULL;
    int i;
    
    fp=fopen(ionPool_filename_, "r");
    /* if the input file is not there, or the user
     * has already specified one particular ion to use, 
     * then don't attempt to set up the ion pool. */
    if (!fp||cfg_ionctrl_) return;
    ion_.pool=1;
    cfg_ionctrl_=1;
    i=0;
    printf("# reading ion pool from %s\n", ionPool_filename_);
    while (fgets(ln, 255, fp))
    {
        sscanf(ln, "%s %lf", dumion, &dumwt);
        pool[i].spec = ion_a2i(dumion);
        pool[i].wt=dumwt;
        printf("# %s(%i) %.5lf\n", dumion, ion_a2i(dumion), dumwt);
        i++;
    }
    nIons=i;
    printf("# %i ions in pool.\n", nIons);
    fclose(fp);
    /* normalize if needed */
    dumwt=0.0;
    for (i=0;i<nIons;i++) dumwt+=pool[i].wt;
    if (dumwt!=1.0) for (i=0;i<nIons;i++) pool[i].wt/=dumwt;
}

int ion_a2i (char * a)
{
    int i=0;
    for (i=0;i<NULL_ION&&strcmp(a,IonIDstr[i]);i++);
    if (i==NULL_ION)
	for (i=0;i<NULL_ION&&strcmp(a,IonIDstr_alt[i]);i++);
    return (i==NULL_ION?-1:i);
}

char * ion_i2a (int i)
{
    return IonIDstr[i];
}

static double tmpMat1[4][4]={{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};
static double tmpMat2[4][4]={{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};
static double rotMat[4][4]={{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};
static double ei,sinTheta,cosTheta,cosT,sinT,Theta,Phi,vel,sinPhi,cosPhi,zmin;
static double zbuf;
static double vx, vy, vz;
static sym_type symmin;
static pt R, V, TMPPT;
static ptPtr r=&R, v=&V, tpt=&TMPPT;
static  atomPtr a=NULL;
static nNodePtr n=NULL;
void ion_newion (void)
{
    void ion_reset(void);
    int i, k;
    
    if (ion_.pool)
    {
        ei=(double)rand()/RAND_MAX;
        vel=0.0;
        for (i=0;pool[i].spec!=-1;i++) 
        {
            vel+=pool[i].wt;
            if (vel>ei) break;
        }
        ion_.spec=pool[i].spec;
        printf("# ion selection from pool: r=%.5lf, ion=%s\n", 
            ei, ion_i2a(ion_.spec)); fflush(stdout);
    }
    if (ion_.spec==-1) {printf("error: Ion is not specified\n");exit(0);}
    if (!ion_.E_i) {printf("error: Ion energy is not specified\n");exit(0);}

    ion_reset();

    srand((unsigned)clock() + (unsigned)time(NULL));

    if (!quiet_) printf("# md series 2 ion creator (c) 1999 cfa\n");
    
    /* Using the data in the IonInternalCoords[] array, create
     * the ion with its base atom at (0, 0, 0). */
    a=cfg_newatom();
    ptPtr_clear(a->pos);
    nNode_push(&(ion_.mem), nNode_pop(&(a->bCards)));
    a->sym=IonInternalCoords[ion_.spec].base;
    a->state=IS_PROJECTILE;
    ion_.mass=per_table[a->sym].mass;
    n=ion_.mem;
    for (i=0;i<IonInternalCoords[ion_.spec].n;i++)
    {
	a=cfg_newatom();
	if (!(ion_.mem)) ion_.mem=nNode_pop(&(a->bCards));
	else 
	{
	    n->next=nNode_pop(&(a->bCards));
	    n=n->next;
	}
	a->sym=IonInternalCoords[ion_.spec].x[i].sym;
	if (paramset_a_)
	    ptPtr_copy(a->pos, 
		(ptPtr)&(IonInternalCoords__A[ion_.spec].x[i].r));
	else
	    ptPtr_copy(a->pos, 
		(ptPtr)&(IonInternalCoords[ion_.spec].x[i].r));
	a->state=IS_PROJECTILE;
	ion_.mass+=per_table[a->sym].mass;
    }
        
    /* Compute scalar ion velocity and cos and sin of incident polar angle */
    if (U_==STILLWEB) ei=ion_.E_i/EV_PER_EPSILON;
    else ei=ion_.E_i;
    vel=sqrt(2*ei/ion_.mass);
    Theta=M_PI/180.0*ion_.Th_i;  /* Theta in radians */
    sinTheta=sin(Theta);
    cosTheta=cos(Theta);
    
    vz=-vel*cosTheta;
    if (ion_.pRand) 
	ion_.P_i=180.0/M_PI*(Phi=(double)rand()/RAND_MAX*2*M_PI);
    else Phi=M_PI/180.0*ion_.P_i;
    cosPhi=cos(Phi);
    sinPhi=sin(Phi);
    vx=vel*sinTheta*cosPhi;
    vy=vel*sinTheta*sinPhi;
    printf("# Theta %.2lf d Phi %.2lf d ionMass %.5le\n", 
	 ion_.Th_i, ion_.P_i, ion_.mass);
    printf("# ion atomic velocity components: %.5lf %.5lf %.5lf\n", 
	vx, vy, vz);
    fflush(stdout);
    
    /* If we are randomizing the orientational angles of the ion,
     * compute these random angles. */
    if (ion_.aRand)
    {
	ion_.a1 = (double)rand()/RAND_MAX*360;
	ion_.a2 = (double)rand()/RAND_MAX*360;
 	ion_.a3 = (double)rand()/RAND_MAX*360;
    }
    printf("# ion xyz rot angles: %.2lf %.2lf %.2lf\n", 
        ion_.a1, ion_.a2, ion_.a3);
    /* Build the rotation matrix from the orientational angles */
    if (ion_.a3)
    {
	dblmat_identity(tmpMat1, 4);
	dblmat_identity(tmpMat2, 4);
	sinT = sin(M_PI/180.0*ion_.a3);
	cosT = cos(M_PI/180.0*ion_.a3);
	tmpMat1[0][0] = cosT;	/* xx */
	tmpMat1[0][1] = -sinT;	/* xy */
	tmpMat1[1][0] = sinT;	/* yx */
	tmpMat1[1][1] = cosT;	/* yy */
	dblmat_copy(tmpMat2, rotMat, 4);
	dblmat_matmult(rotMat, tmpMat2, tmpMat1, 4);
    }
    if (ion_.a2)
    {
	dblmat_identity(tmpMat1, 4);
	dblmat_identity(tmpMat2, 4);
	sinT = sin(M_PI/180.0*ion_.a2);
	cosT = cos(M_PI/180.0*ion_.a2);
	tmpMat1[0][0] = cosT;       /* xx */
	tmpMat1[0][2] = -sinT;      /* xz */
	tmpMat1[2][0] = sinT;       /* zx */
	tmpMat1[2][2] = cosT;       /* zz */
	dblmat_copy(tmpMat2, rotMat, 4);
	dblmat_matmult(rotMat, tmpMat2, tmpMat1, 4);
    }
    if (ion_.a1)
    {
	dblmat_identity(tmpMat1, 4);
	dblmat_identity(tmpMat2, 4);
	sinT = sin(M_PI/180.0*ion_.a1);
	cosT = cos(M_PI/180.0*ion_.a1);
	tmpMat1[1][1] = cosT;	/* yy */
	tmpMat1[1][2] = -sinT;	/* yz */
	tmpMat1[2][1] = sinT;	/* zy */
	tmpMat1[2][2] = cosT;	/* zz */
	dblmat_copy(tmpMat2, rotMat, 4);
	dblmat_matmult(rotMat, tmpMat2, tmpMat1, 4);
    }
    
    /* Assign velocity vector components; rotate around origin
     * by prescribed orientational angles; compute center-of-mass. */
    for (n=ion_.mem;n;n=n->next)
    {
	r=(a=n->addr)->pos;
	v=a->vel;
	v->x=vx; v->y=vy; v->z=vz;
	/* Rotate the ion around the origin using a post-multiply
	 * of the rotation matrix: */
	tpt->x=rotMat[0][0]*r->x+rotMat[1][0]*r->y+rotMat[2][0]*r->z;	
	tpt->y=rotMat[0][1]*r->x+rotMat[1][1]*r->y+rotMat[2][1]*r->z;	
 	tpt->z=rotMat[0][2]*r->x+rotMat[1][2]*r->y+rotMat[2][2]*r->z;
	r->x=tpt->x;
	r->y=tpt->y;
	r->z=tpt->z;
	ion_.r0.x+=per_table[a->sym].mass*r->x;
	ion_.r0.y+=per_table[a->sym].mass*r->y;
	ion_.r0.z+=per_table[a->sym].mass*r->z;
	
	/* if the user explicitly entered an ion COM-velocity,
	 * assign it here. */
	if (ion_.vExplic&&(ion_.v0.x||ion_.v0.y||ion_.v0.z))
	    ptPtr_copy(v, &(ion_.v0));
    }
    ptPtr_scalmult(&(ion_.r0), &(ion_.r0), 1.0/ion_.mass);
    
    /* In preparation for a random placement of the ion,
     * shift the ion's center-of-mass to (0, 0, 0), then
     * raise it in z until it is above the cutoff with
     * respect to the highest surface atom. */
    if (!ion_.r0Set)
    {
	/* Shift the positions of all the atoms in the ion so that
	 * its center-of-mass is at (0, 0, 0); set the ion's c.o.m.
	 * to (0, 0, 0); compute the zpos of the lowest atom in the ion. */
	zmin=1.e9;
	symmin=ion_.mem->addr->sym;
	for (n=ion_.mem;n;n=n->next)
	{
	    r=(a=n->addr)->pos;
	    r->x-=ion_.r0.x;	
	    r->y-=ion_.r0.y;	
	    r->z-=ion_.r0.z;
	    symmin=(r->z<zmin?a->sym:symmin);
	    zmin=(r->z<zmin?r->z:zmin);
	}
	ptPtr_clear(&(ion_.r0));	
    
	/* Based on (1) Zhi_ and (2) the z pos of the lowest atom in the ion and
	 * (3) these two atoms respective cutoff values, add a z-buffer to
	 * each atom's position so that it is beyond cutoff. */
	
	zbuf=atom_rcuij(symZhi_,symmin)+Zhi_-zmin;
	printf("# ion zbuf = %.5lf\n", zbuf);
	ion_.r0.z+=zbuf;
	for (n=ion_.mem;n;n=n->next)
	{
	    r=(a=n->addr)->pos;
	    r->z+=zbuf;
	}
    }
    
    /* Position the ion either by (1) placing it at the user's
     * specified center-of-mass position, or (2) placint it
     * at a random position in x and y.  If the position
     * is random, then the z-position of the ion is already
     * set by the code immediately above this. */
    tpt->x=(ion_.r0Set)?(ion_.r0set.x-ion_.r0.x):
                        (double)rand()/RAND_MAX*Lr_.x-half_Lr_.x;
    tpt->y=(ion_.r0Set)?(ion_.r0set.y-ion_.r0.y):
                        (double)rand()/RAND_MAX*Lr_.y-half_Lr_.y;
    tpt->z=(ion_.r0Set)?(ion_.r0set.z-ion_.r0.z):0.0;
    ion_.r0.x+=tpt->x;
    ion_.r0.y+=tpt->y;
    ion_.r0.z+=tpt->z;
    printf("# ion COMr %.10lf %.10lf %.10lf\n", 
        ion_.r0.x, ion_.r0.y, ion_.r0.z);
    k=0;
    for (n=ion_.mem;n;n=n->next)
    {
	r=(a=n->addr)->pos;
	v=a->vel;
	ptPtr_minimg(r, ptPtr_add(r, r, tpt));
	KE_+=0.5*ptPtr_sqabs(a->vel)*per_table[a->sym].mass;
	/* if Ion atom positions and/or velocities 
	 * were explicitly entered by the user,
	 * assign them here. */
	if (ion_.rExplic) ptPtr_copy(r, &(ion_.r0ex[k]));
	if (ion_.vExplic) ptPtr_copy(v, &(ion_.v0ex[k]));
	k++;
    }
    ion_.Dr.x=ion_.Dr.y=ion_.Dr.z=0.0;
}

void ion_reset (void)
{
    ion_.mem=NULL;
    ion_.ofex=0;
    ptPtr_clear(&(ion_.r0));
    ptPtr_clear(&(ion_.v0));
    ptPtr_clear(&(ion_.r));
    ptPtr_clear(&(ion_.v));
    ptPtr_clear(&(ion_.Dr));
    ion_.impt=0.0;
    ion_.imp_=0;
    ion_.nBonds0=ion_.nBonds=0;
}

static nNodePtr in, inn, jnn;
static atomPtr ia, ina;
static double ik, it, itex, thit;
extern double rt_, tlim0_;
static pt TPT={0, 0, 0}, *tpt2=&TPT;
void ion_update (void)
{
    #if DIAG
    if (!quiet_) printf("# md series 2 ion updater (c) 1999 cfa\n");
    #endif
    ik=it=itex=0.0;
    ion_.nBonds=0;
    ptPtr_copy(tpt2, &(ion_.r));
    ptPtr_clear(&(ion_.r));
    ptPtr_clear(&(ion_.v));
    for (in=ion_.mem;in;in=in->next)
    {
	ia=in->addr;
	ik+=ia->mv2;
	for (inn=ia->nList;inn;inn=inn->next)
	{
	    ina=inn->addr;
	    for (jnn=ion_.mem;jnn&&jnn->addr!=ina;jnn=jnn->next);
	    if (jnn) 
	    {
		thit=atom_peij(ia, jnn->addr);
		it+=thit;
		if (thit<0.0) ion_.nBonds++;
	    }
	    else itex+=0.5*atom_peij(ia, ina);
	}
	ptPtr_add(&(ion_.r), &(ion_.r), 
	    ptPtr_scalmult(tpt, ia->pos, per_table[ia->sym].mass));
	ptPtr_add(&(ion_.v), &(ion_.v), 
	    ptPtr_scalmult(tpt, ia->vel, per_table[ia->sym].mass));
	
    }
    ion_.nBonds/=2;
    if (!ion_.nBonds0) ion_.nBonds0=ion_.nBonds;
    ion_.k=0.5*ik;
    ion_.t=0.5*it;
    ion_.tex=itex;
    ptPtr_scalmult(&(ion_.r), &(ion_.r), 1.0/ion_.mass);
    ptPtr_scalmult(&(ion_.v), &(ion_.v), 1.0/ion_.mass);
    if (!(tpt2->x)&&!(tpt2->y)&&!(tpt2->z)) ptPtr_copy(tpt2, &(ion_.r));
    ptPtr_add(&(ion_.Dr), &(ion_.Dr),
	      ptPtr_minimg(tpt, ptPtr_subtract(tpt, &(ion_.r), tpt2)));
    if (!ion_.imp_ && 
       ((ion_.k/ion_.E_i <= ion_.imp_ef) || (ion_.nBonds < ion_.nBonds0)))
    {
	ion_.imp_=1;
	ion_.impt=rt_-tlim0_;
	ptPtr_copy(&(ion_.rx), &(ion_.r));
	ptPtr_copy(&(ion_.vx), &(ion_.v));
    }
}

static FILE * ifp=NULL;
static char * ifp_fmt =
    "%i %.5lf %.5le %.5le %.5le %.5lf %.5lf %.5lf %.5lf %i %i\n";
extern int t_;
void ion_output (void)
{
    #if DIAG
    printf("# md series 2 ion info outputter (c) 1999 cfa\n");
    printf("# Template output file name: %s\n", 
	ion_.outfile_0[0]?ion_.outfile_0:"stdout");
    printf("# Current output file name: %s\n", 
	ion_.outfile[0]?ion_.outfile:"stdout");
    #endif
    if (ion_.outfile[0]) 
    {
	if (ion_.ofex) 
	{
	    #if DIAG
	    printf("# explicit-append-open %s\n", ion_.outfile);
	    fflush(stdout);
	    #endif
	    ifp=fopen(ion_.outfile, "a");
	}
	else
	{
	    #if DIAG
	    printf("# attempting explicit-read-open %s\n", ion_.outfile);
	    fflush(stdout);
	    #endif
	    if (t_) ifp=fopen(ion_.outfile, "r");
	    if (ifp)
	    {
		fclose(ifp);
		#if DIAG
		printf("# implicit-append-open %s\n", ion_.outfile);
		fflush(stdout);
		#endif
		ifp=fopen(ion_.outfile, "a"); 
	    }
	    else ifp=fopen(ion_.outfile, "w");
	    #if DIAG
	    printf("# explicit-write-open %s<-%x\n", ion_.outfile, ifp);
	    fflush(stdout);
	    #endif
	    if (!ifp) {printf("Error: cannot open ion file %s\n",ion_.outfile);exit(0);}

	    ion_.ofex=1;
	    fprintf(ifp, "# md series 2 impacting ion data (c) 1999 cfa\n");
	    fprintf(ifp, "# Constituent atoms: ");
	    for (n=ion_.mem;n;n=n->next) 
		fprintf(ifp, "%s#%i ", per_table[n->addr->sym].sym, n->addr->id);
	    fprintf(ifp, "\n");
	    fprintf(ifp, "# Constituent atom initial configuration\n");
	    for (n=ion_.mem;n;n=n->next) 
	    {
		fprintf(ifp, "# ");
	        cfg_atomout(n->addr, ifp);
	    }
	    fprintf(ifp, "# Initial COM position: %.5lf %.5lf %.5lf\n", 
		ion_.r0.x, ion_.r0.y, ion_.r0.z);
	    fprintf(ifp, "# Initial energy: %.2lf\n", ion_.E_i);
	    fprintf(ifp, "# Impact @ ke/E_i = %.3lf\n", ion_.imp_ef);
	    fprintf(ifp, "# st, t, ke/E_i, ipe, epe, d.x, d.y, d.z, |d|, #b, imp?\n");
	}
	if (!ifp) {printf("Error: cannot open ion file %s\n",ion_.outfile);exit(0);}
    }
    else ifp=stdout;
    fprintf(ifp, ifp_fmt, 
	    t_, rt_, ion_.k/ion_.E_i, ion_.t, ion_.tex, 
	    ion_.Dr.x, ion_.Dr.y, ion_.Dr.z, sqrt(ptPtr_sqabs(&(ion_.Dr))), 
	    ion_.nBonds, ion_.imp_);
    if (ifp!=stdout) fclose(ifp);
    ifp=NULL;
}

