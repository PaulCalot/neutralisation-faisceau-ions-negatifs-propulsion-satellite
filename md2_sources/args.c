#include "args.h"
/* 
 * md series 2
 * (c) 1999 cameron abrams & regents of uc berkeley
 * dept of chemical engineering
 * university of california, berkeley 
 * 
 *
 * args:  the command-line/setup file/cfg file argument handling module.
 * Globals are assigned based on command-line keyword/value pairs.
 * Arguments are processed in the order in which they are issued.  The 
 * special argument "-s" or "--setup" is followed by the name of a setup 
 * file that contains keyword/value pairs.  These values are overwritten 
 * if any valid arguments pairs appear after the "-s"-pair on the command 
 * line.  Likewise, any valid argument pairs *preceding* the "-s"-pair 
 * are possibly overwritten if the same keyword/value pairs appear in the
 * setup file.
 *
 * PLEASE NOTE:
 * Parameters found in the cfg file header have total override priveledges.
 * This is because the cfg scanner actually builds a second argv[] array
 * from header information in the cfg file and then calls arg_handler()
 * process them.  Certain parameters are only specified in the cfg file, 
 * but any parameters that appear in a setup or on the command line
 * can be put in the header section of the cfg file.  If such a parameter
 * is detected in the cfg and its value was previously set either by 
 * a command line argument or a line in a setup file,  that value is
 * overwritten by the value found in the cfg file.
 */

/* Here is the list of globals accessible by this module: */
extern char cfg_infile_[];  /* name of input cfg file */
char cfg_infile_0[MAXCHAR]="cfg/####.cfg";
extern char cfg_outfile_[]; /* name of output cfg file */
char cfg_outfile_0[MAXCHAR]="";
extern double dt_;	    /* time step */
extern short dtvar_;	    /* turns on/off time step adaptive variation */
extern double Lc_;	    /* max position del allowed per time step */
extern double maxdt_;	    /* maximum allowed time step value */
extern double mindt_;	    /* maximum allowed time step value */
extern double tlim0_;	    /* integration time lower limit */
extern double tlim1_;	    /* integration time upper limit */
extern int Ndt_;	    /* explicit number of time steps */
extern unit_type U_;	    /* unit system denoter */
extern short quiet_;	    /* flag to supress headers and 
			     * messages in the output */
extern pt Lr_, half_Lr_;    /* Box size and 1/2 box size */
extern ipt nCell_;	    /* Box size in number of unit cells */
extern double cfgTime_;	    /* Current cfg creation time */
extern ipt per_;	    /* periodicity flags (def=110) */
extern char cfgName_[];	    /* Name of config -- not the file name */
extern char cfgCreateDate_[]; /* string cfg creation date */
extern char unitCell_[];    /* string unit cell designation */

int impact_num_;	    /* impact number designation */
int impact0_num_=1;	    /* first impact number designation */
int impact1_num_=1;	    /* last impact number designation */

extern int cfg_outint_;     /* number of time steps between in-run cfg
			     * outputs */
extern int data_int_;	    /* number of time steps between integration
			     * data output */
			     
extern char cfg_snapfile_[];/* name (as a template) for output of in-run cfgs */

extern double bhb_Tset_;    /* Berendsen-type heat bath set point, Kelvins */
extern double bhb_tau_;	    /* Berendsen-type heat bath rise time */
extern int bhb_start_;	    /* Time step at which heat bath is applied */
extern double bmb_Pset_;    /* Berendsen-type pressure bath set point, GPa */
extern double bmb_tau_;	    /* Berendsen-type pressure bath rise time */
extern int bmb_start_;	    /* Time step at which pressure bath is applied */

extern short cfg_ionctrl_;  /* control parameter for the ion module */
extern char clust_outfile_[];
extern short clust_report_;  /* allows user to turn off cluster reporting */
char clust_outfile_0[MAXCHAR]="clu/####.clu";

/* flags for controlling choice of parameter sets in
 * the bicubic and tricubic modules for the Tanaka CF
 * potential. */
extern short hcc_opt_;
extern short hcf_opt_;
extern short fcc_opt_;
extern short paramset_a_;   /* default is to use parameter set B */

double ftol_;		    /* tolerance for the conjugate gradient module */
extern double ff_;	    /* Fraction of atoms to fix at bottom
			     * if this run is just to make a crystal */
extern double xl_Re_;	    /* For creation of a crystal, the user may
			     * specify an equilibrium bond length */			     
extern ion_type ion_;	    /* Bombarding ion information */
extern char ionPool_filename_[];

extern int bond_ctrl_;	    /* bond definition control parameter */
extern int des_ctrl_;	    /* desorption scheme control parameter */
extern short USE_FCC;	    /* toggles the Fcc correction in the tbt PEF */
extern short FCC_TESTER;    /* toggles the Fcc tester in the tbt PEF */
extern short HIJ_TESTER;    /* toggles the Hij tester in the tbt PEF */
extern short USE_HIJ;	    /* toggles the Hij correction in the tbt PEF */
extern short USE_RCU;	    /* toggles the use of cutoffs in the tbt PEF */
extern short USE_HIE_CORE;  /* toggles the use of high-energy repuslive
			     * core PEFs a la Beardmore and Smith */
extern double phi_rate_tol_;
extern double fcc_a_t;
extern double hcc_a_t;
extern int SIC_DIAG2_;
extern int SIC_DIAG3_;
extern double td_tau_, td_A_;
			    /* Binding energy parameters for 1st-order
			     * desorption theory */

extern short RDF_;
extern short RDF_34_;
extern int rdf_start_;
extern double rdf_loZlim_;
extern double RDF_dR_;
extern char RDF_filename_[];

short echo_=0;		    /* echos all parameter assignments */

char setup_infile_[MAXCHAR]="";
char setup_outfile_[MAXCHAR]="";

/* Argument handling:  arguments are handled as "keyword/value" pairs.
 * In general, arguments have "values" that are vector quantities.
 * For example, the argument for the simulation domain box size
 * has a value of three doubles, the x, y, and z dimensions of the
 * box.  Additionally, the argument may have a switch variable
 * associated with it.  An argument has two keywords:  one is
 * the abbreviated "-"prefixed command line keyword,  and the second
 * is a longer, more descriptive version that might appear in an
 * input file.  Each argument has a bit field that is set if 
 * that argument is set by the user, and remains unset if the 
 * hardcoded default value is used.  Finally, each argument uses
 * one of several functions to assign its specific value from 
 * the string value that is read in at run time. */
#define MAXARG 100
typedef struct ARGFUNC 
{
    int i;			/* Unique id number of this keyword, 
				 * assigned at runtime, never really used */
    short s;			/* Flag that is set if the runtime 
				 * user sets this keyword */
    int l;			/* Number of elements in value */
    char * key1;		/* Keyword 1; command-line */
    char * key2;		/* Keyword 2; input file */
    void * var;			/* Value; points to a vector with "l" members */
    short * var2;		/* Switch (perhaps a flag) */
    void (*asnf)(void*, char*); /* Assignment function */
    void (*asnf2)(void*, char*);/* Auxiliary assignment function */
    void (*prnf)(FILE*,char*,void*); /* Output function */
} arg_type;

/* Pointers for type conversion */
static int *intp;
static double *dblp;
static char *chrp;
static pt   *ptp;
static ipt  *ipp;
static unit_type *utp;
static short *shtp;
/* Assignment functions */
void dblasn (void *v, char *s){dblp=(double*)v;*dblp=(double)atof(s);}
void intasn (void *v, char *s){intp=(int*)v;*intp=(int)atoi(s);}
void shtasn (void *v, char *s){shtp=(short*)v;*shtp=(short)atoi(s);}
void shton  (void *v, char *s){shtasn(v,"1");}
void shtoff (void *v, char *s){shtasn(v,"0");}
void usasn  (void *v, char *s){utp=(unit_type*)v;*utp=chem_a2ut(s);}
void ionasn (void *v, char *s){intp=(int*)v;*intp=ion_a2i(s);}
void strasn (void *v, char *s){chrp=(char*)v;strcpy(chrp, s);}
void pntasn (void *v, char *s){ptp=(pt*)v;sscanf(s,"%lf %lf %lf",
			     &(ptp->x),&(ptp->y),&(ptp->z));}
void paeasn (void *v, char *s)
{
    int i=0;double x, y, z;
    ptp=(pt*)v;
    sscanf(s,"%i %lf %lf %lf",&i,&(x),&(y),&(z));
    (ptp+i)->x=x;(ptp+i)->y=y;(ptp+i)->z=z;
}
void iptasn (void *v, char *s){ipp=(ipt*)v;sscanf(s,"%i %i %i",
			     &(ipp->i),&(ipp->j),&(ipp->k));}
/* Output functions */
void dblprn (FILE *f, char *s, void *v)
{dblp=(double*)v;fprintf(f, "%s %.5lf", s, *dblp);fflush(f);}
void intprn (FILE *f, char *s, void *v)
{intp=(int*)v;fprintf(f, "%s %i", s, *intp);fflush(f);}
void shtprn (FILE *f, char *s, void *v)
{shtp=(short*)v;fprintf(f, "%s %i", s, *shtp);fflush(f);}
void strprn (FILE *f, char *s, void *v)
{chrp=(char*)v;fprintf(f, "%s %s", s, chrp);fflush(f);}
void pntprn (FILE *f, char *s, void *v)
{ptp=(pt*)v;fprintf(f, "%s %lf %lf %lf",s,ptp->x,ptp->y,ptp->z);fflush(f);}
void iptprn (FILE *f, char *s, void *v)
{ipp=(ipt*)v;fprintf(f, "%s %i %i %i",s,ipp->i,ipp->j,ipp->k);fflush(f);}
void usprn (FILE *f, char *s, void *v)
{utp=(unit_type*)v;fprintf(f, "%s %s",s,chem_ut2a(*utp));fflush(f);}
void ionprn (FILE *f, char *s, void *v)
{intp=(int*)v;fprintf(f, "%s %s",s,ion_i2a(*intp));fflush(f);}


/* Fields of arg_type are:
 * (Unique ID number of this key; set by the fcn arg_initialize)
 * (Set-at-runtime-by-user flag)
 * (# of values anticipated after keyword)
 * (Command-line abbreviated keyword)
 * (More descriptive equivalent keyword as may appear in setup file)
 * (Address of variable to which keyword corresponds)
 * (Address of auxiliary short that corresponds to the keyword, if any)
 * (Pointer to function that assigns the value to the variable)
 * (Pointer to function that turns on/off the aux. short, if any)
 * (Pointer to function that outputs the variable)
 */
static arg_type keys[] = {
/* Strings */
{0,0,1,"-ic", "CfgInputFile",   cfg_infile_0,	0L,strasn,0L,strprn}, 
{0,0,1,"-is", "SetupInputFile", setup_infile_,  0L,strasn,0L,strprn}, 
{0,0,1,"-os", "SetupOutputFile",setup_outfile_, 0L,strasn,0L,strprn}, 
{0,0,1,"-oc", "CfgOutputFile",  cfg_outfile_0,	0L,strasn,0L,strprn}, 
{0,0,1,"-ocl","ClustOutputFile",clust_outfile_0,0L,strasn,0L,strprn}, 
{0,0,1,"-oi", "IonOutputFile",  ion_.outfile_0,	0L,strasn,0L,strprn}, 
{0,0,1,"-csn","CfgSnapshotFile",cfg_snapfile_,  0L,strasn,0L,strprn}, 
{0,0,1,"-cnm","CfgNameDesc",    cfgName_,	0L,strasn,0L,strprn}, 
{0,0,1,"-cdt","CreationDate",   cfgCreateDate_, 0L,strasn,0L,strprn}, 
{0,0,1,"-uc", "UnitCell",	unitCell_,	0L,strasn,0L,strprn}, 
{0,0,1,"-rdff", "RDFFileName",  RDF_filename_,  0L,strasn,0L,strprn}, 
{0,0,1,"-ipf","IonPoolFile",ionPool_filename_,&(cfg_ionctrl_),strasn,shton,strprn}, 
/* Floating-Points */
{0,0,1,"-dt",  "TimeStep",      &dt_,		0L,dblasn,0L,dblprn}, 
{0,0,1,"-mndt","MinTimeStep",   &mindt_,	0L,dblasn,0L,dblprn}, 
{0,0,1,"-mxdt","MaxTimeStep",   &maxdt_,	0L,dblasn,0L,dblprn}, 
{0,0,1,"-t0",  "LowerIntegLim", &tlim0_,	0L,dblasn,0L,dblprn}, 
{0,0,1,"-t1",  "UpperIntegLim", &tlim1_,	0L,dblasn,0L,dblprn}, 
{0,0,1,"-ff",  "FixFrac",	&ff_,		0L,dblasn,0L,dblprn}, 
{0,0,1,"-xlre","CrystRe",	&xl_Re_,	0L,dblasn,0L,dblprn}, 
{0,0,1,"-Lx",  "BoxSize.x",     &(Lr_.x),	0L,dblasn,0L,dblprn}, 
{0,0,1,"-Ly",  "BoxSize.y",     &(Lr_.y),	0L,dblasn,0L,dblprn}, 
{0,0,1,"-Lz",  "BoxSize.z",     &(Lr_.z),	0L,dblasn,0L,dblprn}, 
{0,0,1,"-ctm", "ConfigTime",    &cfgTime_,	0L,dblasn,0L,dblprn}, 
{0,0,1,"-Tset","BerendsenT",    &bhb_Tset_,	0L,dblasn,0L,dblprn}, 
{0,0,1,"-tau", "BerendsenTau",  &bhb_tau_,	0L,dblasn,0L,dblprn}, 
{0,0,1,"-Pset","BerendsenP",    &bmb_Pset_,	0L,dblasn,0L,dblprn}, 
{0,0,1,"-ptau","BerendsenPTau", &bmb_tau_,	0L,dblasn,0L,dblprn}, 
{0,0,1,"-tdt", "ThmlDesThryTau",&td_tau_,	0L,dblasn,0L,dblprn}, 
{0,0,1,"-tdA", "ThmlDesThryA",  &td_A_,		0L,dblasn,0L,dblprn}, 
{0,0,1,"-ftol","ConjGradTol",   &ftol_,		0L,dblasn,0L,dblprn}, 
{0,0,1,"-ionief","IonImpEFrac", &(ion_.imp_ef),	0L,dblasn,0L,dblprn}, 
{0,0,1,"-ionE","IonEnergy",     &(ion_.E_i),	0L,dblasn,0L,dblprn}, 
{0,0,1,"-ionT","IonTheta",      &(ion_.Th_i),	0L,dblasn,0L,dblprn}, 
{0,0,1,"-ionP","IonPhi",	&(ion_.P_i),&(ion_.pRand),dblasn,shtoff,dblprn},
{0,0,1,"-prt","PhiRateTol",     &(phi_rate_tol_),0L,dblasn,0L,dblprn}, 
{0,0,1,"-rdflz", "RDFLowZlimit",&rdf_loZlim_,   0L,dblasn,0L,dblprn}, 
{0,0,1,"-rdfdr", "RDFResolution",&RDF_dR_,      0L,dblasn,0L,dblprn}, 
{0,0,1,"-fccta","fcctestA",     &(fcc_a_t),	0L,dblasn,0L,dblprn}, 
{0,0,1,"-hccta","hcctestA",     &(hcc_a_t),	0L,dblasn,0L,dblprn}, 

/* Integers */
{0,0,1,"-ncx", "nCellx",	&(nCell_.i),    0L,intasn,0L,intprn},
{0,0,1,"-ncy", "nCelly",	&(nCell_.j),    0L,intasn,0L,intprn},
{0,0,1,"-ncz", "nCellz",	&(nCell_.k),    0L,intasn,0L,intprn},
{0,0,1,"-n",   "NumOfTimeSteps",&Ndt_,		0L,intasn,0L,intprn},
{0,0,1,"-i",   "ImpactNum",   &impact_num_,	0L,intasn,0L,intprn},
{0,0,1,"-i0",  "FirstImpactNum",&impact0_num_,  0L,intasn,0L,intprn},
{0,0,1,"-i1",  "LastImpactNum", &impact1_num_,  0L,intasn,0L,intprn},
{0,0,1,"-soi", "SnapshotInt",   &cfg_outint_,   0L,intasn,0L,intprn},
{0,0,1,"-idi", "IntDataInt",    &data_int_,	0L,intasn,0L,intprn},
{0,0,1,"-iui", "IonUpdateInt",  &(ion_.uint),   0L,intasn,0L,intprn},
{0,0,1,"-ioi", "IonOutputInt",  &(ion_.oint),   0L,intasn,0L,intprn},
{0,0,1,"-bc",  "BondDefCtrlPar",&(bond_ctrl_),	0L,intasn,0L,intprn},
{0,0,1,"-dc",  "DesorbCtrlPar", &(des_ctrl_),   0L,intasn,0L,intprn},
{0,0,1,"-hbn", "HeatBathStartStep",&(bhb_start_),0L,intasn,0L,intprn},
{0,0,1,"-mbn", "MomentumBathStartStep",&(bmb_start_),0L,intasn,0L,intprn},
{0,0,1,"-rdfn","RDFStartStep", &(rdf_start_),0L,intasn,0L,intprn},
{0,0,1,"-rdf", "RDFOn", &(RDF_),0L,shton,0L,intprn},
{0,0,1,"-d3",  "Sicf_diag3",	&(SIC_DIAG3_),	0L,intasn,0L,intprn},
{0,0,1,"-d2",  "Sicf_diag2",	&(SIC_DIAG2_),	0L,intasn,0L,intprn},

/* Flags */
{0,0,0,"-oldhcc", "UseTanakaIHcc",&hcc_opt_,	0L,shton, 0L,shtprn}, 
{0,0,0,"-oldhcf", "UseTanakaIHcf",&hcf_opt_,	0L,shton, 0L,shtprn}, 
{0,0,0,"-oldfcc", "UseTanakaIFcc",&fcc_opt_,	0L,shton, 0L,shtprn}, 
{0,0,0,"-paramA", "UseParamSetA",&paramset_a_,	0L,shton, 0L,shtprn}, 
{0,0,0,"-dtv", "VarTimeStep",   &dtvar_,        0L,shton, 0L,shtprn}, 
{0,0,0,"+dtv", "AdaptTimeStep", &dtvar_,	0L,shtoff,0L,shtprn}, 
{0,0,0,"-nocl", "NoClusterReport", &clust_report_, 0L,shtoff,0L,shtprn}, 
{0,0,0,"-q",   "Quiet",         &quiet_,        0L,shton, 0L,shtprn}, 
{0,0,0,"-echo","Echo",		&echo_,		0L,shton, 0L,shtprn}, 
{0,0,0,"-nocore", "NoCore",	&USE_HIE_CORE,	0L,shtoff,0L,shtprn}, 
{0,0,0,"+fcc", "NoFcc",		&USE_FCC,	0L,shtoff,0L,shtprn}, 
{0,0,0,"-fcctest", "fcctest",	&FCC_TESTER,	0L,shton,0L,shtprn}, 
{0,0,0,"-hijtest", "hijtest",	&HIJ_TESTER,	0L,shton,0L,shtprn}, 
{0,0,0,"+hij", "NoHij",		&USE_HIJ,	0L,shtoff,0L,shtprn}, 
{0,0,0,"-rdf34", "RDF34",		&RDF_34_,	0L,shton,0L,shtprn}, 
/*{0,0,0,"+rcu", "NoCutoffs",   &USE_RCU,	0L,shtoff,0L,shtprn}, */
/* Points */
{0,0,3,"-L",   "BoxSize.xyz",   &Lr_,		0L,pntasn,0L,pntprn}, 
{0,0,3,"-ionRA","IonRotAngles.xyz",&(ion_.a1),&(ion_.aRand),pntasn,shtoff,pntprn},
{0,0,3,"-ionCOMr","IonCOMpos.xyz",&(ion_.r0set),&(ion_.r0Set),pntasn,shton,pntprn},
/* I-Points */
{0,0,3,"-nc",  "nCell.xyz",     &nCell_,	0L,iptasn,0L,iptprn}, 
{0,0,3,"-per", "Periodicities", &per_,		0L,iptasn,0L,iptprn}, 
/* Point array elements */
{0,0,4,"-ionr","IonAtomPos",    ion_.r0ex,&(ion_.rExplic),paeasn,shton,NULL}, 
{0,0,4,"-ionv","IonAtomVel",   &(ion_.v0ex),&(ion_.vExplic),paeasn,shton,NULL},
/* Special types */
{0,0,1,"-u",   "UnitSystem",    &U_,	    0L,usasn, 0L,usprn}, 
{0,0,1,"-ion", "Ion", &(ion_.spec),&(cfg_ionctrl_),ionasn,shton,ionprn}, 
/* Terminal */
{0,0,0, "xxx", "XXX", 0L, 0L, 0L, 0L, 0L}
};
int nKws_=0;
static short arg_initialized=0;
void arg_initialize (void)
{
    nKws_=0;
    while (keys[nKws_].var) keys[nKws_].i=nKws_++;
    arg_initialized=1;
}

static char buf0[MAXCHAR]="";
static short setup_=0;
static int nArgCalls_=0;
void arg_handler (int argc, char * argv[], short bwc)
{
    int i=0,a=0,l=0;
    char * p;
    short found=0;
    void help (char *);
    void scan_setup (char *);
    void output_setup (char *, short);

    if (!arg_initialized) arg_initialize();
    
    if (echo_) echo_args(argc, argv);
    
    for (a=0;a<argc;a++)
    {
	i=0;
	found=0;
	while (keys[i].var&&!found)
	{
	    if (!strcmp(argv[a], keys[i].key1)||
		!strcmp(argv[a], keys[i].key2)) 
	    {
		buf0[0]='\0';
		p=buf0;
		for (l=1;l<=keys[i].l;l++) 
		{
		    sprintf(p, "%s%s", (l>1?" ":""), argv[++a]);
		    p+=strlen(argv[a])+(l>1?1:0);
		}
		if (echo_)
		{printf("# %s (%s)\n", keys[i].key2, buf0);fflush(stdout);}
		keys[i].asnf(keys[i].var, (keys[i].l?buf0:NULL));
		if (keys[i].asnf2) keys[i].asnf2(keys[i].var2,NULL);
		keys[i].s=strcmp(keys[i].key1, "-is")?1:0;
		found=1; 
	    }
	    i++;
	}
	if (a&&!found&&!keys[i].var&&bwc==TRAP_BADWORDS) help(argv[a]);
    }
    if (!quiet_)
    {
	printf("# md series 2 argument handler, call %i\n", ++nArgCalls_);
	fflush(stdout);
    }

    /* If a setup input file is named, read it ONCE. */
    if (bwc!=NO_RECURSION&&setup_infile_[0]&&!setup_) 
	{scan_setup(setup_infile_);setup_=1;}
    /* If a setup output file is named, write it ONCE. */
    if (!setup_infile_[0]&&setup_outfile_[0]&&!setup_) 
	{output_setup(setup_outfile_, 0);setup_=1;}

    /* write the setup info to stdout */
    if (!quiet_) output_setup(NULL, 1);

}

void output_setup (char * fn, short cmt)
{
    FILE * fp=NULL;
    int i;
    
    if (!fn) fp=stdout;
    else fp=fopen(fn, "w");
    if (!fp) {printf("Error: cannot create setup file %s\n", fn);exit(0);}
    i=0;
    while (keys[i].var)
    { 
	if (keys[i].s) 
	{
	    if (cmt) printf("# ");fflush(fp);
	    keys[i].prnf(fp, keys[i].key2, keys[i].var);
	    fprintf(fp, "\n");
	}
	i++;
    }
    if (fp!=stdout) fclose(fp);
}

char ln[MAXCHAR];
static char argv_BUF[MAXARG][MAXCHAR];
static char * argv_buf[MAXARG];
static int argc_buf;
void scan_setup (char * fn)
{
    FILE * fp=NULL;
    int i, j;
    char *cp;
    fp=fopen(fn, "r");
    if (!fp) {printf("error: %s setup file not found.\n", fn);return;}
    for (i=0;i<MAXARG;i++) argv_buf[i]=&(argv_BUF[i][0]);
    argc_buf=0;
    while (fgets(ln, MAXCHAR, fp))
    {
	for (j = 0; j < MAXCHAR && isspace(ln[j]); j++);
	if (!isdigit(ln[j]) && ln[j] != '.' && ln[j] != '-')
	{
	    cp=ln;
	    while (cp[0]=='#' || isspace(cp[0])) cp++;
	    if (cp[0]!='%'&&cp[0]!='!'&&cp[0]!=';')
	    {
		#if DIAG
		printf("parsing %s", cp);fflush(stdout);
		#endif
		parseline(cp, argv_buf, &argc_buf);
	    }
	    #if DIAG
	    else {printf("skipping %s", cp);fflush(stdout);}
	    #endif
	}
    }
    arg_handler(argc_buf, argv_buf, NO_RECURSION);
    fclose(fp);
}

void help (char * arg)
{
    printf("# md series 2 help: ");
    printf("keyword %s is not recognized.\n", arg);
    exit(0);
}

void parseline (char * ln, char * argv[], int * argc)
{
    char *c;
    int i=-1;
    c=ln;
    while (*c)
    {
	if (*c&&!isspace(*c))
	    argv[*argc][++i]=*c;
	else 
	{
	    while(*c&&isspace(*c)) c++;
	    if (*c) argv[++*argc][i=0]=*c;
	}
	if (*c) c++;
    }
    ++*argc;
}

void echo_args (int argc, char * argv[])
{
    int i;
    printf("# ");
    for (i=0;i<argc;i++) printf("%s ", argv[i]);
    printf("\n");
}
    
void filenames_initialize (void)
{
    if (!quiet_) printf("# md series 2 filenames initializer (c) 1999 cfa\n");
    if (cfg_infile_0[0]) strcpy(cfg_infile_, cfg_infile_0);
    if (cfg_outfile_0[0]) strcpy(cfg_outfile_, cfg_outfile_0);
    if (clust_outfile_0[0]) strcpy(clust_outfile_, clust_outfile_0);
    if (ion_.outfile_0[0]) strcpy(ion_.outfile, ion_.outfile_0);
}

void filenames_apply (void)
{
    if (!quiet_&&echo_) 
	{printf("# md series 2 filenames applier (c) 1999 cfa\n");fflush(stdout);}
    if (cfg_infile_0[0])
    {
	strcpy(cfg_infile_, cfg_infile_0);
	str_q2n(cfg_infile_, impact_num_-1);
	str_p2n(cfg_infile_, impact_num_-1);
    }
    if (cfg_outfile_0[0])
    {
	strcpy(cfg_outfile_, cfg_outfile_0);
	str_q2n(cfg_outfile_, impact_num_);
	str_p2n(cfg_outfile_, impact_num_);
    }
    if (clust_outfile_0[0])
    {
	strcpy(clust_outfile_, clust_outfile_0);
	str_q2n(clust_outfile_, impact_num_);
	str_p2n(clust_outfile_, impact_num_);
    }
    if (ion_.outfile_0[0])
    {
	strcpy(ion_.outfile, ion_.outfile_0);
	str_q2n(ion_.outfile, impact_num_);
	str_p2n(ion_.outfile, impact_num_);
    }
    if (!quiet_&&echo_) 
	{printf("# md series 2 filenames applier complete\n");fflush(stdout);}
}

void str_q2n (char * f, int n)
/* str_q2n replaces a substring of ?'s in f with a 0-padded 
 * integer n.  e.g., if f is "prefix.?????.suffix" and n is 342, 
 * then f is altered to become "prefix.00342.suffix". */
{
    char *p0, *p1, *p, *q;
    char nstr[10];
    if (!f) return;
    p0=strchr(f, '?');
    p1=strrchr(f, '?');
    if (!p0||!p1) return;
    sprintf(nstr, "%i", n);
    q=p1;
    p=&(nstr[strlen(nstr)-1]);
    while (p>=nstr&&q>=p0) (*(q--))=(*(p--));
    while (q>=p0) (*(q--))='0';
}

void str_p2n (char * f, int n)
/* str_q2n replaces a substring of #'s in f with a 0-padded 
 * integer n.  e.g., if f is "prefix.#####.suffix" and n is 342, 
 * then f is altered to become "prefix.00342.suffix". */
{
    char *p0, *p1, *p, *q;
    char nstr[10];
    if (!f) return;
    p0=strchr(f, '#');
    p1=strrchr(f, '#');
    if (!p0||!p1) return;
    sprintf(nstr, "%i", n);
    q=p1;
    p=&(nstr[strlen(nstr)-1]);
    while (p>=nstr&&q>=p0) (*(q--))=(*(p--));
    while (q>=p0) (*(q--))='0';
}

void filenames_echo (void)
{
    if (cfg_infile_[0]) printf("# input cfg filename: %s\n", 
	cfg_infile_);
    if (cfg_outfile_[0]) printf("# output cfg filename: %s\n", 
	cfg_outfile_);
    if (clust_outfile_[0]) printf("# clust data filename: %s\n", 
	clust_outfile_);
    if (ion_.outfile[0]) printf("# ion data filename: %s\n", 
	ion_.outfile);
    fflush(stdout);
}
