/*
 * MD series 2 -- Molecular Dynamics of Plasma-Surface Chemistry
 * (c) 1999 
 * Cameron Abrams
 * and 
 * The Regents of the University of California, Berkeley
 * 
 * Department of Chemical Engineering
 * University of California, Berkeley	1995-1999
 * 
 * "cfg" is responsible for defining and handling configuration data.
 * 
 */

#include "cfg.h"

extern double version_;
extern int build_;

/* The configuration; the arrays are static because the data in them
 * can only be accessed though the "atom" interface. */
static pt pos_[MAXNUMATOMS];		/* positions */
static pt vel_[MAXNUMATOMS];		/* velocities */
static pt hac_[MAXNUMATOMS];		/* half-accelerations */
static pt frc_[MAXNUMATOMS];		/* forces */

/* Configuration interface */
static atom atomarray_[MAXNUMATOMS];	/* array of atoms */
atomPtr L_=&(atomarray_[0]);		/* address of head of array of atoms,
					 * as atoms are accessed in a linked-list
					 * fashion by some routines */ 
atomPtr TH_=NULL;			/* the "trashHeap": place to keep 
					 * "trashed" atoms, or atoms removed
					 * from the configuration. */
static FILE * cfg_fp=NULL;		/* For input/output */
/* A word about L_ and TH_:  Both of these can be thought of as pointers
 * to heads of linked lists of atom data types.  The actual atom data type
 * instances reside in the static array "atomarray_[]", but their "next"
 * pointers can be assigned to build any combination of ordering into
 * an arbitrary number of "linked lists".  Here, I only care about two
 * such lists: (1) the interface list that reflects the configuration that
 * is the operand for the force routine and integrator, and (2) atoms
 * that have been removed (e.g. because they are sputtered) from the
 * configuration. */

/* Business Cards:  These are neighbor-list members, or neighbor identifiers.
 * Each atom begins owning a supply of cards, each of which bears its owner's
 * address.  When a neighbor of an atom is identified, the atom and this
 * neighbor exchange business cards, so that now each has a record that the
 * other is its neighbor.  These nodes are referenced in two kinds of lists for
 * each atom:  (1) "bCards" points to the head of a linked list containing
 * the remaining supply of business cards owned by the atom; 
 * (2) nList is the owner's list of neighbors. */
static nNode bcardarray_[MAXNUMATOMS][MAXNUMNEIGHBORS];

/* Configuration ID/control information */
short cfg_initialized_=0;
short cfg_established_=0;
short cfg_newcrystal_=0;
int nAtom_=0;			/* Total number of atoms defined in the
				 * static array */
int nFixed_=0, nTrash_=0;
int nElem_[MAXNUMELEMENTS];
char cfgName_[MAXCHAR];
char cfgCreateDate_[MAXCHAR];
double cfgTime_=0.0;
double T_=0.0;
double KE_=0.0;
double sKE_=0.0;
double Zhi_=0.0, Zlo_=0.0;
sym_type symZhi_=XX;
char cfg_infile_[MAXCHAR]="";	/* File name of input cfg data */
char cfg_outfile_[MAXCHAR]="";	/* File name of output cfg data (end-of-run) */
char cfg_snapfile_[MAXCHAR]="";	/* File name of output cfg data (in-run) */
char cfg_snapfile_0[MAXCHAR]="";
short cfg_ionctrl_=0;	    /* Control parameter:
			     * 0 == do nothing (default).
			     * 1 == introduce the bombarding ion.
			     */
int RDF_Hist_[MAXBIN];
int RDF_Hist33_[MAXBIN];
int RDF_Hist34_[MAXBIN];
int RDF_Hist43_[MAXBIN];
int RDF_Hist44_[MAXBIN];
int rdf_start_=0;	    /* time step at which to start rdf computation */
short RDF_=0;		    /* flag */
short RDF_34_=0;	    /* flag: compute 33, 34, 43, and 44 rdfs for C */
double rdf_loZlim_=-99.99;
double RDF_dR_=0.1;
char RDF_filename_[]="rdf.dat";
int RDF_nAtoms_zoned_=0, RDF_nAtoms_zoned3_=0, RDF_nAtoms_zoned4_=0;
int RDF_hist_accesses_=0;
double RDF_Nrho_=0.0, RDF_Nrho3_=0.0, RDF_Nrho4_=0.0;
int RDF_nBins_=0;
double rdf_zmin_, rdf_zmax_;

/* Configuration coordinate space information */
pt Lr_={0, 0, 0};	    /* Box dimensions */
pt half_Lr_={0, 0, 0};	    /* 1/2 box dimensions */
ipt nCell_={0, 0, 0};	    /* # of crystal unit cells in each dimension */
char unitCell_[]="";	    /* name of crystal unit cell; see cryst.c */
ipt per_={1, 1, 0};	    /* Periodic boundary condition switches;
			     * default is X, Y, !Z */

/* Externals:  periodic table owned by chem module */
extern element per_table[];
extern unit_type U_;	    /* specifies unit system.  default is
			     * Angstrom(length)/picosecond(time)/
			     * eV(energy)/Kelvin(temperature) */
extern char * UnitTypeLabels[];
extern short pe_set_;
extern short ke_set_;
extern short quiet_, echo_;
/* Diagnostics */
int CFG_DIAG_=0;

#define JT_LENNORM 1.54
#define JT_TIMENORM 0.020047164 /* ps */

static double zm;
static atomPtr a;
static int i, j;
static int datsize_=0;
void cfg_showsizes (void)
{
    fprintf(stdout, "# Data sizes:\n");
    fprintf(stdout, "#  atoms:\t%i\tx\t%i\t\t\t= %i\n", 
	sizeof(atom), MAXNUMATOMS, MAXNUMATOMS*sizeof(atom));
    datsize_+=MAXNUMATOMS*sizeof(atom);
    fprintf(stdout, "# nNodes:\t%i\tx\t%i\tx\t%i\t= %i\n", 
	sizeof(nNode), MAXNUMATOMS, MAXNUMNEIGHBORS, 
	sizeof(nNode)*MAXNUMATOMS*MAXNUMNEIGHBORS);
    datsize_+=sizeof(nNode)*MAXNUMATOMS*MAXNUMNEIGHBORS;
    fprintf(stdout, "# cNodes:\t%i\tx\t%i\t\t\t= %i\n", 
	sizeof(cNode), MAXNUMCLUSTERS, sizeof(cNode)*MAXNUMCLUSTERS);
    datsize_+=sizeof(cNode)*MAXNUMCLUSTERS;
    fprintf(stdout, "# total: %i\n", datsize_);
}

void cfg_initialize (void)
{
    void cfg_setup (void);
    
    if (!quiet_&&echo_) printf("# md series 2 cfg initializer (c) 1999 cfa\n");
    
    cfg_setup();
    cfg_clear();
}

void cfg_setup (void)
{
    for (i=0;i<MAXNUMATOMS;i++)
    {
	atomarray_[i].pos=&(pos_[i]);
	atomarray_[i].vel=&(vel_[i]);
	atomarray_[i].hac=&(hac_[i]);
	atomarray_[i].frc=&(frc_[i]);
	atomarray_[i].bCards=&(bcardarray_[i][0]);
    }
    cfg_initialized_=1;
}

void cfg_clear (void)
{
    if (!cfg_initialized_) 
	{printf("error: cfg not initialized (clear)\n");exit(0);}
    if (!quiet_&&echo_) 
	{printf("# md series 2 cfg clearer (c) 1999 cfa\n");fflush(stdout);}
    for (i=0;i<MAXNUMATOMS;i++)
    {
	ptPtr_clear(atomarray_[i].pos);
	ptPtr_clear(atomarray_[i].vel);
	ptPtr_clear(atomarray_[i].hac);
	ptPtr_clear(atomarray_[i].frc);
	atomarray_[i].sym=XX;
	atomarray_[i].next=NULL;
	atomarray_[i].c=NULL;
	atomarray_[i].q=atomarray_[i].mv2=0.0;
	for (j=0;j<MAXNUMNEIGHBORS;j++) 
	    bcardarray_[i][j].addr=&(atomarray_[i]);
	for (j=0;j<MAXNUMNEIGHBORS-1;j++) 
	    bcardarray_[i][j].next=&(bcardarray_[i][j+1]);
	bcardarray_[i][MAXNUMNEIGHBORS-1].next=NULL;
	atomarray_[i].bCards=&(bcardarray_[i][0]);
	atomarray_[i].nList=NULL;
	atomarray_[i].nNbrs=0;
    }
    cfg_established_=0;
    pe_set_=0;
    ke_set_=0;
    L_=&(atomarray_[0]);
    TH_=NULL;
    nTrash_=0;
    if (!quiet_&&echo_) 
	{printf("# md series 2 cfg clearer complete\n");fflush(stdout);}
}

double cfg_v (void)
{
    /* Compute the volume of the configuration */
    /* find the z-position of highest atom in cfg */
    zm=-1.e9;
    for (a=L_;a;a=a->next) if (a->pos->z>zm) zm=a->pos->z;
    /* Zspan is (z-pos of highest atom)+(1/2 box height) */
    return (Lr_.x*Lr_.y*(zm+half_Lr_.z));
}

void cfg_nreset (void)
{
    for (i=0;i<MAXNUMATOMS;i++) atom_nreset(&(atomarray_[i]));
}
void cfg_creset (void)
{
    for (i=0;i<MAXNUMATOMS;i++) atomarray_[i].c=NULL;
}

extern short clust_report_;
void cfg_report (void)
{
    void cfg_out (FILE * fp);
    
    if (!quiet_&&echo_) 
    {
	printf("# md series 2 cfg reporter (%s) (c) 1999 cfa\n", 
	    cfg_outfile_[0]?cfg_outfile_:"no report");
	fflush(stdout);
    }
    if (!cfg_newcrystal_&&clust_report_) clust_clean();
    
    if (cfg_outfile_[0])
    {
	if (!strcmp(cfg_outfile_, "stdout")) cfg_fp=stdout;
	else cfg_fp=fopen(cfg_outfile_, "w");
	cfg_out(cfg_fp);
	fclose(cfg_fp);
    }
    
    if (!quiet_&&echo_) 
    {
	printf("# md series 2 cfg reporter (%s) complete.\n", 
	    cfg_outfile_[0]?cfg_outfile_:"no report");
	fflush(stdout);
    }
}

static const char * cfg_line_fmt_= 
    "%i\t%s\t%.16lf\t%.16lf\t%.16lf\t%.16lf\t%.16lf\t%.16lf\t%i";
void cfg_atomout (atomPtr p, FILE *fp)
{
    if (!fp) fp=stdout;
    fprintf(fp, cfg_line_fmt_, p->id, per_table[p->sym].sym, 
	    p->pos->x, p->pos->y, p->pos->z,
	    p->vel->x, p->vel->y, p->vel->z, (p->state==IS_FIXED));
    fprintf(fp,"\n");
}

void cfg_out (FILE *fp)
{
    time_t Time, *tp = &Time;

    a=NULL;
    i=0;
    *tp = time(NULL);
    strftime(cfgCreateDate_, MAXCHAR, "%H:%M%p;%a%d%b%Y", localtime(tp));
    
    if (!quiet_&&echo_) printf("# md series 2 cfg writer (c) 1999 cfa\n");
    if (!fp) fp=stdout;
    
    fprintf(fp,"%% md series 2 atomic configuration file (c) 1999 cfa\n");
    fprintf(fp,"%% Created by md2::cfg v. %.2lf b. %i\n", version_, build_);
    fprintf(fp,"%% '%%' or ';' in first column means line is a comment\n");
    fprintf(fp,"%% '#' in first column means line is a keyword/value pair\n");
    if (cfgName_[0] != '\0') fprintf(fp, "#ConfigName\t%s\n", cfgName_);
    if (cfgCreateDate_[0] != '\0')
	fprintf(fp, "#CreationDate\t%s\n", cfgCreateDate_);
    fprintf(fp, "#UnitSystem\t%s\n", UnitTypeLabels[U_]);
    fprintf(fp,"#BoxSize.xyz\t%.14lf %.14lf %.14lf\n", 
	Lr_.x, Lr_.y, Lr_.z);	
    fprintf(fp, "#nCell.xyz\t%i %i %i\n", nCell_.i, nCell_.j, nCell_.k);
    fprintf(fp, "#cfgTime\t%.3lf\n", cfgTime_);
    #ifdef CRYST_H
    fprintf(fp, "#UnitCell\t%s\n", Cryst_strUC());
    #endif
    fprintf(fp, "#Periodicities\t%i %i %i\n", per_.i, per_.j, per_.k);
    fprintf(fp, "%% Number of Atoms in this cfg = %i\n", atomList_sizeof(L_));
    fprintf(fp, "%% By elements: ");
    for (i=0;i<MAXNUMELEMENTS;i++) nElem_[i]=0;
    for (a=L_;a;a=a->next) nElem_[a->sym]++;
    for (i=0;i<MAXNUMELEMENTS;i++)
	if (nElem_[i]) fprintf(fp, "%s %i ", per_table[i].sym, nElem_[i]);
    fprintf(fp, "\n");
    fprintf(fp, "%% Number of Fixed Atoms = %i\n", nFixed_);
    fprintf(fp, "%% Number of Trashed Atoms (not included) = %i\n", nTrash_);
    for (i=0,a=L_;a;a=a->next) {a->id=i++;cfg_atomout(a,fp);}
}

void cfg_snapout (int idnum)
{
    cfg_fp=NULL;
    if (!cfg_snapfile_0[0]) strcpy(cfg_snapfile_0, cfg_snapfile_);
    strcpy(cfg_snapfile_, cfg_snapfile_0);
    str_q2n(cfg_snapfile_, idnum);
    str_p2n(cfg_snapfile_, idnum);
    cfg_fp=fopen(cfg_snapfile_, "w");
    printf("# Created %s\n", cfg_snapfile_);fflush(stdout);
    if (cfg_fp)
    {
        if (!idnum) fprintf(cfg_fp, "# Initial configuration\n");
	cfg_out(cfg_fp);
	fclose(cfg_fp);
    }
    else {printf("error cfg_snapout\n");exit(0);}
}

extern double bhb_Tset_;
double Eb_, td_tau_=0.001, td_A_=1.e12;
void cfg_establish (void)
{
    void cfg_scan(FILE * fp);
    cfg_fp=NULL;
    
    if (!quiet_&&echo_) 
    {printf("# md series 2 cfg establisher (c) 1999 cfa\n");fflush(stdout);}
    
    if (cfg_established_) {printf("error: cfg already established\n");exit(0);}
    
    /* If the unit cell and number of unit cells and periodicities 
     * have been set, the user wishes to create a new crystal configuration.
     * Do this, and then exit. */
    if (unitCell_[0])
    {
	uc_a2i(unitCell_);
	uc_initialize();
	nAtom_=cryst_cnt();
	if (!quiet_) printf("# Creating new crystal: %s (%i, %i, %i) %i\n", 
	    unitCell_, nCell_.i, nCell_.j, nCell_.k, nAtom_);

	cryst_pos();
	cryst_vel();
	
	cfg_newcrystal_=1;
    }
    
    /* If the input file for the cfg is named, then read from it */
    else if (cfg_infile_[0])
    {
	if (!strcmp(cfg_infile_, "-"))
	    cfg_fp=stdin;
	else cfg_fp=fopen(cfg_infile_, "r");
	if (!cfg_fp) {printf("# Error: Cannot open input cfg %s\n", 
	    cfg_infile_[0]?cfg_infile_:"(null)");exit(0);}
	cfg_scan(cfg_fp);
	if (cfg_fp!=stdin) fclose(cfg_fp);
    }
    
    if (!Lr_.x||!Lr_.y||!Lr_.z) {printf("0001 error: bad boxsize\n");exit(0);}
    ptPtr_scalmult(&half_Lr_, &Lr_, 0.5);
    
    /* Compute thermal binding energy */
    Eb_=KB_EVPERK*bhb_Tset_*log(td_tau_*td_A_);
    if (U_==STILLWEB) Eb_/=EV_PER_EPSILON;
    if (Eb_<0) Eb_=0.0; /* if the log gives a negative binding energy,
			 * we assume the user means that nothing thermally
			 * desorbs */
    /* Note: by convention, we use Eb_ as a positive quantity.
     * If it is set to 0, only clusters with no binding energy
     * (or positive binding energy) are desorbed.  This situation
     * would correspond to physical sputtering if the run time
     * is kept sufficiently short. */
    
    if (!quiet_) printf("# Eb_(%.2lf) = %.5le %s\n", bhb_Tset_, Eb_, 
	(U_==STILLWEB?"ep":"eV"));
    fflush(stdout);
    
    /* Introduce a new incident ion to the cfg */
    switch (cfg_ionctrl_)
    {
	case 0:	    break;
	case 1:	    ion_newion();
		    break;
    }

    cfg_established_=1;
    if (!quiet_&&echo_) 
	{printf("# md series 2 cfg establisher complete\n");fflush(stdout);}
}

static char scr_ln[MAXCHAR], symStr[10], first_word[MAXCHAR], fv[MAXCHAR], *cp;
#define MAXARG 100
static char argv_BUF[MAXARG][MAXCHAR];
static char * argv_buf[MAXARG];
static int argc_buf;
void cfg_scan (FILE * fp)
{
    int fixed = 0, active = 0;
    atomPtr p = NULL;
    char c;
    double lowZ = 1.0e9, hiZ = -1.0e9;
    double lowFZ = 1.0e9, hiFZ = -1.0e9;
    enum { cfa_fmt, bah_fmt, cfa_fmt_pe, bah_fmt_pe, jms_fmt, gamess_fmt};
    #define NFMTS  6
    int dfmt = cfa_fmt;
    #if DIAG
    char * msg_hdr = "cfg::AtomList_ScanCfg";
    char * d_h = msg_hdr;
    #endif
    int idum=0;
    char tStr[30];
    #if DIAG
    const char * cfg_fmtStr[NFMTS] = 
    {
	"cfa", "bah", "cfa-pe", "bah-pe", "jms", "gamess"
    };
    #endif
    void cfg_scanGamessInp(FILE * fp);
    
    i=j=0;
    if (!cfg_initialized_) {printf("error: CFG not initialized.\n");exit(0);}

    if (!quiet_&&echo_) {printf("# md series 2 cfg scanner (c) 1999 cfa\n");fflush(stdout);}

    for (i=0;i<MAXARG;i++) argv_buf[i]=&(argv_BUF[i][0]);

    KE_=T_=0.0;
    for (i=0;i<MAXNUMELEMENTS;i++) nElem_[i]=0;

    first_word[0] = '\0';
    fv[0] = '\0';

    #if DIAG
    if (CFG_DIAG_)
    {
	printf("%s Checking CFG format...\n", d_h);
	fflush(stdout);
    }
    #endif

    /* Determine data file format:
     * -> NOTE: if input stream is stdin, 
     * then the required format is cfa --> bypass format
     * check.
     *
     *	Scan for '@@' and 'jj' and '!'.
     *	If '@@' is found, it is bah_fmt.
     *	If 'jj' is found, it is jms_fmt.
     *	If '!' is found as first character,  it is gamess_fmt.
     *  Otherwise, and by default, it is cfa_fmt.
     */
    if (fp!=stdin)
    {
	i = 1;
	while (i && (c = fgetc(fp)) != 255 
	    #if  defined(_HPUX_SOURCE) || defined(__USE_BSD)
	    && c != EOF 
	    #endif 
	)
	{
	    if (c == '@')
	    {
		c = fgetc(fp);
		if (c == '@')
		{   /* '@@' is detected!  File is Bryan's format. */
		    dfmt = bah_fmt;
		    i = 0;
		}
	    }
	    if (c == 'j')
	    {
		c = fgetc(fp);
		if (c == 'j')
		{   /* 'jj' is detected!  File is Jeanette's format. */
		    dfmt = jms_fmt;
		    i = 0;
		}
	    }
	    if (c == '!')
	    {
		dfmt = gamess_fmt;
		i = 0;
	    }
		// AAL: Added this to stop reading
        if( feof(fp) ) {
           break ;
        }
	}
	#if DIAG
	if (CFG_DIAG_)
	{
	    printf("Finished checking CFG format.\n");
	    printf("# cfgFormat is %s\n", cfg_fmtStr[dfmt]);fflush(stdout);
	}
	#endif
    
	i = 1;
	nAtom_ = nFixed_ = 0;
	lowZ = 1.e3;
	hiZ = -1.e3;
	if (dfmt == bah_fmt)
	{
	    #if DIAG
	    if (CFG_DIAG_)
		printf("%s Determined CFG file is Bryan's format.\n", d_h);
	    #endif
	    rewind(fp);
	    fgets(scr_ln, 255, fp);
	    sscanf(scr_ln, "%s %s", first_word, fv);
	    while (strcmp(first_word, "@@"))
	    {
		if (!strcmp(first_word, "XBOX:"))
		    Lr_.x = JT_LENNORM*atof(fv);
		else if (!strcmp(first_word, "YBOX:"))
		    Lr_.y = JT_LENNORM*atof(fv);
		else if (!strcmp(first_word, "ZBOX:"))
		    Lr_.z = JT_LENNORM*atof(fv);
	    
		fgets(scr_ln, 255, fp);
		sscanf(scr_ln, "%s %s", first_word, fv);
	    }
	    ptPtr_scalmult(&half_Lr_, &Lr_, 0.5);
	    c = fgetc(fp);
	    while (!isdigit(c) && c != '-' && c != '.') c = fgetc(fp);
	    ungetc(c, fp);
	    /* Begin reading in cfg data in Bryan's format */
	    i=0;
	    p=NULL;
	    while (fgets(scr_ln, MAXCHAR, fp))
	    {
		if (p) 
		{
		    p->next=&(atomarray_[i]);
		    p=p->next;
		    p->next=NULL;
		}
		else p=&(atomarray_[0]);
    
		sscanf(scr_ln, "%i\t%i%lf\t%lf\t%lf\t%lf\t%lf\t%lf",
		    &(p->sym), &(active), 
		    &(p->pos->x), &(p->pos->y), &(p->pos->z),
		    &(p->vel->x), &(p->vel->y), &(p->vel->z));
		    
		ptPtr_scalmult(p->pos, p->pos, JT_LENNORM);
		ptPtr_scalmult(p->vel, p->vel, JT_LENNORM/JT_TIMENORM);
		
		p->mv2 = per_table[p->sym].mass*ptPtr_sqabs(p->vel);
    
		KE_ += 0.5*p->mv2;
		
		if (active) p->state = IS_BULK;
		else p->state = IS_FIXED;
    
		switch (p->sym)
		{
		    case 1:	p->sym = C; /* Junichi's '1' is my '0' */
			    break;
		    case 3: p->sym = F;
			    break;
		    case 4: p->sym = Cl;
			    break;
		    case 5: p->sym = Br;
			    break;
		    case 6: p->sym = Cl_S;
			    break;
		    default:break;
		}
		nElem_[p->sym]++;
		p->id = i++;
		if (i==MAXNUMATOMS) {printf("error: too many atoms\n");exit(0);}
		p->nList = NULL;
	    }
	    nAtom_ = i;
	    
	    #if DIAG
	    if (CFG_DIAG_)
	    {
	       printf("AtomList_ScanCfg: returning...\n");
	       fflush(stdout); 
	    }
	    #endif
	    nCell_.i = nCell_.i = nCell_.i = -1;
	    return;
	}
	else if (dfmt == jms_fmt)
	{
	    printf("Config file is Jeanette's format\n");
	    printf("This format not supported yet.  Sorry!\n");
	    exit(0);
	}
	else if (dfmt == gamess_fmt)
	{
	    #if DIAG
	    printf("# Config file is Gamess format\n");
	    #endif
	    cfg_scanGamessInp(fp);
	    return;
	}
	
	/* Rewind the file */
	rewind(fp);
    }
   /* Making it this far necessarily means no '@@' or 'jj' was detected in the
    * first sweep of the input file -- so it's my format. */
    
    /* Read in lines of data until EOF */
    i = 0;
    p=NULL;
    argc_buf=0;
    while (fgets(scr_ln, MAXCHAR, fp))
    {
	/* is this line actual data or a comment?
	 * if first non-whitespace character in scr_ln
	 * is not a number, '.', or '-', then assume
	 * the line is a comment.
	 */
	
	for (j = 0; j < MAXCHAR && isspace(scr_ln[j]); j++);
	if (!isdigit(scr_ln[j]) && scr_ln[j] != '.' && scr_ln[j] != '-')
	{
	    cp=scr_ln;
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
	else
	{ 
	    if (p) 
	    {
		p->next=&(atomarray_[i]);
		p=p->next;
		p->next=NULL;
	    }
	    else p=&(atomarray_[0]);

	    p->state=IS_BULK;	    /* default */
	    ptPtr_clear(p->frc);
	    ptPtr_clear(p->hac); 
	    fixed = 0;
	    
	    /* Here, we determine whether this is an md series 1
	     * style cfg file or an md series 2 style cfg file.
	     * The difference is evident if we examine the second
	     * field on the input line.  If this field is NOT numeric, 
	     * then it must be an atomic symbol, and therefore a 
	     * series 2 type cfg.  */
	    sscanf(scr_ln, "%i %s", &idum, tStr);
	    if (isalpha(tStr[0]))   /* Version 2 format atom line */
	    {
		#if DIAG
		printf("v.2 atom line %i %s", i, scr_ln);fflush(stdout);
		#endif
		sscanf(scr_ln, "%i %s %lf %lf %lf %lf %lf %lf %i", 
		    &idum, symStr,
		    &(p->pos->x), &(p->pos->y), &(p->pos->z),
		    &(p->vel->x), &(p->vel->y), &(p->vel->z),
		    &(fixed));
		if (idum!=i) printf("Warning; cfg atoms not numbered sequentially.\n");
		p->id=i;
		p->sym=Chem_PeriodicTable_SymOfSym(symStr);
	    }
	    else		    /* Version 1 format atom line */
	    {
		#if DIAG
		printf("v.1 atom line %i %s", i, scr_ln);fflush(stdout);
		#endif
		sscanf(scr_ln, "%lf %lf %lf %lf %lf %lf %i %i", 
		    &(p->pos->x), &(p->pos->y), &(p->pos->z),
		    &(p->vel->x), &(p->vel->y), &(p->vel->z),
		    &(p->sym), &(fixed));
		p->id=i;
	    }
	    p->nList = NULL;
	    
	    nElem_[p->sym]++;
    
	    if (fixed)
	    {
		p->state = IS_FIXED;
		ptPtr_clear(p->vel);
		nFixed_++;
		if (p->pos->z < lowFZ) lowFZ = p->pos->z;
		if (p->pos->z > hiFZ) hiFZ = p->pos->z;
	    }
	    p->mv2 = per_table[p->sym].mass*ptPtr_sqabs(p->vel);
	    KE_ += 0.5*p->mv2;

	    if (p->pos->z < lowZ) lowZ = p->pos->z;
	    if (p->pos->z > hiZ) 
	    {
		hiZ = p->pos->z;
		symZhi_=p->sym;
	    }
	    #if DIAG
	    if (CFG_DIAG_ > 1)
	    {
		atom_display(stdout, p);
	    }
	    #endif
	    
	    i++;
	    if (i==MAXNUMATOMS) {printf("error: too many atoms\n");exit(0);}
	}
    }
    nAtom_=i;
    if (!quiet_)
    {
	printf("# %s: %i atoms: ", 
	   (!strcmp(cfg_infile_,"-")?"stdin":cfg_infile_), nAtom_);
	for (i=0;i<MAXNUMELEMENTS;i++) 
	    if (nElem_[i]) printf("%s_%i ", per_table[i].sym, nElem_[i]);
	printf("\n");
	fflush(stdout);
    }
    Zhi_=hiZ;
    Zlo_=lowZ;
    
    /* Call the argument handler to handle parameters listed in
     * the config file. */
    #if DIAG
    printf("#calling in-cfgscan arg handler\n");fflush(stdout);
    echo_args(argc_buf, argv_buf);
    #endif
    arg_handler(argc_buf, argv_buf, IGNORE_BADWORDS);
    
    #ifdef CRYST_H
    Cryst_setUC();
    Cryst_AssignOffsets();
    #endif
    
    
    #if DIAG
    if (CFG_DIAG_)
    {
	if (CFG_DIAG_ > 1) printf("||\n");
	printf("%s nAtom_ = %i.\n", d_h, nAtom_);
    }
    #endif
    /* Use the computed KE_ and nAtom_ to calculate T_ */
    if (U_ == APVK) T_ = KE_/1.5/(nAtom_-nFixed_)/KB_EVPERK;
    else T_ = KE_/1.5/(nAtom_-nFixed_)/KB;

    /* blank the RDF histograms */
    for (i=0;i<MAXBIN;i++) 
	RDF_Hist_[i]=RDF_Hist34_[i]=
	RDF_Hist43_[i]=RDF_Hist44_[i]=0;

    #if DIAG
    if (CFG_DIAG_)
    {
	printf("%s returning...\n", d_h);
	fflush(stdout);
    }
    #endif

}

void cfg_removeatom (atomPtr * a, int * i)
{
    if ((!a||!(*a))&&!i) return;
    if (!(*a)&&i) for ((*a)=L_;*a&&(*a)->id!=*i;(*a)=(*a)->next);
    if (!(*a)) {printf("error 1 in c_ra\n");exit(0);}
    atomList_remove((*a), &L_);
    atomList_push((*a), &TH_);
    nTrash_++;
}

atomPtr cfg_addr (int i)
{
    return (&(atomarray_[0])+i);
}

void cfg_scanGamessInp(FILE * fp)
{
    char atomname[30];
    double z, rx, ry, rz;
    atomPtr p=NULL;
    int i=0;
    
    Lr_.x=Lr_.y=Lr_.z=25.0;
    ptPtr_scalmult(&half_Lr_, &Lr_, 0.5);
    rewind(fp);
    fgets(scr_ln, 255, fp);
    sscanf(scr_ln, "%s", first_word);
    while (strcmp(first_word, "$DATA"))
    {
	fgets(scr_ln, 255, fp);
	sscanf(scr_ln, "%s", first_word);
    }
    fgets(scr_ln, 255, fp);
    fgets(scr_ln, 255, fp);
    fgets(scr_ln, 255, fp);
    while (strcmp(first_word, "$END"))
    {
	fgets(scr_ln, 255, fp);
	sscanf(scr_ln, "%s", first_word);
	sscanf(scr_ln, "%s %lf %lf %lf %lf", atomname, &z, &rx, &ry, &rz);
	if (p) 
	{
	    p->next=&(atomarray_[i]);
	    p=p->next;
	    p->next=NULL;
	}
	else p=&(atomarray_[0]);
	p->state = IS_BULK;
	p->pos->x=rx;
	p->pos->y=ry;
	p->pos->z=rz;
	p->sym=Chem_PeriodicTable_SymOfZ((int)z);
	nElem_[p->sym]++;
	p->id = i++;
	if (i==MAXNUMATOMS) {printf("error: too many atoms\n");exit(0);}
	p->nList = NULL;
    }
    nAtom_ = i;
}    

atomPtr cfg_newatom (void)
{
    L_[nAtom_-1].next=&(L_[nAtom_]);
    #if DIAG
    printf("#new atom %i\n", nAtom_);
    #endif
    /* Introducing a new atom -> new cfg -> current PE info no longer
     * valid -> turn off the pe_set_ lock. */
    pe_set_=0;
    L_[nAtom_].id=nAtom_;
    return &(L_[nAtom_++]);
}

void cfg_returnnewatom (void)
{
    L_[nAtom_-1].next=NULL;
    nAtom_--;
    pe_set_=0;
}

void RDF_normalize_output_hist()
{
    int bin, oldBinContents;
    FILE *fp;
    double cnst, cnst3, cnst4, ri;
    double rlower, rupper, nideal, nideal3, nideal4, scl;
    #if DIAG
    char * d_h = "cfg:RDF_normalize_output_hist:";
    #endif

    #if DIAG
    printf("RDF_nBins_ = %i\n", RDF_nBins_); fflush(stdout);
    printf("RDF_Nrho_ = %.5le\n", RDF_Nrho_); fflush(stdout);
    for (bin = 1; bin <= RDF_nBins_; bin++)
    {
	printf("%i\n", RDF_Hist_[bin]);
    }
    #endif
    
    fp = fopen(RDF_filename_,"w");
    
    RDF_nBins_ = (int)(half_Lr_.x * 1.0/RDF_dR_);
    if (RDF_nBins_>MAXBIN) RDF_nBins_=MAXBIN-1;
    
    fprintf(fp, "# RDF data: nBins = %i\n", RDF_nBins_);
    fprintf(fp, "# number = %i,  zone volume = %.5lf cu.Ang \n",
	RDF_nAtoms_zoned_, (rdf_zmax_-rdf_zmin_)*Lr_.x*Lr_.y);
    fprintf(fp, "# (zrange: %.5lf to %.5lf)\n",rdf_zmax_, rdf_zmin_);
    fprintf(fp, "# density = %.5lf #/cu.Ang.\n", RDF_Nrho_);
    fprintf(fp, "# Lower z limit of RDF zone = %.5lf\n", rdf_loZlim_);
    
    
    cnst=4.0/3.0*M_PI*RDF_Nrho_;
    if (RDF_34_) 
    {
	cnst3=4.0/3.0*M_PI*RDF_Nrho3_;
	cnst4=4.0/3.0*M_PI*RDF_Nrho4_;
    }
    scl=(double)RDF_hist_accesses_*(double)RDF_nAtoms_zoned_;
    #if DIAG
    printf("%s hist accesses: %i\t#atoms considered: %i\n",
	    d_h, RDF_hist_accesses_, RDF_nAtoms_zoned_);
    #endif
    for (bin=1;bin<=RDF_nBins_;bin++)
    {
	ri = (double)(2*bin-1)/2.0*RDF_dR_;
	rlower = (bin-1)*RDF_dR_;
	rupper = rlower+RDF_dR_;
	nideal = cnst*(pow(rupper,3)-pow(rlower,3));
	if (RDF_34_)
	{
	    nideal3=cnst3*(pow(rupper,3)-pow(rlower,3));
	    nideal4=cnst4*(pow(rupper,3)-pow(rlower,3));
	}
	#if DIAG
	printf("%.5le %.5le %.5le %.5le\n", cnst, rupper, rlower, nideal);
	#endif
	fprintf(fp, "%.10lf\t%.10lf\t%i", ri, 
	  ((double)RDF_Hist_[bin])/(nideal*scl),RDF_Hist_[bin]); 
	if (RDF_34_)
	{
	    fprintf(fp, "\t%.10lf\t%i\t%.10lf\t%i\t%.10lf\t%i",
	      ((double)RDF_Hist33_[bin])/(nideal3*scl),RDF_Hist33_[bin], 
	      ((double)RDF_Hist34_[bin])/(nideal4*scl)+
	      ((double)RDF_Hist43_[bin])/(nideal3*scl),
	      RDF_Hist34_[bin]+RDF_Hist43_[bin], 
	      ((double)RDF_Hist44_[bin])/(nideal4*scl),RDF_Hist44_[bin]);
	}
	fprintf(fp, "\n");
    }
    fprintf(fp, "# end of RDF data\n");
    fclose(fp);
}
