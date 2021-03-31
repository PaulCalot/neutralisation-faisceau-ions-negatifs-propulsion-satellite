/*
 * MD series 2 -- Molecular Dynamics of Plasma-Surface Chemistry
 * (c) 1999 Cam Abrams, Dept. of Chemical Engineering
 * Regents of the University of California, Berkeley
 * 
 * Module Name:  cfg2r3d:  creates a Raster3D-format datafile from
 * the supplied config. 
 */

#include <stdio.h>
#include "atom.h"
#include "chem.h"
#include "cfg.h"
#include "dblmat.h"
#include "colormap.h"
#include "point.h"
#include "r3d_utils.h"
#include "genforce.h"

#include "cam_colors.ext"

double version_=1.00;
int build_=1;

enum {NATURAL, DEPTH_RES, KE_RES, PE_RES, COORD, N_COLMODES};
char * color_modeLabels[N_COLMODES] =
{
    "Natural", "Depth-resolved", "KE-resolved", "PE-resolved", "Coord"
};
extern float hue_max_;	/* colormap.c */
extern float hue_min_;	/* colormap.c */

extern char * rot_label[];
extern	int	nAtom_;
extern	atomPtr L_, TH_;
extern	int	Zper_TOG_;
extern	pt	Lr_;
extern	pt	half_Lr_;
extern	unit_type   U_;
extern double T_, KE_,  PE_;
extern	element	per_table[];
extern char cfgName_[];
extern char cfg_infile_[];
double normMax_ = 0.0;
int c_flow = -1;
color_t c0, *cp0 = &c0;
color_t c1, *cp1 = &c1;
int color_mode = NATURAL;
double zmin = 0.0, zmax = 0.0;
double mv2max = 0.0, mv2min = 0.0;
double pemax = 0.0, pemin = 0.0;
int coolTransOn_ = 0;
int ion_id_ = -1;
short darken_fixed_=1;

/* render Header parameters and their default values */
extern double ImageBufferSize_;
extern int xTiles_, yTiles_;
extern int xPixPerTile_, yPixPerTile_;
extern int Scheme_;
extern color_t * bgClr;
color_t * ionClr;
extern int Shadows_;
extern int Phong_;
extern double SecLtContr_;
extern double AmbLtContr_;
extern double SpecRefComp_;
extern double Eyepos_;
extern const ptPtr MainLtPos;
extern short transparency_;
extern short RotateLightWithObject_;

extern short quiet_;

void usage (FILE * fp)
{
    if (!fp) return;
    fprintf(fp, "Usage: cfg2r3d <cfgName> [control params] [> <r3dName>]\n");
    fprintf(fp, "\t<cfgName>   = name of cfg data file\n");
    fprintf(fp, "\t<r3dName>   = name of r3d data file\n");
    fprintf(fp, "\t-df <directions output file>\n");
    fprintf(fp, "Rotation Controls: Values are in degrees\n");
    fprintf(fp, "\t -xr <x-rotation> -yr <y-rotation> -zr <z-rotation>\n");
    fprintf(fp, "Translation Controls: Values are in cfg-unit length\n");
    fprintf(fp, "\t -xtr <x-translation> -ytr <y-translation> -ztr <z-translation>\n");
    fprintf(fp, "Shifting Controls: Values are in cfg-unit length\n");
    fprintf(fp, "\t -xs <x-shift> -ys <y-shift> -zs <z-shift>\n");
    fprintf(fp, "(Shifting means that the actual atom coordinates are\n");
    fprintf(fp, "translated *prior* to rendering the image.)\n");
    fprintf(fp, "Zoom and Normalizing Length Controls:\n");
    fprintf(fp, "\t -zoom <scaleFactor> -norm <normMax>\n");
    fprintf(fp, "Colormode:  Natural, KE-resolved, or Depth-resolved\n");
    fprintf(fp, "\t -colmode <nat|ke|dep|pe>\n");
    fprintf(fp, "\t -kemin <kemin,eV> -kemax <kemax,eV> -cooltrans <YES|NO> for ke colormode only\n");
    fprintf(fp, "\t -pemin <pemin,eV> -pemax <pemax,eV> for pe colormode only\n");
    fprintf(fp, "\t -zmin <zmin> -zmax <zmax> for dep colormode only\n");
    fprintf(fp, "Object controls:\n");
    fprintf(fp, "\t -noatoms does not render atoms\n");
    fprintf(fp, "\t -ionid <ionIndex> colors ion\n");
    fprintf(fp, "\t -ioncolor <r g b>\n");
    fprintf(fp, "\t -pbc renders periodic boundary planes\n");
    fprintf(fp, "\t -pbcHt <float> height extent of vertical pbc planes\n");
    fprintf(fp, "\t -pbcClr <float> clarity of pbc planes\n");
    fprintf(fp, "\t -pbcColor <rgb or color key> color of pbc planes\n");
    fprintf(fp, "\t -pbcCtwy cuts the -y pbc plane off at x for z > 0.0\n");
    fprintf(fp, "\t -velarrfor <id-list> atoms for which velocity arrows are drawn\n");
    fprintf(fp, "\t -velarrlf <float> length of velocity arrows\n");
    fprintf(fp, "\t -velarrshaftrad <float> radius of velocity arrow shafts\n");
    fprintf(fp, "\t -velarrheadw <float> width of velocity arrow heads\n");
    fprintf(fp, "\t -velarrheadl <float> length of velocity arrow heads\n");
    fprintf(fp, "\t -velarrheadl <float> length of velocity arrow heads\n");
    fprintf(fp, "\t -velarrheadnface <int> number of faces of pyramidal arrow heads\n");
    fprintf(fp, "\t -hlids <id-list> list of id's for atoms that are to\n");
    fprintf(fp, "\t                 be highlighted.\n");
    fprintf(fp, "\t -hlclrs <color-namelist> list of colors for highlighted atoms\n");
    fprintf(fp, "\t -onlyelem <elem-str> only atoms of this element are rendered.\n");
    fprintf(fp, "Raster3D render Header Macros:\n");
    fprintf(fp, "\t -xt <xTiles(32)> -yt <yTiles(32)>\n");
    fprintf(fp, "\t -xp <xPixelsPerTile(18)> -yp <yPixelsPerTile(18)>\n");
    fprintf(fp, "\t -scheme <scheme(3)>\n");
    fprintf(fp, "\t -bg (colorName or <r> <g> <b> (0 0 0))\n");
    fprintf(fp, "\t -noshadows -phong <phongPower(25)>\n");
    fprintf(fp, "\t -slc <secondaryLightContributionFraction(0.15)>\n");
    fprintf(fp, "\t -alc <ambientLightContributionFraction(0.05>\n");
    fprintf(fp, "\t -src <specularReflectionComponentFraction(0.25)>\n");
    fprintf(fp, "\t -eyepos <EYEPOS(4.0)>\n");
    fprintf(fp, "\t -mlpos <x> <y> <z> (1 1 1) (Main Light Position)\n");
    fprintf(fp, "\t -orbit or -oc <x> <y> <z> Orbit Center\n");
    fprintf(fp, "Orbit means rotate the light source direction with any\n");
    fprintf(fp, "object rotation to simulate the camera orbiting a stationary\n");
    fprintf(fp, "object.  The (x, y, z) coordinates are the orbit center.\n");
    fprintf(fp, "The -orbit keyword assumes the orbit center is 0, 0, 0.\n");
    fprintf(fp, "The -oc keyword requires exactly 3 numbers to specify\n");
    fprintf(fp, "the orbit center if the user wishes it to be something\n");
    fprintf(fp, "other than (0, 0, 0).  Instead of a floating point value, \n");
    fprintf(fp, "you may specify 'X' for +Box.x/2 and 'x' for -Box.x/2, \n");
    fprintf(fp, "and likewise for <y> and <z> orbit center coords.\n");
}

#define MAXVELARRS 10

short elem_only_=0, belem_only_=0;
#define MAX_ELEMONLYS 5
sym_type elem_[MAX_ELEMONLYS] = {0, 0, 0, 0, 0};
sym_type belem_[MAX_ELEMONLYS] = {0, 0, 0, 0, 0};
int neo_, nbeo_;
int xexplode=0, yexplode=0;
int xexplfac=2, yexplfac=2;
double xexpldel=0.0, yexpldel=0.0;

void main (int argc, char * argv[])
{
    atomPtr L = NULL, a = NULL, hL = NULL, b = NULL;
    FILE * fp = NULL;
    pt Shift ={0.0, 0.0, 0.0}, *shift=&Shift;
    pt  per_shift = {0, 0, 0};
    pt  orbit_center = {0, 0, 0},  * oc = &orbit_center;
    int oc_Lr[3] = {0, 0, 0};
    pt TMP_PT={0, 0, 0}, *tmppt=&TMP_PT;
    char cfn[50];
    char * d_h = "cfg2r3d::main";
    int i = 0,  j = 0;
    double scale = 1.0, span = 0.0, x, y, z;
    void r3d_atomballs (FILE * fp, atomPtr L, double span);
    void r3d_clratomballs (FILE * fp, atomPtr L, color_t **clrs, double span);
    void usage (FILE * fp);
    int rot_order[MAXROTS], nRots = 0;
    double rots[MAXROTS];
    short drawPBCs = 0, frntPbcCt = 0;
    double pbcHt = 0.0, pbcClrty = 0.6;
    ptPtr xyz = NULL, xyZ = NULL, xYz = NULL, xYZ = NULL, 
	  Xyz = NULL, XyZ = NULL, XYz = NULL, XYZ = NULL;
    ptPtr va_base = NULL, va_point = NULL;
    double va_lf = 0.01, va_length = 0.0, va_shaftrad = 0.01, 
	   va_headw = 0.03, va_headl = 0.05, rff2;
    int va_headfaces = 16;
    short noatoms = 0;
    int velArrForID[MAXVELARRS], nVelArr=0;
    color_t * velArrClr = &blue;
    color_t * pbcColor = &peagreen;
    color_t * hlclrs[500];
    int	hlIDs[500], nhl=0, nhlc=0;
    char * dof = NULL;
    short setforce=1;
    
    color_init(cp0, hue_min_, hue_min_, hue_max_);
    color_init(cp1, hue_max_, hue_min_, hue_min_);

    /* Defaults */
    ionClr = &red;
    
    quiet_=1;
    Chem_InitializePeriodicTable();
    pef_initialize();
    cfg_initialize();
    
    
    for (i = 0; i < MAXROTS; i++) rot_order[i] = 0;
    for (i = 0; i < MAXVELARRS; i++) velArrForID[i] = 0;
    for (i = 0; i < 500; i++)
    {
	hlIDs[i] = 0;
	hlclrs[i] = NULL;
    }
    if (argc < 2)
    {
	usage(stdout);
	exit(-1);
    }
    j = 0;
    nRots = 0;
    for (i = 2; i < argc; i++)
    {
	if (!strcmp(argv[i], "-xr"))
	{
	    if (nRots == MAXROTS)
	    {
		printf("Error -- too many rotations requested.\n");exit(0);
	    }
	    else
	    {
		rots[nRots] = atof(argv[++i]);
		rot_order[nRots] = X;
		nRots++;
	    }
	}
	else if (!strcmp(argv[i], "-yr"))
	{
	    if (nRots == MAXROTS)
	    {
		printf("Error -- too many rotations requested.\n");exit(0);
	    }
	    else
	    {
		rots[nRots] = atof(argv[++i]);
		rot_order[nRots] = Y;
		nRots++;
	    }
	}
	else if (!strcmp(argv[i], "-zr"))
	{
	    if (nRots == MAXROTS)
	    {
		printf("Error -- too many rotations requested.\n");exit(0);
	    }
	    else
	    {
		rots[nRots] = atof(argv[++i]);
		rot_order[nRots] = Z;
		nRots++;
	    }
	}
	else if (!strcmp(argv[i], "-velarrfor"))
	{
	    nVelArr = 0;
	    i++;
	    for (; i < argc && isdigit(argv[i][0]); i++)
		velArrForID[nVelArr++] = atoi(argv[i]);
	    i--;
	}
	else if (!strcmp(argv[i], "-hlids"))
	{
	    nhl = 0;
	    i++;
	    while (i<argc&&isdigit(argv[i][0]))
		hlIDs[nhl++] = atoi(argv[i++]);
	    i--;
	}
	else if (!strcmp(argv[i], "-hlclrs"))
	{
	    nhlc = 0;
	    i++;
	    for (; i < argc && argv[i][0] != '-'; i++)
	    {
		if (isdigit(argv[i][0]))
		    color_init(hlclrs[j++], 
				   atof(argv[i]), 
				   atof(argv[++i]), 
				   atof(argv[++i]));
		else
		    hlclrs[nhlc++] = colorPtr(argv[i]);
	    }
	    i--;
	}
	else if (!strcmp(argv[i], "-df")) dof = argv[++i];
	else if (!strcmp(argv[i], "-xtr")) shift->x = atof(argv[++i]);
	else if (!strcmp(argv[i], "-ytr")) shift->y = atof(argv[++i]);
	else if (!strcmp(argv[i], "-ztr")) shift->z = atof(argv[++i]);
	else if (!strcmp(argv[i], "-xs")) per_shift.x = atof(argv[++i]);
	else if (!strcmp(argv[i], "-ys")) per_shift.y = atof(argv[++i]);
	else if (!strcmp(argv[i], "-zs")) per_shift.z = atof(argv[++i]);
 	else if (!strcmp(argv[i], "-zoom")) scale = atof(argv[++i]);
 	else if (!strcmp(argv[i], "-norm")) normMax_ = atof(argv[++i]);
 	else if (!strcmp(argv[i], "-zmin")) zmin = atof(argv[++i]);
 	else if (!strcmp(argv[i], "-zmax")) zmax = atof(argv[++i]);
 	else if (!strcmp(argv[i], "-kemin")) mv2min = 2*atof(argv[++i]);
 	else if (!strcmp(argv[i], "-kemax")) mv2max = 2*atof(argv[++i]);
 	else if (!strcmp(argv[i], "-pemin")) pemin = atof(argv[++i]);
 	else if (!strcmp(argv[i], "-pemax")) pemax = atof(argv[++i]);
 	else if (!strcmp(argv[i], "-xt")) xTiles_ = atoi(argv[++i]);
 	else if (!strcmp(argv[i], "-yt")) yTiles_ = atoi(argv[++i]);
 	else if (!strcmp(argv[i], "-xp")) xPixPerTile_ = atoi(argv[++i]);
 	else if (!strcmp(argv[i], "-yp")) yPixPerTile_ = atoi(argv[++i]);
 	else if (!strcmp(argv[i], "-scheme")) Scheme_ = atoi(argv[++i]);
 	else if (!strcmp(argv[i], "-phong")) Phong_ = atof(argv[++i]);
 	else if (!strcmp(argv[i], "-slc")) SecLtContr_ = atof(argv[++i]);
 	else if (!strcmp(argv[i], "-alc")) AmbLtContr_ = atof(argv[++i]);
 	else if (!strcmp(argv[i], "-src")) SpecRefComp_ = atof(argv[++i]);
 	else if (!strcmp(argv[i], "-eyepos")) Eyepos_ = atof(argv[++i]);
 	else if (!strcmp(argv[i], "-ibs")) ImageBufferSize_ = atof(argv[++i]);
 	else if (!strcmp(argv[i], "-ionid")) ion_id_ = atoi(argv[++i]);
 	else if (!strcmp(argv[i], "-noatoms")) noatoms = 1;
 	else if (!strcmp(argv[i], "-nof")) setforce = 0;
 	else if (!strcmp(argv[i], "+quiet")) quiet_ = 0;
 	else if (!strcmp(argv[i], "-nodarkfix")) darken_fixed_ = 0;
	else if (!strcmp(argv[i], "-xexplode")) xexplode=atoi(argv[++i]);
	else if (!strcmp(argv[i], "-yexplode")) yexplode=atoi(argv[++i]);
	else if (!strcmp(argv[i], "-xexplodefac")) xexplfac=atoi(argv[++i]);
	else if (!strcmp(argv[i], "-yexplodefac")) yexplfac=atoi(argv[++i]);
 	else if (!strcmp(argv[i], "-elemonly"))
	{
	    elem_only_ = 1;
	    neo_ = 0;
	    i++;
	    while (i < argc && argv[i][0] != '-')
		elem_[neo_++] = Chem_PeriodicTable_SymOfSym(argv[i++]);
	    i--;
	}
 	else if (!strcmp(argv[i], "-belemonly"))
	{
	    belem_only_ = 1;
	    nbeo_ = 0;
	    i++;
	    while (i < argc && argv[i][0] != '-')
		belem_[nbeo_++] = Chem_PeriodicTable_SymOfSym(argv[i++]);
	    i--;
	}
 	else if (!strcmp(argv[i], "-pbc")) drawPBCs = 1;
 	else if (!strcmp(argv[i], "-pbcHt")) 
	{
	    drawPBCs = 1;
	    pbcHt = atof(argv[++i]);
 	}
	else if (!strcmp(argv[i], "-pbcClr")) 
	{
	    pbcClrty = atof(argv[++i]);
 	    drawPBCs = 1;
	}
	else if (!strcmp(argv[i], "-pbcCtwy")) frntPbcCt = 1;
 	else if (!strcmp(argv[i], "-velarrlf")) va_lf = atof(argv[++i]);
 	else if (!strcmp(argv[i], "-velarrshaftrad")) va_shaftrad = atof(argv[++i]);
 	else if (!strcmp(argv[i], "-velarrheadw")) va_headw = atof(argv[++i]);
 	else if (!strcmp(argv[i], "-velarrheadl")) va_headl = atof(argv[++i]);
 	else if (!strcmp(argv[i], "-velarrheadnface")) va_headfaces = atoi(argv[++i]);
 	else if (!strcmp(argv[i], "-pbcColor"))
	{
	    if (isdigit(argv[++i][0]))
	    {
		pbcColor = color_init(pbcColor, atof(argv[i]), atof(argv[++i]), atof(argv[++i]));
	    }
	    else pbcColor = colorPtr(argv[i]);
	}
 	else if (!strcmp(argv[i], "-ioncolor"))
	{
	    if (isdigit(argv[++i][0]))
	    {
		ionClr = color_init(ionClr, atof(argv[i]), atof(argv[++i]), atof(argv[++i]));
	    }
	    else ionClr = colorPtr(argv[i]);
	}
 	else if (!strcmp(argv[i], "-velarrcolor"))
	{
	    if (isdigit(argv[++i][0]))
	    {
		velArrClr = color_init(velArrClr, atof(argv[i]), atof(argv[++i]), atof(argv[++i]));
	    }
	    else velArrClr = colorPtr(argv[i]);
	}
 	else if (!strcmp(argv[i], "-bg"))
	{
	    if (isdigit(argv[++i][0]))
	    {
		bgClr = color_init(bgClr, atof(argv[i]), atof(argv[++i]), atof(argv[++i]));
	    }
	    else bgClr = colorPtr(argv[i]);
	}
 	else if (!strcmp(argv[i], "-mlpos"))
	{
	    x = atof(argv[++i]);
	    y = atof(argv[++i]);
	    z = atof(argv[++i]);
	    MainLtPos->x = x;
	    MainLtPos->y = y;
	    MainLtPos->z = z;
	}
 	else if (!strcmp(argv[i], "-noshadows")) 
	{
	    Shadows_ = 0;
	}
 	else if (!strcmp(argv[i], "-orbit"))
	{
	    RotateLightWithObject_ = 1;
	}
 	else if (!strcmp(argv[i], "-oc"))
	{
	    RotateLightWithObject_ = 1;
	    i++;
	    if (argv[i][0] == 'X')
	    {
		oc_Lr[0] = 1;
		sscanf(&(argv[i][1]), "%lf", &(oc->x));
	    }
	    else if (argv[i][0] == 'x')
	    {
		oc_Lr[0] = -1;
		sscanf(&(argv[i][1]), "%lf", &(oc->x));
	    }
	    else orbit_center.x = atof(argv[i]);
	    i++;
	    if (argv[i][0] == 'Y')
	    {
		oc_Lr[1] = 1;
		sscanf(&(argv[i][1]), "%lf", &(oc->y));
	    }
	    else if (argv[i][0] == 'y')
	    {
		oc_Lr[1] = -1;
		sscanf(&(argv[i][1]), "%lf", &(oc->y));
	    }
	    else orbit_center.y = atof(argv[i]);
	    i++;
	    if (argv[i][0] == 'Z')
	    {
		oc_Lr[2] = 1;
		sscanf(&(argv[i][1]), "%lf", &(oc->z));
	    }
	    else if (argv[i][0] == 'z')
	    {
		oc_Lr[2] = -1;
		sscanf(&(argv[i][1]), "%lf", &(oc->z));
	    }
	    else orbit_center.z = atof(argv[i]);
	}
 	else if (!strcmp(argv[i], "-cooltrans"))
	{
	    coolTransOn_ = 1;
	}
 	else if (!strcmp(argv[i], "-colmode"))
	{
	    i++;
	    if (!strcmp(argv[i], "ke")) color_mode = KE_RES;
	    else if (!strcmp(argv[i], "dep")) color_mode = DEPTH_RES;
	    else if (!strcmp(argv[i], "pe")) color_mode = PE_RES;
	    else if (!strcmp(argv[i], "coord")) color_mode = COORD;
	    else color_mode = NATURAL;
	}
	else
	{
	    printf("Error:  Keyword %s not recognized.\n", argv[i]);
	    usage(stdout);
	    exit(-1);
	}
    }
    
    i = 0;
    while (i < 500 && hlclrs[i]) i++;
    if (i == 0) hlclrs[i++] = &white;
    for (j = i; j < 500; j++) hlclrs[j] = hlclrs[i-1];
	
    strcpy(cfg_infile_, argv[1]);
    fp=NULL;
    if ((!strcmp(cfg_infile_, "-")) || (fp=fopen(cfg_infile_, "r")))
    {
	if (fp) fclose(fp);
	cfg_establish();
	if (setforce) 
	{
	    if (!quiet_) printf("# force calculation...\n");
	    cfg_setforce();
	    if (!quiet_) printf("# total E=%.5lf\n", PE_+KE_);
	}
	if (oc_Lr[0]) orbit_center.x += oc_Lr[0]*half_Lr_.x;
	if (oc_Lr[1]) orbit_center.y += oc_Lr[1]*half_Lr_.y;
	if (oc_Lr[2]) orbit_center.z += oc_Lr[2]*half_Lr_.z;
	if (drawPBCs && pbcHt == 0.0) pbcHt = Lr_.z;
	for (a = L_; a; a = a->next)
	    ptPtr_minimg(a->pos, ptPtr_add(a->pos, a->pos, &per_shift));
	for (a = L_; a; a = a->next) 
	    ptPtr_add(a->pos, ptPtr_scalmult(&TMP_PT, oc, -1.0), a->pos);

	if (!normMax_)
	{
	    normMax_ = sqrt(ptPtr_sqabs(&Lr_));
	}
	span = normMax_;
	/* add the default buffer to the image span */
	span += ImageBufferSize_;
	/* find the max and min z-values */
	if ((!zmin || !zmax) && color_mode == DEPTH_RES) 
	    for (a = L_; a; a = a->next)
	    {
		if (a->pos->z < zmin) zmin = a->pos->z;
		if (a->pos->z > zmax) zmax = a->pos->z;
	    }

	/* find the max and min mv2's */
	if (!(mv2min && mv2max) && color_mode == KE_RES) 
	    for (a = L_; a; a = a->next)
	    {
		if (a->mv2 > mv2max) mv2max = a->mv2;
		if (a->mv2 < mv2min) mv2min = a->mv2;
	    }
	
	/* Pull the atoms to be highlighted out of L and push them
	 * onto TH_. */
	i = 0;
	while (hlIDs[i]) 
	{
	    for (a=L_;a&&a->id!=hlIDs[i];a=a->next);
	    if (a) cfg_removeatom(&a, &(hlIDs[i++]));
	}
	if (dof)
	{
	    fp = fopen(dof, "w");
	    fprintf(fp, "# %s direction file %s\n", argv[0], dof);
	    fprintf(fp, "cfg %s\n", argv[1]);
	    fprintf(fp, "norm %.10lf\n", normMax_);
	    for (i = 0; i < MAXROTS && rot_order[i] != 0; i++)
		fprintf(fp, "%s %.10lf\n", rot_label[rot_order[i]], rots[i]);
	    fprintf(fp, "zoom %.10lf\n", scale);
	    fprintf(fp, "oc.x %.10lf\n", oc->x);
	    fprintf(fp, "oc.y %.10lf\n", oc->y);
	    fprintf(fp, "oc.z %.10lf\n", oc->z);
	    fprintf(fp, "tr.x %.10lf\n", shift->x);
	    fprintf(fp, "tr.y %.10lf\n", shift->y);
	    fprintf(fp, "tr.z %.10lf\n", shift->z);
	    if (drawPBCs) fprintf(fp, "pbc\n");
	    fprintf(fp, "mlp.x %.10lf\n", MainLtPos->x);
	    fprintf(fp, "mlp.y %.10lf\n", MainLtPos->y);
	    fprintf(fp, "mlp.z %.10lf\n", MainLtPos->z);
	    fprintf(fp, "#end of direction file %s\n", dof);
	    fclose(fp);
	}
	
	
	r3d_header(stdout, rots, rot_order, scale, shift, span);
	printf("# cfg2r3d version %.2lf by cam abrams\n", version_);
	printf("#<cfgName> = %s\n", cfn);
        printf("#Input mass-velocities imply total KE of %.5e and T of %.3f\n", 
            KE_, T_);
        printf("#Box size is %.5le %.5le %.5le\n", Lr_.x, Lr_.y, Lr_.z);
	fflush(stdout);
	if (rot_order[0] != 0)
	{
	    printf("#Rotation order is: ");
	    for (i = 0; i < MAXROTS && rot_order[i] != 0; i++)
		printf("%s[%.2lf] ", rot_label[rot_order[i]], rots[i]);
	    printf("\n");
	}
	fflush(stdout);
	if (scale) printf("#Scale is %.3lf.\n", scale);
	fflush(stdout);
	printf("#nAtom = %i\n", nAtom_);
	fflush(stdout);
	printf("#Color mode is %s\n", color_modeLabels[color_mode]);
	switch(color_mode)
	{
	    case KE_RES:    printf("# ke-min is %.5lf\n", 0.5*mv2min);
			    printf("# ke-max is %.5lf\n", 0.5*mv2max);
			    break;
	    case DEPTH_RES: printf("# z-min is %.5lf\n", zmin);
			    printf("# z-max is %.5lf\n", zmax);
			    break;			    
	}
	printf("#span is %.5lf\n", span);
	fflush(stdout);
	if (xexplode) xexpldel=half_Lr_.x/xexplode;
	if (yexplode) yexpldel=half_Lr_.y/yexplode;
	if (!noatoms) 
	{
	    r3d_atomballs(stdout, L_, span);
	    r3d_clratomballs(stdout, TH_, hlclrs, span);
	}
    }
    
    if (drawPBCs&&!(xexplode||yexplode))
    {
	xyz = ptPtr_create((half_Lr_.x-oc->x)/span, 
			 (half_Lr_.y-oc->y)/span, 
			 (-half_Lr_.z+pbcHt-oc->z)/span);
	Xyz = ptPtr_create((-half_Lr_.x-oc->x)/span, 
			 (half_Lr_.y-oc->y)/span, 
			 (-half_Lr_.z+pbcHt-oc->z)/span);
	XYz = ptPtr_create((-half_Lr_.x-oc->x)/span, 
			(-half_Lr_.y-oc->y)/span, 
			((-half_Lr_.z+pbcHt-oc->z))/span);
	XyZ = ptPtr_create((-half_Lr_.x-oc->x)/span, 
			(half_Lr_.y-oc->y)/span, 
			(-half_Lr_.z-oc->z)/span);
	XYZ = ptPtr_create((-half_Lr_.x-oc->x)/span,
			(-half_Lr_.y-oc->y)/span, 
			(-half_Lr_.z-oc->z)/span);
	xYz = ptPtr_create((half_Lr_.x-oc->x)/span, 
			(-half_Lr_.y-oc->y)/span, 
			((-half_Lr_.z+pbcHt-oc->z))/span);
	xYZ = ptPtr_create((half_Lr_.x-oc->x)/span, 
			(-half_Lr_.y-oc->y)/span, 
			(-half_Lr_.z-oc->z)/span);
	xyZ = ptPtr_create((half_Lr_.x-oc->x)/span, 
			(half_Lr_.y-oc->y)/span, 
			(-half_Lr_.z-oc->z)/span);

	printf("# pbcs:\n");
	r3d_transparency_on(stdout, pbcClrty);
	r3d_tri(stdout, XYz, XYZ, XyZ, *pbcColor);
	r3d_tri(stdout, XYz, Xyz, XyZ, *pbcColor);
	r3d_tri(stdout, Xyz, xyz, xyZ, *pbcColor);
	r3d_tri(stdout, Xyz, XyZ, xyZ, *pbcColor);
	r3d_tri(stdout, xyz, xyZ, xYZ, *pbcColor);
	r3d_tri(stdout, xyz, xYz, xYZ, *pbcColor);
	if (frntPbcCt)
	{
	    xyz->z = 0.0;
	    Xyz->z = 0.0;
	    XYz->z = 0.0;
	    xYz->z = 0.0;
	}
	r3d_tri(stdout, xYz, xYZ, XYZ, *pbcColor);
	r3d_tri(stdout, xYz, XYz, XYZ, *pbcColor);
	
	r3d_transparency_off(stdout);
	r3d_tri(stdout, XYZ, XyZ, xyZ, *pbcColor);
	r3d_tri(stdout, XYZ, xYZ, xyZ, *pbcColor);
    }
    if (velArrForID[0])
    {
	for (i = 0; i < MAXVELARRS && velArrForID[i]; i++)
	{
	    a = cfg_addr(velArrForID[i]);
	    if (a)
	    {
		va_base = ptPtr_create(a->pos->x/span, a->pos->y/span, a->pos->z/span);
		va_length = sqrt(ptPtr_sqabs(a->vel))/span;
		printf("#va_length = %.5le\n", va_length);
		printf("#va_lf = %.5le\n", va_lf);
		va_point = ptPtr_create(a->vel->x, a->vel->y, a->vel->z);
		ptPtr_scalmult(va_point, va_point, va_lf);
		ptPtr_add(va_point, va_point, va_base);
		r3d_3darrow(stdout, va_base, va_point, va_shaftrad,
			    va_headl, va_headw, va_headfaces, *velArrClr);
	    }
	}
    }
    
    printf("#Program Ends.\n");
}


void r3d_atomballs (FILE * fp, atomPtr L, double span)
{
    void r3d_ball_recursive (FILE * fp, atomPtr L, double Max);
    if (!fp || !L) return;
    r3d_ball_recursive (fp, L, span);
}

void r3d_clratomballs (FILE * fp, atomPtr L, color_t **clrs, double span)
{
    void r3d_clrball_recursive (FILE * fp, atomPtr L, color_t **clrs, 
	double Max, int * cnt);
    int cnt = 0;
    if (!fp || !L) return;
    
    r3d_transparency_off(fp);
    r3d_clrball_recursive (fp, L, clrs, span, &cnt);
}

color_t c, *cp = &c;
double pe;
nNodePtr np_=NULL;
short renderme;
int xexplshift=0.0, yexplshift=0.0;
void r3d_ball_recursive (FILE * fp, atomPtr L, double Max)
{
    int nn;
    int k;
    if (!fp || !L) return;
    color_init(cp, per_table[L->sym].color.r, 
		        per_table[L->sym].color.g, 
			per_table[L->sym].color.b );
    if (color_mode == DEPTH_RES)
	color_mapFloat(zmin, cp0, zmax, cp1, L->pos->z, cp, c_flow);
    else if (color_mode == KE_RES)
    {
	color_mapFloat(mv2min, cp0, mv2max, cp1, L->mv2, cp, c_flow);
	if (coolTransOn_ && (L->mv2 < mv2min) && !(L->state == IS_FIXED))
	    r3d_transparency_on(fp, 0.8);
	else r3d_transparency_off(fp);
    }
    else if (color_mode == PE_RES)
    {
	pe=0.0;
	for (np_=L->nList;np_;np_=np_->next) pe+=0.5*atom_peij(L, np_->addr);
	color_mapFloat(pemin, cp0, pemax, cp1, pe, cp, c_flow);
    }
    else if (color_mode == COORD)
    {
	nn=nNode_n(L->nList);
	printf("# number of neighbors is %i\n", nn);fflush(stdout);
	if (nn==1) color_scalMultColor(cp, 0.01);
	if (nn==2) color_scalMultColor(cp, 0.1);
	if (nn==3) color_scalMultColor(cp, 0.3);
    }
    if (L->state == IS_FIXED)
    {
	if (darken_fixed_) color_scalMultColor(cp, 0.2);
//	printf("fixed-c\n");fflush(stdout);
    }
    if (L->id == ion_id_)
    {
	cp->r = ionClr->r;
	cp->g = ionClr->g;
	cp->b = ionClr->b;
    }
    for (k=0;k<MAX_ELEMONLYS&&L->sym!=elem_[k];k++);
    renderme=elem_only_&&k<MAX_ELEMONLYS;
    if (renderme&&elem_only_&&belem_only_)
    {
//	printf("? %s(%x)#%i -> %s\n", per_table[L->sym].sym, 
//	    (L->nList), L->id, per_table[belem_[k]].sym);
//	fflush(stdout);
	if (belem_[k])
	{
	    for (np_=L->nList;np_&&np_->addr->sym!=belem_[k];np_=np_->next)
	    {
//		printf("n%s#%i ", per_table[np_->addr->sym].sym, np_->addr->id);
//		fflush(stdout);
	    }
	    renderme=(np_&&np_->addr->sym==belem_[k]);
	}
    }
//    printf("? %s#%i %i\n", per_table[L->sym].sym, L->id, renderme);fflush(stdout);
    renderme=renderme||(!belem_only_&&!elem_only_);
    xexplshift=0.0;
    if (xexplode) 
	xexplshift=xexplfac*(((int)(L->pos->x/xexpldel))+0.5*(L->pos->x<0?-1:1))
		    *xexpldel;
    yexplshift=0.0;
    if (yexplode) 
	yexplshift=yexplfac*((int)(L->pos->y/yexpldel))+0.5*(L->pos->y<0?-1:1))
		    *yexpldel;
    if (renderme)
    {
	pe=0.0;
	for (np_=L->nList;np_;np_=np_->next) pe+=0.5*atom_peij(L, np_->addr);
	fprintf(fp, "# atom %s_{%i}, r = %.4lf, mv2 = %.5lf, pe = %.5lf, rgb=(%.3lf %.3lf %.3lf)\n", 
	    per_table[L->sym].sym, L->id, per_table[L->sym].radius, 
	    L->mv2, pe, cp->r, cp->g, cp->b);
	fprintf(fp, "2\n");
	fprintf(fp, "%.3lf %.3lf %.3lf %.3lf %.3lf %.3lf %.3lf\n", 
	    (L->pos->x+xexplshift)/Max, (L->pos->y+yexplshift)/Max, 
	    L->pos->z/Max, 
	    per_table[L->sym].radius/Max, 
	    cp->r, cp->g, cp->b);
    }
    r3d_ball_recursive (fp, L->next, Max);
}

void r3d_clrball_recursive (FILE * fp, atomPtr L, color_t **clrs, 
			    double Max, int * cnt)
{
    color_t * cp;
    if (!fp || !L) return;
    cp = clrs[*cnt];
    (*cnt)++;
    if (!elem_only_ || (elem_only_ && 
		       (L->sym == elem_[0] || L->sym == elem_[1] || 
			L->sym == elem_[2] || L->sym == elem_[3] || 
			L->sym == elem_[4])))
    {
	fprintf(fp, "# highlight atom %s_{%i}, r = %.4lf, mv2 = %.5lf\n", 
	    per_table[L->sym].sym, L->id, per_table[L->sym].radius, L->mv2);
	fprintf(fp, "2\n");
	fprintf(fp, "%.3lf %.3lf %.3lf %.3lf ", 
	    L->pos->x/Max, L->pos->y/Max, L->pos->z/Max, 
	    per_table[L->sym].radius/Max);
	fflush(fp);
	fprintf(fp, "%.3lf %.3lf %.3lf\n", 
	    cp->r, cp->g, cp->b);
	fflush(fp);
    }
    r3d_clrball_recursive (fp, L->next, clrs, Max, cnt);
}

