#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "dblmat.h"
#include "point.h"
#include "colormap.h"
#include "r3d_utils.h"
extern color_t black, white, red, blue, green, darkgold1, peagreen, 
               violet, violet2, bronze, darkpeagreen, 
		grey1, grey2, yellow;
extern char * rot_label[];
extern double ImageBufferSize_;
extern int xTiles_, yTiles_;
extern int xPixPerTile_, yPixPerTile_;
extern int Scheme_;
extern color_t bgClr;
extern int Shadows_;
extern int Phong_;
extern double SecLtContr_;
extern double AmbLtContr_;
extern double SpecRefComp_;
extern double Eyepos_;
extern ptPtr MainLtPos;

void usage (FILE * fp)
{
    if (!fp) return;
    fprintf(fp, "Usage: r3d_test [control params] [> <r3dName>]\n");
    fprintf(fp, "Raster3D render control parameters:\n");
    fprintf(fp, "\t [-xt <xTiles(32)>] [-yt <yTiles(32)>]\n");
    fprintf(fp, "\t [-xp <xPixelsPerTile(18)>] [-yp <yPixelsPerTile(18)>]\n");
    fprintf(fp, "\t [-scheme <scheme(3)>]\n");
    fprintf(fp, "\t [-bg <r> <g> <b> (0 0 0)]\n");
    fprintf(fp, "\t [-noshadows] [-phong <phongPower(25)>]\n");
    fprintf(fp, "\t [-slc <secondaryLightContributionFraction(0.15)>]\n");
    fprintf(fp, "\t [-alc <ambientLightContributionFraction(0.05>]\n");
    fprintf(fp, "\t [-src <specularReflectionComponentFraction(0.25)>]\n");
    fprintf(fp, "\t [-eyepos <EYEPOS(4.0)>]\n");
    fprintf(fp, "\t [-mlpos <x> <y> <z> (1 1 1)] (Main Light Position)\n");
}

void main (int argc, char * argv[])
{
    char * d_h = "ptd::main", cmt = '#';
    int i = 0,  j = 0;
    double scale = 1.0, span = 1.0;
    ptPtr shift = ptPtr_Init(0.0, 0.0, 0.0);
    int rot_order[MAXROTS];
    double rots[MAXROTS];
    ptPtr origin = ptPtr_Init(0.0, 0.0, 0.0);
    ptPtr rtp1 = ptPtr_Init(0.2, 45.0, 45.0);
    int arcsegs = 20;
    double liner = 0.01;
    double x = 0.0, y = 0.0, z = 0.0;

    for (i = 0; i < MAXROTS; i++) rots[i] = 0.0;
    for (i = 0; i < MAXROTS; i++) rot_order[i] = 0;
    for (i = 1; i < argc; i++)
    {
	if (!strcmp(argv[i], "-xr"))
	{
	    if (j == MAXROTS)
	    {
		printf("Error -- too many rotations requested.\n");
	    }
	    else
	    {
		rots[j] = atof(argv[++i]);
		rot_order[j] = X;
		j++;
	    }
	}
	else if (!strcmp(argv[i], "-yr"))
	{
	    if (j == MAXROTS)
	    {
		printf("Error -- too many rotations requested.\n");
	    }
	    else
	    {
		rots[j] = atof(argv[++i]);
		rot_order[j] = Y;
		j++;
	    }
	}
	else if (!strcmp(argv[i], "-zr"))
	{
	    if (j == MAXROTS)
	    {
		printf("Error -- too many rotations requested.\n");
	    }
	    else
	    {
		rots[j] = atof(argv[++i]);
		rot_order[j] = Z;
		j++;
	    }
	}
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
 	else if (!strcmp(argv[i], "-zoom")) scale = atof(argv[++i]);
 	else if (!strcmp(argv[i], "-arcsegs")) arcsegs = atoi(argv[++i]);
	else if (!strcmp(argv[i], "-liner")) liner = atof(argv[++i]);
 	else if (!strcmp(argv[i], "-bg"))
	{
	    bgClr.r = atof(argv[++i]);
	    bgClr.g = atof(argv[++i]);
	    bgClr.b = atof(argv[++i]);
	}
 	else if (!strcmp(argv[i], "-mlpos"))
	{
	    x = atof(argv[++i]);
	    y = atof(argv[++i]);
	    z = atof(argv[++i]);
	    MainLtPos = ptPtr_Init(x, y, z);
	}
 	else if (!strcmp(argv[i], "-trans"))
	{
	    shift->x = atof(argv[++i]);
	    shift->y = atof(argv[++i]);
	    shift->z = atof(argv[++i]);
	}
 	else if (!strcmp(argv[i], "-noshadows")) 
	{
	    Shadows_ = 0;
	}
	else
	{
	    printf("Error:  Keyword %s not recognized.\n", argv[i]);
	    usage(stdout);
	    exit(-1);
	}
    }
    
    r3d_header(stdout, rots, rot_order, scale, shift, span);
    for (i = 0; i < 10; i++)
    {
	rtp1->z += 10;
	r3d_volumeElement_polar(stdout, rtp1, 0.04, 5.0, 5.0, grey1, 0.0);
    }
    printf("#Program Ends.\n");  
}
