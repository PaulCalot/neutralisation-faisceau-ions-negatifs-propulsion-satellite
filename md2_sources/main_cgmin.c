#include <stdio.h>
#include "args.h"
#include "atom.h"
#include "cfg.h"
#include "chem.h"
#include "point.h"
#include "cfg_cgmin.h"

double version_=1.00;
int build_=1;

extern char cfg_outfile_[];
extern short quiet_;

int SIC_DIAG2_=1;
int SIC_DIAG3_=2;
void main (int argc, char * argv[])
{
    FILE * fp = NULL;
    double fret;
    int iter;
    
    /* Handle all command line arguments */
    arg_handler(argc, argv, TRAP_BADWORDS);
    
    /* Perform initializations */
    Chem_InitializePeriodicTable();
    pef_initialize();
    filenames_initialize();
    cfg_initialize();
        
    /* Read in the cfg data */
    cfg_establish();
    
    /* Perform the minimization */
    cfg_cgmin(&iter, &fret, NULL);

    /* Output the configuration data */
    if (cfg_outfile_[0]) fp=fopen(cfg_outfile_, "w");
    else fp=stdout;
    fprintf(fp, "# FRPR ConjGrad Minimimized Energy: %.5le\n", fret);
    cfg_out(fp);

    fclose(fp);
    
    if (!quiet_) printf("# Thank you for using md series 2.\n");
}
