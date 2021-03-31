#include <stdio.h>
#include "args.h"
#include "atom.h"
#include "cfg.h"
#include "chem.h"
#include "point.h"
#include "cfg_sat.h"

double version_=1.00;
int build_=1;

extern char cfg_outfile_[];
extern short quiet_;

void main (int argc, char * argv[])
{
    FILE * fp = NULL;
    int succ, nTrials=1;
    
    /* Handle all command line arguments */
    arg_handler(argc, argv, TRAP_BADWORDS);
    
    /* Perform initializations */
    Chem_InitializePeriodicTable();
    pef_initialize();
    filenames_initialize();
    cfg_initialize();
        
    /* Read in the cfg data */
    cfg_establish();
    
    /* Perform the saturation */
    cfg_sat(nTrials, &succ);

    /* Output the configuration data */
    if (cfg_outfile_[0]) fp=fopen(cfg_outfile_, "w");
    else fp=stdout;
    fprintf(fp, "# Trials %i  Successes %i  Ratio %.2lf\n", 
	nTrials, succ, ((double)nTrials/succ));
    cfg_out(fp);

    fclose(fp);
    
    if (!quiet_) printf("# Thank you for using md series 2.\n");
}
