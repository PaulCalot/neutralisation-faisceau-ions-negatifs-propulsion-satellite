#include <stdio.h>
#include "args.h"
#include "atom.h"
#include "cfg.h"
#include "chem.h"
#include "point.h"
#include "push.h"
#include "genforce.h"

double version_=1.00;
int build_=1;

extern short quiet_;

int SIC_DIAG2_=1;
int SIC_DIAG3_=2;
void main (int argc, char * argv[])
{
    /* Handle all command line arguments */
    arg_handler(argc, argv);
    
    /* Perform initializations */
    Chem_InitializePeriodicTable();
    pef_initialize();
    cfg_initialize();
    
    if (!quiet_) cfg_showsizes();
    
    /* Establish the cfg data */
    cfg_establish();
    
    /* Perform the integration */
    cfg_push();

    /* Output the configuration data */
    cfg_report();
    
    if (!quiet_) printf("# Thank you for using md series 2.\n");
}
