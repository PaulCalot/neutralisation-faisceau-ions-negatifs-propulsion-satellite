/*
 * md series 2 (c) 1999 cameron abrams
 * 
 * domd.c: single-run driver
 */
#include "domd.h"
extern short quiet_, cfg_newcrystal_;
void do_md (int runId)
{
    if (!quiet_) 
    {
	printf("# md series 2 do_md (%i) driver (c) 1999 cfa\n", runId);
	fflush(stdout);
    }

    /* Insert the run number into any "????" or "####" strings in the 
     * file name templates. */
    filenames_apply();
    if (!quiet_) filenames_echo();
        
    /* Establish the configuration data: 
     * Read-in from input file, maybe introduce ion. */
    cfg_establish();
    
    /* Perform the integration */
    cfg_push();

    /* Output the configuration data, including info on clusters */
    cfg_report();
    
    /* Clear the configuration and cluster information for the next run */
    cfg_clear();
    if (!cfg_newcrystal_) clust_clear();
}
