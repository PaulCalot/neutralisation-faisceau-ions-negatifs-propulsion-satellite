#include <stdio.h>

/*
 * This program simply reports the PE of the cluster
 *     F(2)
 *    /
 *  Si(1)
 *    \
 *     F(3)
 * 
 * as a function of the angle subtended at the Silicon
 * where r12 and r13 are fixed at 1.6 A.
 * 
 */
 
#include "tbt_sicf.h"
#include "atom.h"
#include "cfg.h"
#include "chem.h"
#include "math.h"

double version_;
int build_;

extern char cfg_infile_[];
extern short pe_set_;
extern double PE_;

int main (int argc,  char * argv[])
{
    pt R01,R02,R12;
    ptPtr r0, r1, r2;
    ptPtr r01=&R01, r02=&R02, r12=&R12;
    double rr01, rr02, rr12, theta102;
    atomPtr si=NULL, f1=NULL, f2=NULL;
    Chem_InitializePeriodicTable();
    pef_initialize();
    cfg_initialize();

    strcpy(cfg_infile_, argv[1]);
    printf("%s\n", cfg_infile_);
    
    cfg_establish();
    
    si=cfg_addr(0);
    f1=cfg_addr(1);
    f2=cfg_addr(2);

    r0=si->pos;
    r1=f1->pos;
    r2=f2->pos;
    
    r0->x=r0->y=r0->z=0.0;
    r1->x=1.6; r1->y=r1->z=0.0;
    r2->x=r2->z=0.0; r2->y=-1.6;
    
    
    ptPtr_subtract(r01, r0, r1);
    ptPtr_subtract(r02, r0, r2);
    ptPtr_subtract(r12, r1, r2);
    
    rr01=ptPtr_abs(r01);
    rr02=ptPtr_abs(r02);
    rr12=ptPtr_abs(r12);

    theta102=180.0/M_PI*acos((rr01*rr01+rr02*rr02-rr12*rr12)/2*rr01*rr02);

    for (theta102=90.0;theta102<180;theta102+=1.0)
    {
	r1->x=rr01*cos(M_PI/180.0*(theta102-90.0));
	r1->y=rr01*sin(M_PI/180.0*(theta102-90.0));
	pe_set_=0;
	cfg_setforce();
	printf("%.5lf %.5lf\n", theta102, PE_);
	fflush(stdout);
    }
    
    cfg_clear();
    
}
