/*
 * MD series2 -- Molecular Dynamics of Plasma-Surface Chemistry
 * (c) 1999 
 * Cameron Abrams
 * and 
 * The Regents of the University of California, Berkeley
 * 
 * Department of Chemical Engineering
 * University of California, Berkeley	1995-1999
 * 
 * Module Name:	clust
 * 
 */

#include "clust.h"

extern double version_;
extern int build_;

extern atomPtr L_;		/* The configuration & trashHeap*/
extern element per_table[];	/* The periodic table */
extern short quiet_;		/* Output supression control */
extern short echo_;		/* Output supression control */
extern pt Lr_;
static atomPtr ip;
/* Clusters:  a cluster is nothing more than a way to identify and keep
 * track of a collection of atoms that are "bound".  Atom "i" and "j" are
 * considered part of the same cluster ("bound") if one there exists a path
 * from i to j along atom-atom bonds (e.g., i->m->n->j). */					 
static cNode clustarray_[MAXNUMCLUSTERS];
					/* array of clusters */
cNodePtr C_=&(clustarray_[0]);		/* head of linked list representation
					 * of the static array */
int nClust_=0;

char clust_outfile_[MAXCHAR]="";	/* name for output file for cluster info */
char cname[MAXCHAR]="";			/* string for clust_name() */

short clust_report_=1;

static const char * clust_statelabel[CS_NULLSTATE] = 
{
    "Base", "Ion", "Ok", "Trash", "Unused"
};

extern ion_type ion_;

int des_ctrl_=2;	/* Thermal desorption model */

static nNodePtr n=NULL, m=NULL, l=NULL;
static int ci=0;
static cNodePtr dc;
void clust_clear (void)
{
    int i=0;
    #if DIAG
    if (!quiet_&&echo_) 
	{printf("# md series 2 clust clearer (c) 1999 cfa\n");fflush(stdout);}
    #endif
    C_=&(clustarray_[0]);
    for (i=0;i<MAXNUMCLUSTERS;i++)
    {
	while (n=nNode_pop(&(C_[i].mem))) nNode_push(&(n->addr->bCards), n);
	C_[i].next=NULL;
	C_[i].id=i;
	C_[i].state=CS_UNUSED;
    }
    nClust_=0;
    cfg_creset();
    #if DIAG
    if (!quiet_&&echo_) 
	{printf("# md series 2 clust clearer complete\n");fflush(stdout);}
    #endif
}

int clust_iontest(cNodePtr c)
{
    if (!c) return 0;
    if (!ion_.mem) return 0;
    for (n=ion_.mem;n;n=n->next)
    {
	for (m=c->mem;m&&m->addr!=n->addr;m=m->next);
	if (!m) return 0;
    }
    for (n=c->mem;n;n=n->next)
    {
	for (m=ion_.mem;m&&m->addr!=n->addr;m=m->next);
	if (!m) return 0;
    }
    return 1;
}

void clust_clean (void)
{
    void clust_trash (void);

    if (!quiet_) printf("# md series 2 cluster cleaner (c) 1999 cfa\n");
    if (!quiet_) printf("# desorption control is %i: %s\n", 
	des_ctrl_, (des_ctrl_==2?"Thrml on":"Thrml off"));

    clust_trash();	/* move member atoms to the trashHeap so they
			 * are not output as part of the simulation cfg */

    clust_cfgreport();	/* output trashed clusters in cfg format */
 
    if (!quiet_) {printf("# md series 2 cluster cleaner complete\n");fflush(stdout);}
}

static int ndb;
int clust_ndb (cNodePtr c)
{
    if (!c||!(c->mem)) return 0;
    if (chem_isInert(c->mem->addr->sym)) return 0;
    ndb=0;
    for (n=c->mem;n;n=n->next)
    {
	ndb+=chem_valence(n->addr->sym)-n->addr->nNbrs;
    }
    if (ndb<0) return 0;  /* overcoordination is possible */
    return ndb;
}

extern double Zlo_;
extern double Eb_;
double tbe=0.0;
double mom=0.0;
static short rule[7];
void clust_trash (void)
/* Applies the desorption model to decide which clusters are kept and
 * which are trashed. */
{
    int i=0;
    cNodePtr c=C_;
    short final=0;
    for(c=C_;c;c=c->next) 
    {
	/* Cluster trash rules defined here */
	for (i=0;i<7;i++) rule[i]=0;
	/* Rule 0: Cluster must not be the base cluster */
	rule[0]=(c->state!=CS_BASE);
	tbe=mom=0.0;
	if (rule[0])
	{
	    /* Rule 1: Cluster is a leftover inert ion */
	    rule[1]=(c->state==CS_ION)&&(chem_isInert(c->mem->addr->sym));
	    /* Rule 2: Cluster is moving up */
	    rule[2]=((mom=clust_comv(c)->z)>0.0);
	    /* Rule 3: Cluster must either be bound with less than the
	     * thermal desorption energy. */
	    if (!rule[1])
	    {
		tbe=clust_b(c);
		/* Note: Eb_ is by convention used as a positive quantity.
		 * A cluster's binding energy is reported as a negative quantity
		 * if that cluster is interacting favorably with it's neighbors.
		 * So, the opposite of the binding energy must be smaller than
		 * the thermal desorption cutoff in order for a cluster to desorb.
		 * Note that if a cluster's binding energy is 0 (or less likely, 
		 * positive) we assume that desorption for such a cluster is 
		 * guaranteed.  In this case, the opposite of the binding energy
		 * is guaranteed to be less than or equal to the specified cutoff
		 * energy, so rule[4] is accepted. */
		rule[3]=((-tbe)<=Eb_);
		/* Rule 4: if cluster is moving down, delete it if it is below
		 * the surface (a push-through). */
		if (mom<0.0) rule[4]=(clust_com(c)->z<Zlo_);
	    }
	}
	final=rule[1]||rule[4];
	/* Desorption control modes:
	 *	Always remove: inert ions and pushthroughs
	 *
	 * 1 == desorb any isolated cluster that
	 *      is moving up. (Physical sputtering)
	 * 2 == desorb any isolated cluster that
	 *	is bound with less than the thermal binding energy Eb_.
	 *	(Physical sputtering and 1st order thermal desorption theory)
	 */
	if (des_ctrl_==1) final=(final||rule[2]);
	if (des_ctrl_==2) final=(final||rule[3]);
	if (final)
	{
	    c->state=CS_TRASH;
	    for (n=c->mem;n;n=n->next) 
	    {
		/* remove this atom from the atomList and put it on
		 * the trashHeap */
		cfg_removeatom(&(n->addr), NULL);
	    }
	    printf("# Clean %s %s: %s %.5le (ctrl=[%s])\n",
		(rule[1]?"Inert":(rule[4]?"Pushthru":"Product")), 
		clust_name(cname, c), 
		(des_ctrl_==1?"mom ":"be "), (des_ctrl_==1?mom:tbe), 
		(des_ctrl_==1?"PS":"PS+TDT"));
	    fflush(stdout);
	}
	else 
	    printf("# Keep %s %s: %s %.5le (ctrl=[%s])\n", 
		(!rule[0]?"Base":(rule[1]?"Inert":(rule[4]?"Pushthru":"Product"))),
		clust_name(cname, c), 
		(des_ctrl_==1?"mom ":"be "), (des_ctrl_==1?mom:tbe), 
		(des_ctrl_==1?"PS":"PS+TDT"));
	fflush(stdout);
    }
}

static int eTally[MAXNUMELEMENTS];
char * clust_name (char * space, cNodePtr c)
{
    int i=0;
    char * p=space;
    if (!space||!c) return NULL;
    for (i=0;i<MAXNUMELEMENTS;i++) eTally[i]=0;
    for (n=c->mem;n;n=n->next) eTally[n->addr->sym]++;
    
    sprintf(p, "%i(%s): ", c->id, clust_statelabel[c->state]);
    p+=strlen(p);
    for (i=0;i<MAXNUMELEMENTS;i++)
    {
	if (eTally[i])
	{
	    sprintf(p, "%s_%i ", per_table[i].sym, eTally[i]);
	    p+=strlen(p);
	}
    }
    return space;
}

static pt COM={0, 0, 0};
static pt COMV={0, 0, 0};
static pt TPT={0, 0, 0};
static ptPtr tpt=&TPT;
static double c_mass=0.0;
ptPtr clust_com (cNodePtr c)
{
    c_mass=0.0;
    if (!c) return NULL;
    ptPtr_clear(&COM);
    for (n=c->mem;n;n=n->next)
    {
	c_mass+=per_table[n->addr->sym].mass;
	ptPtr_add(&COM, &COM, 
	    ptPtr_scalmult(tpt, n->addr->pos, per_table[n->addr->sym].mass));
    }
    return ptPtr_scalmult(&COM, &COM, 1.0/c_mass);
}
ptPtr clust_comv (cNodePtr c)
{
    c_mass=0.0;
    if (!c) return NULL;
    ptPtr_clear(&COMV);
    for (n=c->mem;n;n=n->next)
    {
	c_mass+=per_table[n->addr->sym].mass;
	ptPtr_add(&COMV, &COMV, 
	    ptPtr_scalmult(tpt, n->addr->vel, per_table[n->addr->sym].mass));
    }
    return ptPtr_scalmult(&COMV, &COMV, 1.0/c_mass);
}
double c_k_=0.0;
double clust_k (cNodePtr c)
{
    c_k_=0.0;
    if (!c) return 0.0;
    for (n=c->mem;n;n=n->next) c_k_+=n->addr->mv2;
    return 0.5*c_k_;
}

double c_b_=0.0;
double clust_b (cNodePtr c)
/* Returns the potential energy of interaction between 
 * atoms of this cluster *and* any atoms outside this
 * cluster:  This is how a cluster's "binding energy" 
 * is computed. */
{
    c_b_=0.0;
    if (!c) return 0.0;
    for (n=c->mem;n;n=n->next) 
    /* For each member of the cluster */
    {
	for (m=n->addr->nList;m;m=m->next)
	/* For each neighbor of each atom in the cluster */
	{
	    /* Is this neighbor in the cluster? */
	    for (l=c->mem;l&&l->addr!=m->addr;l=l->next);
	    /* This neighbor isn't in the cluster */
	    if (!l) c_b_+=array_peij(n->addr->id, m->index);
	}
    }
    return c_b_;
}    

static FILE * ofp_=NULL;
void clust_cfgreport (void)
{
    cNodePtr c;
    
    printf("# Cluster report to %s\n", clust_outfile_);
    fflush(stdout);
    
    if (clust_outfile_[0])
    {
	ofp_=fopen(clust_outfile_, "w");
	fprintf(ofp_, "%% Cluster report\n");
	fprintf(ofp_,"#BoxSize.xyz\t%.14lf %.14lf %.14lf\n", 
	    Lr_.x, Lr_.y, Lr_.x);	
	ci=0;
	for (c=C_;c;c=c->next)
	{
	    if (c->state==CS_TRASH)
	    {
		ci++;
		fprintf(ofp_, "%% Cluster %s %.5lf %.5lf\n", 
		    clust_name(cname, c), clust_k(c), clust_b(c));
		fflush(ofp_);
		for (n=c->mem;n;n=n->next) cfg_atomout(n->addr, ofp_);
	    }
	}
	
	fprintf(ofp_, "%% Cluster count %i\n", ci);
	fprintf(ofp_, "%% End of cluster report. (c) 1999 cfa\n");
	fclose(ofp_);
    }
}

cNodePtr new_clust (void)
{
    C_[nClust_].id=nClust_;
    C_[nClust_].state=CS_OK;
    C_[nClust_].mem=NULL;
    #if DIAG
    printf("# new cluster %i\n", C_[nClust_].id);fflush(stdout);
    #endif
    if (nClust_+1==MAXNUMCLUSTERS) {printf("error: too many clusters\n");exit(0);}
    return &(C_[nClust_++]);
}

static cNodePtr c=NULL, oc=NULL, ac, bc;
void clust_enclust (atomPtr a, atomPtr b)
{
    c=oc=ac=bc=NULL;
    n=NULL;
    /* Four possibilities: */
    if (!a->c&&!b->c)	/* neither a nor b belong to a cluster */
    {
	c=new_clust();
	a->c=b->c=(void*)c;
	nNode_push(&(c->mem), nNode_pop(&(a->bCards)));
	nNode_push(&(c->mem), nNode_pop(&(b->bCards)));
	#if DIAG
	printf("# atom %s_%i & %s_%i joined cluster %i\n", 
	    per_table[a->sym].sym, a->id, per_table[b->sym].sym, b->id, c->id);
	fflush(stdout);
	#endif
    }
    else if (a->c&&!(b->c)) /* a has a cluster; b does not */
    {
	c=(cNodePtr)(a->c);
	b->c=(void*)c;
	nNode_push(&(c->mem), nNode_pop(&(b->bCards)));
	#if DIAG
	printf("# atom %s_%i joined atoms %s_%i's cluster %i\n", 
	    per_table[b->sym].sym, b->id, per_table[a->sym].sym, a->id, c->id);
	fflush(stdout);
	#endif
    }
    else if (!(a->c)&&b->c) /* b has a cluster; a does not */
    {
	c=(cNodePtr)(b->c);
	a->c=(void*)c;
	nNode_push(&(c->mem), nNode_pop(&(a->bCards)));
	#if DIAG
	printf("# atom %s_%i joined atoms %s_%i's cluster %i\n", 
	    per_table[a->sym].sym, a->id, per_table[b->sym].sym, b->id, c->id);
	fflush(stdout);
	#endif
    }
    else		    /* both a & b have a cluster */
    {
	if (a->c!=b->c)	    /* a & b have different clusters */
	{
	    ac=(cNodePtr)(a->c);
	    bc=(cNodePtr)(b->c);
	    c=(cNodePtr)(ac->id>bc->id?a->c:b->c);
	    oc=(cNodePtr)(ac->id>bc->id?b->c:a->c);
	    /* preserve the Base cluster designation */
	    c->state=(oc->state==CS_BASE?CS_BASE:c->state);
	    while (n=nNode_pop(&(oc->mem)))
	    {
		n->addr->c=(void*)c;
		nNode_push(&(c->mem), n);
	    }
	    #if DIAG
	    printf("# atoms of cluster %i joined cluster %i\n", 
		oc->id, c->id);
	    fflush(stdout);
	    #endif
	}
	else c=(cNodePtr)a->c;
    }
    if (c->state!=CS_BASE&&(a->state==IS_FIXED||b->state==IS_FIXED))
	c->state=CS_BASE;
}

void clust_sweepsingles (void)
{
    for (ip=L_;ip;ip=ip->next) if (!(ip->c)) 
    {
	ip->c=(void*)new_clust();
	((cNode*)ip->c)->state=CS_OK;
	((cNode*)ip->c)->mem=nNode_pop(&(ip->bCards));
    }
}

int nAic_=0;
void clust_fixlist (void)
{
    #if DIAG
    printf("# md series 2 cluster list fixer (c) 1999 cfa\n");
    fflush(stdout);
    #endif

    C_=&(clustarray_[ci=0]);
    while (ci<MAXNUMCLUSTERS&&!(C_->mem)) C_=&(clustarray_[++ci]);
    if (ci==MAXNUMCLUSTERS) 
    {
	#if DIAG
	printf("# no atoms in clusters\n");
	#endif
	return;
    }
    c=C_;
    nClust_=1;
    nAic_=0;
    while (c)
    {
	if (clust_iontest(c)) c->state=CS_ION;
	if (c->state==CS_OK) 
	{
	    for (n=c->mem;n;n=n->next)
		n->addr->state=IS_UNCLEARED;
	}
	for (n=c->mem;n;n=n->next) nAic_++;
	#if DIAG
	printf("# cluster %s\n", clust_name(cname, c));
	fflush(stdout);
	#endif
	ci++;
	if (ci==MAXNUMCLUSTERS) dc=NULL;
	else dc=&(clustarray_[ci]);
	if (dc) while (ci<MAXNUMCLUSTERS&&!(dc->mem)) dc=&(clustarray_[++ci]);
	if (ci==MAXNUMCLUSTERS) dc=NULL;
	else nClust_++;
	c->next=dc;
	c=c->next;
    }
    #if DIAG
    {
	int nF;
	nF=0;
	for (ip=L_;ip;ip=ip->next) 
	    if (ip->state==IS_FIXED){nF++;}
	printf("clusr_fixlist: number of fixed atoms %i\n", nF);
    }
    #endif
    #if DIAG
    printf("# md series 2 clust fixer: nC=%i nAiC=%i\n", nClust_, nAic_);
    fflush(stdout);
    #endif
}
