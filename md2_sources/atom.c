/*
 * MD series 2 -- Molecular Dynamics of Plasma-Surface Chemistry
 * (c) 1999 Cam Abrams, Dept. of Chemical Engineering
 * Regents of the University of California, Berkeley
 * 
 * Module Name:		atom 
 * Module Class:	base
 * 
 */

#include "atom.h"

extern element per_table[];		    /* periodic table */
#if DIAG
int ATOM_DIAG_;
#endif

void atom_clear (atomPtr a)
{
    if (!a) return;
    a->sym	= 0;
    a->pbc	= 0;
    a->nNbrs	= 0;
    a->state	= IS_UNDETERMINED;
    a->nList	= NULL;
    a->bCards	= NULL;
    a->next	= NULL;
    a->mv2	= 0.0;
    a->q	= 0.0;
    a->c	= NULL;
    ptPtr_clear(a->pos);
    ptPtr_clear(a->vel);
    ptPtr_clear(a->hac);
    ptPtr_clear(a->frc);
    a->flags.is_RDF_a	= 0;
    a->flags.is_BAD_a	= 0;
    a->flags.is_PEF_a	= 1;
}

void atom_copy (atomPtr a, atomPtr b)
/* Copy all non-array-based info in *a to *b */
{
    if (!a||!b) return;
    a->sym	= b->sym;
    a->pbc	= b->pbc;
    a->nNbrs	= b->nNbrs;
    a->state	= b->state;
    a->next	= b->next;
    a->mv2	= b->mv2;
    a->q	= b->q;
    a->flags.is_RDF_a	= b->flags.is_RDF_a;
    a->flags.is_BAD_a	= b->flags.is_BAD_a;
    a->flags.is_PEF_a	= b->flags.is_PEF_a;
}

void atom_out (FILE * fp, atomPtr a)
{
    if (!a) return;
    atom_display(fp, a);
    fprintf(fp, "End info on atom [%i].\n\n", a->id);
    atom_out(fp, a->next);
}

pt displayAtom_tm;
ptPtr displayAtom_tmp = &displayAtom_tm;
void atom_display (FILE * fp, atomPtr a)
{
    if (!a || !fp) return;
    fprintf(fp, "ATOM ID [%i]:\n", a->id);
    fprintf(fp, "\tSymbol:\t{%s}\n", per_table[a->sym].sym);
    fprintf(fp, "\tKinetic E:\t%.4le\n", 0.5*a->mv2);
    fprintf(fp, "\tPosition:\t");
    ptPtr_out(fp, a->pos);
    fprintf(fp, "\n");
    if (a->flags.is_RDF_a) fprintf(fp, "\tRDF-active.\n");
    if (a->flags.is_BAD_a) fprintf(fp, "\tBAD-active.\n");
    fprintf(fp, "\tVelocity:\t");
    ptPtr_out(fp, a->vel);
    fprintf(fp, "\n");
    fprintf(fp, "\tAcceleration:\t");
    ptPtr_out(fp, ptPtr_scalmult(displayAtom_tmp, a->hac, 2.0));
    fprintf(fp, "\n");
    fprintf(fp, "\tForce:\t\t");
    ptPtr_out(fp, a->frc);
    fprintf(fp, "\n");
    fprintf(fp, "\tState = [%s].\n", atom_stateStr(a));
}

static nNodePtr n;
void atom_becomeNeighbors (atomPtr a, atomPtr b, int * na, int * nb)
{
    if (!a||!b) return;
    
   #if DIAG
    printf("bc0: %s_%i & %s_%i\n", 
	per_table[a->sym].sym, a->id, per_table[b->sym].sym, b->id);
    fflush(stdout);
    #endif 
    n=nNode_pop(&(b->bCards));
    if (!n) {printf("error 1 in a_bc, bobc %s_%i\n",
    per_table[b->sym].sym,b->id);exit(0);}
    nNode_push(&(a->nList), n);

    *na = a->nList->index = (a->nList->next?a->nList->next->index+1:0);
    a->nNbrs++;
    
    n=nNode_pop(&(a->bCards));
    if (!n) {printf("error 1 in a_bc, aobc %s_%i\n",
    per_table[a->sym].sym,a->id);exit(0);}
    nNode_push(&(b->nList), n);

    *nb = b->nList->index = (b->nList->next?b->nList->next->index+1:0);
    b->nNbrs++;
    #if DIAG
    printf("bc1: %s_%i & %s_%i, %i, %i\n", 
	per_table[a->sym].sym, a->id, per_table[b->sym].sym, b->id,
	*na, *nb);
    fflush(stdout);
    #endif 
}

void atom_ungreet (atomPtr a)
{
    nNodePtr n, m;
    if (!a) return;
    #if DIAG
    printf("# atom %s_%i ungreets:\n", per_table[a->sym].sym, a->id);
    #endif
    fflush(stdout);
    while (n=nNode_pop(&(a->nList)))
    {
	a->nNbrs--;
	m=atom_bCardrecall(n->addr, a);
	nNode_push(&(a->bCards), m);
	nNode_push(&(n->addr->bCards), n);
	#if DIAG
	printf("#  -> atom %s_%i.\n", per_table[n->addr->sym].sym, n->addr->id);
	#endif
	fflush(stdout);
    }
}

int nNode_n (nNodePtr L)
{
    if (!L) return 0;
    return (1+nNode_n(L->next));
}

void nNode_push (nNodePtr * L, nNodePtr a)
{
    if (!a||!L) return;
    a->next = *L;
    *L = a;
}

nNodePtr nNode_pop (nNodePtr * L)
{
    nNodePtr a = NULL;
    if (!L||!(*L)) return NULL;
    a=*L;
    *L=(*L)->next;
    a->next=NULL;
    return a;
}

nNodePtr nNode_remove (nNodePtr * L, nNodePtr a)
{
    nNodePtr b=NULL;
    if (!L||!(*L)) return NULL;
    if (a==*L)
    {
	*L=(*L)->next;
	a->next=NULL;
	return a;
    }
    else
    {
	for (b=*L;b&&b->next!=a;b=b->next);
	if (b)
	{
	    b->next=a->next;
	    a->next=NULL;
	    return a;
	}
    }
    return NULL;
}

nNodePtr atom_bCardrecall (atomPtr v, atomPtr a)
{
    nNodePtr n;
    if (!v||!a) return NULL;
    for (n=v->nList;n&&n->addr!=a;n=n->next);
    if (n) 
    {
	v->nNbrs--;
	return nNode_remove(&(v->nList), n);
    }
    return NULL;
}

void atom_nreset (atomPtr a)
{
    nNodePtr n=NULL;
    while (n=nNode_pop(&(a->nList)))
    {
	nNode_push(&(n->addr->bCards), n);
	n->index=0;
    }
    a->nNbrs=0;
}

void atomList_push (atomPtr a, atomPtr * L)
{
    if (!L||!a) return;
    a->next=(*L);
    (*L)=a;
}

static atomPtr tp=NULL;
void atomList_remove (atomPtr a, atomPtr * L)
{
    if (!a||!L||!(*L)) return;
    if ((*L)==a)
    {
	(*L)=a->next;
	a->next=NULL;
    }
    else
    {
	for (tp=(*L);tp&&tp->next!=a;tp=tp->next);
	if (tp)
	{
	    tp->next=a->next;
	    a->next=NULL;
	}
    }
}

int atomList_sizeof (atomPtr L)
{
    if (!L) return 0;
    else return 1+atomList_sizeof(L->next);
}

atomPtr atomList_tailof (atomPtr L)
{
    if (!L) return NULL;
    if (!(L->next)) 
    {printf("tail: %s_%i\n", per_table[L->sym].sym, L->id);fflush(stdout);return L;}
    else return atomList_tailof(L->next);
}

atomPtr atomList_pop (atomPtr * L)
{
    atomPtr a=NULL;
    if (!L||!(*L)) return NULL;
    a=(*L);
    (*L)=(*L)->next;
    a->next=NULL;
    return a;
}

void atomList_indexsort (atomPtr * L)
{
    atomPtr a=NULL, RL=NULL;
    atomPtr atomList_removehiid (atomPtr * L);
    if (!L||!(*L)) return;
    while (a=atomList_removehiid(L)) atomList_push(a, &RL);
    (*L)=RL;
}
static int hiid=0;
static atomPtr tap=NULL, hitap=NULL;
atomPtr atomList_removehiid (atomPtr * L)
{
    hiid=-1;
    if (!L||!(*L)) return NULL;
    for (tap=*L;tap;tap=tap->next) if (tap->id > hiid) {hiid=tap->id;hitap=tap;}
    atomList_remove(hitap, L);
    return hitap;
}

void atomList_zpossort (atomPtr * L)
{
    atomPtr a=NULL, RL=NULL;
    atomPtr atomList_removehi(atomPtr * L);
    if (!L||!(*L)) return;
    while (a=atomList_removehi(L)) atomList_push(a, &RL);
    (*L)=RL;
}
static double hi=0.0;
atomPtr atomList_removehi (atomPtr * L)
{
    hi=-1.e9;
    if (!L||!(*L)) return NULL;
    hitap=NULL;
    for (tap=*L;tap;tap=tap->next) 
	if (tap->pos->z>hi) {hi=tap->pos->z;hitap=tap;}
    atomList_remove(hitap, L);
    return hitap;
}

void atomList_reverse (atomPtr * L)
{
    tap=NULL;
    if (!L||!(*L)) return;
    while (hitap=atomList_pop(L)) atomList_push(hitap, &tap);
    *L=tap;
}


int atom_getState (atomPtr a)
{
    return a->state;
}

void atom_setState (atomPtr a, int State)
{
    a->state = State;
}

char * Atom_StateLabels[ATOM_NSTATES] = 
{
    "Undetermined", "Fixed", "Surface", "Sputtered", 
    "Pushthru", "Bulk", "Falling", "Scattered", "Projectile"
};
char * Atom_StateLabelSyms[ATOM_NSTATES] = 
{
    "?", "F", "S", "$", "P", "B", ">", "+", "*"
};
char * atom_stateStr (atomPtr a)
{
    return Atom_StateLabels[a->state];    
}
char * atom_stateSymStr (atomPtr a)
{
    return Atom_StateLabelSyms[a->state];    
}

