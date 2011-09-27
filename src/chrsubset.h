/* $Id: chrsubset.h 40271 2011-05-28 02:29:18Z twu $ */
#ifndef CHRSUBSET_INCLUDED
#define CHRSUBSET_INCLUDED
#include "bool.h"
#include "genomicpos.h"
#include "iit-read.h"
#include "chrnum.h"

#define T Chrsubset_T
typedef struct T *T;

extern void
Chrsubset_free (T *old);

extern T
Chrsubset_make (Chrnum_T chrnum, IIT_T chromosome_iit);
extern void
Chrsubset_print (T this);
extern void
Chrsubset_print_chromosomes (T this, IIT_T chromosome_iit);
extern char *
Chrsubset_name (T this);
extern int
Chrsubset_nincluded (T this, IIT_T chromosome_iit);
bool
Chrsubset_includep (T this, Genomicpos_T position, IIT_T chromosome_iit);
extern int *
Chrsubset_newindices (T this);
extern int
Chrsubset_newindex (T this, int index);
extern int
Chrsubset_oldindex (T this, int index);
extern unsigned int *
Chrsubset_transitions (int **signs, int *nedges, T this, IIT_T chromosome_iit);
extern T
Chrsubset_new_single (Chrnum_T chrnum, IIT_T chromosome_iit);
extern T
Chrsubset_read (char *user_chrsubsetfile, char *genomesubdir, char *fileroot, 
		char *user_chrsubsetname, IIT_T chromosome_iit);

#undef T
#endif
