/* $Id: chrsubset.h 99737 2013-06-27 19:33:03Z twu $ */
#ifndef CHRSUBSET_INCLUDED
#define CHRSUBSET_INCLUDED
#include "bool.h"
#include "types.h"
#include "genomicpos.h"
#include "iit-read-univ.h"
#include "chrnum.h"

#define T Chrsubset_T
typedef struct T *T;

extern void
Chrsubset_free (T *old);

extern T
Chrsubset_make (Chrnum_T chrnum, Univ_IIT_T chromosome_iit);
extern void
Chrsubset_print (T this);
extern void
Chrsubset_print_chromosomes (T this, Univ_IIT_T chromosome_iit);
extern char *
Chrsubset_name (T this);
extern int
Chrsubset_nincluded (T this, Univ_IIT_T chromosome_iit);
bool
Chrsubset_includep (T this, Univcoord_T position, Univ_IIT_T chromosome_iit);
extern int *
Chrsubset_newindices (T this);
extern int
Chrsubset_newindex (T this, int index);
extern int
Chrsubset_oldindex (T this, int index);
extern unsigned int *
Chrsubset_transitions (int **signs, int *nedges, T this, Univ_IIT_T chromosome_iit);
extern T
Chrsubset_new_single (Chrnum_T chrnum, Univ_IIT_T chromosome_iit);
extern T
Chrsubset_read (char *user_chrsubsetfile, char *genomesubdir, char *fileroot, 
		char *user_chrsubsetname, Univ_IIT_T chromosome_iit);

#undef T
#endif
