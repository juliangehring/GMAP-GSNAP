/* $Id: chrsubset.h,v 1.4 2005/06/23 22:46:24 twu Exp $ */
#ifndef CHRSUBSET_INCLUDED
#define CHRSUBSET_INCLUDED
#include "bool.h"
#include "genomicpos.h"
#include "iit-read.h"

#define T Chrsubset_T
typedef struct T *T;

extern void
Chrsubset_print (T this);
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
extern T
Chrsubset_read (char *user_chrsubsetfile, char *genomesubdir, char *fileroot, 
		char *user_chrsubsetname, IIT_T chromosome_iit);

#undef T
#endif
