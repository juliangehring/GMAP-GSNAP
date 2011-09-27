/* $Id: match.h,v 1.36 2005/02/15 01:58:50 twu Exp $ */
#ifndef MATCH_INCLUDED
#define MATCH_INCLUDED
#include "bool.h"
#include "genomicpos.h"
#include "iit-read.h"

#define T Match_T
typedef struct T *T;

extern int
Match_querypos (T this);
extern bool
Match_forwardp (T this);
extern bool
Match_fivep (T this);
extern Genomicpos_T
Match_position (T this);
extern int
Match_chrnum (T this);
extern Genomicpos_T
Match_chrpos (T this);
extern void
Match_set_pairedp (T this);
extern bool
Match_pairedp (T this);

extern int
Match_cmp (const void *a, const void *b);

extern T
Match_new (int querypos, bool forwardp, bool fivep, 
	   Genomicpos_T position, IIT_T chromosome_iit);
extern T
Match_copy (T this);

extern void
Match_free (T *old);

extern void
Match_print_two (int pathnum, T start, T end, IIT_T chromosome_iit, IIT_T contig_iit, 
		 char *dbroot, bool zerobasedp);

#undef T
#endif

