/* $Id: matchpool.h,v 1.2 2006/10/09 16:59:23 twu Exp $ */
#ifndef MATCHPOOL_INCLUDED
#define MATCHPOOL_INCLUDED
#include "bool.h"
#include "iit-read.h"
#include "genomicpos.h"
#include "match.h"
#include "list.h"

#define T Matchpool_T
typedef struct T *T;

extern void
Matchpool_free (T *old);
extern T
Matchpool_new (void);
extern void
Matchpool_reset (T this);
extern void
Matchpool_save (T this);
extern void
Matchpool_restore (T this);
extern List_T
Matchpool_push (List_T list, T this, int querypos, bool forwardp, bool fivep,
		Genomicpos_T position, IIT_T chromosome_iit);
extern List_T
Matchpool_push_existing (List_T list, T this, Match_T match);
extern List_T
Matchpool_pop (List_T list, Match_T *x);
extern List_T
Matchpool_transfer (List_T dest, List_T source);

#undef T
#endif
