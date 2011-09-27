/* $Id: diagpool.h,v 1.4 2010-07-10 19:16:41 twu Exp $ */
#ifndef DIAGPOOL_INCLUDED
#define DIAGPOOL_INCLUDED
#include "diag.h"
#include "list.h"

#define USE_DIAGPOOL 1

#define T Diagpool_T
typedef struct T *T;

extern void
Diagpool_free (T *old);
extern void
Diagpool_free_memory (T this);
extern T
Diagpool_new (void);
extern void
Diagpool_reset (T this);
extern List_T
Diagpool_push (List_T list, T this, int diagonal, int querystart, int queryend, int nconsecutive);
extern List_T
Diagpool_pop (List_T list, Diag_T *x);

#undef T
#endif
