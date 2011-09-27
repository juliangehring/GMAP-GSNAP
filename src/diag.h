/* $Id: diag.h,v 1.5 2007/02/07 08:35:27 twu Exp $ */
#ifndef DIAG_INCLUDED
#define DIAG_INCLUDED
#include "bool.h"

#define T Diag_T
typedef struct T *T;

extern int
Diag_diagonal (T this);
extern int
Diag_querystart (T this);
extern int
Diag_queryend (T this);
extern int
Diag_nconsecutive (T this);
extern bool
Diag_dominatedp (T this);
extern void
Diag_set_dominatedp (T this);
extern int
Diag_compare_nconsecutive (const void *x, const void *y);
extern int
Diag_compare_diagonal (const void *x, const void *y);

#undef T
#endif

