/* $Id: intpool.h,v 1.1 2007/02/05 07:45:55 twu Exp $ */
#ifndef INTPOOL_INCLUDED
#define INTPOOL_INCLUDED
#include "intlist.h"

#define T Intpool_T
typedef struct T *T;

extern void
Intpool_free (T *old);
extern T
Intpool_new (void);
extern void
Intpool_reset (T this);
extern Intlist_T
Intpool_push (Intlist_T intlist, T this, int integer);
extern Intlist_T
Intpool_pop (Intlist_T intlist, int *integer);

#undef T
#endif
