/* $Id: smooth.h,v 1.7 2005/02/15 01:58:51 twu Exp $ */
#ifndef SMOOTH_INCLUDED
#define SMOOTH_INCLUDED
#include "list.h"
#include "pairpool.h"

extern List_T
Smooth_pairs (int *nshortexons, List_T pairs, Pairpool_T pairpool);
extern void
Smooth_reset (List_T pairs);

#endif


