/* $Id: smooth.h,v 1.10 2006/03/21 01:46:25 twu Exp $ */
#ifndef SMOOTH_INCLUDED
#define SMOOTH_INCLUDED
#include "bool.h"
#include "list.h"
#include "pairpool.h"

extern List_T
Smooth_pairs (int *nshortexons, int *ndeleteexons, List_T pairs, Pairpool_T pairpool, int indexsize);
extern void
Smooth_reset (List_T pairs);

#endif


