/* $Id: smooth.h,v 1.14 2008/09/04 20:32:13 twu Exp $ */
#ifndef SMOOTH_INCLUDED
#define SMOOTH_INCLUDED
#include "list.h"
#include "pairpool.h"

extern List_T
Smooth_pairs_by_netgap (bool *deletep, List_T pairs, Pairpool_T pairpool);
extern List_T
Smooth_pairs_by_size (bool *shortp, bool *deletep, List_T pairs, Pairpool_T pairpool, int stage2_indexsize);

#endif


