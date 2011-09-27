/* $Id: smooth.h,v 1.13 2006/12/13 20:34:13 twu Exp $ */
#ifndef SMOOTH_INCLUDED
#define SMOOTH_INCLUDED
#include "list.h"
#include "pairpool.h"

extern List_T
Smooth_pairs_by_netgap (bool *deletep, List_T pairs, Pairpool_T pairpool, int stage2_indexsize);
extern List_T
Smooth_pairs_by_size (bool *shortp, bool *deletep, List_T pairs, Pairpool_T pairpool, int stage2_indexsize);

#endif


