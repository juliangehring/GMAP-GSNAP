/* $Id: smooth.h 90983 2013-04-01 19:42:40Z twu $ */
#ifndef SMOOTH_INCLUDED
#define SMOOTH_INCLUDED
#include "list.h"
#include "pairpool.h"

extern List_T
Smooth_pairs_by_netgap (bool *deletep, List_T pairs, Pairpool_T pairpool);
extern List_T
Smooth_pairs_by_size (bool *shortp, bool *deletep, List_T pairs, Pairpool_T pairpool, int stage2_indexsize);
extern List_T
Smooth_mark_short_exons (List_T pairs, Pairpool_T pairpool, int stage2_indexsize);

#endif


