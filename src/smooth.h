/* $Id: smooth.h 40271 2011-05-28 02:29:18Z twu $ */
#ifndef SMOOTH_INCLUDED
#define SMOOTH_INCLUDED
#include "list.h"
#include "pairpool.h"

extern List_T
Smooth_pairs_by_netgap (bool *deletep, List_T pairs, Pairpool_T pairpool);
extern List_T
Smooth_pairs_by_size (bool *shortp, bool *deletep, List_T pairs, Pairpool_T pairpool, int stage2_indexsize);

#endif


