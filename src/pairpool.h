/* $Id: pairpool.h 52068 2011-11-09 19:32:06Z twu $ */
#ifndef PAIRPOOL_INCLUDED
#define PAIRPOOL_INCLUDED

typedef struct Pairpool_T *Pairpool_T;

#include "bool.h"
#include "pair.h"
#include "list.h"

#define T Pairpool_T

extern void
Pairpool_free (T *old);
extern void
Pairpool_free_memory (T this);
extern void
Pairpool_report_memory (T this);
extern T
Pairpool_new (void);
extern void
Pairpool_reset (T this);
extern List_T
Pairpool_push (List_T list, T this, int querypos, int genomepos, char cdna, char comp, char genome, int dynprogindex);
extern List_T
Pairpool_push_gapalign (List_T list, T this, int querypos, int genomepos, char cdna, char comp, char genome, bool extraexonp);
extern List_T
Pairpool_push_gapholder (List_T list, T this, int queryjump, int genomejump);
extern List_T
Pairpool_push_existing (List_T list, T this, Pair_T pair);
extern List_T
Pairpool_pop (List_T list, Pair_T *x);
extern List_T
Pairpool_transfer (List_T dest, List_T source);
extern List_T
Pairpool_transfer_n (List_T dest, List_T source, int n);
extern int
Pairpool_count_bounded (int *nstart, List_T source, int minpos, int maxpos);
extern List_T
Pairpool_transfer_bounded (List_T dest, List_T source, int minpos, int maxpos);
extern List_T
Pairpool_copy (List_T source, T this);
extern struct Pair_T *
Pairpool_copy_array (struct Pair_T *source, int npairs);

#undef T
#endif
