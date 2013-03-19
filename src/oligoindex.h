/* $Id: oligoindex.h 79521 2012-11-19 22:11:24Z twu $ */
#ifndef OLIGOINDEX_INCLUDED
#define OLIGOINDEX_INCLUDED
#include "bool.h"
#include "types.h"
#include "genomicpos.h"
#include "list.h"
#include "intlist.h"
#include "diagpool.h"

#define OVERABUNDANCE_CHECK 50
#define OVERABUNDANCE_PCT 0.97
#define OVERABUNDANCE_MIN 200

typedef UINT4 Shortoligomer_T;

#define T Oligoindex_T
typedef struct T *T;
struct T {

#ifdef PMAP
  int indexsize_aa;
  Shortoligomer_T msb;
#else
  int indexsize;
  Shortoligomer_T mask;
#endif

  int diag_lookback;
  int suffnconsecutive;

  bool query_evaluated_p;
  Oligospace_T oligospace;
  bool *overabundant;
  bool *inquery;
  int *counts;
  int *relevant_counts;
  Genomicpos_T **positions;
  Genomicpos_T **pointers;
};

extern int
Oligoindex_indexsize (T this);

extern int *
Oligoindex_counts_copy (T this);

extern void
Oligoindex_dump (T this);

extern void
Oligoindex_counts_dump (T this, int *counts);

extern bool
Oligoindex_counts_equal (T this, int *counts);

extern T *
Oligoindex_new_major (int *noligoindices);

extern T *
Oligoindex_new_minor (int *noligoindices);

extern int
Oligoindex_indexsize (T this);

extern double
Oligoindex_set_inquery (int *badoligos, int *repoligos, int *trimoligos, int *trim_start, int *trim_end,
			T this, char *queryuc_ptr, int querylength, bool trimp);
extern void
Oligoindex_tally (T this, char *genomicuc_trimptr, int genomicuc_trimlength,
		  char *queryuc_ptr, int querylength, int sequencepos);
extern void
Oligoindex_untally (T this);
extern void
Oligoindex_clear_inquery (T this);
extern void
Oligoindex_free_array (T **oligoindices, int noligoindices);

extern List_T
Oligoindex_get_mappings (List_T diagonals, bool *coveredp, Genomicpos_T **mappings, int *npositions,
			 int *totalpositions, bool *oned_matrix_p, int *maxnconsecutive, 
			 T this, char *queryuc_ptr, int querylength,
			 Genomicpos_T chrstart, Genomicpos_T chrend,
			 Genomicpos_T chroffset, Genomicpos_T chrhigh, bool plusp,
			 Diagpool_T diagpool);

#undef T
#endif

