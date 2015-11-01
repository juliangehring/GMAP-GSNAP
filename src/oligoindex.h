/* $Id: oligoindex.h 100638 2013-07-05 20:20:50Z twu $ */
#ifndef OLIGOINDEX_INCLUDED
#define OLIGOINDEX_INCLUDED
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "bool.h"
#include "types.h"
#include "genomicpos.h"
#include "list.h"
#include "intlist.h"
#include "diagpool.h"
#ifdef HAVE_SSE2
#include <emmintrin.h>
#endif

#define OVERABUNDANCE_CHECK 50
#define OVERABUNDANCE_PCT 0.97
#define OVERABUNDANCE_MIN 200

typedef UINT4 Shortoligomer_T;
typedef unsigned char Count_T;

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
#ifdef HAVE_SSE2
  __m128i *inquery_allocated;
  __m128i *counts_allocated;
  Count_T *inquery;
#else
  bool *inquery;
#endif
  Count_T *counts;
#ifdef PMAP
  int *relevant_counts;
  bool *overabundant;
#endif
  Chrpos_T **positions;
  Chrpos_T **pointers;
};

extern int
Oligoindex_indexsize (T this);

extern Count_T *
Oligoindex_counts_copy (T this);

extern void
Oligoindex_dump (T this);

extern void
Oligoindex_counts_dump (T this, Count_T *counts);

extern bool
Oligoindex_counts_equal (T this, Count_T *counts);

extern T *
Oligoindex_new_major (int *noligoindices);

extern T *
Oligoindex_new_minor (int *noligoindices);

extern int
Oligoindex_indexsize (T this);

extern double
Oligoindex_set_inquery (int *badoligos, int *repoligos, int *trimoligos, int *trim_start, int *trim_end,
			T this, char *queryuc_ptr, int querylength, bool trimp);
#ifdef PMAP
extern void
Oligoindex_tally (T this, char *genomicuc_trimptr, int genomicuc_trimlength,
		  char *queryuc_ptr, int querylength, int sequencepos);
#endif
extern void
Oligoindex_untally (T this, char *queryuc_ptr, int querylength);
extern void
Oligoindex_clear_inquery (T this, char *queryuc_ptr, int querylength);
extern void
Oligoindex_free_array (T **oligoindices, int noligoindices);

extern List_T
Oligoindex_get_mappings (List_T diagonals, bool *coveredp, Chrpos_T **mappings, int *npositions,
			 int *totalpositions, bool *oned_matrix_p, int *maxnconsecutive, 
			 T this, char *queryuc_ptr, int querylength,
			 Chrpos_T chrstart, Chrpos_T chrend,
			 Univcoord_T chroffset, Univcoord_T chrhigh, bool plusp,
			 Diagpool_T diagpool);

#undef T
#endif

