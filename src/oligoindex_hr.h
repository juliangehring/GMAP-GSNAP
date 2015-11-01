/* $Id: oligoindex_hr.h 121509 2013-12-13 21:56:56Z twu $ */
#ifndef OLIGOINDEX_HR_INCLUDED
#define OLIGOINDEX_HR_INCLUDED

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "bool.h"
#include "types.h"
#include "mode.h"
#include "genomicpos.h"
#include "list.h"
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

  int indexsize;
  Shortoligomer_T mask;

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

  Chrpos_T **positions;
  Chrpos_T **pointers;
};


extern void
Oligoindex_hr_setup (Genomecomp_T *ref_blocks_in, Mode_T mode_in);

extern int
Oligoindex_indexsize (T this);

extern T *
Oligoindex_new_major (int *noligoindices);

extern T *
Oligoindex_new_minor (int *noligoindices);

extern double
Oligoindex_set_inquery (int *badoligos, int *repoligos, int *trimoligos, int *trim_start, int *trim_end,
			T this, char *queryuc_ptr, int querylength, bool trimp);
extern void
Oligoindex_hr_tally (T this, Univcoord_T mappingstart, Univcoord_T mappingend, bool plusp,
		     char *queryuc_ptr, int querylength, Chrpos_T chrpos, int genestrand);
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

