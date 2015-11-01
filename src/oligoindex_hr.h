/* $Id: oligoindex_hr.h 132475 2014-04-06 04:14:11Z twu $ */
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

typedef struct Oligoindex_array_T *Oligoindex_array_T;

extern void
Oligoindex_hr_setup (Genomecomp_T *ref_blocks_in, Mode_T mode_in);

extern int
Oligoindex_indexsize (T this);

extern int
Oligoindex_array_length (Oligoindex_array_T oligoindices);
extern T
Oligoindex_array_elt (Oligoindex_array_T oligoindices, int source);

extern Oligoindex_array_T
Oligoindex_array_new_major (int max_querylength, int max_genomiclength);

extern Oligoindex_array_T
Oligoindex_array_new_minor (int max_querylength, int max_genomiclength);

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
Oligoindex_array_free(Oligoindex_array_T *old);

extern List_T
Oligoindex_get_mappings (List_T diagonals, bool *coveredp, Chrpos_T **mappings, int *npositions,
			 int *totalpositions, bool *oned_matrix_p, int *maxnconsecutive, 
			 Oligoindex_array_T array, T this, char *queryuc_ptr, int querylength,
			 Chrpos_T chrstart, Chrpos_T chrend,
			 Univcoord_T chroffset, Univcoord_T chrhigh, bool plusp,
			 Diagpool_T diagpool);

#undef T
#endif

