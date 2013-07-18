/* $Id: spanningelt.h 99756 2013-06-27 21:20:19Z twu $ */
#ifndef SPANNINGELT_INCLUDED
#define SPANNINGELT_INCLUDED
#include "bool.h"
#include "genomicpos.h"
#include "types.h"
#include "indexdb_hr.h"
#include "indexdb.h"
#include "list.h"

#define T Spanningelt_T
typedef struct T *T;
struct T {
  bool partnerp;
  int querylength;
  int partner_querypos;		/* for debugging */
  int querypos;			/* for debugging */

  /* Intersectionr results are in native format, not littleendian */
  Univcoord_T *intersection_diagonals;
  int intersection_ndiagonals;

  Univcoord_T *partner_positions;
  int partner_diagterm;
  int partner_npositions;

  Compoundpos_T compoundpos;
  int compoundpos_diagterm;

  Univcoord_T *positions_allocated;
  Univcoord_T *positions;
  int diagterm;
  int npositions;
  
  int candidates_score;		/* score used for generating candidates */
  int pruning_score;		 /* score used for pruning */
  int miss_querypos5; /* If partnerp is true, this is the overlap of the two partners */
  int miss_querypos3;

  /* Reset values */
  Univcoord_T *intersection_diagonals_reset;
  int intersection_ndiagonals_reset;
  Univcoord_T *partner_positions_reset;
  int partner_npositions_reset;
  Univcoord_T *positions_reset;
  int npositions_reset;
};

extern Width_T
Spanningelt_setup (Width_T index1part_in, Width_T index1interval_in);
extern void
Spanningelt_init_positions_free (bool positions_fileio_p);
extern void
Spanningelt_free (T *old);
extern T
Spanningelt_reset (T this);
extern void
Spanningelt_print (T this);
extern void
Spanningelt_print_set (List_T spanningset);

extern List_T
Spanningelt_set (int *minscore, Storedoligomer_T *stage1_oligos, bool **stage1_retrievedp,
		 Univcoord_T ***stage1_positions, int **stage1_npositions,
		 Indexdb_T indexdb, int query_lastpos, int querylength,
		 int mod, bool plusp);

extern int
Spanningelt_candidates_cmp (const void *a, const void *b);
extern int
Spanningelt_pruning_cmp (const void *a, const void *b);

extern Univcoord_T *
Spanningelt_diagonals (int *ndiagonals, T this, int *miss_querypos5, int *miss_querypos3);

#undef T
#endif

