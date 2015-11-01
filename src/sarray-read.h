/* $Id: sarray-read.h 109412 2013-09-27 19:27:59Z twu $ */
#ifndef SARRAY_READ_INCLUDED
#define SARRAY_READ_INCLUDED
#include "access.h"
#include "genome.h"
#include "compress.h"
#include "genomicpos.h"
#include "splicetrie.h"
#include "iit-read-univ.h"


#define T Sarray_T
typedef struct T *T;

/* For benchmarking */
extern UINT4 *
Sarray_lcpptrs (Sarray_T this);
extern UINT4 *
Sarray_lcpcomp (Sarray_T this);
Univcoord_T
Sarray_size (Sarray_T this);

extern void
Sarray_setup (T sarray_in, Genome_T genome_in, Univ_IIT_T chromosome_iit_in, int circular_typeint_in,
	      Chrpos_T shortsplicedist_in, int splicing_penalty_in,
	      int max_deletionlength, int max_end_deletions,
	      int max_middle_insertions, int max_end_insertions,
	      Univcoord_T *splicesites_in, Splicetype_T *splicetypes_in,
	      Chrpos_T *splicedists_in, int nsplicesites_in);

extern T
Sarray_new (char *directory, char *fileroot, char *snps_root, Access_mode_T access);
extern void
Sarray_free (T *old);

extern void
Sarray_search_greedy (int *found_score, List_T *subs, List_T *indels, List_T *singlesplicing,
		      List_T *doublesplicing, char *queryuc_ptr, char *queryrc, int querylength,
		      Compress_T query_compress_fwd, Compress_T query_compress_rev, int nmisses_allowed,
		      bool first_read_p);

#undef T
#endif


