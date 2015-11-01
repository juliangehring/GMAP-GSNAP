/* $Id: indel.h 133760 2014-04-20 05:16:56Z twu $ */
#ifndef INDEL_INCLUDED
#define INDEL_INCLUDED
#include "bool.h"
#include "list.h"
#include "chrnum.h"
#include "genomicpos.h"
#include "compress.h"

extern void
Indel_setup (int min_indel_end_matches_in, int indel_penalty_middle_in);

extern List_T
Indel_solve_middle_insertion (bool *foundp, int *found_score, int *nhits, List_T hits,
			      Univcoord_T left, Chrnum_T chrnum, Univcoord_T chroffset,
			      Univcoord_T chrhigh, Chrpos_T chrlength,
			      int indels, Compress_T query_compress,
			      int querylength, int max_mismatches_allowed,
			      bool plusp, int genestrand, bool first_read_p, bool sarrayp);

extern List_T
Indel_solve_middle_deletion (bool *foundp, int *found_score, int *nhits, List_T hits,
			     Univcoord_T left, Chrnum_T chrnum, Univcoord_T chroffset,
			     Univcoord_T chrhigh, Chrpos_T chrlength,
			     int indels, Compress_T query_compress,
			     int querylength, int max_mismatches_allowed,
			     bool plusp, int genestrand, bool first_read_p, bool sarrayp);

#endif

