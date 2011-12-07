#ifndef SPLICETRIE_INCLUDED
#define SPLICETRIE_INCLUDED

#include "bool.h"
#include "types.h"
#include "genomicpos.h"
#include "list.h"
#include "intlist.h"
#include "uintlist.h"
#include "splicetrie_build.h"	/* For Splicetype_T */
#include "compress.h"

#include "dynprog.h"
#include "pairpool.h"

extern void
Splicetrie_setup (
#ifdef GSNAP
		  UINT4 *splicecomp_in,
#endif
		  Genomicpos_T *splicesites_in, UINT4 *splicefrags_ref_in, UINT4 *splicefrags_alt_in,
		  unsigned int *trieoffsets_obs_in, unsigned int *triecontents_obs_in,
		  unsigned int *trieoffsets_max_in, unsigned int *triecontents_max_in,
		  bool snpp_in, bool amb_closest_p_in, bool amb_clip_p_in, int min_shortend_in);

#ifdef GSNAP
extern bool
Splicetrie_splicesite_p (Genomicpos_T left, int pos5, int pos3);
#endif

extern List_T
Splicetrie_solve_end5 (List_T best_pairs, unsigned int *triestart,
		       Genomicpos_T knownsplice_limit_low, Genomicpos_T knownsplice_limit_high,

		       int *finalscore, int *nmatches, int *nmismatches,
		       int *nopens, int *nindels, bool *knownsplicep, int *ambig_end_length,
		       int *threshold_miss_score, int obsmax_penalty, int perfect_score,

		       Genomicpos_T anchor_splicesite, char *splicejunction, int splicelength,
		       Splicetype_T far_splicetype,
		       Genomicpos_T chroffset, Genomicpos_T chrpos, int genomiclength,
		       int *dynprogindex, Dynprog_T dynprog, 
		       char *revsequence1, char *revsequenceuc1,
		       int length1, int length2, int revoffset1, int revoffset2,
#ifdef PMAP
		       char *queryaaseq,
#endif
		       int cdna_direction, bool watsonp, bool jump_late_p, Pairpool_T pairpool,
		       int extraband_end, double defect_rate);

extern List_T
Splicetrie_solve_end3 (List_T best_pairs, unsigned int *triestart,
		       Genomicpos_T knownsplice_limit_low, Genomicpos_T knownsplice_limit_high,

		       int *finalscore, int *nmatches, int *nmismatches,
		       int *nopens, int *nindels, bool *knownsplicep, int *ambig_end_length,
		       int *threshold_miss_score, int obsmax_penalty, int perfect_score,

		       Genomicpos_T anchor_splicesite, char *splicejunction, int splicelength, int contlength,
		       Splicetype_T far_splicetype,
		       Genomicpos_T chroffset, Genomicpos_T chrpos, int genomiclength,
		       int *dynprogindex, Dynprog_T dynprog, 
		       char *sequence1, char *sequenceuc1,
		       int length1, int length2, int offset1, int offset2,
#ifdef PMAP
		       char *queryaaseq,
#endif
		       int cdna_direction, bool watsonp, bool jump_late_p, Pairpool_T pairpool,
		       int extraband_end, double defect_rate);

#ifdef GSNAP
#ifdef END_KNOWNSPLICING_SHORTCUT

extern Uintlist_T
Splicetrie_dump_coords_left (int *best_nmismatches, unsigned int *triestart, int pos5, int pos3,
			     Compress_T query_compress, char *queryptr, bool plusp,
			     Genomicpos_T knownsplice_limit_low, Genomicpos_T knownsplice_limit_high);

extern Uintlist_T
Splicetrie_dump_coords_right (int *best_nmismatches, unsigned int *triestart, int pos5, int pos3,
			      Compress_T query_compress, char *queryptr, bool plusp,
			      Genomicpos_T knownsplice_limit_low, Genomicpos_T knownsplice_limit_high);
#endif
#endif

#ifdef GSNAP
extern Intlist_T
Splicetrie_find_left (int *best_nmismatches, Intlist_T *nmismatches_list, int i,
		      Genomicpos_T origleft, int pos5, int pos3, Genomicpos_T chroffset,
		      Compress_T query_compress, char *queryptr, int querylength,
		      int max_mismatches_allowed, bool plusp, int genestrand,
		      bool collect_all_p);

extern Intlist_T
Splicetrie_find_right (int *best_nmismatches, Intlist_T *nmismatches_list, int i,
		       Genomicpos_T origleft, int pos5, int pos3, Genomicpos_T chrhigh,
		       Compress_T query_compress, char *queryptr, int max_mismatches_allowed,
		       bool plusp, int genestrand, bool collect_all_p);
#endif

#endif
