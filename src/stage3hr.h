/* $Id: stage3hr.h 46162 2011-08-31 16:28:11Z twu $ */
#ifndef STAGE3HR_INCLUDED
#define STAGE3HR_INCLUDED

#include <stdio.h>
#include "bool.h"
#include "sense.h"
#include "chrnum.h"
#include "genomicpos.h"
#include "intlist.h"
#include "iit-read.h"
#include "shortread.h"
#include "genome.h"
#include "compress.h"
#include "resulthr.h"
#include "substring.h"
#include "pair.h"


/* Should arrange in order of goodness, best to worst */
typedef enum {EXACT, SUB, INSERTION, DELETION,
	      HALFSPLICE_DONOR, HALFSPLICE_ACCEPTOR, SPLICE, DISTANT_SPLICE,
	      ONE_THIRD_SHORTEXON, TWO_THIRDS_SHORTEXON, SHORTEXON,
	      TERMINAL, GMAP} Hittype_T;


#define T Stage3end_T
typedef struct T *T;

typedef struct Stage3pair_T *Stage3pair_T;


extern void
Stage3hr_setup (bool invert_first_p_in, bool invert_second_p_in,
		IIT_T genes_iit_in, int *genes_divint_crosstable_in,
		IIT_T tally_iit_in, int *tally_divint_crosstable_in,
		IIT_T runlength_iit_in, int *runlength_divint_crosstable_in,
		int mapq_unique_score_in, bool distances_observed_p,
		int pairmax_in, int antistranded_penalty_in, bool favor_multiexon_p_in);

extern Hittype_T
Stage3end_hittype (T this);
extern char *
Stage3end_hittype_string (T this);
extern Chrnum_T
Stage3end_chrnum (T this);
extern Chrnum_T
Stage3end_effective_chrnum (T this);
extern Genomicpos_T
Stage3end_chroffset (T this);
extern Genomicpos_T
Stage3end_chrhigh (T this);
extern Genomicpos_T
Stage3end_genomicstart (T this);
extern Genomicpos_T
Stage3end_genomicend (T this);
extern int
Stage3end_query_alignment_length (T this);
extern Genomicpos_T
Stage3end_genomic_alignment_length (T this);
extern Genomicpos_T
Stage3end_chrpos_low_trim (T this);
extern int
Stage3end_mapq_score (T this);
extern int
Stage3end_score (T this);
extern int
Stage3end_best_score (List_T hits);
extern int
Stage3end_best_score_paired (List_T hits);
extern int
Stage3end_nmatches (T this);
extern int
Stage3end_nmismatches_whole (T this);
extern int
Stage3end_nmismatches_bothdiff (T this);
extern int
Stage3end_nmismatches_refdiff (T this);
extern Endtype_T
Stage3end_gmap_start_endtype (T this);
extern Endtype_T
Stage3end_gmap_end_endtype (T this);
extern int
Stage3end_nindels (T this);
extern int
Stage3end_indel_pos (T this);
extern bool
Stage3end_plusp (T this);
extern bool
Stage3end_paired_usedp (T this);

extern Substring_T
Stage3end_substring1 (T this);
extern Substring_T
Stage3end_substring2 (T this);
extern char *
Stage3end_deletion_string (T this);
extern Substring_T
Stage3end_substring_donor (T this);
extern Substring_T
Stage3end_substring_acceptor (T this);
extern Substring_T
Stage3end_substringD (T this);
extern Substring_T
Stage3end_substringA (T this);
extern Substring_T
Stage3end_substring_low (T this);
extern Substring_T
Stage3end_substring_high (T this);
extern Substring_T
Stage3end_substring_containing (T this, int querypos);
extern struct Pair_T *
Stage3end_pairarray (T this);
extern int
Stage3end_npairs (T this);
extern Genomicpos_T
Stage3end_distance (T this);
extern Genomicpos_T
Stage3end_shortexon_acceptor_distance (T this);
extern Genomicpos_T
Stage3end_shortexon_donor_distance (T this);
extern int
Stage3end_sensedir (T this);
extern int
Stage3end_sensedir_nonamb (T this);
extern int
Stage3end_cdna_direction (T this);
extern bool
Stage3end_gmap_triedp (T this);
extern void
Stage3end_set_gmap_triedp (T this);
extern int
Stage3end_gmap_querystart (T this);
extern int
Stage3end_gmap_queryend (T this);
extern int
Stage3end_terminal_trim (T this);
extern bool
Stage3end_contains_known_splicesite (T this);
extern bool
Stage3end_bad_stretch_p (T this, Compress_T query_compress_fwd, Compress_T query_compress_rev);

extern bool
Stage3end_genomicbound_from_start (Genomicpos_T *genomicbound, T this, int overlap, Genomicpos_T chroffset);
extern bool
Stage3end_genomicbound_from_end (Genomicpos_T *genomicbound, T this, int overlap, Genomicpos_T chroffset);

extern void
Stage3end_free (T *old);


extern Stage3end_T
Stage3pair_hit5 (Stage3pair_T this);
extern Stage3end_T
Stage3pair_hit3 (Stage3pair_T this);
extern int
Stage3pair_mapq_score (Stage3pair_T this);
extern Genomicpos_T
Stage3pair_pairlength (Stage3pair_T this);
extern int
Stage3pair_nmatches (Stage3pair_T this);
extern int
Stage3pair_overlap (Stage3pair_T this);
extern void
Stage3pair_set_private5p (Stage3pair_T this);
extern void
Stage3pair_clear_private5p (Stage3pair_T this);
extern void
Stage3pair_set_private3p (Stage3pair_T this);
extern void
Stage3pair_clear_private3p (Stage3pair_T this);

extern void
Stage3pair_free (Stage3pair_T *old);

extern T
Stage3end_copy (T old);

extern T
Stage3end_new_exact (int *found_score, Genomicpos_T left, int genomiclength, Compress_T query_compress,
		     bool plusp, Chrnum_T chrnum, Genomicpos_T chroffset, Genomicpos_T chrhigh);
extern T
Stage3end_new_substitution (int *found_score, int nmismatches, Genomicpos_T left,
			    int genomiclength, Compress_T query_compress,
			    bool plusp, Chrnum_T chrnum, Genomicpos_T chroffset, Genomicpos_T chrhigh);
extern T
Stage3end_new_insertion (int *found_score, int nindels, int indel_pos, int nmismatches1, int nmismatches2,
			 Genomicpos_T left, int genomiclength, Compress_T query_compress,
			 int querylength, bool plusp, Chrnum_T chrnum, Genomicpos_T chroffset,
			 Genomicpos_T chrhigh, int indel_penalty, bool end1_indel_p, bool end2_indel_p);
extern T
Stage3end_new_deletion (int *found_score, int nindels, int indel_pos, int nmismatches1, int nmismatches2,
			Genomicpos_T left, int genomiclength, Compress_T query_compress,
			int querylength, bool plusp, Chrnum_T chrnum, Genomicpos_T chroffset,
			Genomicpos_T chrhigh, int indel_penalty, bool end1_indel_p, bool end2_indel_p);

extern T
Stage3end_new_terminal (int querystart, int queryend, int nmismatches,
			Genomicpos_T left, Compress_T query_compress,
			int querylength, bool plusp, Endtype_T left_endtype, Endtype_T right_endtype,
			Chrnum_T chrnum, Genomicpos_T chroffset, Genomicpos_T chrhigh,
			int max_mismatches_allowed);
extern T
Stage3end_new_splice (int *found_score, int donor_nmismatches, int acceptor_nmismatches,
		      Substring_T donor, Substring_T acceptor, Genomicpos_T distance,
		      bool shortdistancep, int splicing_penalty, int querylength,
		      int amb_nmatches, Intlist_T ambi_left, Intlist_T ambi_right,
		      Intlist_T amb_nmismatches_left, Intlist_T amb_nmismatches_right,
		      bool copy_donor_p, bool copy_acceptor_p,
		      bool first_read_p, int sensedir);
extern T
Stage3end_new_shortexon (int *found_score, Substring_T donor, Substring_T acceptor, Substring_T shortexon,
			 Genomicpos_T acceptor_distance, Genomicpos_T donor_distance,
			 int amb_nmatches_donor, int amb_nmatches_acceptor,
			 Intlist_T ambi_left, Intlist_T ambi_right,
			 Intlist_T amb_nmismatches_left, Intlist_T amb_nmismatches_right,
			 bool copy_donor_p, bool copy_acceptor_p, bool copy_shortexon_p,
			 int splicing_penalty, int querylength, int sensedir);

extern T
Stage3end_new_gmap (int nmismatches_whole, int nmatches_pretrim,
		    int ambig_end_length_5, int ambig_end_length_3,
		    Splicetype_T ambig_splicetype_5, Splicetype_T ambig_splicetype_3,
		    struct Pair_T *pairarray, int npairs, int nsegments, Genomicpos_T left, int genomiclength,
		    bool plusp, int querylength,
		    Chrnum_T chrnum, Genomicpos_T chroffset, Genomicpos_T chrhigh, int sensedir);

extern List_T
Stage3end_sort_bymatches (List_T hits);

extern Stage3end_T *
Stage3end_eval_and_sort (int *npaths, Stage3end_T *stage3array, int maxpaths, Shortread_T queryseq,
			 Compress_T query_compress_fwd, Compress_T query_compress_rev,
			 Genome_T genome, char *quality_string, bool displayp);
extern List_T
Stage3pair_remove_excess_terminals (List_T hitpairlist);
extern List_T
Stage3end_optimal_score (List_T hitlist, int cutoff_level, int suboptimal_mismatches);
extern int
Stage3end_noptimal (List_T hitlist);
extern List_T
Stage3end_remove_duplicates (List_T hitlist, Shortread_T queryseq1, Shortread_T queryseq2);
extern List_T
Stage3end_remove_overlaps (List_T hitlist);
extern List_T
Stage3end_resolve_multimapping (List_T hitlist);
extern Pairtype_T
Stage3_determine_pairtype (T hit5, T hit3);
extern Pairtype_T
Stage3pair_pairtype (Stage3pair_T this);


/* If hit5 and hit3 are not NULL, then we know this is part of a pair */
extern void
Stage3end_print (FILE *fp, T this, int score,
		 IIT_T chromosome_iit, Shortread_T queryseq,
		 bool invertp, T hit5, T hit3, int pairedlength, int pairscore,
		 Pairtype_T pairtype, int mapq_score);

extern void
Stage3pair_print (Result_T result, Resulttype_T resulttype,
		  IIT_T chromosome_iit, Shortread_T queryseq1, Shortread_T queryseq2,
		  int maxpaths, bool quiet_if_excessive_p,
		  bool nofailsp, bool failsonlyp,
		  bool fails_as_input_p, bool fastq_format_p, int quality_shift,
		  FILE *fp_nomapping_1, FILE *fp_nomapping_2,
		  FILE *fp_unpaired_uniq, FILE *fp_unpaired_transloc, FILE *fp_unpaired_mult,
		  FILE *fp_halfmapping_uniq, FILE *fp_halfmapping_transloc, FILE *fp_halfmapping_mult,
		  FILE *fp_paired_uniq_inv, FILE *fp_paired_uniq_scr,
		  FILE *fp_paired_uniq_long, FILE *fp_paired_mult,
		  FILE *fp_concordant_uniq, FILE *fp_concordant_transloc, FILE *fp_concordant_mult);

extern Stage3pair_T
Stage3pair_new (T hit5, T hit3, Genomicpos_T *splicesites,
		Compress_T query5_compress_fwd, Compress_T query5_compress_rev,
		Compress_T query3_compress_fwd, Compress_T query3_compress_rev,
		Genomicpos_T expected_pairlength,
		Pairtype_T pairtype, int splicing_penalty,
		bool private5p, bool private3p, bool expect_concordant_p);
extern void
Stage3pair_privatize (Stage3pair_T *array, int npairs);

extern List_T
Stage3pair_sort_bymatches (List_T hits);

extern List_T
Stage3pair_remove_overlaps (List_T hitpairlist);

extern List_T
Stage3pair_resolve_multimapping (List_T hitpairs);

extern Stage3pair_T *
Stage3pair_eval_and_sort (int *npaths, Stage3pair_T *stage3pairarray, int maxpaths,
			  Shortread_T queryseq1, Shortread_T queryseq2,
			  Compress_T query5_compress_fwd, Compress_T query5_compress_rev, 
			  Compress_T query3_compress_fwd, Compress_T query3_compress_rev, 
			  Genome_T genome, char *quality_string_5, char *quality_string_3);

extern List_T
Stage3pair_optimal_score (List_T hitpairlist, int cutoff_level, int suboptimal_mismatches);

extern List_T
Stage3_pair_up_concordant (bool *abort_pairing_p, int *found_score, int *nconcordant,
			   List_T *samechr, List_T *conc_transloc, List_T hitpairs,
			   List_T *hitarray5, int narray5, List_T *hitarray3, int narray3,
			   int cutoff_level_5, int cutoff_level_3, int subopt_levels,
			   Genomicpos_T *splicesites,
			   Compress_T query5_compress_fwd, Compress_T query5_compress_rev,
			   Compress_T query3_compress_fwd, Compress_T query3_compress_rev,
			   Genomicpos_T expected_pairlength,
			   int querylength5, int querylength3, int maxpairedpaths,
			   int splicing_penalty);

#undef T
#endif

