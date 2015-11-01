/* $Id: stage3hr.h 148359 2014-09-19 22:09:34Z twu $ */
#ifndef STAGE3HR_INCLUDED
#define STAGE3HR_INCLUDED

#include <stdio.h>
#include "bool.h"
#include "sense.h"
#include "chrnum.h"
#include "genomicpos.h"
#include "intlist.h"
#include "iit-read-univ.h"
#include "iit-read.h"
#include "shortread.h"
#include "genome.h"
#include "compress.h"
#include "resulthr.h"
#include "substring.h"
#include "pair.h"


/* Should arrange in order of goodness, best to worst */
typedef enum {EXACT, SUB, INSERTION, DELETION,
	      HALFSPLICE_DONOR, HALFSPLICE_ACCEPTOR, SPLICE, SAMECHR_SPLICE, TRANSLOC_SPLICE,
	      ONE_THIRD_SHORTEXON, TWO_THIRDS_SHORTEXON, SHORTEXON,
	      GMAP, TERMINAL} Hittype_T;


#define T Stage3end_T
typedef struct T *T;

typedef struct Stage3pair_T *Stage3pair_T;


extern void
Stage3hr_setup (bool invert_first_p_in, bool invert_second_p_in,
		IIT_T genes_iit_in, int *genes_divint_crosstable_in,
		IIT_T tally_iit_in, int *tally_divint_crosstable_in,
		IIT_T runlength_iit_in, int *runlength_divint_crosstable_in,
		int terminal_output_minlength_in, bool distances_observed_p, int pairmax_in,
		Chrpos_T expected_pairlength, Chrpos_T pairlength_deviation,
		int localsplicing_penalty_in, int indel_penalty_middle_in,
		int antistranded_penalty_in, bool favor_multiexon_p_in,
		int gmap_min_nconsecutive_in, int index1part, int index1interval,
		bool novelsplicingp_in, bool merge_samechr_p_in,
		bool *circularp_in, char *failedinput_root_in, bool fastq_format_p_in,
		bool print_m8_p_in, bool want_random_p_in);

extern void
Stage3hr_file_setup_single (FILE *fp_failedinput_in, FILE *fp_nomapping_in,
			    FILE *fp_unpaired_uniq_in, FILE *fp_unpaired_circular_in, FILE *fp_unpaired_transloc_in,
			    FILE *fp_unpaired_mult_in, FILE *fp_unpaired_mult_xs_1_in);

extern void
Stage3hr_file_setup_paired (FILE *fp_failedinput_1_in, FILE *fp_failedinput_2_in, FILE *fp_nomapping_in,
			    FILE *fp_halfmapping_uniq_in, FILE *fp_halfmapping_circular_in, FILE *fp_halfmapping_transloc_in,
			    FILE *fp_halfmapping_mult_in, FILE *fp_halfmapping_mult_xs_1_in, FILE *fp_halfmapping_mult_xs_2_in,
			    FILE *fp_paired_uniq_circular_in, FILE *fp_paired_uniq_inv_in, FILE *fp_paired_uniq_scr_in,
			    FILE *fp_paired_uniq_long_in, FILE *fp_paired_mult_in, FILE *fp_paired_mult_xs_1_in, FILE *fp_paired_mult_xs_2_in,
			    FILE *fp_concordant_uniq_in, FILE *fp_concordant_circular_in, FILE *fp_concordant_transloc_in, 
			    FILE *fp_concordant_mult_in, FILE *fp_concordant_mult_xs_1_in, FILE *fp_concordant_mult_xs_2_in);

extern void
Stage3hr_file_setup_all (FILE *fp_failedinput_1_in, FILE *fp_failedinput_2_in, FILE *fp_nomapping_in,
			 FILE *fp_unpaired_uniq_in, FILE *fp_unpaired_circular_in, FILE *fp_unpaired_transloc_in,
			 FILE *fp_unpaired_mult_in, FILE *fp_unpaired_mult_xs_1_in, FILE *fp_unpaired_mult_xs_2_in,
			 FILE *fp_halfmapping_uniq_in, FILE *fp_halfmapping_circular_in, FILE *fp_halfmapping_transloc_in,
			 FILE *fp_halfmapping_mult_in, FILE *fp_halfmapping_mult_xs_1_in, FILE *fp_halfmapping_mult_xs_2_in,
			 FILE *fp_paired_uniq_circular_in, FILE *fp_paired_uniq_inv_in, FILE *fp_paired_uniq_scr_in,
			 FILE *fp_paired_uniq_long_in, FILE *fp_paired_mult_in, FILE *fp_paired_mult_xs_1_in, FILE *fp_paired_mult_xs_2_in,
			 FILE *fp_concordant_uniq_in, FILE *fp_concordant_circular_in, FILE *fp_concordant_transloc_in, 
			 FILE *fp_concordant_mult_in, FILE *fp_concordant_mult_xs_1_in, FILE *fp_concordant_mult_xs_2_in);

extern Hittype_T
Stage3end_hittype (T this);
extern char *
Stage3end_hittype_string (T this);
extern int
Stage3end_genestrand (T this);
extern bool
Stage3end_sarrayp (T this);
extern bool
Stage3end_improved_by_gmap_p (T this);
extern void
Stage3end_set_improved_by_gmap (T this);
extern bool
Stage3end_anomalous_splice_p (T this);
extern Chrnum_T
Stage3end_chrnum (T this);
extern Chrnum_T
Stage3end_effective_chrnum (T this);
extern Univcoord_T
Stage3end_chroffset (T this);
extern Univcoord_T
Stage3end_chrhigh (T this);
extern Chrpos_T
Stage3end_chrlength (T this);
extern Univcoord_T
Stage3end_genomicstart (T this);
extern Univcoord_T
Stage3end_genomicend (T this);
extern int
Stage3end_query_alignment_length (T this);
extern Chrpos_T
Stage3end_genomic_alignment_length (T this);
extern Chrpos_T
Stage3end_chrpos_low_trim (T this);
extern int
Stage3end_mapq_score (T this);
extern int
Stage3end_absmq_score (T this);
extern int
Stage3end_score (T this);
extern int
Stage3end_gmap_max_match_length (T this);
extern double
Stage3end_gmap_min_splice_prob (T this);
extern int
Stage3end_best_score (List_T hits);
extern bool
Stage3end_equiv_score_unpaired_p (List_T hits, int best_score);
extern int
Stage3end_best_score_paired (List_T hits);
extern int
Stage3end_nmatches_posttrim (T this);
extern int
Stage3end_nmismatches_whole (T this);
extern int
Stage3end_nmismatches_bothdiff (T this);
extern int
Stage3end_nmismatches_refdiff (T this);
extern Endtype_T
Stage3end_start_endtype (T this);
extern Endtype_T
Stage3end_end_endtype (T this);
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
extern int
Stage3end_trim_left (T this);
extern int
Stage3end_trim_right (T this);
extern int
Stage3end_trim_left_raw (T this);
extern int
Stage3end_trim_right_raw (T this);
extern int
Stage3end_circularpos (T this);


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
extern Chrpos_T
Stage3end_distance (T this);
extern Chrpos_T
Stage3end_shortexonA_distance (T this);
extern Chrpos_T
Stage3end_shortexonD_distance (T this);
extern double
Stage3end_chimera_prob (T this);
extern double
Stage3end_shortexon_prob (T this);
extern Univcoord_T
Stage3end_chimera_segmenti_left (T this);
extern Univcoord_T
Stage3end_chimera_segmentj_left (T this);
extern int
Stage3end_chimera_segmenti_cmp (const void *a, const void *b);
extern int
Stage3end_chimera_segmentj_cmp (const void *a, const void *b);
extern int
Stage3end_shortexon_substringD_cmp (const void *a, const void *b);
extern int
Stage3end_shortexon_substringA_cmp (const void *a, const void *b);

extern int
Stage3end_sensedir (T this);
extern int
Stage3end_sensedir_nonamb (T this);
extern int
Stage3end_cdna_direction (T this);
extern int
Stage3end_nintrons (T this);
extern bool
Stage3end_start_ambiguous_p (T this);
extern bool
Stage3end_end_ambiguous_p (T this);
extern int
Stage3end_amb_nmatches_start (T this);
extern int
Stage3end_amb_nmatches_end (T this);

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
extern int
Stage3end_terminal_length (T this);
extern bool
Stage3end_contains_known_splicesite (T this);
extern bool
Stage3end_indel_contains_known_splicesite (bool *leftp, bool *rightp, T this);
extern bool
Stage3end_bad_stretch_p (T this, Compress_T query_compress_fwd, Compress_T query_compress_rev);

extern bool
Stage3end_genomicbound_from_start (Univcoord_T *genomicbound, T this, int overlap, Univcoord_T chroffset);
extern bool
Stage3end_genomicbound_from_end (Univcoord_T *genomicbound, T this, int overlap, Univcoord_T chroffset);

extern void
Stage3end_free (T *old);
extern void
Stage3end_list_free (List_T *values);


extern bool
Stage3pair_anomalous_splice_p (Stage3pair_T this);
extern int
Stage3pair_genestrand (Stage3pair_T this);
extern Stage3end_T
Stage3pair_hit5 (Stage3pair_T this);
extern Stage3end_T
Stage3pair_hit3 (Stage3pair_T this);
extern int
Stage3pair_mapq_score (Stage3pair_T this);
extern int
Stage3pair_absmq_score (Stage3pair_T this);


extern Chrpos_T
Stage3pair_pairlength (Stage3pair_T this);
extern int
Stage3pair_nmatches_posttrim (int *nmatches5, int *nmatches3, Stage3pair_T this);
extern bool
Stage3pair_concordantp (List_T hitpairs);
extern List_T
Stage3pair_filter_nonconcordant (List_T hitpairs);
extern int
Stage3pair_overlap (int *hardclip5, int *hardclip3, Stage3pair_T this);
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
Stage3end_new_exact (int *found_score, Univcoord_T left, int genomiclength, Compress_T query_compress,
		     bool plusp, int genestrand, bool first_read_p,
		     Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh, Chrpos_T chrlength,
		     bool sarrayp);
extern T
Stage3end_new_substitution (int *found_score, int nmismatches, Univcoord_T left,
			    int genomiclength, Compress_T query_compress,
			    bool plusp, int genestrand, bool first_read_p,
			    Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh, Chrpos_T chrlength,
			    bool sarrayp);
extern T
Stage3end_new_insertion (int *found_score, int nindels, int indel_pos, int nmismatches1, int nmismatches2,
			 Univcoord_T left, int genomiclength, Compress_T query_compress,
			 int querylength, bool plusp, int genestrand, bool first_read_p,
			 Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh, Chrpos_T chrlength,
			 int indel_penalty, bool sarrayp);
extern T
Stage3end_new_deletion (int *found_score, int nindels, int indel_pos, int nmismatches1, int nmismatches2,
			Univcoord_T left, int genomiclength, Compress_T query_compress,
			int querylength, bool plusp, int genestrand, bool first_read_p,
			Chrnum_T chrnum, Univcoord_T chroffset,	Univcoord_T chrhigh, Chrpos_T chrlength,
			int indel_penalty, bool sarrayp);

extern T
Stage3end_new_terminal (int querystart, int queryend, Univcoord_T left, Compress_T query_compress,
			int querylength, bool plusp, int genestrand, bool first_read_p,
			Endtype_T start_endtype, Endtype_T end_endtype,
			Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh, Chrpos_T chrlength,
			int max_mismatches_allowed, bool sarrayp);
extern T
Stage3end_new_splice (int *found_score, int donor_nmismatches, int acceptor_nmismatches,
		      Substring_T donor, Substring_T acceptor, Chrpos_T distance,
		      bool shortdistancep, int splicing_penalty, int querylength, int amb_nmatches,
#ifdef LARGE_GENOMES
		      Uint8list_T ambcoords_donor, Uint8list_T ambcoords_acceptor,
#else
		      Uintlist_T ambcoords_donor, Uintlist_T ambcoords_acceptor,
#endif
		      Intlist_T amb_knowni_donor, Intlist_T amb_knowni_acceptor,
		      Intlist_T amb_nmismatches_donor, Intlist_T amb_nmismatches_acceptor,
		      bool copy_donor_p, bool copy_acceptor_p,
		      bool first_read_p, int sensedir, bool sarrayp);
extern T
Stage3end_new_shortexon (int *found_score, Substring_T donor, Substring_T acceptor, Substring_T shortexon,
			 int amb_nmatches_donor, int amb_nmatches_acceptor,
#ifdef LARGE_GENOMES
			 Uint8list_T ambcoords_donor, Uint8list_T ambcoords_acceptor,
#else
			 Uintlist_T ambcoords_donor, Uintlist_T ambcoords_acceptor,
#endif
			 Intlist_T amb_knowni_donor, Intlist_T amb_knowni_acceptor,
			 Intlist_T amb_nmismatches_donor, Intlist_T amb_nmismatches_acceptor,
			 bool copy_donor_p, bool copy_acceptor_p, bool copy_shortexon_p,
			 int splicing_penalty, int querylength, int sensedir, bool sarrayp);

extern T
Stage3end_new_gmap (int nmismatches_whole, int nmatches_posttrim, int max_match_length,
		    int ambig_end_length_5, int ambig_end_length_3,
		    Splicetype_T ambig_splicetype_5, Splicetype_T ambig_splicetype_3,
		    double ambig_prob_5, double ambig_prob_3, double min_splice_prob,
		    struct Pair_T *pairarray, int npairs, int nsegments, int nintrons, int nindelbreaks,
		    Univcoord_T left, int genomiclength, bool plusp, int genestrand, bool first_read_p, int querylength,
		    Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh, Chrpos_T chrlength,
		    int cdna_direction, int sensedir);

extern List_T
Stage3end_sort_bymatches (List_T hits);
extern List_T
Stage3end_sort_by_paired_seenp (List_T hits);

extern Stage3end_T *
Stage3end_eval_and_sort (int *npaths, int *first_absmq, int *second_absmq,
			 Stage3end_T *stage3array, int maxpaths, Shortread_T queryseq,
			 Compress_T query_compress_fwd, Compress_T query_compress_rev,
			 Genome_T genome, char *quality_string, bool displayp);
extern Stage3end_T *
Stage3end_eval_and_sort_guided (int *npaths, int *first_absmq, int *second_absmq, Stage3end_T guide,
				Stage3end_T *stage3array, int maxpaths, Shortread_T queryseq,
				Compress_T query_compress_fwd, Compress_T query_compress_rev,
				Genome_T genome, char *quality_string, bool displayp);
extern List_T
Stage3pair_remove_excess_terminals (List_T hitpairlist);
extern List_T
Stage3end_optimal_score (List_T hitlist, int cutoff_level, int suboptimal_mismatches,
			 Compress_T query_compress_fwd, Compress_T query_compress_rev,
			 int querylength, bool keep_gmap_p, bool finalp);
extern bool
Stage3pair_sense_consistent_p (List_T hitpairlist);
extern List_T
Stage3end_linearize_5 (List_T hitlist);
extern List_T
Stage3end_linearize_3 (List_T hitlist);
extern List_T
Stage3end_remove_circular_alias (List_T hitlist);
extern List_T
Stage3end_remove_duplicates (List_T hitlist);
extern List_T
Stage3end_filter_terminals (List_T hits, int querylength);
extern List_T
Stage3end_remove_overlaps (List_T hitlist, bool finalp);
extern List_T
Stage3end_resolve_multimapping (List_T hitlist);
extern Pairtype_T
Stage3_determine_pairtype (T hit5, T hit3);



/* If hit5 and hit3 are not NULL, then we know this is part of a pair */
extern void
Stage3end_print (FILE *fp, T this, int score, Univ_IIT_T chromosome_iit, Shortread_T queryseq,
		 Shortread_T headerseq, char *acc_suffix, bool invertp,
		 T hit5, T hit3, int pairedlength, int pairscore,
		 Pairtype_T pairtype, int mapq_score);

extern Pairtype_T
Stage3pair_pairtype (Stage3pair_T this);
extern bool
Stage3pair_circularp (Stage3pair_T this);

extern void
Stage3pair_print (Result_T result, Resulttype_T resulttype,
		  Univ_IIT_T chromosome_iit, Shortread_T queryseq1, Shortread_T queryseq2,
		  int maxpaths, bool quiet_if_excessive_p,
		  bool nofailsp, bool failsonlyp, bool fastq_format_p, int quality_shift);

extern Stage3pair_T
Stage3pair_new (T hit5, T hit3, Univcoord_T *splicesites,
		Compress_T query5_compress_fwd, Compress_T query5_compress_rev,
		Compress_T query3_compress_fwd, Compress_T query3_compress_rev,
		int genestrand, Pairtype_T pairtype, int splicing_penalty,
		bool private5p, bool private3p, bool expect_concordant_p);

struct Pair_T *
Stage3pair_merge (int *npairs, int *querylength_merged, char **queryseq_merged, char **quality_merged,
		  Stage3pair_T this, Shortread_T queryseq5, Shortread_T queryseq3,
		  int querylength5, int querylength3, int clipdir, int hardclip5, int hardclip3);

extern void
Stage3pair_privatize (Stage3pair_T *array, int npairs);

extern List_T
Stage3pair_sort_bymatches (List_T hits);

extern List_T
Stage3pair_remove_overlaps (List_T hitpairlist, bool translocp, bool finalp);

extern List_T
Stage3pair_resolve_multimapping (List_T hitpairs);

extern Stage3pair_T *
Stage3pair_eval_and_sort (int *npaths, int *first_absmq, int *second_absmq,
			  Stage3pair_T *stage3pairarray, int maxpaths,
			  Shortread_T queryseq1, Shortread_T queryseq2,
			  Compress_T query5_compress_fwd, Compress_T query5_compress_rev, 
			  Compress_T query3_compress_fwd, Compress_T query3_compress_rev, 
			  Genome_T genome, char *quality_string_5, char *quality_string_3);

extern List_T
Stage3pair_optimal_score (List_T hitpairlist, int cutoff_level, int suboptimal_mismatches,
			  Compress_T query5_compress_fwd, Compress_T query5_compress_rev,
			  Compress_T query3_compress_fwd, Compress_T query3_compress_rev,
			  int querylength5, int querylength3, bool keep_gmap_p, bool finalp);

#if 0
extern List_T
Stage3end_unalias_circular (List_T hitlist);
#endif

extern List_T
Stage3pair_remove_circular_alias (List_T hitpairlist);

extern List_T
Stage3_pair_up_concordant (bool *abort_pairing_p, int *found_score, int *nconcordant, int *nsamechr,
			   List_T *samechr, List_T *conc_transloc, List_T *with_terminal,
			   List_T hitpairs, List_T *hitarray5, int narray5, List_T *hitarray3, int narray3,
			   List_T terminals5, List_T terminals3,
			   int cutoff_level_5, int cutoff_level_3, int subopt_levels,
			   Univcoord_T *splicesites,
			   Compress_T query5_compress_fwd, Compress_T query5_compress_rev,
			   Compress_T query3_compress_fwd, Compress_T query3_compress_rev,
			   int querylength5, int querylength3, int maxpairedpaths,
			   int splicing_penalty, int genestrand);

#undef T
#endif

