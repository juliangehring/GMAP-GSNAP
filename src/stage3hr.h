/* $Id: stage3hr.h 37257 2011-03-28 16:39:14Z twu $ */
#ifndef STAGE3HR_INCLUDED
#define STAGE3HR_INCLUDED

#include <stdio.h>
#include "bool.h"
#include "chrnum.h"
#include "genomicpos.h"
#include "intlist.h"
#include "iit-read.h"
#include "shortread.h"
#include "genome.h"
#include "compress.h"
#include "resulthr.h"
#include "substring.h"


/* Should arrange in order of goodness, best to worst */
typedef enum {EXACT, SUB, INSERTION, DELETION, HALFSPLICE_DONOR, HALFSPLICE_ACCEPTOR,
	      SPLICE, ONE_THIRD_SHORTEXON, TWO_THIRDS_SHORTEXON, SHORTEXON, TERMINAL} Hittype_T;

typedef enum {CONCORDANT, PAIRED_INVERSION, PAIRED_SCRAMBLE, PAIRED_TOOLONG, UNPAIRED} Pairtype_T;


#define SENSE_NULL 0x0
#define SENSE_ANTI 0x1
#define SENSE_FORWARD 0x2

#define SENSE_CONSISTENT_P(x,y) ((x | y) != 0x3)
#define SENSE_INCONSISTENT_P(x,y) ((x | y) == 0x3)

#define SENSE_CONSISTENT_FOR_INVERSION_P(x,y) ((x & y) == 0x0)
#define SENSE_INCONSISTENT_FOR_INVERSION_P(x,y) ((x & y) != 0x0)


#define T Stage3_T
typedef struct T *T;

typedef struct Stage3pair_T *Stage3pair_T;


extern void
Stage3hr_setup (bool invert_first_p_in, bool invert_second_p_in);

extern Hittype_T
Stage3_hittype (T this);
extern Chrnum_T
Stage3_chrnum (T this);
extern Genomicpos_T
Stage3_chroffset (T this);
extern Genomicpos_T
Stage3_genomicstart (T this);
extern Genomicpos_T
Stage3_genomicend (T this);
extern int
Stage3_query_alignment_length (T this);
extern Genomicpos_T
Stage3_genomic_alignment_length (T this);
extern Genomicpos_T
Stage3_chrpos_low_trim (T this);
extern int
Stage3_mapq_score (T this);
extern int
Stage3_score (T this);
extern int
Stage3_nmismatches_whole (T this);
extern int
Stage3_nmismatches_bothdiff (T this);
extern int
Stage3_nmismatches_refdiff (T this);
extern int
Stage3_nindels (T this);
extern int
Stage3_indel_pos (T this);
extern bool
Stage3_plusp (T this);

extern Substring_T
Stage3_substring1 (T this);
extern Substring_T
Stage3_substring2 (T this);
extern char *
Stage3_deletion_string (T this);
extern Substring_T
Stage3_substringD (T this);
extern Substring_T
Stage3_substringA (T this);
extern Substring_T
Stage3_substring_low (T this);
extern Genomicpos_T
Stage3_distance (T this);
extern Genomicpos_T
Stage3_shortexon_acceptor_distance (T this);
extern Genomicpos_T
Stage3_shortexon_donor_distance (T this);
extern int
Stage3_sensedir (T this);

extern void
Stage3_free (T *old);



extern Stage3_T
Stage3pair_hit5 (Stage3pair_T this);
extern Stage3_T
Stage3pair_hit3 (Stage3pair_T this);
extern int
Stage3pair_mapq_score (Stage3pair_T this);
extern Genomicpos_T
Stage3pair_pairlength (Stage3pair_T this);

extern void
Stage3pair_free (Stage3pair_T *old);

extern T
Stage3_new_exact (int *found_score, Genomicpos_T left, int genomiclength,
		  Compress_T query_compress, UINT4 *genome_blocks, UINT4 *snp_blocks,
		  bool plusp, Chrnum_T chrnum, Genomicpos_T chroffset);
extern T
Stage3_new_substitution (int *found_score, int nmismatches, int ncolordiffs, Genomicpos_T left,
			 int genomiclength, Compress_T query_compress, UINT4 *genome_blocks, UINT4 *snp_blocks,
			 bool plusp, char *query, Chrnum_T chrnum, Genomicpos_T chroffset,
			 bool dibasep, bool cmetp);
extern T
Stage3_new_insertion (int *found_score, int nindels, int indel_pos, int nmismatches1, int nmismatches2,
		      int ncolordiffs1, int ncolordiffs2, Genomicpos_T left, int genomiclength,
		      Compress_T query_compress, UINT4 *genome_blocks, UINT4 *snp_blocks,
		      int querylength, bool plusp, char *query, Chrnum_T chrnum, Genomicpos_T chroffset,
		      int indel_penalty, bool dibasep, bool cmetp);
extern T
Stage3_new_deletion (int *found_score, int nindels, int indel_pos, int nmismatches1, int nmismatches2,
		     int ncolordiffs1, int ncolordiffs2, Genomicpos_T left, int genomiclength,
		     Compress_T query_compress, UINT4 *genome_blocks, UINT4 *snp_blocks,
		     int querylength, bool plusp, char *query, Chrnum_T chrnum, Genomicpos_T chroffset,
		     int indel_penalty, bool dibasep, bool cmetp);

extern T
Stage3_new_terminal (int querystart, int queryend, int nmismatches, int ncolordiffs,
		     Genomicpos_T left, Compress_T query_compress, UINT4 *genome_blocks, UINT4 *snp_blocks,
		     int querylength, bool plusp, Endtype_T left_endtype, Endtype_T right_endtype, char *query, 
		     Chrnum_T chrnum, Genomicpos_T chroffset,
		     int terminal_penalty, int max_mismatches_allowed,
		     bool dibasep, bool cmetp);
extern T
Stage3_new_splice (int *found_score, int donor_nmismatches, int acceptor_nmismatches,
		   Substring_T donor, Substring_T acceptor, Genomicpos_T distance,
		   bool shortdistancep, int splicing_penalty, int querylength,
		   Intlist_T ambi_left, Intlist_T ambi_right,
		   Intlist_T amb_nmismatches_left, Intlist_T amb_nmismatches_right,
		   bool copy_donor_p, bool copy_acceptor_p, bool first_read_p, int sensedir);
extern T
Stage3_new_shortexon (int *found_score, Substring_T donor, Substring_T acceptor, Substring_T shortexon,
		      Genomicpos_T acceptor_distance, Genomicpos_T donor_distance,
		      Intlist_T ambi_left, Intlist_T ambi_right,
		      Intlist_T amb_nmismatches_left, Intlist_T amb_nmismatches_right,
		      bool copy_donor_p, bool copy_acceptor_p, bool copy_shortexon_p,
		      int splicing_penalty, int querylength, bool first_read_p, int sensedir);

extern Stage3_T *
Stage3_eval_and_sort (Stage3_T *stage3array, int npaths, int maxpaths, Shortread_T queryseq,
		      Compress_T query_compress_fwd, Compress_T query_compress_rev,
		      UINT4 *genome_blocks, UINT4 *snp_blocks, Genome_T genome,
		      char *quality_string, bool dibasep, bool cmetp);
extern List_T
Stage3pair_remove_excess_terminals (List_T hitpairlist);
extern List_T
Stage3_optimal_score (List_T hitlist, int cutoff_level, int suboptimal_mismatches);
extern int
Stage3_noptimal (List_T hitlist);
extern List_T
Stage3_remove_duplicates (List_T hitlist);
extern Pairtype_T
Stage3_determine_pairtype (T hit5, T hit3);
extern Pairtype_T
Stage3pair_pairtype (Stage3pair_T this);


/* If hit5 and hit3 are not NULL, then we know this is part of a pair */
extern void
Stage3_print (FILE *fp, T this, int score, UINT4 *genome_blocks,
	      IIT_T chromosome_iit, Shortread_T queryseq,
	      bool invertp, T hit5, T hit3, int pairedlength, int pairscore,
	      Pairtype_T pairtype, int mapq_score);

extern void
Stage3_print_paired (Result_T result, Resulttype_T resulttype, bool translocationp, UINT4 *genome_blocks,
		     IIT_T chromosome_iit, Shortread_T queryseq1, Shortread_T queryseq2,
		     Genomicpos_T pairmax, int maxpaths, bool quiet_if_excessive_p,
		     bool nofailsp, bool failsonlyp,
		     bool fails_as_input_p, bool fastq_format_p, int quality_shift,
		     FILE *fp_nomapping_1, FILE *fp_nomapping_2,
		     FILE *fp_unpaired_uniq, FILE *fp_unpaired_mult,
		     FILE *fp_halfmapping_uniq, FILE *fp_halfmapping_mult,
		     FILE *fp_paired_uniq_inv, FILE *fp_paired_uniq_scr,
		     FILE *fp_paired_uniq_long, FILE *fp_paired_mult,
		     FILE *fp_concordant_uniq, FILE *fp_concordant_mult);

extern Stage3pair_T
Stage3pair_new (T hit5, T hit3, char *query5, char *query3,
		Genomicpos_T *splicesites, Compress_T query5_compress_fwd, Compress_T query5_compress_rev,
		Compress_T query3_compress_fwd, Compress_T query3_compress_rev,
		UINT4 *genome_blocks, UINT4 *snp_blocks,
		Genomicpos_T expected_pairlength, Genomicpos_T pairlength_deviation,
		Pairtype_T pairtype, int splicing_penalty, bool dibasep, bool cmetp);

extern List_T
Stage3pair_remove_duplicates (List_T hitpairlist);

extern void
Stage3pair_eval (Stage3pair_T *stage3pairarray, int npaths, int maxpaths,
		 Shortread_T queryseq1, Shortread_T queryseq2,
		 Compress_T query5_compress_fwd, Compress_T query5_compress_rev, 
		 Compress_T query3_compress_fwd, Compress_T query3_compress_rev, 
		 UINT4 *genome_blocks, UINT4 *snp_blocks, Genome_T genome,
		 char *quality_string_5, char *quality_string_3, bool dibasep, bool cmetp);

extern List_T
Stage3pair_optimal_score (List_T hitpairlist, int cutoff_level, int suboptimal_mismatches);

extern List_T
Stage3_pair_up_concordant (bool *abort_pairing_p, int *found_score, int *nconcordant,
			   List_T *samechr, List_T hitpairs,
			   List_T *hitarray5, int narray5, List_T *hitarray3, int narray3,
			   int cutoff_level_5, int cutoff_level_3, int subopt_levels,
			   Genomicpos_T *splicesites, char *query5, char *query3,
			   Compress_T query5_compress_fwd, Compress_T query5_compress_rev,
			   Compress_T query3_compress_fwd, Compress_T query3_compress_rev,
			   UINT4 *genome_blocks, UINT4 *snp_blocks, Genomicpos_T pairmax,
			   Genomicpos_T expected_pairlength, Genomicpos_T pairlength_deviation,
			   int querylength5, int querylength3, int maxpairedpaths,
			   bool allow_concordant_translocations_p,
			   int splicing_penalty, bool dibasep, bool cmetp);

#undef T
#endif

