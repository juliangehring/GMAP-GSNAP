/* $Id: stage3hr.h,v 1.64 2010-07-27 00:02:23 twu Exp $ */
#ifndef STAGE3HR_INCLUDED
#define STAGE3HR_INCLUDED

#include "bool.h"
#include "chrnum.h"
#include "genomicpos.h"
#include "iit-read.h"
#include "sequence.h"
#include "genome.h"
#include "resulthr.h"
#include "substring.h"

#define SENSE_NULL 0x0
#define SENSE_ANTI 0x1
#define SENSE_FORWARD 0x2

#define SENSE_CONSISTENT_P(x,y) ((x | y) < 0x3)
#define SENSE_INCONSISTENT_P(x,y) ((x | y) == 0x3)


typedef enum {UNPAIRED, SAMECHR, CONCORDANT} Pairtype_T;

#define T Stage3_T
typedef struct T *T;

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
extern Genomicpos_T
Stage3_genomiclow (T this);
extern int
Stage3_score (T this);
extern int
Stage3_geneprob (T this);
extern int
Stage3_nmismatches (T this);
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
extern Substring_T
Stage3_substring_low (T this);
extern Genomicpos_T
Stage3_distance (T this);
extern int
Stage3_sensedir (T this);


typedef struct Stage3pair_T *Stage3pair_T;

extern void
Stage3hr_print_nsnpdiffs (bool labelsp);
extern void
Stage3hr_print_ncolordiffs ();

extern void
Stage3_free (T *old);

extern Pairtype_T
Stage3pair_pairtype (Stage3pair_T this);
extern Stage3_T
Stage3pair_hit5 (Stage3pair_T this);
extern Stage3_T
Stage3pair_hit3 (Stage3pair_T this);
extern int
Stage3pair_hit5_score (Stage3pair_T this);
extern int
Stage3pair_hit3_score (Stage3pair_T this);
extern Genomicpos_T
Stage3pair_pairlength (Stage3pair_T this);

extern void
Stage3pair_free (Stage3pair_T *old);

extern T
Stage3_new_exact (int *found_score, Genomicpos_T left, int genomiclength, bool plusp, Chrnum_T chrnum, Genomicpos_T chroffset);
extern T
Stage3_new_substitution (int *found_score, int nmismatches, int ncolordiffs, Genomicpos_T left,
			 int genomiclength, bool plusp, char *genomicfwd, char *query, 
			 Chrnum_T chrnum, Genomicpos_T chroffset, int trim_maxlength, bool dibasep, bool cmetp);
extern T
Stage3_new_insertion (int *found_score, int nindels, int indel_pos, int nmismatches1, int nmismatches2,
		      int ncolordiffs1, int ncolordiffs2,
		      Genomicpos_T left, int genomiclength, int querylength, bool plusp,
		      char *genomicfwd, char *query, Chrnum_T chrnum, Genomicpos_T chroffset,
		      int indel_penalty, int trim_maxlength, bool dibasep, bool cmetp);
extern T
Stage3_new_deletion (int *found_score, int nindels, int indel_pos, int nmismatches1, int nmismatches2,
		     int ncolordiffs1, int ncolordiffs2,
		     Genomicpos_T left, int genomiclength, int querylength, bool plusp,
		     char *genomicfwd, char *query, Chrnum_T chrnum, Genomicpos_T chroffset,
		     int indel_penalty, int trim_maxlength, bool dibasep, bool cmetp);
extern T
Stage3_new_terminal (int querystart, int queryend, int nmismatches, int ncolordiffs,
		     Genomicpos_T left, int querylength,
		     bool plusp, char *genomicseg, char *query, 
		     Chrnum_T chrnum, Genomicpos_T chroffset,
		     int terminal_penalty, int trim_maxlength, bool dibasep, bool cmetp);

extern int
Stage3_output_cmp (const void *a, const void *b);
extern List_T
Stage3_optimal_score (List_T hitlist, int cutoff_level, int suboptimal_mismatches);
extern int
Stage3_noptimal (List_T hitlist);
extern List_T
Stage3pair_sort_distance (List_T hitpairlist);
extern List_T
Stage3_mark_ambiguous_splices (bool *ambiguousp, List_T hitlist);
extern List_T
Stage3_remove_duplicates (List_T hitlist);


extern void
Stage3_geneprob_eval (List_T list, IIT_T geneprob_iit, IIT_T chromosome_iit, Sequence_T queryseq);
extern void
Stage3pair_geneprob_eval (List_T list, IIT_T geneprob_iit, IIT_T chromosome_iit,
			  Sequence_T queryseq5, Sequence_T queryseq3);


extern void
Stage3_print_nomapping_sam (Sequence_T queryseq, Resulttype_T resulttype, bool first_read_p,
			    int nhits_mate, Sequence_T queryseq_mate);
extern void
Stage3_print_sam (T this, T mate, int pathnum, int npaths,
		  int score, Genome_T genome, IIT_T chromosome_iit, Sequence_T queryseq,
		  Sequence_T queryseq2, IIT_T snps_iit, int *snps_divint_crosstable,
		  IIT_T splicesites_iit, int donor_typeint, int acceptor_typeint,
		  int *splicesites_divint_crosstable, int pairedlength, Resulttype_T resulttype,
		  bool first_read_p, int npaths_mate);

extern void
Stage3_print (T this, int score, int geneprob, Genome_T genome, IIT_T chromosome_iit, Sequence_T queryseq,
	      IIT_T snps_iit, int *snps_divint_crosstable,
	      IIT_T splicesites_iit, int donor_typeint, int acceptor_typeint,
	      int *splicesites_divint_crosstable, bool invertp,
	      T hit5, T hit3, Pairtype_T pairtype, int pairedlength, int pairscore);

extern void
Stage3_print_paired (Result_T result, Genome_T genome, IIT_T chromosome_iit, Sequence_T queryseq1, Sequence_T queryseq2,
		     Genomicpos_T pairmax, IIT_T snps_iit, int *snps_divint_crosstable,
		     IIT_T splicesites_iit, int donor_typeint, int acceptor_typeint,
		     int *splicesites_divint_crosstable, int maxpaths, bool quiet_if_excessive_p,
		     bool circularp, bool invertp, bool nofailsp, bool failsonlyp);

extern T
Stage3_new_splice (int *found_score, Substring_T donor, Substring_T acceptor, Genomicpos_T distance,
		   bool shortdistancep, int splicing_penalty, int querylength,
		   bool copy_donor_p, bool copy_acceptor_p, int sensedir);
extern List_T
Stage3_filter_bymatch (List_T hitlist);

extern List_T
Stage3pair_remove_samechr (List_T hitpairlist);
extern List_T
Stage3pair_remove_duplicates (List_T hitpairlist);
extern List_T
Stage3pair_optimal_score (List_T hitpairlist, int cutoff_level, int suboptimal_mismatches);

extern List_T
Stage3_pair_up_concordant (bool *abort_pairing_p, int *found_score, int *nconcordant, List_T hitpairs,
			   List_T *hitarray5, int narray5, List_T *hitarray3, int narray3,
			   int cutoff_level_5, int cutoff_level_3, int subopt_levels,
			   Genomicpos_T pairmax, Genomicpos_T expected_pairlength, int querylength5, int querylength3,
			   int maxpairedpaths);

extern List_T
Stage3_pair_up_samechr (int *found_score, List_T hitpairs,
			List_T *hitarray5, int narray5, List_T *hitarray3, int narray3,
			int cutoff_level_5, int cutoff_level_3, int subopt_levels,
			Genomicpos_T expected_pairlength, int querylength5, int querylength3,
			int maxpairedpaths);

#undef T
#endif

