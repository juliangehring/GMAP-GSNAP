/* $Id: substring.h 52221 2011-11-11 00:46:16Z twu $ */
#ifndef SUBSTRING_INCLUDED
#define SUBSTRING_INCLUDED

#include <stdio.h>
#include "mode.h"
#include "genomicpos.h"
#include "chrnum.h"
#include "shortread.h"
#include "genome.h"
#include "compress.h"
#include "iit-read.h"
#include "bool.h"
#include "pairdef.h"

typedef enum {END, INS, DEL, DON, ACC, AMB_DON, AMB_ACC, TERM} Endtype_T;

extern void
Substring_setup (bool print_nsnpdiffs_p_in, bool print_snplabels_p_in,
		 bool show_refdiff_p_in, IIT_T snps_iit_in, int *snps_divint_crosstable_in,
		 IIT_T genes_iit_in, int *genes_divint_crosstable_in,
		 IIT_T splicesites_iit_in, int *splicesites_divint_crosstable_in,
		 int donor_typeint_in, int acceptor_typeint_in, int trim_mismatch_score_in,
		 bool output_sam_p_in, Mode_T mode_in);

#define T Substring_T
typedef struct T *T;

extern T
Substring_new (int nmismatches_whole, Chrnum_T chrnum, Genomicpos_T chroffset,
	       Genomicpos_T chrhigh, Genomicpos_T left, Genomicpos_T genomicstart, Genomicpos_T genomicend,
	       Compress_T query_compress, Endtype_T start_endtype, Endtype_T end_endtype,
	       int querystart, int queryend, int querylength,
	       Genomicpos_T alignstart, Genomicpos_T alignend, int genomiclength,
	       int extraleft, int extraright, bool exactp,
	       bool plusp, int genestrand, bool trim_left_p, bool trim_right_p, int minlength);

extern double
Substring_compute_mapq (T this, Compress_T query_compress, char *quality_string);

extern int
Substring_display_prep (char **deletion, T this, char *query, Compress_T query_compress_fwd, Compress_T query_compress_rev,
			Genome_T genome, int deletion_pos, int deletion_length);

extern bool
Substring_bad_stretch_p (T this, Compress_T query_compress_fwd, Compress_T query_compress_rev);

extern void
Substring_free (T *old);

extern bool
Substring_contains_p (T this, int querypos);
extern void
Substring_print_ends (T this, int chrnum);
extern int
Substring_compare (T substring1, T substring2);
extern bool
Substring_overlap_p (T substring1, T substring2);
extern Genomicpos_T
Substring_insert_length (T substring5, T substring3);

extern int
Substring_splicesites_i (T this);
extern int
Substring_splicesites_i_A (T this);
extern int
Substring_splicesites_i_D (T this);

extern bool
Substring_plusp (T this);
extern char *
Substring_genomic_bothdiff (T this);
extern char *
Substring_genomic_refdiff (T this);
extern char *
Substring_genomic_querydir (T this);
extern int
Substring_nmismatches_whole (T this);
extern int
Substring_nmismatches_bothdiff (T this);
extern int
Substring_nmismatches_refdiff (T this);
extern int
Substring_nmatches (T this);
extern int
Substring_nmatches_posttrim (T this);
extern void
Substring_set_nmismatches_terminal (T this, int nmismatches_whole);
extern Endtype_T
Substring_start_endtype (T this);
extern Endtype_T
Substring_end_endtype (T this);
extern double
Substring_mapq_loglik (T this);
extern int
Substring_trim_left (T this);
extern int
Substring_trim_right (T this);
extern int
Substring_querystart (T this);
extern int
Substring_querystart_orig (T this);
extern int
Substring_queryend (T this);
extern int
Substring_querylength (T this);
extern int
Substring_match_length (T this);
extern int
Substring_match_length_orig (T this);
extern Genomicpos_T
Substring_genomic_alignment_length (T this);

extern Chrnum_T
Substring_chrnum (T this);
extern Genomicpos_T
Substring_chroffset (T this);
extern Genomicpos_T
Substring_chrhigh (T this);
extern Genomicpos_T
Substring_alignstart (T this);
extern Genomicpos_T
Substring_alignend (T this);
extern Genomicpos_T
Substring_alignstart_trim (T this);
extern Genomicpos_T
Substring_alignend_trim (T this);
extern Genomicpos_T
Substring_genomicstart (T this);
extern Genomicpos_T
Substring_genomicend (T this);
extern Genomicpos_T
Substring_genomiclength (T this);

extern Genomicpos_T
Substring_chrstart (T this);
extern Genomicpos_T
Substring_chrend (T this);

extern double
Substring_chimera_prob (T this);
extern int
Substring_chimera_pos (T this);
extern int
Substring_chimera_pos_A (T this);
extern int
Substring_chimera_pos_D (T this);
extern bool
Substring_chimera_knownp (T this);
extern int
Substring_nchimera_known (T this);
extern int
Substring_nchimera_novel (T this);
extern bool
Substring_chimera_sensep (T this);


extern T
Substring_copy (T old);

extern T
Substring_new_donor (int splicesites_i, int splicesites_offset, int donor_pos, int donor_nmismatches,
		     double donor_prob, Genomicpos_T left, Compress_T query_compress,
		     int querylength, bool plusp, int genestrand, bool sensep,
		     Chrnum_T chrnum, Genomicpos_T chroffset, Genomicpos_T chrhigh);
extern T
Substring_new_acceptor (int splicesites_i, int splicesites_offset, int acceptor_pos, int acceptor_nmismatches,
			double acceptor_prob, Genomicpos_T left, Compress_T query_compress,
			int querylength, bool plusp, int genestrand, bool sensep,
			Chrnum_T chrnum, Genomicpos_T chroffset, Genomicpos_T chrhigh);
extern T
Substring_new_shortexon (int acceptor_splicesites_i, int donor_splicesites_i, int splicesites_offset,
			 int acceptor_pos, int donor_pos, int nmismatches,
			 double acceptor_prob, double donor_prob, Genomicpos_T left,
			 Compress_T query_compress, int querylength,
			 bool plusp, int genestrand, bool sensep,
			 bool acceptor_ambp, bool donor_ambp,
			 Chrnum_T chrnum, Genomicpos_T chroffset, Genomicpos_T chrhigh);

extern List_T
Substring_sort_chimera_halves (List_T hitlist, bool ascendingp);


extern void
Substring_print_single (FILE *fp, T substring, Shortread_T queryseq,
			char *chr, bool invertp);
extern void
Substring_print_insertion_1 (FILE *fp, T substring1, T substring2, int nindels, 
			     Shortread_T queryseq, char *chr, bool invertp);
extern void
Substring_print_insertion_2 (FILE *fp, T substring1, T substring2, int nindels,
			     Shortread_T queryseq, char *chr, bool invertp);
extern void
Substring_print_deletion_1 (FILE *fp, T substring1, T substring2, int nindels, 
			    char *deletion, Shortread_T queryseq, char *chr,
			    bool invertp);
extern void
Substring_print_deletion_2 (FILE *fp, T substring1, T substring2, int nindels, 
			    Shortread_T queryseq, char *chr, bool invertp);
extern void
Substring_print_donor (FILE *fp, T donor, bool sensep, bool invertp, Shortread_T queryseq,
		       IIT_T chromosome_iit, T acceptor, Genomicpos_T chimera_distance);
extern void 
Substring_print_acceptor (FILE *fp, T acceptor, bool sensep, bool invertp, Shortread_T queryseq,
			  IIT_T chromosome_iit, T donor, Genomicpos_T chimera_distance);
extern void
Substring_print_shortexon (FILE *fp, T shortexon, bool sensep, bool invertp, Shortread_T queryseq,
			   IIT_T chromosome_iit, Genomicpos_T distance1, Genomicpos_T distance2);

extern void
Substring_print_gmap (FILE *fp, struct Pair_T *pairs, int npairs, int nsegments, bool invertedp,
		      Endtype_T start_endtype, Endtype_T end_endtype,
		      Chrnum_T chrnum, Genomicpos_T chroffset, Genomicpos_T chrhigh,
		      int querylength, bool watsonp, int cdna_direction, int score,
		      int insertlength, int pairscore, int mapq_score, IIT_T chromosome_iit);

extern bool
Substring_contains_known_splicesite (T this);

extern Overlap_T
Substring_gene_overlap (T this, bool favor_multiexon_p);

extern long int
Substring_tally (T this, IIT_T tally_iit, int *tally_divint_crosstable);

extern bool
Substring_runlength_p (T this, IIT_T runlength_iit, int *runlength_divint_crosstable);


#ifdef USE_OLD_MAXENT
extern void
Substring_assign_donor_prob (T donor, Genome_T genome, IIT_T chromosome_iit);
extern void
Substring_assign_acceptor_prob (T acceptor, Genome_T genome, IIT_T chromosome_iit);
extern void
Substring_assign_shortexon_prob (T shortexon, Genome_T genome, IIT_T chromosome_iit);
#else
extern void
Substring_assign_donor_prob (T donor);
extern void
Substring_assign_acceptor_prob (T acceptor);
extern void
Substring_assign_shortexon_prob (T shortexon);
#endif

#undef T
#endif


