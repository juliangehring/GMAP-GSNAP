/* $Id: pair.h 53169 2011-11-28 02:56:15Z twu $ */
#ifndef PAIR_INCLUDED
#define PAIR_INCLUDED

typedef struct Pair_T *Pair_T;

#include "bool.h"
#include "genomicpos.h"
#include "chrnum.h"
#include "list.h"
#include "iit-read.h"
#include "sequence.h"
#include "reader.h"		/* For cDNAEnd_T */
#include "uintlist.h"
#include "genome.h"
#include "chimera.h"
#include "substring.h"		/* For Endtype_T */

#define MATCHESPERGAP 3


#define T Pair_T

extern void
Pair_setup (int trim_mismatch_score_in, int trim_indel_score_in);
extern int
Pair_querypos (T this);
extern int
Pair_genomepos (T this);
extern char
Pair_cdna (T this);
extern char
Pair_comp (T this);
extern char
Pair_genome (T this);
extern bool
Pair_gapp (T this);
extern bool
Pair_shortexonp (T this);
extern void
Pair_set_shortexonp (T this);
extern void
Pair_set_genomepos (struct Pair_T *pairarray, int npairs, Genomicpos_T chrpos,
		    Genomicpos_T genomiclength, bool watsonp);
extern List_T
Pair_protect (List_T pairs);

extern T
Pair_new (int querypos, int genomepos, char cdna, char comp, char genome);
extern void
Pair_free (T *old);

extern int
Pair_translation_length (struct T *pairs, int npairs);
extern void
Pair_print_continuous (FILE *fp, struct T *pairs, int npairs, bool watsonp,
		       bool diagnosticp, bool genomefirstp, int invertmode,
		       bool nointronlenp);

extern void
Pair_print_continuous_byexon (FILE *fp, struct T *pairs, int npairs, bool watsonp, bool diagnosticp, int invertmode);
extern void
Pair_print_alignment (FILE *fp, struct T *pairs, int npairs, Chrnum_T chrnum,
		      Genomicpos_T chroffset, IIT_T chromosome_iit, bool watsonp,
		      bool diagnosticp, int invertmode, bool nointronlenp,
		      int wraplength);

extern void
Pair_print_pathsummary (FILE *fp, int pathnum, T start, T end, Chrnum_T chrnum,
			Genomicpos_T chroffset, IIT_T chromosome_iit, bool referencealignp, 
			IIT_T altstrain_iit, char *strain, IIT_T contig_iit, char *dbversion,
			int querylength_given, int skiplength, int trim_start, int trim_end,
			int nexons, int matches, int unknowns, int mismatches, 
			int qopens, int qindels, int topens, int tindels, int goodness,
			bool watsonp, int cdna_direction,
			int translation_start, int translation_end, int translation_length,
			int relaastart, int relaaend, bool maponlyp,
			bool diagnosticp, int stage1_genomicstart, int stage1_genomiclength,
			int stage2_source, int stage2_indexsize, double stage3_defectrate);

extern void
Pair_print_coordinates (FILE *fp, struct T *pairs, int npairs, Chrnum_T chrnum,
			Genomicpos_T chroffset, IIT_T chromosome_iit,
			bool watsonp, int invertmode);

extern void
Pair_dump_one (T this, bool zerobasedp);
extern void
Pair_dump_list (List_T pairs, bool zerobasedp);
extern void
Pair_dump_array (struct T *pairs, int npairs, bool zerobasedp);
extern Genomicpos_T
Pair_genomicpos (struct T *pairs, int npairs, int querypos, bool headp);
extern int
Pair_codon_changepos (struct T *pairs, int npairs, int aapos, int cdna_direction);


extern void
Pair_check_list (List_T pairs);
extern bool
Pair_check_array (struct T *pairs, int npairs);

extern void
Pair_print_exonsummary (FILE *fp, struct T *pairs, int npairs, Chrnum_T chrnum,
			Genomicpos_T chroffset, Genome_T genome, IIT_T chromosome_iit,
			bool watsonp, int cdna_direction, bool genomefirstp, int invertmode);
extern void
Pair_print_gff3 (FILE *fp, struct T *pairs, int npairs, int pathnum, char *accession, 
		 T start, T end, Chrnum_T chrnum, IIT_T chromosome_iit, Sequence_T usersegment,
		 int translation_end, 
		 int querylength_given, int skiplength, int matches, int mismatches, 
		 int qindels, int tindels, bool watsonp, int cdna_direction,
		 bool gff_gene_format_p, bool gff_estmatch_format_p, char *sourcename);

extern void
Pair_print_gsnap (FILE *fp, struct T *pairs, int npairs, int nsegments, bool invertedp,
		  Endtype_T start_endtype, Endtype_T end_endtype,
		  Chrnum_T chrnum, Genomicpos_T chroffset, Genomicpos_T chrhigh,
		  int querylength, bool watsonp, int cdna_direction, int score,
		  int insertlength, int pairscore, int mapq_score,
		  IIT_T chromosome_iit, IIT_T splicesites_iit,
		  int *splicesites_divint_crosstable, int donor_typeint, int acceptor_typeint);
extern int
Pair_gsnap_nsegments (int *total_nmismatches, int *total_nindels, int *nintrons,
		      int *nindelbreaks, struct T *pairs, int npairs);


extern void
Pair_print_sam (FILE *fp, struct T *pairs, int npairs,
		char *accession, Chrnum_T chrnum, IIT_T chromosome_iit, Sequence_T usersegment,
		char *queryseq_ptr, char *quality_string,
		int hardclip5, int hardclip3, int querylength_given,
		bool watsonp, int cdna_direction, int chimera_part, Chimera_T chimera,
		int quality_shift, bool firstp, int pathnum, int npaths, int absmq_score, int second_absmq,
#ifdef GSNAP
		unsigned int flag, int pair_mapq_score, int end_mapq_score,
		Genomicpos_T chrpos, Genomicpos_T mate_chrpos, int pairedlength,
#else
		int mapq_score, bool sam_paired_p,
#endif
		char *sam_read_group_id);

extern void
Pair_print_sam_nomapping (FILE *fp, char *accession, char *queryseq_ptr,
			  char *quality_string, int querylength, int quality_shift,
			  bool firstp, bool sam_paired_p, char *sam_read_group_id);

extern Uintlist_T
Pair_exonbounds (struct T *pairs, int npairs, Genomicpos_T chroffset);
extern void
Pair_print_pslformat_nt (FILE *fp, struct T *pairs, int npairs, T start, T end,
			 Sequence_T queryseq, Chrnum_T chrnum,
			 IIT_T chromosome_iit, Sequence_T usersegment,
			 int matches, int unknowns, int mismatches, 
			 bool watsonp);


extern void
Pair_print_pslformat_pro (FILE *fp, struct T *pairs, int npairs, T start, T end,
			  Sequence_T queryseq, Chrnum_T chrnum,
			  IIT_T chromosome_iit, Sequence_T usersegment,
			  bool watsonp, int cdna_direction);

extern void
Pair_print_exons (FILE *fp, struct T *pairs, int npairs, int wraplength, int ngap, bool cdnap);

extern void
Pair_print_protein_genomic (FILE *fp, struct T *ptr, int npairs, int wraplength, bool forwardp);
#ifdef PMAP
extern void
Pair_print_nucleotide_cdna (FILE *fp, struct T *ptr, int npairs, int wraplength);
#else
extern void
Pair_print_protein_cdna (FILE *fp, struct T *ptr, int npairs, int wraplength, bool forwardp);
#endif

extern void
Pair_print_compressed (FILE *fp, int pathnum, int npaths, T start, T end, Sequence_T queryseq, char *dbversion, 
		       Sequence_T usersegment, int nexons, double fracidentity,
		       struct T *pairs, int npairs, Chrnum_T chrnum,
		       Genomicpos_T chroffset, IIT_T chromosome_iit, int querylength_given,
		       int skiplength, int trim_start, int trim_end, bool checksump,
		       int chimerapos, int chimeraequivpos, double donor_prob, double acceptor_prob,
		       int chimera_cdna_direction, char *strain, bool watsonp, int cdna_direction);

extern void
Pair_print_iit_map (FILE *fp, Sequence_T queryseq, char *accession,
		    T start, T end, Chrnum_T chrnum, IIT_T chromosome_iit);
extern void
Pair_print_iit_exon_map (FILE *fp, struct T *pairs, int npairs, Sequence_T queryseq, char *accession,
			 T start, T end, Chrnum_T chrnum, IIT_T chromosome_iit);
extern void
Pair_print_splicesites (FILE *fp, struct T *pairs, int npairs, char *accession,
			int nexons, Chrnum_T chrnum, IIT_T chromosome_iit, bool watsonp);
extern void
Pair_print_introns (FILE *fp, struct T *pairs, int npairs, char *accession,
		    int nexons, Chrnum_T chrnum, IIT_T chromosome_iit);

extern int
Pair_nmatches (List_T pairs);
extern void
Pair_fracidentity_simple (int *matches, int *unknowns, int *mismatches, List_T pairs);
extern void
Pair_fracidentity (int *matches, int *unknowns, int *mismatches, 
		   int *qopens, int *qindels, int *topens, int *tindels,
		   int *ncanonical, int *nsemicanonical, int *nnoncanonical,
		   List_T pairs, int cdna_direction);
extern int
Pair_fracidentity_score (List_T pairs, int cdna_direction);

extern double
Pair_frac_error (List_T pairs, int cdna_direction);

extern void
Pair_fracidentity_bounded (int *matches, int *unknowns, int *mismatches, 
			   int *qopens, int *qindels, int *topens, int *tindels,
			   int *ncanonical, int *nsemicanonical, int *nnoncanonical,
			   struct T *pairs, int npairs, 
			   int cdna_direction, int minpos, int maxpos);
extern int *
Pair_matchscores (struct T *ptr, int npairs, int querylength);
extern int *
Pair_matchscores_list (int *nmatches, int *ntotal, int *length, List_T pairs);

extern void
Pair_pathscores (bool *gapp, int *pathscores, struct T *ptr, int npairs, 
		 int cdna_direction, int querylength, cDNAEnd_T cdnaend);

extern int
Pair_cdna_direction (List_T pairs);
extern int
Pair_nexons_approx (List_T pairs);
extern int
Pair_nexons (struct T *pairs, int npairs);
extern bool
Pair_consistentp (int *ncanonical, struct T *pairs, int npairs, int cdna_direction);

extern Genomicpos_T
Pair_binary_search_ascending (int *querypos, int lowi, int highi, struct T *pairarray,
			      Genomicpos_T goal_start, Genomicpos_T goal_end);
extern Genomicpos_T
Pair_binary_search_descending (int *querypos, int lowi, int highi, struct T *pairarray,
			       Genomicpos_T goal_start, Genomicpos_T goal_end);
extern Genomicpos_T
Pair_genomicpos_low (int *hardclip_low, int *hardclip_high, struct T *pairarray, int npairs,
		     int querylength, bool watsonp);

extern Genomicpos_T
Pairarray_genomicbound_from_start (struct T *pairarray, int npairs, int overlap);
extern Genomicpos_T
Pairarray_genomicbound_from_end (struct T *pairarray, int npairs, int overlap);

extern List_T
Pair_trim_ends (List_T pairs);

#ifdef GSNAP
extern double
Pair_compute_mapq (struct T *pairarray, int npairs, char *quality_string);
extern Overlap_T
Pair_gene_overlap (struct T *pairarray, int npairs, IIT_T genes_iit, int divno,
		   bool favor_multiexon_p);
extern void
Pair_init (int quality_score_adj_in);
#endif


#undef T
#endif
