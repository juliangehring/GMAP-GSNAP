/* $Id: pair.h,v 1.122 2010/02/03 18:12:54 twu Exp $ */
#ifndef PAIR_INCLUDED
#define PAIR_INCLUDED
#include "bool.h"
#include "genomicpos.h"
#include "chrnum.h"
#include "list.h"
#include "iit-read.h"
#include "sequence.h"
#include "reader.h"		/* For cDNAEnd_T */
#include "uintlist.h"
#include "genome.h"


#define MATCHESPERGAP 3

#define T Pair_T
typedef struct T *T;

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

extern int
Pair_find_genomicpos (Genomicpos_T position, struct T *ptr, int npairs, Genomicpos_T chrpos,
		      Genomicpos_T chroffset, Genomicpos_T genomiclength, bool watsonp);
extern int
Pair_find_chrpos (Genomicpos_T position, struct T *ptr, int npairs, Genomicpos_T chrpos,
		  Genomicpos_T genomiclength, bool watsonp);

extern T
Pair_new (int querypos, int genomepos, char cdna, char comp, char genome);
extern void
Pair_free (T *old);

extern void
Pair_print_continuous (struct T *pairs, int npairs, Chrnum_T chrnum, Genomicpos_T chrpos,
		       Genomicpos_T chroffset, Genomicpos_T genomiclength,
		       bool watsonp, int cdna_direction, bool universalp, bool zerobasedp, bool diagnosticp, 
		       bool genomefirstp, int invertmode, bool nointronlenp);

extern void
Pair_print_continuous_byexon (struct T *pairs, int npairs, bool watsonp, bool diagnosticp, int invertmode);
extern void
Pair_print_alignment (struct T *pairs, int npairs, Chrnum_T chrnum, Genomicpos_T chrpos,
		      Genomicpos_T chroffset, IIT_T chromosome_iit, Genomicpos_T genomiclength, 
		      bool watsonp, int cdna_direction, bool universalp, bool zerobasedp, bool diagnosticp, 
		      bool genomicprop, int invertmode, bool nointronlenp, int wraplength);

extern void
Pair_print_pathsummary (int pathnum, T start, T end, Chrnum_T chrnum, Genomicpos_T chrpos,
			Genomicpos_T chroffset, IIT_T chromosome_iit, bool referencealignp, 
			IIT_T altstrain_iit, char *strain, IIT_T contig_iit, char *dbversion,
			int querylength_given, int skiplength, int trim_start, int trim_end, Genomicpos_T genomiclength,
			int nexons, int matches, int unknowns, int mismatches, 
			int qopens, int qindels, int topens, int tindels, int goodness,
			bool watsonp, int cdna_direction, double defect_rate, 
			int translation_start, int translation_end, int translation_length,
			int relaastart, int relaaend, bool zerobasedp, bool maponlyp,
			bool diagnosticp, int stage1_genomicstart, int stage1_genomiclength,
			double stage2_diag_runtime, double stage2_align_runtime, int stage2_indexsize,
			double stage3_runtime, double stage3_defectrate);

extern void
Pair_print_mutation (struct T *pairs, int npairs, int cdna_direction,
		     int translation_frame, int translation_start, int translation_end,
		     Genomicpos_T chrpos, Genomicpos_T genomiclength, bool watsonp, 
		     bool zerobasedp, Genomicpos_T mutposition, char *change);

extern void
Pair_print_coordinates (struct T *pairs, int npairs, Chrnum_T chrnum, Genomicpos_T chrpos,
			Genomicpos_T chroffset, IIT_T chromosome_iit, Genomicpos_T genomiclength, 
			bool watsonp, bool zerobasedp, int invertmode);

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
Pair_print_exonsummary (struct T *pairs, int npairs, Chrnum_T chrnum, Genomicpos_T chrpos,
			Genomicpos_T chroffset, Genome_T genome, IIT_T chromosome_iit, Genomicpos_T genomiclength,
			bool watsonp, bool universalp, bool zerobasedp, bool genomefirstp, 
			int invertmode);
extern void
Pair_print_gff3 (struct T *pairs, int npairs, int pathnum, char *dbversion, char *accession, 
		 T start, T end, Chrnum_T chrnum, Genomicpos_T chrpos,
		 IIT_T chromosome_iit, Genomicpos_T genomiclength,
		 int translation_start, int translation_end, 
		 int querylength_given, int skiplength, int matches, int unknowns, int mismatches, 
		 int qopens, int qindels, int topens, int tindels, 
		 bool watsonp, int cdna_direction,
		 bool gff_gene_format_p, char *user_genomicseg);

extern Uintlist_T
Pair_exonbounds (struct T *pairs, int npairs, Genomicpos_T chrpos,
		 Genomicpos_T chroffset, Genomicpos_T genomiclength, bool watsonp);
extern void
Pair_print_pslformat_nt (struct T *pairs, int npairs, T start, T end,
			 Sequence_T queryseq, Chrnum_T chrnum, Genomicpos_T chrpos,
			 IIT_T chromosome_iit, Sequence_T usersegment, Genomicpos_T genomiclength,
			 int nexons, int matches, int unknowns, int mismatches, 
			 int qopens, int qindels, int topens, int tindels,
			 bool watsonp, int cdna_direction);


extern void
Pair_print_pslformat_pro (struct T *pairs, int npairs, T start, T end,
			  Sequence_T queryseq, Chrnum_T chrnum, Genomicpos_T chrpos,
			  IIT_T chromosome_iit, Sequence_T usersegment, Genomicpos_T genomiclength, int nexons, 
			  int qopens, int qindels, int topens, int tindels,
			  bool watsonp, int cdna_direction);

extern void
Pair_print_exons (struct T *pairs, int npairs, int wraplength, int ngap, bool cdnap);

extern void
Pair_print_protein_genomic (struct T *ptr, int npairs, int wraplength, bool forwardp);
#ifdef PMAP
extern void
Pair_print_nucleotide_cdna (struct T *ptr, int npairs, int wraplength);
#else
extern void
Pair_print_protein_cdna (struct T *ptr, int npairs, int wraplength, bool forwardp);
#endif

extern void
Pair_print_compressed (int pathnum, int npaths, T start, T end, Sequence_T queryseq, char *dbversion, 
		       int nexons, double fracidentity,
		       struct T *pairs, int npairs, Chrnum_T chrnum, Genomicpos_T chrpos,
		       Genomicpos_T chroffset, IIT_T chromosome_iit, 
		       int querylength_given, int skiplength, int trim_start, int trim_end,
		       Genomicpos_T genomiclength, bool checksump,
		       int chimerapos, int chimeraequivpos, double donor_prob, double acceptor_prob,
		       int chimera_cdna_direction, char *strain, bool watsonp, int cdna_direction,
		       bool zerobasedp, int worker_id);

extern void
Pair_print_iit_map (Sequence_T queryseq, int pathnum, char *accession,
		    T start, T end, Chrnum_T chrnum, Genomicpos_T chrpos,
		    IIT_T chromosome_iit, Genomicpos_T genomiclength, bool watsonp);
extern void
Pair_print_iit_exon_map (struct T *pairs, int npairs, Sequence_T queryseq, int pathnum, char *accession,
			 T start, T end, Chrnum_T chrnum, Genomicpos_T chrpos,
			 IIT_T chromosome_iit, Genomicpos_T genomiclength, bool watsonp);
extern void
Pair_print_splicesites (struct T *pairs, int npairs, char *accession,
			int nexons, Chrnum_T chrnum, Genomicpos_T chrpos,
			IIT_T chromosome_iit, Genomicpos_T genomiclength, bool watsonp);

extern void
Pair_fracidentity_simple (int *matches, int *unknowns, int *mismatches, List_T pairs);
extern void
Pair_fracidentity (int *matches, int *unknowns, int *mismatches, 
		   int *qopens, int *qindels, int *topens, int *tindels,
		   int *ncanonical, int *nsemicanonical, int *nnoncanonical,
		   List_T pairs, int cdna_direction);
extern int
Pair_fracidentity_max (int *changepoint, List_T pairs, int cdna_direction);

extern double
Pair_frac_error (List_T pairs, int cdna_direction);

extern void
Pair_fracidentity_bounded (int *matches, int *unknowns, int *mismatches, 
			   int *qopens, int *qindels, int *topens, int *tindels,
			   int *ncanonical, int *nsemicanonical, int *nnoncanonical,
			   struct T *pairs, int npairs, 
			   int cdna_direction, int minpos, int maxpos);
extern int *
Pair_matchscores (int *nmatches, int *ntotal, int *length, struct T *ptr, int npairs);
extern int *
Pair_matchscores_list (int *nmatches, int *ntotal, int *length, List_T pairs);

extern void
Pair_pathscores (int *pathscores, struct T *ptr, int npairs, 
		 int cdna_direction, int querylength, cDNAEnd_T cdnaend);

extern int
Pair_cdna_direction (List_T pairs);
extern int
Pair_nexons_approx (List_T pairs);
extern int
Pair_nexons (struct T *pairs, int npairs);
extern bool
Pair_consistentp (int *ncanonical, struct T *pairs, int npairs, int cdna_direction);

#undef T
#endif
