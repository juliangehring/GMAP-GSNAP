/* $Id: pair.h,v 1.91 2006/03/05 03:16:30 twu Exp $ */
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
extern void
Pair_flip (T this);

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
		       bool genomefirstp, int invertmode, bool nointronlenp, int ngap);

extern void
Pair_print_alignment (struct T *pairs, int npairs, Chrnum_T chrnum, Genomicpos_T chrpos,
		      Genomicpos_T chroffset, IIT_T chromosome_iit, Genomicpos_T genomiclength, 
		      bool watsonp, int cdna_direction, bool universalp, bool zerobasedp, bool diagnosticp, 
		      bool genomicprop, int invertmode, bool nointronlenp, int wraplength, int ngap);
extern void
Pair_debug_alignment (List_T list, int ngap);

extern void
Pair_print_pathsummary (int pathnum, T start, T end, Chrnum_T chrnum, Genomicpos_T chrpos,
			Genomicpos_T chroffset, IIT_T chromosome_iit, bool referencealignp, 
			IIT_T altstrain_iit, char *strain, 
			IIT_T contig_iit, char *dbversion, Genomicpos_T genomiclength,
			int nexons, double coverage, int matches, int unknowns, int mismatches, 
			int qopens, int qindels, int topens, int tindels, int goodness,
			bool watsonp, int cdna_direction, double defect_rate, 
			int translation_start, int translation_end, int translation_length,
			int relaastart, int relaaend, bool zerobasedp, bool maponlyp);
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
extern void
Pair_check_array (struct T *pairs, int npairs);

extern void
Pair_print_exonsummary (struct T *pairs, int npairs, Chrnum_T chrnum, Genomicpos_T chrpos,
			Genomicpos_T chroffset, IIT_T chromosome_iit, Genomicpos_T genomiclength,
			bool watsonp, bool universalp, bool zerobasedp, bool genomefirstp, 
			int invertmode);
extern Uintlist_T
Pair_exonbounds (struct T *pairs, int npairs, Genomicpos_T chrpos,
		 Genomicpos_T chroffset, Genomicpos_T genomiclength, bool watsonp);
extern void
Pair_print_pslformat_nt (struct T *pairs, int npairs, T start, T end,
			 Sequence_T queryseq, Chrnum_T chrnum, Genomicpos_T chrpos,
			 IIT_T chromosome_iit, Genomicpos_T genomiclength,
			 int nexons, int matches, int unknowns, int mismatches, 
			 int qopens, int qindels, int topens, int tindels,
			 bool watsonp, int cdna_direction);


extern void
Pair_print_pslformat_pro (struct T *pairs, int npairs, T start, T end,
			  Sequence_T queryseq, Chrnum_T chrnum, Genomicpos_T chrpos,
			  IIT_T chromosome_iit, Genomicpos_T genomiclength, int nexons, 
			  int qopens, int qindels, int topens, int tindels,
			  bool watsonp, int cdna_direction);

extern void
Pair_print_exons (struct T *pairs, int npairs, int wraplength, int ngap, bool cdnap);

extern void
Pair_print_protein_genomic (struct T *ptr, int npairs, int wraplength);
#ifdef PMAP
extern void
Pair_print_nucleotide_cdna (struct T *ptr, int npairs, int wraplength);
#else
extern void
Pair_print_protein_cdna (struct T *ptr, int npairs, int wraplength);
#endif

extern void
Pair_print_compressed (Sequence_T queryseq, char *version, int pathnum, int npaths, 
		       int nexons, double coverage, double fracidentity,
		       struct T *pairs, int npairs, Chrnum_T chrnum, Genomicpos_T chrpos,
		       Genomicpos_T chroffset, IIT_T chromosome_iit,
		       Genomicpos_T genomiclength, bool checksump, 
		       int chimerapos, int chimeraequivpos, double donor_prob, double acceptor_prob,
		       int chimera_cdna_direction, char *strain, bool watsonp, int cdna_direction, 
		       bool zerobasedp, int worker_id);

extern void
Pair_fracidentity (int *matches, int *unknowns, int *mismatches, 
		   int *qopens, int *qindels, int *topens, int *tindels,
		   int *ncanonical, int *nsemicanonical, int *nnoncanonical,
		   List_T pairs, int cdna_direction);
extern void
Pair_fracidentity_bounded (int *matches, int *unknowns, int *mismatches, 
			   int *qopens, int *qindels, int *topens, int *tindels,
			   int *ncanonical, int *nsemicanonical, int *nnoncanonical,
			   struct T *pairs, int npairs, 
			   int cdna_direction, int minpos, int maxpos);
extern int *
Pair_matchscores (struct T *ptr, int npairs, 
		  int cdna_direction, int querylength);
extern void
Pair_pathscores (int *pathscores, struct T *ptr, int npairs, 
		 int cdna_direction, int querylength, cDNAEnd_T cdnaend);

extern int
Pair_cdna_direction (List_T pairs);
extern int
Pair_nexons (List_T pairs);
extern bool
Pair_consistentp (int *ncanonical, struct T *pairs, int npairs, int cdna_direction);

extern struct T *
Pair_block_copy (List_T list, int *npairs, int ngap);

#undef T
#endif
