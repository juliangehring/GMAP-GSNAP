/* $Id: stage3.h,v 1.143 2010/02/03 18:18:26 twu Exp $ */
#ifndef STAGE3_INCLUDED
#define STAGE3_INCLUDED
#include "bool.h"
#include "chrnum.h"
#include "genomicpos.h"
#include "list.h"
#include "sequence.h"
#include "genome.h"
#include "stage2.h"
#include "pairdef.h"
#include "pairpool.h"
#include "match.h"
#include "dynprog.h"
#include "iit-read.h"
#include "reader.h"		/* For cDNAEnd_T */
#include "gbuffer.h"
#include "gregion.h"
#include "stopwatch.h"

#define EXTRAQUERYGAP 10

typedef enum {NO_STAGE3DEBUG, POST_STAGE2, POST_SMOOTHING, POST_SINGLES, 
	      POST_INTRONS, POST_CANONICAL, POST_CHANGEPOINT, POST_DISTAL_MEDIAL,
              POST_TRIM_MIDDLE} Stage3debug_T;

#define T Stage3_T
typedef struct T *T;

extern Genomicpos_T
Stage3_genomicstart (T this);
extern Genomicpos_T
Stage3_genomicend (T this);
Genomicpos_T
Stage3_genomiclength (T this);
extern bool
Stage3_watsonp (T this);
extern int
Stage3_cdna_direction (T this);
extern int
Stage3_straintype (T this);
extern int
Stage3_goodness (T this);
extern struct Pair_T *
Stage3_pairarray (int *npairs, T this);
extern int
Stage3_matches (T this);
extern int
Stage3_mismatches (T this);
extern int
Stage3_indels (T this);
extern int
Stage3_querystart (T this);
extern int
Stage3_queryend (T this);
extern int
Stage3_translation_start (T this);
extern int
Stage3_translation_end (T this);
extern int
Stage3_domain (T this);
extern int
Stage3_largemargin (int *newstart, int *newend, T this, int queryntlength);

extern double
Stage3_fracidentity (T this);
extern Genomicpos_T
Stage3_genomicpos (T this, int querypos, bool headp);
extern void
Stage3_pathscores (int *pathscores, T this, int querylength, cDNAEnd_T cdnaend);
extern int
Stage3_chimeric_goodness (int *matches1, int *matches2, T part1, T part2, int breakpoint, int querylength);

extern int
Stage3_cmp (const void *a, const void *b);
extern bool
Stage3_overlap (T x, T y);

extern void
Stage3_recompute_goodness (List_T stage3list);
extern void
Stage3_free (T *old);

extern void
Stage3_genomicbounds (Genomicpos_T *genomicstart, Genomicpos_T *genomiclength, T this);

extern bool
Stage3_test_bounds (T this, int minpos, int maxpos);
extern T
Stage3_apply_bounds (T this, int minpos, int maxpos, bool revertp);

#ifdef PMAP
extern T
Stage3_translate_cdna (T this, Sequence_T queryaaseq, bool strictp);
extern T
Stage3_backtranslate_cdna (T this, bool diagnosticp);
#else
extern T
Stage3_translate_genomic (T this, bool fulllengthp, bool truncatep, bool strictp);
#endif
extern void
Stage3_translate_cdna_via_reference (T this, T reference, bool literalrefp);
extern void
Stage3_fix_cdna_direction (T this, T reference);
extern void
Stage3_print_pathsummary (T this, int pathnum, IIT_T chromosome_iit, IIT_T contig_iit, 
			  IIT_T altstrain_iit, Sequence_T queryseq, bool fulllengthp, 
			  bool truncatep, bool strictp, char *dbversion, int maxmutations, 
			  bool zerobasedp, bool diagnosticp, bool maponlyp);
extern void
Stage3_print_pslformat_nt (T this, int pathnum, IIT_T chromosome_iit, Sequence_T usersegment, Sequence_T queryseq);
#ifdef PMAP
extern void
Stage3_print_pslformat_pro (T this, int pathnum, IIT_T chromosome_iit, Sequence_T usersegment, Sequence_T queryseq, bool strictp);
#endif
extern void
Stage3_print_gff3 (T this, int pathnum, IIT_T chromosome_iit, Sequence_T queryseq, char *dbversion,
		   bool diagnosticp, bool fulllengthp, bool truncatep, bool strictp,
		   bool gff_gene_format_p, char *user_genomicseg);
extern void
Stage3_print_iit_map (T this, int pathnum, IIT_T chromosome_iit, Sequence_T queryseq);
extern void
Stage3_print_iit_exon_map (T this, int pathnum, IIT_T chromosome_iit, Sequence_T queryseq);
extern void
Stage3_print_splicesites (T this, IIT_T chromosome_iit, Sequence_T queryseq);

extern void
Stage3_print_mutations (T this, T reference, IIT_T chromosome_iit, Sequence_T queryseq,
			char *dbversion, bool showalignp, bool zerobasedp, 
			bool continuousp, bool diagnosticp, int proteinmode,
			int invertmode, bool nointronlenp, int wraplength,
			int maxmutations);
extern void
Stage3_print_map (T this, IIT_T map_iit, int *map_divint_crosstable,
		  bool map_iit_universal_p, int map_iit_forward_type,
		  int map_iit_reverse_type, IIT_T chromosome_iit,
		  int pathnum, bool map_exons_p, bool map_bothstrands_p, int nflanking);
extern void
Stage3_print_alignment (T this, Sequence_T queryaaseq, Genome_T genome,
			IIT_T chromosome_iit, bool alignsummaryonlyp, bool universalp, bool zerobasedp,
			bool continuousp, bool continuous_by_exon_p, bool diagnosticp, bool strictp, bool genomefirstp, int proteinmode,
			int invertmode, bool nointronlenp, int wraplength, bool maponlyp);

extern void
Stage3_print_coordinates (T this, Sequence_T queryaaseq, IIT_T chromosome_iit, bool zerobasedp,
			  int invertmode, bool fulllengthp, bool truncatep, bool strictp, bool maponlyp);
extern void
Stage3_print_cdna_exons (T this, int wraplength, int ngap);

extern void
Stage3_print_genomic_exons (T this, int wraplength, int ngap);

extern void
Stage3_print_cdna (T this, Sequence_T queryaaseq, bool fulllengthp, bool truncatep, bool strictp, int wraplength);

extern void
Stage3_print_protein_genomic (T this, Sequence_T queryaaseq, bool fulllengthp, bool truncatep, bool strictp,
			      int wraplength);

extern void
Stage3_print_compressed (T this, Sequence_T queryseq, IIT_T chromosome_iit,
			 char *dbversion, int pathnum, int npaths,
			 bool checksump, int chimerapos, int chimeraequivpos,
			 double donor_prob, double acceptor_prob, 
			 int chimera_cdna_direction, bool zerobasedp, bool truncatep, bool strictp,
			 int worker_id);

extern T
Stage3_compute (List_T path, Genomicpos_T genomicstart, Genomicpos_T genomiclength,
#ifdef PMAP
		Sequence_T queryaaseq,
#endif
		Sequence_T queryseq, Sequence_T queryuc,
		Sequence_T genomicseg, Sequence_T genomicuc,
		int straintype, char *strain,
		Genome_T genome, Genome_T genomealt, Chrnum_T chrnum,
		Genomicpos_T chroffset, Genomicpos_T chrpos, Genomicpos_T chrlength,
		bool watsonp, int maxpeelback, int maxlookback, 
		int extramaterial_end, int extramaterial_paired,
		int extraband_single, int extraband_end, int extraband_paired,
		bool end_microexons_p, int minendexon,
		Pairpool_T pairpool, Gbuffer_T gbuffer, Dynprog_T dynprogL,
		Dynprog_T dynprogM, Dynprog_T dynprogR,
		IIT_T altstrain_iit, int ngap, Stage3debug_T stage3debug,
		bool diagnosticp, bool checkp, Stopwatch_T stopwatch, bool do_final_p, int sense,
		Oligoindex_T *oligoindices, int noligoindices, Diagpool_T diagpool,
		int sufflookback, int nsufflookback, int maxintronlen);

extern T
Stage3_direct (Gregion_T gregion,
#ifdef PMAP
	       Sequence_T queryaaseq,
#endif
	       Sequence_T queryseq, Sequence_T queryuc, Pairpool_T pairpool, Genome_T genome,
	       Chrnum_T chrnum,  Genomicpos_T chroffset, Genomicpos_T chrpos, bool watsonp,
	       int ngap, Gbuffer_T gbuffer, Dynprog_T dynprogL, Dynprog_T dynprogR,
	       int extramaterial_end, int extraband_end);

extern bool
Stage3_mergeable (char *comp, Stage3_T firstpart, Stage3_T secondpart, int chimera_direction, 
		  double donor_prob, double acceptor_prob);

extern void
Stage3_merge (T firstpart, T secondpart, char comp, Pairpool_T pairpool, Genome_T genome, int ngap);

#undef T
#endif
