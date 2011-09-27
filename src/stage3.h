/* $Id: stage3.h 36802 2011-03-18 20:42:40Z twu $ */
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

typedef enum {SIMPLE, SUMMARY, ALIGNMENT, COMPRESSED, CONTINUOUS, CONTINUOUS_BY_EXON,
	      EXONS_CDNA, EXONS_GENOMIC, CDNA, PROTEIN_GENOMIC,
	      PSL_NT, PSL_PRO, GFF3_GENE, GFF3_MATCH_CDNA, GFF3_MATCH_EST,
	      SAM, COORDS, SPLICESITES, MAP_GENES, MAP_EXONS} Printtype_T;

typedef enum {NO_STAGE3DEBUG, POST_STAGE2, POST_SMOOTHING, POST_SINGLES, 
	      POST_INTRONS, POST_HMM, POST_DUAL_BREAKS, POST_CYCLES,
	      POST_CANONICAL, POST_CHANGEPOINT, POST_DISTAL_MEDIAL} Stage3debug_T;

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
Stage3_pairarray (T this);
extern int
Stage3_npairs (T this);
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
extern int
Stage3_identity_cmp (const void *a, const void *b);
extern bool
Stage3_overlap (T x, T y);

extern void
Stage3_recompute_goodness (List_T stage3list);
extern void
Stage3_free (T *old, bool free_pairarray_p);

extern void
Stage3_genomicbounds (Genomicpos_T *genomicstart, Genomicpos_T *genomiclength, T this);

extern bool
Stage3_test_bounds (T this, int minpos, int maxpos);
extern T
Stage3_apply_bounds (T this, int minpos, int maxpos, bool revertp);
extern void
Stage3_merge_chimera (T this_left, T this_right, int minpos1, int maxpos1, int minpos2, int maxpos2);

#ifdef PMAP
extern void
Stage3_translate_cdna (T this, Sequence_T queryaaseq, bool strictp);
extern void
Stage3_backtranslate_cdna (T this, bool diagnosticp);
#else
extern void
Stage3_translate_genomic (T this, int npairs, bool fulllengthp, int cds_startpos, int querylength,
			  bool truncatep, bool strictp);
#endif
extern void
Stage3_translate_cdna_via_reference (T this, T reference, bool literalrefp);
extern void
Stage3_fix_cdna_direction (T this, T reference);
extern void
Stage3_translate (T this, Sequence_T queryseq, int querylength, bool fulllengthp,
		  int cds_startpos, bool truncatep, bool strictp,
		  bool diagnosticp, bool maponlyp);
extern void
Stage3_translate_chimera (T this, T mate, Sequence_T queryseq, int querylength, bool fulllengthp,
			  int cds_startpos, bool truncatep, bool strictp,
			  bool diagnosticp, bool maponlyp);
extern void
Stage3_print_pathsummary (FILE *fp, T this, int pathnum, IIT_T chromosome_iit, IIT_T contig_iit, 
			  IIT_T altstrain_iit, Sequence_T queryseq, int querylength, 
			  char *dbversion, int maxmutations, bool diagnosticp, bool maponlyp);
extern void
Stage3_print_pslformat_nt (FILE *fp, T this, int pathnum, IIT_T chromosome_iit, Sequence_T usersegment, Sequence_T queryseq);
#ifdef PMAP
extern void
Stage3_print_pslformat_pro (FILE *fp, T this, int pathnum, IIT_T chromosome_iit, Sequence_T usersegment, Sequence_T queryseq, bool strictp);
#endif
extern void
Stage3_print_gff3 (FILE *fp, T this, int pathnum, IIT_T chromosome_iit, Sequence_T queryseq, int querylength,
		   Printtype_T printtype, char *sourcename);
#ifndef PMAP
extern void
Stage3_print_sam (FILE *fp, T this, int pathnum, int npaths, bool translocationp,
		  IIT_T chromosome_iit, Sequence_T queryseq, int chimera_part,
		  int quality_shift, bool sam_paired_p, bool cigar_noncanonical_splices_p,
		  char *sam_read_group_id);
#endif
extern void
Stage3_print_iit_map (FILE *fp, T this, int pathnum, IIT_T chromosome_iit, Sequence_T queryseq);
extern void
Stage3_print_iit_exon_map (FILE *fp, T this, int pathnum, IIT_T chromosome_iit, Sequence_T queryseq);
extern void
Stage3_print_splicesites (FILE *fp, T this, IIT_T chromosome_iit, Sequence_T queryseq);

extern void
Stage3_print_mutations (FILE *fp, T this, T reference, IIT_T chromosome_iit, Sequence_T queryseq,
			char *dbversion, bool showalignp, bool diagnosticp, int proteinmode,
			int invertmode, bool nointronlenp, int wraplength,
			int maxmutations);
extern void
Stage3_print_map (FILE *fp, T this, IIT_T map_iit, int *map_divint_crosstable, IIT_T chromosome_iit,
		  int pathnum, bool map_exons_p, bool map_bothstrands_p, int nflanking, bool print_comment_p);
extern void
Stage3_print_alignment (FILE *fp, T this, Sequence_T queryaaseq, Genome_T genome,
			IIT_T chromosome_iit, Printtype_T printtype,
			bool continuousp, bool continuous_by_exon_p, bool diagnosticp, bool strictp, bool genomefirstp, int proteinmode,
			int invertmode, bool nointronlenp, int wraplength, bool maponlyp);

extern void
Stage3_print_coordinates (FILE *fp, T this, Sequence_T queryaaseq, IIT_T chromosome_iit,
			  int invertmode, int querylength);
extern void
Stage3_print_cdna (FILE *fp, T this, Sequence_T queryaaseq, int wraplength);

extern void
Stage3_print_protein_genomic (FILE *fp, T this, Sequence_T queryaaseq, int wraplength);

extern void
Stage3_print_compressed (FILE *fp, T this, Sequence_T queryseq, IIT_T chromosome_iit,
			 char *dbversion, int pathnum, int npaths,
			 bool checksump, int chimerapos, int chimeraequivpos,
			 double donor_prob, double acceptor_prob, 
			 int chimera_cdna_direction, bool truncatep, bool strictp);

extern T
Stage3_compute (List_T path, int stage2_source, int stage2_indexsize,
		Genomicpos_T genomicstart, Genomicpos_T genomiclength,
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
		bool diagnosticp, bool checkp, Stopwatch_T stopwatch,
		bool do_final_p, int sense_try, int sense_filter,
		Oligoindex_T *oligoindices_minor, int noligoindices_minor, Diagpool_T diagpool,
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
