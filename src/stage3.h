/* $Id: stage3.h 52863 2011-11-20 23:53:35Z twu $ */
#ifndef STAGE3_INCLUDED
#define STAGE3_INCLUDED

typedef struct Stage3_T *Stage3_T;

#include "bool.h"
#include "sense.h"
#include "chrnum.h"
#include "genomicpos.h"
#include "list.h"
#include "sequence.h"
#include "genome.h"
#include "stage2.h"
#include "pairdef.h"
#include "pairpool.h"
#include "splicetrie.h"
#include "splicetrie_build.h"	/* For Splicetype_T */
#include "dynprog.h"
#include "iit-read.h"
#include "reader.h"		/* For cDNAEnd_T */
#include "chimera.h"
#include "stopwatch.h"

#ifndef GSNAP
#include "gregion.h"
#endif

#define EXTRAQUERYGAP 10

typedef enum {SIMPLE, SUMMARY, ALIGNMENT, COMPRESSED, CONTINUOUS, CONTINUOUS_BY_EXON,
	      EXONS_CDNA, EXONS_GENOMIC, CDNA, PROTEIN_GENOMIC,
	      PSL_NT, PSL_PRO, GFF3_GENE, GFF3_MATCH_CDNA, GFF3_MATCH_EST,
	      SAM, COORDS, SPLICESITES, INTRONS, MAP_GENES, MAP_EXONS} Printtype_T;

typedef enum {NO_STAGE3DEBUG, POST_STAGE2, POST_SMOOTHING, POST_SINGLES, 
	      POST_INTRONS, POST_HMM, POST_DUAL_BREAKS, POST_CYCLES,
	      POST_CANONICAL, POST_CHANGEPOINT, POST_DISTAL_MEDIAL} Stage3debug_T;

#define T Stage3_T

extern void
Stage3_setup (bool splicingp_in,
	      IIT_T splicesites_iit_in, int *splicesites_divint_crosstable_in,
	      int donor_typeint_in, int acceptor_typeint_in,
	      Genomicpos_T *splicesites_in,
	      int min_intronlength_in, int max_deletionlength_in,
	      int expected_pairlength_in, int pairlength_deviation_in);

extern bool
Stage3_watsonp (T this);
extern int
Stage3_cdna_direction (T this);
extern int
Stage3_straintype (T this);
extern int
Stage3_goodness (T this);
extern int
Stage3_absmq_score (T this);
extern int
Stage3_mapq_score (T this);
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
extern Chrnum_T
Stage3_chrnum (T this);
extern Genomicpos_T
Stage3_chrstart (T this);
extern Genomicpos_T
Stage3_chrend (T this);
extern Genomicpos_T
Stage3_genomicstart (T this);
extern Genomicpos_T
Stage3_genomicend (T this);

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
Stage3_pathscores (bool *gapp, int *pathscores, T this, int querylength, cDNAEnd_T cdnaend);
extern int
Stage3_chimeric_goodness (int *matches1, int *matches2, T part1, T part2, int breakpoint);

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
			  IIT_T altstrain_iit, Sequence_T queryseq,
			  char *dbversion, int maxmutations, bool diagnosticp, bool maponlyp);
extern void
Stage3_print_pslformat_nt (FILE *fp, T this, IIT_T chromosome_iit, Sequence_T usersegment, Sequence_T queryseq);
#ifdef PMAP
extern void
Stage3_print_pslformat_pro (FILE *fp, T this, IIT_T chromosome_iit, Sequence_T usersegment, Sequence_T queryseq, bool strictp);
#endif
extern void
Stage3_print_gff3 (FILE *fp, T this, int pathnum, IIT_T chromosome_iit, Sequence_T usersegment,
		   Sequence_T queryseq, int querylength, Printtype_T printtype, char *sourcename);
#ifndef PMAP
extern void
Stage3_print_sam (FILE *fp, T this, int pathnum, int npaths,
		  int absmq_score, int second_mapq, int mapq_score,
		  IIT_T chromosome_iit, Sequence_T usersegment,
		  Sequence_T queryseq, int chimera_part, Chimera_T chimera,
		  int quality_shift, bool sam_paired_p, char *sam_read_group_id);
#endif
extern void
Stage3_print_iit_map (FILE *fp, T this, IIT_T chromosome_iit, Sequence_T queryseq);
extern void
Stage3_print_iit_exon_map (FILE *fp, T this, IIT_T chromosome_iit, Sequence_T queryseq);
extern void
Stage3_print_splicesites (FILE *fp, T this, IIT_T chromosome_iit, Sequence_T queryseq);
extern void
Stage3_print_introns (FILE *fp, T this, IIT_T chromosome_iit, Sequence_T queryseq);

extern void
Stage3_print_mutations (FILE *fp, T this, T reference, IIT_T chromosome_iit, Sequence_T queryseq,
			char *dbversion, bool showalignp, bool diagnosticp,
			int invertmode, bool nointronlenp, int wraplength,
			int maxmutations);
extern void
Stage3_print_map (FILE *fp, T this, IIT_T map_iit, int *map_divint_crosstable, IIT_T chromosome_iit,
		  int pathnum, bool map_exons_p, bool map_bothstrands_p, int nflanking, bool print_comment_p);
extern void
Stage3_print_alignment (FILE *fp, T this, Sequence_T queryaaseq, Genome_T genome,
			IIT_T chromosome_iit, Printtype_T printtype,
			bool continuousp, bool continuous_by_exon_p, bool diagnosticp, bool strictp, bool genomefirstp,
			int invertmode, bool nointronlenp, int wraplength, bool maponlyp);

extern void
Stage3_print_coordinates (FILE *fp, T this, Sequence_T queryaaseq, IIT_T chromosome_iit,
			  int invertmode);
extern void
Stage3_print_cdna (FILE *fp, T this, Sequence_T queryaaseq, int wraplength);

extern void
Stage3_print_protein_genomic (FILE *fp, T this, Sequence_T queryaaseq, int wraplength);

extern void
Stage3_print_compressed (FILE *fp, T this, Sequence_T queryseq, IIT_T chromosome_iit,
			 char *dbversion, Sequence_T usersegment, int pathnum, int npaths,
			 bool checksump, int chimerapos, int chimeraequivpos,
			 double donor_prob, double acceptor_prob, 
			 int chimera_cdna_direction, bool truncatep, bool strictp);
extern T
Stage3_new (struct Pair_T *pairarray, List_T pairs, int npairs, int cdna_direction,
	    Genomicpos_T stage1_genomicstart, Genomicpos_T stage1_genomiclength, 
	    int stage2_source, int stage2_indexsize,
	    int matches, int unknowns, int mismatches, int qopens, int qindels,
	    int topens, int tindels, int ncanonical, int nsemicanonical,
	    int nnoncanonical, double defect_rate,
	    Chrnum_T chrnum, Genomicpos_T chroffset, Genomicpos_T chrpos, bool watsonp,
	    int skiplength, int trimlength, double stage3_runtime,
	    int straintype, char *strain, IIT_T altstrain_iit);

extern struct Pair_T *
Stage3_compute (List_T *pairs, int *npairs, int *cdna_direction, int *sensedir, int *matches,
		int *nmatches_pretrim, int *nmatches_posttrim,
		int *ambig_end_length_5, int *ambig_end_length_3,
		Splicetype_T *ambig_splicetype_5, Splicetype_T *ambig_splicetype_3,
		int *unknowns, int *mismatches, int *qopens, int *qindels, int *topens, int *tindels,
		int *ncanonical, int *nsemicanonical, int *nnoncanonical, 
		double *defect_rate, List_T path, int genomiclength,
#ifdef GSNAP
#ifdef END_KNOWNSPLICING_SHORTCUT
		int cutoff_level, char *queryptr, Compress_T query_compress,
		UINT4 *splicefrags_ref, UINT4 *splicefrags_alt,
#endif
#endif
#ifdef PMAP
		char *queryaaseq_ptr,
#endif
		char *queryseq_ptr, char *queryuc_ptr, int querylength,
		int skiplength, int query_subseq_offset,
		char *genomicseg_ptr, char *genomicuc_ptr,
		Chrnum_T chrnum, Genomicpos_T chroffset, Genomicpos_T chrpos,
		Genomicpos_T knownsplice_limit_low, Genomicpos_T knownsplice_limit_high,
		Genome_T genome, bool usersegment_p, bool watsonp, bool jump_late_p,
		int maxpeelback, int maxpeelback_distalmedial, int nullgap,
		int extramaterial_end, int extramaterial_paired,
		int extraband_single, int extraband_end, int extraband_paired, int minendexon,
		Pairpool_T pairpool, Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,
		int ngap, Stage3debug_T stage3debug, bool diagnosticp, bool checkp,
		bool do_final_p, int sense_try, int sense_filter,
		Oligoindex_T *oligoindices_minor, int noligoindices_minor, Diagpool_T diagpool,
		int sufflookback, int nsufflookback, int maxintronlen, int close_indels_mode,
		int paired_favor_mode, int zero_offset);

#ifndef GSNAP
extern T
Stage3_direct (Gregion_T gregion,
#ifdef PMAP
	       Sequence_T queryaaseq,
#endif
	       Sequence_T queryseq, Sequence_T queryuc, Pairpool_T pairpool, Genome_T genome,
	       Chrnum_T chrnum,  Genomicpos_T chroffset, Genomicpos_T chrpos, bool watsonp,
	       int ngap, Dynprog_T dynprogL, Dynprog_T dynprogR,
	       int extramaterial_end, int extraband_end);
#endif

extern bool
Stage3_mergeable (char *comp, Genomicpos_T *genomegap, Stage3_T firstpart, Stage3_T secondpart,
		  int exonexonpos, int queryntlength, int chimera_direction,
		  double donor_prob, double acceptor_prob);

extern void
Stage3_merge_readthrough (T this_left, T this_right, char comp, Genomicpos_T genomegap,
			  int minpos1, int maxpos1, int minpos2, int maxpos2,
			  Pairpool_T pairpool, Genome_T genome, int ngap);

#undef T
#endif
