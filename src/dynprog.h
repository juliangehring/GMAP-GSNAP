/* $Id: dynprog.h 106269 2013-08-29 18:30:43Z twu $ */
#ifndef DYNPROG_INCLUDED
#define DYNPROG_INCLUDED

/* BEST_LOCAL is a local alignment, whereas QUERYEND_INDELS and
   QUERYEND_NOGAPS are global.  QUERYEND_GAP allows an intron at the
   end */
typedef enum {QUERYEND_GAP, QUERYEND_INDELS, QUERYEND_NOGAPS, BEST_LOCAL} Endalign_T;
typedef struct Dynprog_T *Dynprog_T;

#include "bool.h"
#include "list.h"
#include "pairpool.h"
#include "chrnum.h"
#include "iit-read.h"
#include "splicetrie_build.h"	/* For splicetype */
#include "genome.h"
#include "mode.h"

#ifdef GSNAP
#include "compress.h"
#endif


/* Note: HIGHQ, MEDQ, and LOWQ indicate parameters for high, medium,
   and low sequence quality, respectively */
#define DEFECT_HIGHQ 0.003
#define DEFECT_MEDQ 0.014

#define UNKNOWNJUMP -1000000

#define T Dynprog_T

extern char *
Dynprog_endalign_string (Endalign_T endalign);

extern void
Dynprog_setup (bool novelsplicingp_in,
	       IIT_T splicing_iit_in, int *splicing_divint_crosstable_in,
	       int donor_typeint_in, int acceptor_typeint_in,
	       Univcoord_T *splicesites_in, Splicetype_T *splicetypes_in,
	       Chrpos_T *splicedists_in, int nsplicesites_in,
	       Trieoffset_T *trieoffsets_obs_in, Triecontent_T *triecontents_obs_in,
	       Trieoffset_T *trieoffsets_max_in, Triecontent_T *triecontents_max_in,
	       bool homopolymerp_in);

extern int
Dynprog_score (int matches, int mismatches, int qopens, int qindels, int topens, int tindels,
	       double defect_rate);

extern T
Dynprog_new (int maxlookback, int extraquerygap, int maxpeelback,
	     int extramaterial_end, int extramaterial_paired);
extern void
Dynprog_free (T *old);

#ifdef PMAP
extern char
Dynprog_codon_char (char aa, int codonpos);
#endif

extern int
Dynprog_pairdistance (int c1, int c2);

extern void
Dynprog_term (void);
extern void
Dynprog_init (int maxlookback, int extraquerygap, int maxpeelback,
	      int extramaterial_end, int extramaterial_paired, Mode_T mode);

extern List_T
Dynprog_single_gap (int *dynprogindex, int *finalscore,
		    int *nmatches, int *nmismatches, int *nopens, int *nindels,
		    T dynprog, char *sequence1, char *sequenceuc1,
		    int length1, int length2, int offset1, int offset2,
		    Univcoord_T chroffset, Univcoord_T chrhigh,
#ifdef PMAP
		    char *queryaaseq,
#endif
		    bool watsonp, bool jump_late_p, Pairpool_T pairpool,
		    int extraband_single, double defect_rate, bool widebandp);

extern List_T
Dynprog_cdna_gap (int *dynprogindex, int *finalscore, bool *incompletep,
		  T dynprogL, T dynprogR, char *sequence1L, char *sequenceuc1L,
		  char *revsequence1R, char *revsequenceuc1R,
#if 0
		  char *sequence2, char *sequenceuc2,
#endif
		  int length1L, int length1R, int length2,
		  int offset1L, int revoffset1R, int offset2,
		  Univcoord_T chroffset, Univcoord_T chrhigh,
#ifdef PMAP
		  char *queryaaseq,
#endif
		  int cdna_direction, bool watsonp, bool jump_late_p, Pairpool_T pairpool,
		  int extraband_paired, double defect_rate);

extern List_T
Dynprog_genome_gap (int *dynprogindex, int *finalscore, int *new_leftgenomepos, int *new_rightgenomepos,
		    double *left_prob, double *right_prob,
		    int *nmatches, int *nmismatches, int *nopens, int *nindels,
		    int *exonhead, int *introntype, T dynprogL, T dynprogR, 
		    char *sequence1, char *sequenceuc1,
		    int length1, int length2L, int length2R, 
		    int offset1, int offset2L, int revoffset2R, 
		    Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
#ifdef PMAP
		    char *queryaaseq,
#endif
		    int cdna_direction, bool watsonp, bool jump_late_p, Pairpool_T pairpool, int extraband_paired,
		    double defect_rate, int maxpeelback, bool halfp, bool finalp, bool use_probabilities_p,
		    int score_threshold, bool splicingp);

extern List_T
Dynprog_end5_gap (int *dynprogindex, int *finalscore, int *nmatches, int *nmismatches,
		  int *nopens, int *nindels, T dynprog,
		  char *revsequence1, char *revsequenceuc1,
		  int length1, int length2, int revoffset1, int revoffset2, 
		  Univcoord_T chroffset, Univcoord_T chrhigh,
#ifdef PMAP
		  char *queryaaseq,
#endif
		  int cdna_direction, bool watsonp, bool jump_late_p, Pairpool_T pairpool,
		  int extraband_end, double defect_rate, Endalign_T endalign);

extern List_T
Dynprog_end5_splicejunction (int *dynprogindex, int *finalscore, int *missscore,
			     int *nmatches, int *nmismatches, int *nopens, int *nindels, T dynprog, 
			     char *revsequence1, char *revsequenceuc1,
			     char *revsequence2, char *revsequenceuc2, char *revsequencealt2,
			     int length1, int length2, int revoffset1, int revoffset2_anchor, int revoffset2_far,
			     Univcoord_T chroffset, Univcoord_T chrhigh,
#ifdef PMAP
			     char *queryaaseq,
#endif
			     int cdna_direction, bool watsonp, bool jump_late_p, Pairpool_T pairpool,
			     int extraband_end, double defect_rate, int contlength);

extern List_T
Dynprog_end3_gap (int *dynprogindex, int *finalscore, int *nmatches, int *nmismatches,
		  int *nopens, int *nindels, T dynprog,
		  char *sequence1, char *sequenceuc1,
		  int length1, int length2, int offset1, int offset2, 
		  Univcoord_T chroffset, Univcoord_T chrhigh,
#ifdef PMAP
		  char *queryaaseq,
#endif
		  int cdna_direction, bool watsonp, bool jump_late_p, Pairpool_T pairpool,
		  int extraband_end, double defect_rate, Endalign_T endalign);

extern List_T
Dynprog_end3_splicejunction (int *dynprogindex, int *finalscore, int *missscore,
			     int *nmatches, int *nmismatches, int *nopens, int *nindels, T dynprog,
			     char *sequence1, char *sequenceuc1,
			     char *sequence2, char *sequenceuc2, char *sequencealt2,
			     int length1, int length2, int offset1, int offset2_anchor, int offset2_far,
			     Univcoord_T chroffset, Univcoord_T chrhigh,
#ifdef PMAP
			     char *queryaaseq,
#endif
			     int cdna_direction, bool watsonp, bool jump_late_p, Pairpool_T pairpool,
			     int extraband_end, double defect_rate, int contlength);

extern bool
Dynprog_make_splicejunction_5 (char *splicejunction, char *splicejunction_alt, Univcoord_T splicecoord,
			       int splicelength, int contlength, Splicetype_T far_splicetype,
			       bool watsonp);

extern bool
Dynprog_make_splicejunction_3 (char *splicejunction, char *splicejunction_alt, Univcoord_T splicecoord,
			       int splicelength, int contlength, Splicetype_T far_splicetype,
			       bool watsonp);

#if 0
extern List_T
Dynprog_add_known_splice_5 (int *nmatches_distal, List_T pairs, Univcoord_T anchor_splicesite, Univcoord_T far_splicesite,
			    Univcoord_T chroffset, Univcoord_T chrhigh, bool watsonp, Pairpool_T pairpool);

extern List_T
Dynprog_add_known_splice_3 (int *nmatches_distal, List_T pairs, Univcoord_T anchor_splicesite, Univcoord_T far_splicesite,
			    Univcoord_T chroffset, Univcoord_T chrhigh, bool watsonp, Pairpool_T pairpool);
#endif


extern List_T
Dynprog_end5_known (bool *knownsplicep, int *dynprogindex, int *finalscore,
		    int *ambig_end_length, Splicetype_T *ambig_splicetype,
		    int *nmatches, int *nmismatches, int *nopens, int *nindels, T dynprog, 
		    char *revsequence1, char *revsequenceuc1,
		    int length1, int length2, int revoffset1, int revoffset2, 
		    Univcoord_T chroffset, Univcoord_T chrhigh,
		    Univcoord_T knownsplice_limit_low, Univcoord_T knownsplice_limit_high,
#ifdef PMAP
		    char *queryaaseq,
#endif
#ifdef GSNAP
#ifdef END_KNOWNSPLICING_SHORTCUT
		    int cutoff_level, char *queryptr, int querylength, Compress_T query_compress,
#endif
#endif
		    int cdna_direction, bool watsonp, bool jump_late_p,
		    Pairpool_T pairpool, int extraband_end, double defect_rate);

extern List_T
Dynprog_end3_known (bool *knownsplicep, int *dynprogindex, int *finalscore,
		    int *ambig_end_length, Splicetype_T *ambig_splicetype,
		    int *nmatches, int *nmismatches, int *nopens, int *nindels, T dynprog, 
		    char *sequence1, char *sequenceuc1,
		    int length1, int length2, int offset1, int offset2, int querylength,
		    Univcoord_T chroffset, Univcoord_T chrhigh,
		    Univcoord_T knownsplice_limit_low, Univcoord_T knownsplice_limit_high,
#ifdef PMAP
		    char *queryaaseq,
#endif
#ifdef GSNAP
#ifdef END_KNOWNSPLICING_SHORTCUT
		    int cutoff_level, char *queryptr, int querylength, Compress_T query_compress,
#endif
#endif
		    int cdna_direction, bool watsonp, bool jump_late_p,
		    Pairpool_T pairpool, int extraband_end, double defect_rate);


#if 0
extern int
Dynprog_internal_gap_stats (T dynprog, char *sequenceuc1, char *sequenceuc2,
			    int length1, int length2, int offset1, int offset2, 
#ifdef PMAP
			    char *queryaaseq,
#endif
			    int cdna_direction, int extraband_end, double defect_rate);
#endif

#if 0
extern List_T
Dynprog_dual_break (int *dynprogindex, int *finalscore, T dynprogL, T dynprogR, 
		    char *sequence1L, char *sequenceuc1L,
		    char *sequence2L, char *sequenceuc2L,
		    char *revsequence1R, char *revsequenceuc1R,
		    char *revsequence2R, char *revsequenceuc2R,
		    int length1, int length2, int offset1L, int offset2L, 
		    int revoffset1R, int revoffset2R, 
#ifdef PMAP
		    char *queryaaseq,
#endif
		    int cdna_direction, bool watsonp, bool jump_late_p, Pairpool_T pairpool,
		    int extraband_end, double defect_rate);
#endif

extern List_T
Dynprog_microexon_int (double *bestprob2, double *bestprob3, int *dynprogindex, int *microintrontype,
		       char *sequence1, char *sequenceuc1,
		       int length1, int length2L, int length2R,
		       int offset1, int offset2L, int revoffset2R, int cdna_direction,
#ifdef PMAP
		       char *queryaaseq, char *genomicuc,
#endif
		       char *queryseq, char *queryuc,
		       Univcoord_T chroffset, Univcoord_T chrhigh, bool watsonp,
		       Pairpool_T pairpool, double defect_rate);

extern List_T
Dynprog_microexon_5 (int *dynprogindex, int *microintrontype, int *microexonlength,
		     char *revsequence1, char *revsequenceuc1,
		     char *revsequence2, char *revsequenceuc2,
		     int length1, int length2, int revoffset1, int revoffset2, int cdna_direction,
#ifdef PMAP
		     char *queryaaseq,
#endif
		     char *queryseq, char *queryuc, char *genomicseg, char *genomicuc,
		     Pairpool_T pairpool, bool end_microexons_p);

extern List_T
Dynprog_microexon_3 (int *dynprogindex, int *microintrontype, int *microexonlength,
		     char *sequence1, char *sequenceuc1,
		     char *sequence2, char *sequenceuc2,
		     int length1, int length2, int offset1, int offset2, int cdna_direction,
#ifdef PMAP
		     char *queryaaseq,
#endif
		     char *queryseq, char *queryuc, char *genomicseg, char *genomicuc,
		     int genomiclength, Pairpool_T pairpool, bool end_microexons_p);

#undef T
#endif
