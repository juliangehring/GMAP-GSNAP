/* $Id: dynprog.h 52836 2011-11-19 00:20:07Z twu $ */
#ifndef DYNPROG_INCLUDED
#define DYNPROG_INCLUDED

typedef enum {BEST_LOCAL, QUERYEND_INDELS, QUERYEND_NOGAPS} Endalign_T;
typedef struct Dynprog_T *Dynprog_T;

#include "bool.h"
#include "list.h"
#include "pairpool.h"
#include "chrnum.h"
#include "iit-read.h"
#include "splicetrie_build.h"	/* For splicetype */

#ifdef GSNAP
#include "compress.h"
#endif


/* Note: HIGHQ, MEDQ, and LOWQ indicate parameters for high, medium,
   and low sequence quality, respectively */
#define DEFECT_HIGHQ 0.003
#define DEFECT_MEDQ 0.014

#define UNKNOWNJUMP -1000000

#define T Dynprog_T

extern void
Dynprog_setup (IIT_T splicesites_iit_in, int *splicesites_divint_crosstable_in,
	       int donor_typeint_in, int acceptor_typeint_in,
	       Genomicpos_T *splicesites_in, Splicetype_T *splicetypes_in,
	       Genomicpos_T *splicedists_in, int nsplicesites_in,
	       unsigned int *trieoffsets_obs_in, unsigned int *triecontents_obs_in,
	       unsigned int *trieoffsets_max_in, unsigned int *triecontents_max_in);

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
	      int extramaterial_end, int extramaterial_paired);

extern List_T
Dynprog_single_gap (int *dynprogindex, int *finalscore,
		    int *nmatches, int *nmismatches, int *nopens, int *nindels,
		    T dynprog, char *sequence1, char *sequenceuc1, char *sequence2, char *sequenceuc2,
		    int length1, int length2, int offset1, int offset2,
#ifdef PMAP
		    char *queryaaseq,
#endif
		    int cdna_direction, bool jump_late_p, Pairpool_T pairpool,
		    int extraband_single, double defect_rate, int close_indels_mode, bool widebandp);

extern List_T
Dynprog_cdna_gap (int *dynprogindex, int *finalscore, bool *incompletep,
		  T dynprogL, T dynprogR, char *sequence1L, char *sequenceuc1L,
		  char *revsequence1R, char *revsequenceuc1R,
		  char *sequence2, char *sequenceuc2,
		  int length1L, int length1R, int length2,
		  int offset1L, int revoffset1R, int offset2,
#ifdef PMAP
		  char *queryaaseq,
#endif
		  int cdna_direction, bool jump_late_p, Pairpool_T pairpool,
		  int extraband_paired, double defect_rate);

extern List_T
Dynprog_genome_gap (int *dynprogindex, int *finalscore, int *new_leftgenomepos, int *new_rightgenomepos,
		    double *left_prob, double *right_prob,
		    int *nmatches, int *nmismatches, int *nopens, int *nindels,
		    int *exonhead, int *introntype, T dynprogL, T dynprogR, 
		    char *sequence1, char *sequenceuc1,
		    char *sequence2L, char *sequenceuc2L,
		    char *revsequence2R, char *revsequenceuc2R,
		    int length1, int length2L, int length2R, 
		    int offset1, int offset2L, int revoffset2R, 
		    Chrnum_T chrnum, Genomicpos_T chroffset, Genomicpos_T chrpos, Genomicpos_T genomiclength,
		    char *genomicuc_ptr, bool use_genomicseg_p,
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
		  char *revsequence2, char *revsequenceuc2,
		  int length1, int length2, int revoffset1, int revoffset2, 
#ifdef PMAP
		  char *queryaaseq,
#endif
		  int cdna_direction, bool jump_late_p, Pairpool_T pairpool,
		  int extraband_end, double defect_rate, Endalign_T endalign);

extern List_T
Dynprog_end3_gap (int *dynprogindex, int *finalscore, int *nmatches, int *nmismatches,
		  int *nopens, int *nindels, T dynprog,
		  char *sequence1, char *sequenceuc1,
		  char *sequence2, char *sequenceuc2,
		  int length1, int length2, int offset1, int offset2, 
#ifdef PMAP
		  char *queryaaseq,
#endif
		  int cdna_direction, bool jump_late_p, Pairpool_T pairpool,
		  int extraband_end, double defect_rate, Endalign_T endalign);

extern void
Dynprog_make_splicejunction_5 (char *splicejunction, Genomicpos_T splicecoord,
			       int splicelength, Splicetype_T far_splicetype,
			       bool watsonp);

extern void
Dynprog_make_splicejunction_3 (char *splicejunction, Genomicpos_T splicecoord,
			       int splicelength, int contlength, Splicetype_T far_splicetype,
			       bool watsonp);

extern List_T
Dynprog_add_known_splice_5 (int *nmatches_distal, List_T pairs, Genomicpos_T anchor_splicesite, Genomicpos_T far_splicesite,
			    Genomicpos_T chroffset, Genomicpos_T chrpos, int genomiclength,
			    bool watsonp, Pairpool_T pairpool);

extern List_T
Dynprog_add_known_splice_3 (int *nmatches_distal, List_T pairs, Genomicpos_T anchor_splicesite, Genomicpos_T far_splicesite,
			    Genomicpos_T chroffset, Genomicpos_T chrpos, int genomiclength,
			    bool watsonp, Pairpool_T pairpool);


extern List_T
Dynprog_end5_known (bool *knownsplicep, int *dynprogindex, int *finalscore,
		    int *ambig_end_length, Splicetype_T *ambig_splicetype,
		    int *nmatches, int *nmismatches, int *nopens, int *nindels, T dynprog, 
		    char *revsequence1, char *revsequenceuc1,
		    char *revsequence2, char *revsequenceuc2,
		    int length1, int length2, int revoffset1, int revoffset2, 
		    Genomicpos_T chroffset, Genomicpos_T chrpos, int genomiclength,
		    Genomicpos_T knownsplice_limit_low, Genomicpos_T knownsplice_limit_high,
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
		    char *sequence2, char *sequenceuc2,
		    int length1, int length2, int offset1, int offset2, int querylength,
		    Genomicpos_T chroffset, Genomicpos_T chrpos, int genomiclength,
		    Genomicpos_T knownsplice_limit_low, Genomicpos_T knownsplice_limit_high,
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
		       char *sequence2L, char *sequenceuc2L,
		       char *revsequence2R, char *revsequenceuc2R,
		       int length1, int length2L, int length2R,
		       int offset1, int offset2L, int revoffset2R, int cdna_direction,
#ifdef PMAP
		       char *queryaaseq,
#endif
		       char *queryseq, char *queryuc, char *genomicseg, char *genomicuc,
		       Genomicpos_T chroffset, Genomicpos_T chrpos, Genomicpos_T genomiclength, bool watsonp,
		       bool use_genomicseg_p, Pairpool_T pairpool, double defect_rate);

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
