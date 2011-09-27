/* $Id: dynprog.h,v 1.53 2007/04/23 15:52:16 twu Exp $ */
#ifndef DYNPROG_INCLUDED
#define DYNPROG_INCLUDED
#include "bool.h"
#include "list.h"
#include "pairpool.h"

/* Note: HIGHQ, MEDQ, and LOWQ indicate parameters for high, medium,
   and low sequence quality, respectively */
#define DEFECT_HIGHQ 0.003
#define DEFECT_MEDQ 0.014

#define T Dynprog_T
typedef struct T *T;

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
		    int cdna_direction, Pairpool_T pairpool, int extraband_single, double defect_rate);

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
		  int cdna_direction, Pairpool_T pairpool, int extraband_paired, double defect_rate,
		  bool forcep);

extern List_T
Dynprog_genome_gap (int *dynprogindex, int *finalscore, int *new_leftgenomepos, int *new_rightgenomepos,
		    int *nmatches, int *nmismatches, int *nopens, int *nindels,
		    int *exonhead, int *introntype, T dynprogL, T dynprogR, 
		    char *sequence1, char *sequenceuc1,
		    char *sequence2L, char *sequenceuc2L,
		    char *revsequence2R, char *revsequenceuc2R,
		    int length1, int length2L, int length2R, 
		    int offset1, int offset2L, int revoffset2R, 
#ifdef PMAP
		    char *queryaaseq,
#endif
		    int cdna_direction, Pairpool_T pairpool, int extraband_paired,
		    double defect_rate, int maxpeelback, bool halfp, bool finalp);

extern List_T
Dynprog_end5_gap (int *dynprogindex, int *finalscore, int *nmatches, int *nmismatches,
		  int *nopens, int *nindels, T dynprog,
		  char *revsequence1, char *revsequenceuc1,
		  char *revsequence2, char *revsequenceuc2,
		  int length1, int length2, int revoffset1, int revoffset2, 
#ifdef PMAP
		  char *queryaaseq,
#endif
		  int cdna_direction, Pairpool_T pairpool, int extraband_end, double defect_rate);

extern List_T
Dynprog_end3_gap (int *dynprogindex, int *finalscore, int *nmatches, int *nmismatches,
		  int *nopens, int *nindels, T dynprog,
		  char *sequence1, char *sequenceuc1,
		  char *sequence2, char *sequenceuc2,
		  int length1, int length2, int offset1, int offset2, 
#ifdef PMAP
		  char *queryaaseq,
#endif
		  int cdna_direction, Pairpool_T pairpool, int extraband_end, double defect_rate);

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
		    int cdna_direction, Pairpool_T pairpool, int extraband_end, double defect_rate);

extern List_T
Dynprog_microexon_int (int *dynprogindex, int *microintrontype,
		       char *sequence1, char *sequenceuc1,
		       char *sequence2L, char *sequenceuc2L,
		       char *revsequence2R, char *revsequenceuc2R,
		       int length1, int length2L, int length2R,
		       int offset1, int offset2L, int revoffset2R, int cdna_direction,
#ifdef PMAP
		       char *queryaaseq,
#endif
		       char *queryseq, char *queryuc, char *genomicseg, char *genomicuc,
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
