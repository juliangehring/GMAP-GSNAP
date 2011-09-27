/* $Id: dynprog.h,v 1.44 2005/10/01 05:35:08 twu Exp $ */
#ifndef DYNPROG_INCLUDED
#define DYNPROG_INCLUDED
#include "bool.h"
#include "list.h"
#include "pairpool.h"

/* Note: HIGHQ, MEDQ, and LOWQ indicate parameters for high, medium,
   and low sequence quality, respectively */
#define DEFECT_HIGHQ 0.003
#define DEFECT_MEDQ 0.014

/* Intron types.  Results of logical AND of dinucleotide pairs.  */
#define GTAG_FWD 0x20		/* 100000 GT-AG */
#define GCAG_FWD 0x10		/* 010000 GC-AG */
#define ATAC_FWD 0x08		/* 001000 AT-AC */
#define GTAG_REV 0x04		/* 000100 CT-AC */
#define GCAG_REV 0x02		/* 000010 CT-GC */
#define ATAC_REV 0x01		/* 000001 GT-AT */
#define NONINTRON 0x00


#define T Dynprog_T
typedef struct T *T;

extern T
Dynprog_new (int maxlookback, int extraquerygap, int maxpeelback,
	     int extramaterial_end, int extramaterial_paired);
extern void
Dynprog_free (T *old);

extern int
Dynprog_pairdistance (int c1, int c2);

extern void
Dynprog_term (void);
extern void
Dynprog_init (int maxlookback, int extraquerygap, int maxpeelback,
	      int extramaterial_end, int extramaterial_paired);

extern List_T
Dynprog_single_gap (int *finalscore, int *nmatches, int *nmismatches, int *nopens, int *nindels,
		    T dynprog, char *sequence1, char *sequenceuc1, char *sequence2, char *sequenceuc2,
		    int length1, int length2, int offset1, int offset2,
#ifdef PMAP
		    char *queryaaseq,
#endif
		    Pairpool_T pairpool, int extraband_single, double defect_rate);

extern List_T
Dynprog_cdna_gap (int *finalscore, T dynprogL, T dynprogR, 
		  char *sequence1L, char *sequenceuc1L,
		  char *revsequence1R, char *revsequenceuc1R,
		  char *sequence2, char *sequenceuc2,
		  int length1L, int length1R, int length2,
		  int offset1L, int revoffset1R, int offset2,
#ifdef PMAP
		  char *queryaaseq,
#endif
		  Pairpool_T pairpool, int extraband_paired, double defect_rate);

extern List_T
Dynprog_genome_gap (int *finalscore, int *nmatches, int *nmismatches, int *nopens, int *nindels,
		    int *exonhead, int *introntype, T dynprogL, T dynprogR, 
		    char *sequence1, char *sequenceuc1,
		    char *sequence2L, char *sequenceuc2L,
		    char *revsequence2R, char *revsequenceuc2R,
		    int length1, int length2L, int length2R, 
		    int offset1, int offset2L, int revoffset2R, 
#ifdef PMAP
		    char *queryaaseq,
#endif
		    int cdna_direction, int ngap, Pairpool_T pairpool, int extraband_paired,
		    bool endp, double defect_rate, bool returnpairsp, bool addgapp);

extern List_T
Dynprog_end5_gap (int *finalscore, int *nmatches, int *nmismatches,
		  int *nopens, int *nindels, T dynprog,
		  char *revsequence1, char *revsequenceuc1,
		  char *revsequence2, char *revsequenceuc2,
		  int length1, int length2, int revoffset1, int revoffset2, 
#ifdef PMAP
		  char *queryaaseq,
#endif
		  Pairpool_T pairpool, int extraband_end, double defect_rate,
		  int cdna_direction, int ngap, bool extend_mismatch_p);

extern List_T
Dynprog_end3_gap (int *finalscore, int *nmatches, int *nmismatches,
		  int *nopens, int *nindels, T dynprog,
		  char *sequence1, char *sequenceuc1,
		  char *sequence2, char *sequenceuc2,
		  int length1, int length2, int offset1, int offset2, 
#ifdef PMAP
		  char *queryaaseq,
#endif
		  Pairpool_T pairpool, int extraband_end, double defect_rate,
		  int cdna_direction, int ngap, bool extend_mismatch_p);

extern List_T
Dynprog_microexon_int (int *microintrontype,
		       char *sequence1, char *sequenceuc1,
		       char *sequence2L, char *sequenceuc2L,
		       char *revsequence2R, char *revsequenceuc2R,
		       int length1, int length2L, int length2R,
		       int offset1, int offset2L, int revoffset2R, int cdna_direction,
#ifdef PMAP
		       char *queryaaseq,
#endif
		       char *queryseq, char *queryuc, char *genomicseg, char *genomicuc,
		       Pairpool_T pairpool, int ngap);

extern List_T
Dynprog_microexon_5 (int *microintrontype, int *microexonlength,
		     char *revsequence1, char *revsequenceuc1,
		     char *revsequence2, char *revsequenceuc2,
		     int length1, int length2, int revoffset1, int revoffset2, int cdna_direction,
#ifdef PMAP
		     char *queryaaseq,
#endif
		     char *queryseq, char *queryuc, char *genomicseg, char *genomicuc,
		     Pairpool_T pairpool, int ngap, bool end_microexons_p);

extern List_T
Dynprog_microexon_3 (int *microintrontype, int *microexonlength,
		     char *sequence1, char *sequenceuc1,
		     char *sequence2, char *sequenceuc2,
		     int length1, int length2, int offset1, int offset2, int cdna_direction,
#ifdef PMAP
		     char *queryaaseq,
#endif
		     char *queryseq, char *queryuc, char *genomicseg, char *genomicuc,
		     int genomiclength, Pairpool_T pairpool, int ngap, bool end_microexons_p);

#undef T
#endif
