static char rcsid[] = "$Id: dynprog.c 112262 2013-10-22 22:07:51Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "dynprog.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>		/* For ceil, log, pow */
#include <ctype.h>		/* For tolower */
#ifdef HAVE_SSE2
#include <emmintrin.h>
#endif
#ifdef HAVE_SSE4_1
#include <smmintrin.h>
#endif


#include "bool.h"
#include "except.h"
#include "assert.h"
#include "mem.h"
#include "comp.h"
#include "pair.h"
#include "pairdef.h"
#include "listdef.h"
#include "boyer-moore.h"
#include "intron.h"
#include "complement.h"
#include "splicetrie.h"
#include "maxent.h"
#include "maxent_hr.h"
#include "fastlog.h"


/* Tests whether get_genomic_nt == genomicseg in compute_scores procedures */
/* #define EXTRACT_GENOMICSEG 1 */


/* Prints parameters and results */
#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/* Prints matrices */
#ifdef DEBUG2
#define debug2(x) x
#else
#define debug2(x)
#endif

/* Prints all winning bridge scores */
#ifdef DEBUG3
#define debug3(x) x
#else
#define debug3(x)
#endif

/* Prints all losing bridge scores */
#ifdef DEBUG3A
#define debug3a(x) x
#else
#define debug3a(x)
#endif

/* Codon instantiation */
#ifdef DEBUG4
#define debug4(x) x
#else
#define debug4(x)
#endif

/* Known splicing */
#ifdef DEBUG5
#define debug5(x) x
#else
#define debug5(x)
#endif

/* Ends */
#ifdef DEBUG6
#define debug6(x) x
#else
#define debug6(x)
#endif

/* Ends, known splicing.  May want to turn on DEBUG3 in splicetrie.c  */
#ifdef DEBUG7
#define debug7(x) x
#else
#define debug7(x)
#endif

/* Getting genomic nt */
#ifdef DEBUG8
#define debug8(x) x
#else
#define debug8(x)
#endif

/* Splice site probabilities */
#ifdef DEBUG9
#define debug9(x) x
#else
#define debug9(x)
#endif


/* Binary search */
#ifdef DEBUG10
#define debug10(x) x
#else
#define debug10(x)
#endif

/* traceback_nogaps */
#ifdef DEBUG11
#define debug11(x) x
#else
#define debug11(x)
#endif

/* Old matrix computations */
#ifdef DEBUG12
#define debug12(x) x
#else
#define debug12(x)
#endif

/* Old matrix computations, details */
#ifdef DEBUG12A
#define debug12a(x) x
#else
#define debug12a(x)
#endif

/* SIMD F loop */
#ifdef DEBUG13
#define debug13(x) x
#else
#define debug13(x)
#endif

/* Comparing SIMD with standard code */
#ifdef DEBUG14
#define debug14(x) x
#else
#define debug14(x)
#endif

/* print_vector */
#ifdef DEBUG15
#define debug15(x) x
#else
#define debug15(x)
#endif

/* Homopolymer (e.g., PacBio) */
#ifdef DEBUG16
#define debug16(x) x
#else
#define debug16(x)
#endif

/* Homopolymer details */
#ifdef DEBUG16A
#define debug16a(x) x
#else
#define debug16a(x)
#endif




#if defined(DEBUG2) || defined(DEBUG14)
#define NEG_INFINITY_DISPLAY -99
#endif

#define NEG_INFINITY_8 -128
#define MAX_CHAR 127

#define NEG_INFINITY_16 -32768
#define MAX_SHORT 32767


/* We can allow -128 and -32768 for NEG_INFINITY in SIMD procedures,
   because we are using saturation */
#if defined(HAVE_SSE4_1) || defined(HAVE_SSE2)
#define ONE_CHAR 1
#define LAST_CHAR 15
#define SIMD_NCHARS 16		/* 16 8-bit chars in 128 bits */

#define ONE_SHORT 2
#define LAST_SHORT 14		/* (8 - 1) * 2 */
#define SIMD_NSHORTS 8		/* 8 16-bit shorts in 128 bits */
#endif

/* Can allow -32768 in non-SIMD procedures, because we are using ints */
#define NEG_INFINITY_32 -32768


#define ONESIDEGAP 1

/*
#define RIGHTANGLE 1
*/

#define PROB_CEILING 0.85
#define PROB_FLOOR 0.75

#define MICROEXON_PVALUE_HIGHQ 0.01
#define MICROEXON_PVALUE_MEDQ 0.001
#define MICROEXON_PVALUE_LOWQ 0.0001
#define ENDSEQUENCE_PVALUE 0.001 /* Have stricter threshold for making end exons */

#define MIN_MICROEXON_LENGTH 3
#ifdef PMAP
#define MAX_MICROEXON_LENGTH 17	/* Should be oligomer length - 1 plus peelback */
#else
#define MAX_MICROEXON_LENGTH 12	/* Should be oligomer length - 1 plus peelback */
#endif
#define MICROINTRON_LENGTH 9
#define INSERT_PAIRS 9

#define SIMD_MAXLENGTH_EPI8 40  /* 128/3 */
#define FULLMATCH 3
#define HALFMATCH 1
#ifdef PMAP
#define AMBIGUOUS 0
#else
#define AMBIGUOUS -1
#endif

typedef enum {HIGHQ, MEDQ, LOWQ, ENDQ} Mismatchtype_T;
#define NMISMATCHTYPES 4

/* Mismatch penalty of 10 or less needed to find

   CCTGTA...CAG   
   |  >>>   >>>
   CAG

   Mismatch penalty of 6 or less needed to find

   TCTGTA...CAG   
      >>>   >>>
   CAG

*/

/* These values were set to -5, -4, -3, but this led to chopped ends
   in GMAP alignments, and failure to find chimeras */
#define MISMATCH_HIGHQ -3
#define MISMATCH_MEDQ -3
#define MISMATCH_LOWQ -3

/* Previously allowed lower mismatch scores on end to allow more
   complete alignments to the end, and because ends are typically of
   lower quality.  Previously made equal to FULLMATCH, because
   criterion is nmatches >= nmismatches.  However, extensions at ends
   appear to defeat purpose of trimming, so increase mismatch at end
   from -3 to -4. */
#define MISMATCH_ENDQ -5


/* Note: In definitions below, extensions don't include the first base */

/* Make PAIRED the same as SINGLE.  Need to avoid the problem where
   gap-match-gap is preferred over two mismatches, so

   OPEN+FULLMATCH+OPEN < MISMATCH+MISMATCH, or

   OPEN < (MISMATCH+MISMATCH-FULLMATCH)/2

   In middle:
            (mismatch+mismatch-fullmatch)/2   open
   HIGHQ     (-20-3)/2=-11.5                  -17
   MEDQ      (-18-3)/2=-10.5                  -15
   LOWQ      (-16-3)/2= -9.5                  -13

   At ends:
            (mismatch+mismatch-fullmatch)/2   open
   HIGHQ     (-14-3)/2= -8.5                  -15
   MEDQ      (-12-3)/2= -7.5                  -13
   LOWQ      (-10-3)/2= -6.5                  -11 
*/

typedef enum {PAIRED_HIGHQ, PAIRED_MEDQ, PAIRED_LOWQ, 
	      SINGLE_HIGHQ, SINGLE_MEDQ, SINGLE_LOWQ,
	      END_HIGHQ, END_MEDQ, END_LOWQ,
	      CDNA_HIGHQ, CDNA_MEDQ, CDNA_LOWQ} Jumptype_T;
#define NJUMPTYPES 12

#define PAIRED_OPEN_HIGHQ -12	/* was -18 */
#define PAIRED_OPEN_MEDQ -8
#define PAIRED_OPEN_LOWQ -4

#define PAIRED_EXTEND_HIGHQ -3	/* was -3 */
#define PAIRED_EXTEND_MEDQ -2
#define PAIRED_EXTEND_LOWQ -1

#define SINGLE_OPEN_HIGHQ -12	/* was -10 */
#define SINGLE_OPEN_MEDQ -8
#define SINGLE_OPEN_LOWQ -4

#define SINGLE_EXTEND_HIGHQ -3	/* was -3 */
#define SINGLE_EXTEND_MEDQ -2
#define SINGLE_EXTEND_LOWQ -1


/* cDNA insertions are biologically not meaningful, so look for a good
   gap opening somewhere */
#define CDNA_OPEN_HIGHQ -10
#define CDNA_OPEN_MEDQ -10
#define CDNA_OPEN_LOWQ -10

#define CDNA_EXTEND_HIGHQ -7
#define CDNA_EXTEND_MEDQ -7
#define CDNA_EXTEND_LOWQ -7

/* Ends tend to be of lower quality, so we don't want to introduce gaps.
   Also, we make then indifferent to the quality of the rest of the
   sequence. */
/* was -10 open and -3 extend */
#define END_OPEN_HIGHQ -10
#define END_OPEN_MEDQ -8
#define END_OPEN_LOWQ -6

#define END_EXTEND_HIGHQ -2
#define END_EXTEND_MEDQ -2
#define END_EXTEND_LOWQ -2


/* To reward one mismatch, but not two, should make

   FULLMATCH < INTRON+MISMATCH, and
   FULLMATCH+FULLMATCH > INTRON+MISMATCH+MISMATCH, or

   FULLMATCH-MISMATCH < INTRON < FULLMATCH+FULLMATCH-MISMATCH-MISMATCH

             1 mismatch    2 mismatches  3 mismatches  intron
   HIGHQ     3-(-10)=13 ** 6-(-20)=26    9-(-30)=39      22
   MEDQ      3-(-9)= 12    6-(-18)=24 ** 9-(-27)=36      25
   LOWQ      3-(-8)= 11    6-(-16)=22 ** 9-(-24)=33      28 */

/* To reward one gap, but not two, in preference to matching part of
   the dinucleotide, 

   FULLMATCH < INTRON+OPEN, and
   FULLMATCH < INTRON+OPEN+EXTEND, or

   FULLMATCH-OPEN < INTRON < FULLMATCH-OPEN-EXTEND

             1 gap         gap+extend    gap+2extend   intron
   HIGHQ     3-(-17)=20 ** 3-(-24)=26    3-(-31)=34      22
   MEDQ      3-(-15)=18    3-(-21)=24 ** 3-(-27)=30      25
   LOWQ      3-(-13)=16    3-(-18)=21    3-(-23)=26 **   28 */

/* Don't want to make too high, otherwise we will harm evaluation of
   dual introns vs. single intron */
#define CANONICAL_INTRON_HIGHQ 10 /* GT-AG */
#define CANONICAL_INTRON_MEDQ  16
#define CANONICAL_INTRON_LOWQ  22

#define FINAL_CANONICAL_INTRON_HIGHQ 30 /* GT-AG */
#define FINAL_CANONICAL_INTRON_MEDQ  36
#define FINAL_CANONICAL_INTRON_LOWQ  42

#define KNOWN_SPLICESITE_REWARD 20

/* Prefer alternate intron to other non-canonicals, but don't
   introduce mismatches or gaps to identify */
#define GCAG_INTRON 15
#define ATAC_INTRON 12
#define FINAL_GCAG_INTRON 20    /* Amount above regular should approximately
				   match FINAL_CANONICAL_INTRON - CANONICAL_INTRON */
#define FINAL_ATAC_INTRON 12

/* .01 = Prob(noncanonical) > Prob(sequence gap) = 0.003*(.20) */
#define MAXHORIZJUMP_HIGHQ 1
#define MAXVERTJUMP_HIGHQ 1

/* .01 = Prob(noncanonical) > Prob(sequence gap) = 0.014*(.20) */
#define MAXHORIZJUMP_MEDQ 1
#define MAXVERTJUMP_MEDQ 1

/* .01 = Prob(noncanonical) < Prob(sequence gap) */
#define MAXHORIZJUMP_LOWQ 1
#define MAXVERTJUMP_LOWQ 1


#if defined(HAVE_SSE4_1) || defined(HAVE_SSE2)
typedef char Score8_T;
typedef char Direction8_T;

typedef short Score16_T;
typedef short Direction16_T;
#endif

typedef short Pairdistance_T;

/* For standard dynamic programming.  Use ints, so NEG_INFINITY_32 works. */
typedef int Score32_T;
typedef int Direction32_T;

/* Genome is on the horizontal axis.  Query sequence is on the vertical axis.  Dynamic programming fills in matrices column by column */
/* The following values are for directions_nogap.  For directions_Egap, the choices are DIAG or not DIAG (meaning HORIZ). */
/*  For directions_Fgap, the choices are DIAG or not DIAG (meaning VERT) */
#define VERT -2			/* or VERT == -3 in SIMD code.  Don't check for dir == VERT.  Check instead if dir == DIAG or dir == HORIZ */
#define HORIZ -1
#define DIAG 0			/* Pre-dominant case.  Directions_alloc clears to this value. */


static IIT_T splicing_iit;
static int *splicing_divint_crosstable;
static int donor_typeint;
static int acceptor_typeint;

static Univcoord_T *splicesites;
static Splicetype_T *splicetypes;
static Chrpos_T *splicedists;
static int nsplicesites;
static Trieoffset_T *trieoffsets_obs;
static Triecontent_T *triecontents_obs;
static Trieoffset_T *trieoffsets_max;
static Triecontent_T *triecontents_max;

static bool novelsplicingp;
static bool homopolymerp;


char *
Dynprog_endalign_string (Endalign_T endalign) {
  switch (endalign) {
  case QUERYEND_GAP: return "queryend_gap";
  case QUERYEND_INDELS: return "queryend_indels";
  case QUERYEND_NOGAPS: return "queryend_nogaps";
  case BEST_LOCAL: return "best_local";
  default: 
    printf("endalign %d not recognized\n",endalign);
    return "";
  }
}



void
Dynprog_setup (bool novelsplicingp_in,
	       IIT_T splicing_iit_in, int *splicing_divint_crosstable_in,
	       int donor_typeint_in, int acceptor_typeint_in,
	       Univcoord_T *splicesites_in, Splicetype_T *splicetypes_in,
	       Chrpos_T *splicedists_in, int nsplicesites_in,
	       Trieoffset_T *trieoffsets_obs_in, Triecontent_T *triecontents_obs_in,
	       Trieoffset_T *trieoffsets_max_in, Triecontent_T *triecontents_max_in,
	       bool homopolymerp_in) {
  novelsplicingp = novelsplicingp_in;

  splicing_iit = splicing_iit_in;
  splicing_divint_crosstable = splicing_divint_crosstable_in;
  donor_typeint = donor_typeint_in;
  acceptor_typeint = acceptor_typeint_in;

  splicesites = splicesites_in;
  splicetypes = splicetypes_in;
  splicedists = splicedists_in;
  nsplicesites = nsplicesites_in;
  trieoffsets_obs = trieoffsets_obs_in;
  triecontents_obs = triecontents_obs_in;
  trieoffsets_max = trieoffsets_max_in;
  triecontents_max = triecontents_max_in;

  homopolymerp = homopolymerp_in;

  return;
}


int
Dynprog_score (int matches, int mismatches, int qopens, int qindels, int topens, int tindels,
	       double defect_rate) {

  if (defect_rate < DEFECT_HIGHQ) {
    return FULLMATCH*matches + MISMATCH_HIGHQ*mismatches + SINGLE_OPEN_HIGHQ*qopens + SINGLE_EXTEND_HIGHQ*qindels
      + SINGLE_OPEN_HIGHQ*topens + SINGLE_EXTEND_HIGHQ*tindels;
  } else if (defect_rate < DEFECT_MEDQ) {
    return FULLMATCH*matches + MISMATCH_MEDQ*mismatches + SINGLE_OPEN_MEDQ*qopens + SINGLE_EXTEND_MEDQ*qindels
      + SINGLE_OPEN_MEDQ*topens + SINGLE_EXTEND_MEDQ*tindels;
  } else {
    return FULLMATCH*matches + MISMATCH_LOWQ*mismatches + SINGLE_OPEN_LOWQ*qopens + SINGLE_EXTEND_LOWQ*qindels
      + SINGLE_OPEN_LOWQ*topens + SINGLE_EXTEND_LOWQ*tindels;
  }
}


/************************************************************************
 * get_genomic_nt
 ************************************************************************/

static char complCode[128] = COMPLEMENT_LC;

static char
get_genomic_nt (char *g_alt, int genomicpos, Univcoord_T chroffset, Univcoord_T chrhigh,
		bool watsonp) {
  char c2, c2_alt;
  Univcoord_T pos;

#if 0
  /* If the read has a deletion, then we will extend beyond 0 or genomiclength, so do not restrict. */
  if (genomicpos < 0) {
    return '*';

  } else if (genomicpos >= genomiclength) {
    return '*';

  }
#endif

  if (watsonp) {
    if ((pos = chroffset + genomicpos) < chroffset) { /* Must be <, and not <=, or dynamic programming will fail */
      *g_alt = '*';
      return '*';

    } else if (pos >= chrhigh) {
      *g_alt = '*';
      return '*';

#if 0
    } else if (genome) {
      /* Not necessary, because Genome_get_char_blocks should work */
      debug8(printf("At %u, genomicnt is %c\n",
		    genomicpos,Genome_get_char(genome,pos)));
      return Genome_get_char(genome,pos);
#endif

    } else {
      /* GMAP with user-supplied genomic segment */
      debug8(printf("At %u, genomicnt is %c\n",
		    genomicpos,Genome_get_char_blocks(pos)));
      return Genome_get_char_blocks(&(*g_alt),pos);
    }

  } else {
    if ((pos = chrhigh - genomicpos) < chroffset) { /* Must be <, and not <=, or dynamic programming will fail */
      *g_alt = '*';
      return '*';

    } else if (pos >= chrhigh) {
      *g_alt = '*';
      return '*';

#if 0
    } else if (genome) {
      /* Not necessary, because Genome_get_char_blocks should work */
      c2 = Genome_get_char(genome,pos);
#endif

    } else {
      /* GMAP with user-supplied genomic segment */
      c2 = Genome_get_char_blocks(&c2_alt,pos);
    }
    debug8(printf("At %u, genomicnt is %c\n",genomicpos,complCode[(int) c2]));
    *g_alt = complCode[(int) c2_alt];
    return complCode[(int) c2];
  }
}


/************************************************************************
 * Matrix
 ************************************************************************/

#if !defined(HAVE_SSE2) || defined(DEBUG14)
/* Makes a matrix of dimensions 0..glength x 0..rlength inclusive */
static Score32_T **
Matrix32_alloc (int rlength, int glength, Score32_T **ptrs, Score32_T *space) {
  Score32_T **matrix;
  int i;

  if (glength <= 0 || rlength <= 0) {
    fprintf(stderr,"dynprog: lengths are negative: %d %d\n",rlength,glength);
    abort();
  }

  matrix = ptrs;
  matrix[0] = space;
  for (i = 1; i <= glength; i++) {
    matrix[i] = &(matrix[i-1][rlength + 1]);
  }

  memset((void *) space,0,(glength+1)*(rlength+1)*sizeof(Score32_T));

  return matrix;
}

#endif

#if defined(DEBUG2) || defined(DEBUG14)
static void
Matrix8_print (Score8_T **matrix, int rlength, int glength, char *rsequence,
	       char *gsequence, char *gsequencealt,
	       int goffset, Univcoord_T chroffset, Univcoord_T chrhigh,
	       bool watsonp, bool revp) {
  int i, j;
  char g_alt;

  _mm_lfence();

  if (gsequence) {
    printf("  ");
    for (j = 0; j <= glength; ++j) {
      if (j == 0) {
	printf("    ");
      } else {
	printf("  %c ",revp ? gsequence[-j+1] : gsequence[j-1]);
      }
    }
    printf("\n");
  }

  printf("  ");
  for (j = 0; j <= glength; ++j) {
    if (j == 0) {
      printf("    ");
    } else {
      if (revp == false) {
	printf("  %c ",get_genomic_nt(&g_alt,goffset+j-1,chroffset,chrhigh,watsonp));
      } else {
	printf("  %c ",get_genomic_nt(&g_alt,goffset+1-j,chroffset,chrhigh,watsonp));
      }
    }
  }
  printf("\n");

  for (i = 0; i <= rlength; ++i) {
    if (i == 0) {
      printf("  ");
    } else {
      printf("%c ",revp ? rsequence[-i+1] : rsequence[i-1]);
    }
    for (j = 0; j <= glength; ++j) {
      if (matrix[j][i] < NEG_INFINITY_DISPLAY) {
	printf("%3d ",NEG_INFINITY_DISPLAY);
      } else {
	printf("%3d ",matrix[j][i]);
      }
    }
    printf("\n");
  }
  printf("\n");

  return;
}

static void
Matrix16_print (Score16_T **matrix, int rlength, int glength, char *rsequence,
		char *gsequence, char *gsequencealt,
		int goffset, Univcoord_T chroffset, Univcoord_T chrhigh,
		bool watsonp, bool revp) {
  int i, j;
  char g_alt;

  _mm_lfence();

  if (gsequence) {
    printf("  ");
    for (j = 0; j <= glength; ++j) {
      if (j == 0) {
	printf("    ");
      } else {
	printf("  %c ",revp ? gsequence[-j+1] : gsequence[j-1]);
      }
    }
    printf("\n");
  }

  printf("  ");
  for (j = 0; j <= glength; ++j) {
    if (j == 0) {
      printf("    ");
    } else {
      if (revp == false) {
	printf("  %c ",get_genomic_nt(&g_alt,goffset+j-1,chroffset,chrhigh,watsonp));
      } else {
	printf("  %c ",get_genomic_nt(&g_alt,goffset+1-j,chroffset,chrhigh,watsonp));
      }
    }
  }
  printf("\n");

  for (i = 0; i <= rlength; ++i) {
    if (i == 0) {
      printf("  ");
    } else {
      printf("%c ",revp ? rsequence[-i+1] : rsequence[i-1]);
    }
    for (j = 0; j <= glength; ++j) {
      if (matrix[j][i] < NEG_INFINITY_DISPLAY) {
	printf("%3d ",NEG_INFINITY_DISPLAY);
      } else {
	printf("%3d ",matrix[j][i]);
      }
    }
    printf("\n");
  }
  printf("\n");

  return;
}

static void
Matrix32_print (Score32_T **matrix, int rlength, int glength, char *rsequence,
		char *gsequence, char *gsequencealt,
		int goffset, Univcoord_T chroffset, Univcoord_T chrhigh,
		bool watsonp, bool revp) {
  int i, j;
  char g_alt;

  if (gsequence) {
    printf("  ");
    for (j = 0; j <= glength; ++j) {
      if (j == 0) {
	printf("    ");
      } else {
	printf("  %c ",revp ? gsequence[-j+1] : gsequence[j-1]);
      }
    }
    printf("\n");
  }

  printf("  ");
  for (j = 0; j <= glength; ++j) {
    if (j == 0) {
      printf("    ");
    } else {
      if (revp == false) {
	printf("  %c ",get_genomic_nt(&g_alt,goffset+j-1,chroffset,chrhigh,watsonp));
      } else {
	printf("  %c ",get_genomic_nt(&g_alt,goffset+1-j,chroffset,chrhigh,watsonp));
      }
    }
  }
  printf("\n");

  for (i = 0; i <= rlength; ++i) {
    if (i == 0) {
      printf("  ");
    } else {
      printf("%c ",revp ? rsequence[-i+1] : rsequence[i-1]);
    }
    for (j = 0; j <= glength; ++j) {
      if (matrix[j][i] < NEG_INFINITY_DISPLAY) {
	printf("%3d ",NEG_INFINITY_DISPLAY);
      } else {
	printf("%3d ",matrix[j][i]);
      }
    }
    printf("\n");
  }
  printf("\n");

  return;
}
#endif



#ifdef DEBUG12
struct Int3_T {
  Score16_T Egap;
  Score16_T Fgap;
  Score16_T nogap;
};
#endif


#ifdef DEBUG12
/* Makes a matrix of dimensions 0..glength x 0..rlength inclusive */
static struct Int3_T **
Matrix3_alloc (int rlength, int glength, struct Int3_T **ptrs, struct Int3_T *space) {
  struct Int3_T **matrix;
  int i;

  if (glength <= 0 || rlength <= 0) {
    fprintf(stderr,"dynprog: lengths are negative: %d %d\n",rlength,glength);
    abort();
  }

  matrix = ptrs;
  matrix[0] = space;
  for (i = 1; i <= glength; i++) {
    matrix[i] = &(matrix[i-1][rlength + 1]);
  }

  memset((void *) space,0,(glength+1)*(rlength+1)*sizeof(struct Int3_T));

  return matrix;
}
#endif


#ifdef DEBUG12A
static void
Matrix3_print (struct Int3_T **matrix, int rlength, int glength, char *rsequence,
	       char *gsequence, char *gsequencealt,
	       int goffset, Univcoord_T chroffset, Univcoord_T chrhigh,
	       bool watsonp, bool revp) {
  int i, j;
  char g_alt;

  printf("G1");
  if (gsequence) {
    for (j = 0; j <= glength; ++j) {
      if (j == 0) {
	printf("    ");
      } else {
	printf("  %c ",revp ? gsequence[-j+1] : gsequence[j-1]);
      }
    }
    printf("\n");
  }

  for (j = 0; j <= glength; ++j) {
    if (j == 0) {
      printf("    ");
    } else {
      if (revp == false) {
	printf("  %c ",get_genomic_nt(&g_alt,goffset+j-1,chroffset,chrhigh,watsonp));
      } else {
	printf("  %c ",get_genomic_nt(&g_alt,goffset+1-j,chroffset,chrhigh,watsonp));
      }
    }
  }
  printf("\n");

  for (i = 0; i <= rlength; ++i) {
    if (i == 0) {
      printf("  ");
    } else {
      printf("%c ",revp ? rsequence[-i+1] : rsequence[i-1]);
    }
    for (j = 0; j <= glength; ++j) {
      if (matrix[j][i].Egap < NEG_INFINITY) {
	printf("%3d ",NEG_INFINITY);
      } else {
	printf("%3d ",matrix[j][i].Egap);
      }
    }
    printf("\n");
  }
  printf("\n");


  printf("NG");
  if (gsequence) {
    for (j = 0; j <= glength; ++j) {
      if (j == 0) {
	printf("    ");
      } else {
	printf("  %c ",revp ? gsequence[-j+1] : gsequence[j-1]);
      }
    }
    printf("\n");
  }

  for (j = 0; j <= glength; ++j) {
    if (j == 0) {
      printf("    ");
    } else {
      if (revp == false) {
	printf("  %c ",get_genomic_nt(&g_alt,goffset+j-1,chroffset,chrhigh,watsonp));
      } else {
	printf("  %c ",get_genomic_nt(&g_alt,goffset+1-j,chroffset,chrhigh,watsonp));
      }
    }
  }
  printf("\n");

  for (i = 0; i <= rlength; ++i) {
    if (i == 0) {
      printf("  ");
    } else {
      printf("%c ",revp ? rsequence[-i+1] : rsequence[i-1]);
    }
    for (j = 0; j <= glength; ++j) {
      if (matrix[j][i].nogap < NEG_INFINITY) {
	printf("%3d ",NEG_INFINITY);
      } else {
	printf("%3d ",matrix[j][i].nogap);
      }
    }
    printf("\n");
  }
  printf("\n");


  printf("G2");
  if (gsequence) {
    for (j = 0; j <= glength; ++j) {
      if (j == 0) {
	printf("    ");
      } else {
	printf("  %c ",revp ? gsequence[-j+1] : gsequence[j-1]);
      }
    }
    printf("\n");
  }

  for (j = 0; j <= glength; ++j) {
    if (j == 0) {
      printf("    ");
    } else {
      if (revp == false) {
	printf("  %c ",get_genomic_nt(&g_alt,goffset+j-1,chroffset,chrhigh,watsonp));
      } else {
	printf("  %c ",get_genomic_nt(&g_alt,goffset+1-j,chroffset,chrhigh,watsonp));
      }
    }
  }
  printf("\n");

  for (i = 0; i <= rlength; ++i) {
    if (i == 0) {
      printf("  ");
    } else {
      printf("%c ",revp ? rsequence[-i+1] : rsequence[i-1]);
    }
    for (j = 0; j <= glength; ++j) {
      if (matrix[j][i].Fgap < NEG_INFINITY) {
	printf("%3d ",NEG_INFINITY);
      } else {
	printf("%3d ",matrix[j][i].Fgap);
      }
    }
    printf("\n");
  }
  printf("\n");

  return;
}
#endif


/************************************************************************/
/*  Directions  */
/************************************************************************/

#if !defined(HAVE_SSE2) || defined(DEBUG14)
/* Makes a matrix of dimensions 0..glength x 0..rlength inclusive */
static Direction32_T **
Directions32_alloc (int rlength, int glength, Direction32_T **ptrs, Direction32_T *space) {
  Direction32_T **directions;
  int i;

  directions = ptrs;
  directions[0] = space;
  for (i = 1; i <= glength; i++) {
    directions[i] = &(directions[i-1][rlength + 1]);
  }

  memset((void *) space,/*DIAG*/0,(glength+1)*(rlength+1)*sizeof(Direction32_T));

  return directions;
}
#endif

#ifdef DEBUG2
static void
Directions8_print (Direction8_T **directions_nogap, Direction8_T **directions_Egap, Direction8_T **directions_Fgap,
		   int rlength, int glength, char *rsequence, char *gsequence, char *gsequence_alt,
		   int goffset, Univcoord_T chroffset, Univcoord_T chrhigh,
		   bool watsonp, bool revp) {
  int i, j;
  char g_alt;

  _mm_lfence();

  if (gsequence) {
    printf("  ");
    for (j = 0; j <= glength; ++j) {
      if (j == 0) {
	printf("      ");
      } else {
	printf("  %c   ",revp ? gsequence[-j+1] : gsequence[j-1]);
      }
    }
    printf("\n");
  }

  printf("  ");
  for (j = 0; j <= glength; ++j) {
    if (j == 0) {
      printf("      ");
    } else {
      if (revp == false) {
	printf("  %c   ",get_genomic_nt(&g_alt,goffset+j-1,chroffset,chrhigh,watsonp));
      } else {
	printf("  %c   ",get_genomic_nt(&g_alt,goffset+1-j,chroffset,chrhigh,watsonp));
      }
    }
  }
  printf("\n");

  for (i = 0; i <= rlength; ++i) {
    if (i == 0) {
      printf("  ");
    } else {
      printf("%c ",revp ? rsequence[-i+1] : rsequence[i-1]);
    }
    for (j = 0; j <= glength; ++j) {
      if (directions_Egap[j][i] == DIAG) {
	printf("D");
      } else {
	/* Must be HORIZ */
	printf("H");
      }
      printf("|");
      if (directions_nogap[j][i] == DIAG) {
	printf("D");
      } else if (directions_nogap[j][i] == HORIZ) {
	printf("H");
      } else {
	/* Must be VERT */
	printf("V");
      }
      printf("|");
      if (directions_Fgap[j][i] == DIAG) {
	printf("D");
      } else {
	/* Must be VERT */
	printf("V");
      }
      printf(" ");
    }
    printf("\n");
  }
  printf("\n");
  return;
}


static void
Directions16_print (Direction16_T **directions_nogap, Direction16_T **directions_Egap, Direction16_T **directions_Fgap,
		    int rlength, int glength, char *rsequence, char *gsequence, char *gsequence_alt,
		    int goffset, Univcoord_T chroffset, Univcoord_T chrhigh,
		    bool watsonp, bool revp) {
  int i, j;
  char g_alt;

  _mm_lfence();

  if (gsequence) {
    printf("  ");
    for (j = 0; j <= glength; ++j) {
      if (j == 0) {
	printf("      ");
      } else {
	printf("  %c   ",revp ? gsequence[-j+1] : gsequence[j-1]);
      }
    }
    printf("\n");
  }

  printf("  ");
  for (j = 0; j <= glength; ++j) {
    if (j == 0) {
      printf("      ");
    } else {
      if (revp == false) {
	printf("  %c   ",get_genomic_nt(&g_alt,goffset+j-1,chroffset,chrhigh,watsonp));
      } else {
	printf("  %c   ",get_genomic_nt(&g_alt,goffset+1-j,chroffset,chrhigh,watsonp));
      }
    }
  }
  printf("\n");

  for (i = 0; i <= rlength; ++i) {
    if (i == 0) {
      printf("  ");
    } else {
      printf("%c ",revp ? rsequence[-i+1] : rsequence[i-1]);
    }
    for (j = 0; j <= glength; ++j) {
      if (directions_Egap[j][i] == DIAG) {
	printf("D");
      } else {
	/* Must be HORIZ */
	printf("H");
      }
      printf("|");
      if (directions_nogap[j][i] == DIAG) {
	printf("D");
      } else if (directions_nogap[j][i] == HORIZ) {
	printf("H");
      } else {
	/* Must be VERT */
	printf("V");
      }
      printf("|");
      if (directions_Fgap[j][i] == DIAG) {
	printf("D");
      } else {
	/* Must be VERT */
	printf("V");
      }
      printf(" ");
    }
    printf("\n");
  }
  printf("\n");
  return;
}


static void
Directions32_print (Direction32_T **directions_nogap, Direction32_T **directions_Egap, Direction32_T **directions_Fgap,
		    int rlength, int glength, char *rsequence, char *gsequence, char *gsequence_alt,
		    int goffset, Univcoord_T chroffset, Univcoord_T chrhigh,
		    bool watsonp, bool revp) {
  int i, j;
  char g_alt;

  if (gsequence) {
    printf("  ");
    for (j = 0; j <= glength; ++j) {
      if (j == 0) {
	printf("      ");
      } else {
	printf("  %c   ",revp ? gsequence[-j+1] : gsequence[j-1]);
      }
    }
    printf("\n");
  }

  printf("  ");
  for (j = 0; j <= glength; ++j) {
    if (j == 0) {
      printf("      ");
    } else {
      if (revp == false) {
	printf("  %c   ",get_genomic_nt(&g_alt,goffset+j-1,chroffset,chrhigh,watsonp));
      } else {
	printf("  %c   ",get_genomic_nt(&g_alt,goffset+1-j,chroffset,chrhigh,watsonp));
      }
    }
  }
  printf("\n");

  for (i = 0; i <= rlength; ++i) {
    if (i == 0) {
      printf("  ");
    } else {
      printf("%c ",revp ? rsequence[-i+1] : rsequence[i-1]);
    }
    for (j = 0; j <= glength; ++j) {
      if (directions_Egap[j][i] == DIAG) {
	printf("D");
      } else {
	/* Must be HORIZ */
	printf("H");
      }
      printf("|");
      if (directions_nogap[j][i] == DIAG) {
	printf("D");
      } else if (directions_nogap[j][i] == HORIZ) {
	printf("H");
      } else {
	/* Must be VERT */
	printf("V");
      }
      printf("|");
      if (directions_Fgap[j][i] == DIAG) {
	printf("D");
      } else {
	/* Must be VERT */
	printf("V");
      }
      printf(" ");
    }
    printf("\n");
  }
  printf("\n");
  return;
}
#endif




#define QUERY_MAXLENGTH 500
#define GENOMIC_MAXLENGTH 2000


#define T Dynprog_T
struct T {
  int max_rlength;
  int max_glength;

#ifdef DEBUG12
  struct Int3_T **matrix3_ptrs, *matrix3_space;
#endif

#if !defined(HAVE_SSE2) || defined(DEBUG14)
  Score32_T **matrix_ptrs, *matrix_space;
  Direction32_T **directions_ptrs_0, *directions_space_0;
  Direction32_T **directions_ptrs_1, *directions_space_1;
  Direction32_T **directions_ptrs_2, *directions_space_2;
#endif
#ifdef HAVE_SSE2
  void **aligned_matrix_ptrs, *aligned_matrix_space;
  void **aligned_directions_ptrs_0, *aligned_directions_space_0;
  void **aligned_directions_ptrs_1, *aligned_directions_space_1;
  void **aligned_directions_ptrs_2, *aligned_directions_space_2;
#endif
};

static void
compute_maxlengths (int *max_rlength, int *max_glength,
		    int maxlookback, int extraquerygap, int maxpeelback,
		    int extramaterial_end, int extramaterial_paired) {
  *max_rlength = maxlookback + maxpeelback;
  if (*max_rlength < QUERY_MAXLENGTH) {
    *max_rlength = QUERY_MAXLENGTH;
  }

  *max_glength = *max_rlength + extraquerygap;
  if (extramaterial_end > extramaterial_paired) {
    *max_glength += extramaterial_end;
  } else {
    *max_glength += extramaterial_paired;
  }

  if (*max_glength < GENOMIC_MAXLENGTH) {
    *max_glength = GENOMIC_MAXLENGTH;
  }

  return;
}


T
Dynprog_new (int maxlookback, int extraquerygap, int maxpeelback,
	     int extramaterial_end, int extramaterial_paired) {
  T new = (T) MALLOC(sizeof(*new));
  int max_rlength, max_glength;

  compute_maxlengths(&max_rlength,&max_glength,
		     maxlookback,extraquerygap,maxpeelback,
		     extramaterial_end,extramaterial_paired);
  new->max_rlength = max_rlength;
  new->max_glength = max_glength;

#ifdef DEBUG12
  new->matrix3_ptrs = (struct Int3_T **) CALLOC(max_glength+1,sizeof(struct Int3_T *));
  new->matrix3_space = (struct Int3_T *) CALLOC((max_glength+1)*(max_rlength+1),sizeof(struct Int3_T));
#endif
#if !defined(HAVE_SSE2) || defined(DEBUG14)
  new->matrix_ptrs = (Score32_T **) CALLOC(max_glength+1,sizeof(Score32_T *));
  new->matrix_space = (Score32_T *) CALLOC((max_glength+1)*(max_glength+1),sizeof(Score32_T));
  new->directions_ptrs_0 = (Direction32_T **) CALLOC(max_glength+1,sizeof(Direction32_T *));
  new->directions_space_0 = (Direction32_T *) CALLOC((max_glength+1)*(max_rlength+1),sizeof(Direction32_T));
  new->directions_ptrs_1 = (Direction32_T **) CALLOC(max_glength+1,sizeof(Direction32_T *));
  new->directions_space_1 = (Direction32_T *) CALLOC((max_glength+1)*(max_rlength+1),sizeof(Direction32_T));
  new->directions_ptrs_2 = (Direction32_T **) CALLOC(max_glength+1,sizeof(Direction32_T *));
  new->directions_space_2 = (Direction32_T *) CALLOC((max_glength+1)*(max_rlength+1),sizeof(Direction32_T));
#endif
#if defined(HAVE_SSE4_1) || defined(HAVE_SSE2)
  /* Use SIMD_NCHARS > SIMD_NSHORTS and sizeof(Score16_T) > sizeof(Score8_T) */
  new->aligned_matrix_ptrs = (void **) CALLOC(max_glength+1,sizeof(void *));
  new->aligned_matrix_space = (void *) _mm_malloc((max_glength+1)*(max_rlength+SIMD_NCHARS+SIMD_NCHARS)*sizeof(Score16_T),16);
  new->aligned_directions_ptrs_0 = (void **) CALLOC(max_glength+1,sizeof(void *));
  new->aligned_directions_space_0 = (void *) _mm_malloc((max_glength+1)*(max_rlength+SIMD_NCHARS+SIMD_NCHARS)*sizeof(Score16_T),16);
  new->aligned_directions_ptrs_1 = (void **) CALLOC(max_glength+1,sizeof(void *));
  new->aligned_directions_space_1 = (void *) _mm_malloc((max_glength+1)*(max_rlength+SIMD_NCHARS+SIMD_NCHARS)*sizeof(Score16_T),16);
  new->aligned_directions_ptrs_2 = (void **) CALLOC(max_glength+1,sizeof(void *));
  new->aligned_directions_space_2 = (void *) _mm_malloc((max_glength+1)*(max_rlength+SIMD_NCHARS+SIMD_NCHARS)*sizeof(Score16_T),16);
#endif
  return new;
}


void
Dynprog_free (T *old) {
  if (*old) {
#ifdef DEBUG12
    FREE((*old)->matrix3_ptrs);
    FREE((*old)->matrix3_space);
#endif
#if !defined(HAVE_SSE2) || defined(DEBUG14)
    FREE((*old)->matrix_ptrs);
    FREE((*old)->matrix_space);
    FREE((*old)->directions_ptrs_2);
    FREE((*old)->directions_space_2);
    FREE((*old)->directions_ptrs_1);
    FREE((*old)->directions_space_1);
    FREE((*old)->directions_ptrs_0);
    FREE((*old)->directions_space_0);
#endif
#if defined(HAVE_SSE4_1) || defined(HAVE_SSE2)
    FREE((*old)->aligned_matrix_ptrs);
    _mm_free((*old)->aligned_matrix_space);
    FREE((*old)->aligned_directions_ptrs_2);
    _mm_free((*old)->aligned_directions_space_2);
    FREE((*old)->aligned_directions_ptrs_1);
    _mm_free((*old)->aligned_directions_space_1);
    FREE((*old)->aligned_directions_ptrs_0);
    _mm_free((*old)->aligned_directions_space_0);
#endif

    FREE(*old);
  }
  return;
}

/************************************************************************/

#ifdef PMAP

static int aa_index_table[128] =
  { -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,

  /*     *                                  */
    -1, 20, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1,

  /* A,  B,  C,  D,  E,  F,  G,  H,  I,  J, */
     0, -1,  1,  2,  3,  4,  5,  6,  7, -1,

  /* K,  L,  M,  N,  O,  P,  Q,  R,  S,  T  */
     8,  9, 10, 11, -1, 12, 13, 14, 15, 16,

  /* U,  V,  W,  X,  Y,  Z  */
    21, 17, 18, -1, 19, -1,

    -1, -1, -1, -1, -1, -1,

  /* a,  b,  c,  d,  e,  f,  g,  h,  i,  j, */
     0, -1,  1,  2,  3,  4,  5,  6,  7, -1,

  /* k,  l,  m,  n,  o,  p,  q,  r,  s,  t  */
     8,  9, 10, 11, -1, 12, 13, 14, 15, 16,

  /* u,  v,  w,  x,  y,  z  */
    21, 17, 18, -1, 19, -1,

    -1, -1, -1, -1, -1};


static char *iupac_table[21+1] =
  {"GCN",			/* A */
   "TGY",			/* C */
   "GAY",			/* D */
   "GAR",			/* E */
   "TTY",			/* F */
   "GGN",			/* G */
   "CAY",			/* H */
   "ATH",			/* I */
   "AAR",			/* K */
   "YTN",			/* L */
   "ATG",			/* M */
   "AAY",			/* N */
   "CCN",			/* P */
   "CAR",			/* Q */
   "MGN",			/* R */
   "WSN",			/* S */
   "ACN",			/* T */
   "GTN",			/* V */
   "TGG",			/* W */
   "TAY",			/* Y */
   "TRR",			/* STOP */
   "TGA"};			/* U */


static char aa_table[21+1] = "ACDEFGHIKLMNPQRSTVWY*U";


char
Dynprog_codon_char (char aa, int codonpos) {
  int index;
  char *codon;

  if ((index = aa_index_table[(int) aa]) < 0) {
    return 'N';
  } else {
    codon = iupac_table[index];
    return codon[codonpos];
  }
}



#ifdef PMAP
/* Same as in boyer-moore.c */
/* Handle only cases in iupac table in dynprog.c */
static bool matchtable[26][26] = 
/*  A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y,Z */
  {{1,0,0,0,0,0,0,1,0,0,0,0,1,1,0,0,0,1,0,0,0,0,1,0,0,0}, /* A */
   {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, /* B */
   {0,0,1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,0,1,0}, /* C */
   {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, /* D */
   {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, /* E */
   {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, /* F */
   {0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,1,1,0,0,0,0,0,0,0}, /* G */
   {1,0,1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,1,1,1,0,0,1,0,1,0}, /* H = [ACT] */
   {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, /* I */
   {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, /* J */
   {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, /* K */
   {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, /* L */
   {1,0,1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,1,1,0,0,0,1,0,1,0}, /* M = [AC] */
   {1,0,1,0,0,0,1,1,0,0,0,0,1,1,0,0,0,1,1,1,0,0,1,0,1,0}, /* N = [ACGT] */
   {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, /* O */
   {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, /* P */
   {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, /* Q */
   {1,0,0,0,0,0,1,1,0,0,0,0,1,1,0,0,0,1,1,0,0,0,1,0,0,0}, /* R = [AG] */
   {0,0,1,0,0,0,1,1,0,0,0,0,1,1,0,0,0,1,1,0,0,0,0,0,1,0}, /* S = [CG] */
   {0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,1,0,1,0}, /* T */
   {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, /* U */
   {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, /* V */
   {1,0,0,0,0,0,0,1,0,0,0,0,1,1,0,0,0,1,0,1,0,0,1,0,1,0}, /* W = [AT] */
   {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, /* X */
   {0,0,1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,1,1,0,0,1,0,1,0}, /* Y = [CT] */
   {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}}; /* Z */
#endif


/* both newntsequence and oldntsequence are offset, but aasequence starts from position 0 */
static char *
instantiate_codons (char *oldntsequence, char *aasequence, int offset, int ntlength) {
  char *newntsequence;
  int aapos, ntpos;
  int index, frame;
  char *codon, aa;

  newntsequence = (char *) CALLOC(ntlength+1,sizeof(char));
  strncpy(newntsequence,oldntsequence,ntlength);

  for (ntpos = 0; ntpos < ntlength; ntpos++) {
    if (oldntsequence[ntpos] == BACKTRANSLATE_CHAR) {
      aapos = (offset+ntpos)/3;
	
      aa = aasequence[aapos];
      index = aa_index_table[(int) aa];
      if (index < 0) {
	newntsequence[ntpos] = 'N';
      } else {
	codon = iupac_table[index];

	frame = (offset+ntpos) % 3;
	newntsequence[ntpos] = codon[frame];
      }
    }
  }
  debug4(printf("New sequence: %s\n",newntsequence));
  return newntsequence;
}
#endif


/************************************************************************/


/* This is still needed, because sequences passed into compute_scores
   might be lower-case */
#define PREUC 1			

static Pairdistance_T **pairdistance_array[NMISMATCHTYPES];
#ifndef HAVE_SSE4_1
static Pairdistance_T **pairdistance_array_plus_128[NMISMATCHTYPES];
#endif

static bool **consistent_array;

int
Dynprog_pairdistance (int c1, int c2) {
  return pairdistance_array[HIGHQ][c1][c2];
}

static void
permute_cases (int NA1, int NA2, Pairdistance_T score) {
  int i;

#ifdef PREUC
  int na1, na2;

  na1 = tolower(NA1);
  na2 = tolower(NA2);
#endif

#ifdef PREUC
  consistent_array[na1][na2] = true;
  consistent_array[na1][NA2] = true;
  consistent_array[NA1][na2] = true;
#endif
  consistent_array[NA1][NA2] = true;

#ifdef PREUC
  consistent_array[na2][na1] = true;
  consistent_array[na2][NA1] = true;
  consistent_array[NA2][na1] = true;
#endif
  consistent_array[NA2][NA1] = true;
  for (i = 0; i < NMISMATCHTYPES; i++) {
#ifdef PREUC
    pairdistance_array[i][na1][na2] = score;
    pairdistance_array[i][na1][NA2] = score;
    pairdistance_array[i][NA1][na2] = score;
#endif
    pairdistance_array[i][NA1][NA2] = score;

#ifdef PREUC
    pairdistance_array[i][na2][na1] = score;
    pairdistance_array[i][na2][NA1] = score;
    pairdistance_array[i][NA2][na1] = score;
#endif
    pairdistance_array[i][NA2][NA1] = score;
  }

  return;
}

static void
permute_cases_oneway (int NA1, int NA2, Pairdistance_T score) {
  int i;

#ifdef PREUC
  int na1, na2;

  na1 = tolower(NA1);
  na2 = tolower(NA2);
#endif

#ifdef PREUC
  consistent_array[na1][na2] = true;
  consistent_array[na1][NA2] = true;
  consistent_array[NA1][na2] = true;
#endif
  consistent_array[NA1][NA2] = true;

  for (i = 0; i < NMISMATCHTYPES; i++) {
#ifdef PREUC
    pairdistance_array[i][na1][na2] = score;
    pairdistance_array[i][na1][NA2] = score;
    pairdistance_array[i][NA1][na2] = score;
#endif
    pairdistance_array[i][NA1][NA2] = score;
  }

  return;
}


static void
pairdistance_init (Mode_T mode) {
  int i, j, ptr;
  int c, c1, c2;

  consistent_array = (bool **) CALLOC(128,sizeof(bool *));
  consistent_array[0] = (bool *) CALLOC(128*128,sizeof(bool));
  ptr = 0;
  for (j = 1; j < 128; j++) {
    ptr += 128;
    consistent_array[j] = &(consistent_array[0][ptr]);
  }
  for (i = 0; i < NMISMATCHTYPES; i++) {
    pairdistance_array[i] = (Pairdistance_T **) CALLOC(128,sizeof(Pairdistance_T *));
    pairdistance_array[i][0] = (Pairdistance_T *) CALLOC(128*128,sizeof(Pairdistance_T));
#ifndef HAVE_SSE4_1
    pairdistance_array_plus_128[i] = (Pairdistance_T **) CALLOC(128,sizeof(Pairdistance_T *));
    pairdistance_array_plus_128[i][0] = (Pairdistance_T *) CALLOC(128*128,sizeof(Pairdistance_T));
#endif
    ptr = 0;
    for (j = 1; j < 128; j++) {
      ptr += 128;
      pairdistance_array[i][j] = &(pairdistance_array[i][0][ptr]);
#ifndef HAVE_SSE4_1
      pairdistance_array_plus_128[i][j] = &(pairdistance_array_plus_128[i][0][ptr]);
#endif
    }
  }

#ifdef PREUC
  for (c1 = 'A'; c1 <= 'z'; c1++) {
    for (c2 = 'A'; c2 < 'z'; c2++) {
      pairdistance_array[HIGHQ][c1][c2] = MISMATCH_HIGHQ;
      pairdistance_array[MEDQ][c1][c2] = MISMATCH_MEDQ;
      pairdistance_array[LOWQ][c1][c2] = MISMATCH_LOWQ;
      pairdistance_array[ENDQ][c1][c2] = MISMATCH_ENDQ;
    }
  }
#else
  for (c1 = 'A'; c1 <= 'Z'; c1++) {
    for (c2 = 'A'; c2 < 'Z'; c2++) {
      pairdistance_array[HIGHQ][c1][c2] = MISMATCH_HIGHQ;
      pairdistance_array[MEDQ][c1][c2] = MISMATCH_MEDQ;
      pairdistance_array[LOWQ][c1][c2] = MISMATCH_LOWQ;
      pairdistance_array[ENDQ][c1][c2] = MISMATCH_ENDQ;
    }
  }
#endif

  permute_cases('U','T',FULLMATCH);

  permute_cases('R','A',HALFMATCH);
  permute_cases('R','G',HALFMATCH);

  permute_cases('Y','T',HALFMATCH);
  permute_cases('Y','C',HALFMATCH);

  permute_cases('W','A',HALFMATCH);
  permute_cases('W','T',HALFMATCH);

  permute_cases('S','G',HALFMATCH);
  permute_cases('S','C',HALFMATCH);

  permute_cases('M','A',HALFMATCH);
  permute_cases('M','C',HALFMATCH);

  permute_cases('K','G',HALFMATCH);
  permute_cases('K','T',HALFMATCH);

  permute_cases('H','A',AMBIGUOUS);
  permute_cases('H','T',AMBIGUOUS);
  permute_cases('H','C',AMBIGUOUS);

  permute_cases('B','G',AMBIGUOUS);
  permute_cases('B','C',AMBIGUOUS);
  permute_cases('B','T',AMBIGUOUS);

  permute_cases('V','G',AMBIGUOUS);
  permute_cases('V','A',AMBIGUOUS);
  permute_cases('V','C',AMBIGUOUS);

  permute_cases('D','G',AMBIGUOUS);
  permute_cases('D','A',AMBIGUOUS);
  permute_cases('D','T',AMBIGUOUS);

  permute_cases('N','T',AMBIGUOUS);
  permute_cases('N','C',AMBIGUOUS);
  permute_cases('N','A',AMBIGUOUS);
  permute_cases('N','G',AMBIGUOUS);

  permute_cases('X','T',AMBIGUOUS);
  permute_cases('X','C',AMBIGUOUS);
  permute_cases('X','A',AMBIGUOUS);
  permute_cases('X','G',AMBIGUOUS);

  if (mode == CMET_STRANDED || mode == CMET_NONSTRANDED) {
    /* Query-T can match Genomic-C */
    permute_cases_oneway('T','C',FULLMATCH);
    permute_cases_oneway('A','G',FULLMATCH);
  }

  for (c = 'A'; c < 'Z'; c++) {
    permute_cases(c,c,FULLMATCH);
  }

#ifndef HAVE_SSE4_1
#ifdef PREUC
  for (i = 0; i < NMISMATCHTYPES; i++) {
    for (c1 = 'A'; c1 <= 'z'; c1++) {
      for (c2 = 'A'; c2 < 'z'; c2++) {
	pairdistance_array_plus_128[i][c1][c2] = 128 + pairdistance_array[i][c1][c2];
      }
    }
  }
#else
  for (i = 0; i < NMISMATCHTYPES; i++) {
    for (c1 = 'A'; c1 <= 'Z'; c1++) {
      for (c2 = 'A'; c2 < 'Z'; c2++) {
	pairdistance_array_plus_128[i][c1][c2] = 128 + pairdistance_array[i][c1][c2];
      }
    }
  }
#endif
#endif

  return;
}


/************************************************************************/


#if 0
/* static int max_jump_penalty_lookup; */
static Pairdistance_T *jump_penalty_array[NJUMPTYPES];

/* For lengths of 1,2,3, returns open.  Then, for lengths of 4,5,6,
   returns open + 3*extend.  Add extra extend to penalty, so ---- is
   preferred over ---|-. */
static int
jump_penalty (int length, int open, int extend) {
#ifdef CODONPENALTY
  int ncodons;
#endif

#ifdef CODONPENALTY
  ncodons = (length - 1)/3;
  return open + extend + ncodons*3*extend;
#else
  return open + extend*length;	/* was extend*(length-1), but now matches result of jump_penalty_init  */
#endif
}
#endif


#if 0
static void
jump_penalty_init (int maxlookback, int extraquerygap, int maxpeelback,
		   int extramaterial_end, int extramaterial_paired) {
  int max_lookup, max_rlength, max_glength, i;
  int length, remainder, phase;
  int paired_highq_penalty = PAIRED_OPEN_HIGHQ + PAIRED_EXTEND_HIGHQ,
    paired_medq_penalty = PAIRED_OPEN_MEDQ + PAIRED_EXTEND_MEDQ,
    paired_lowq_penalty = PAIRED_OPEN_LOWQ + PAIRED_EXTEND_LOWQ,
    single_highq_penalty = SINGLE_OPEN_HIGHQ + SINGLE_EXTEND_HIGHQ,
    single_medq_penalty = SINGLE_OPEN_MEDQ + SINGLE_EXTEND_MEDQ,
    single_lowq_penalty = SINGLE_OPEN_LOWQ + SINGLE_EXTEND_LOWQ,
    end_highq_penalty = END_OPEN_HIGHQ + END_EXTEND_HIGHQ,
    end_medq_penalty = END_OPEN_MEDQ + END_EXTEND_MEDQ,
    end_lowq_penalty = END_OPEN_LOWQ + END_EXTEND_LOWQ,
    cdna_highq_penalty = CDNA_OPEN_HIGHQ + CDNA_EXTEND_HIGHQ,
    cdna_medq_penalty = CDNA_OPEN_MEDQ + CDNA_EXTEND_MEDQ,
    cdna_lowq_penalty = CDNA_OPEN_LOWQ + CDNA_EXTEND_LOWQ;
  int paired_highq_delta = 3*PAIRED_EXTEND_HIGHQ,
    paired_medq_delta = 3*PAIRED_EXTEND_MEDQ,
    paired_lowq_delta = 3*PAIRED_EXTEND_LOWQ,
    single_highq_delta = 3*SINGLE_EXTEND_HIGHQ,
    single_medq_delta = 3*SINGLE_EXTEND_MEDQ,
    single_lowq_delta = 3*SINGLE_EXTEND_LOWQ,
    end_highq_delta = 3*END_EXTEND_HIGHQ,
    end_medq_delta = 3*END_EXTEND_MEDQ,
    end_lowq_delta = 3*END_EXTEND_LOWQ,
    cdna_highq_delta = 3*CDNA_EXTEND_HIGHQ,
    cdna_medq_delta = 3*CDNA_EXTEND_MEDQ,
    cdna_lowq_delta = 3*CDNA_EXTEND_LOWQ;

  compute_maxlengths(&max_rlength,&max_glength,
		     maxlookback,extraquerygap,maxpeelback,
		     extramaterial_end,extramaterial_paired);

  if (max_rlength > max_glength) {
    max_lookup = max_rlength;
  } else {
    max_lookup = max_glength;
  }
  if ((remainder = max_lookup % 3) != 0) {
    max_lookup += (3 - remainder);
  }

  /* Set global */
  /* max_jump_penalty_lookup = max_lookup; */

  for (i = 0; i < NJUMPTYPES; i++) {
    jump_penalty_array[i] = (Pairdistance_T *) CALLOC(max_lookup+1,sizeof(Pairdistance_T));
  }

  length = 1;
  while (length+2 <= max_lookup) {
    for (phase = 0; phase < 3; phase++) {
      jump_penalty_array[PAIRED_HIGHQ][length] = paired_highq_penalty;
      jump_penalty_array[PAIRED_MEDQ][length] = paired_medq_penalty;
      jump_penalty_array[PAIRED_LOWQ][length] = paired_lowq_penalty;
      jump_penalty_array[SINGLE_HIGHQ][length] = single_highq_penalty;
      jump_penalty_array[SINGLE_MEDQ][length] = single_medq_penalty;
      jump_penalty_array[SINGLE_LOWQ][length] = single_lowq_penalty;
      jump_penalty_array[END_HIGHQ][length] = end_highq_penalty;
      jump_penalty_array[END_MEDQ][length] = end_medq_penalty;
      jump_penalty_array[END_LOWQ][length] = end_lowq_penalty;
      jump_penalty_array[CDNA_HIGHQ][length] = cdna_highq_penalty;
      jump_penalty_array[CDNA_MEDQ][length] = cdna_medq_penalty;
      jump_penalty_array[CDNA_LOWQ][length] = cdna_lowq_penalty;
      length++;
    }
    paired_highq_penalty += paired_highq_delta;
    paired_medq_penalty += paired_medq_delta;
    paired_lowq_penalty += paired_lowq_delta;
    single_highq_penalty += single_highq_delta;
    single_medq_penalty += single_medq_delta;
    single_lowq_penalty += single_lowq_delta;
    end_highq_penalty += end_highq_delta;
    end_medq_penalty += end_medq_delta;
    end_lowq_penalty += end_lowq_delta;
    cdna_highq_penalty += cdna_highq_delta;
    cdna_medq_penalty += cdna_medq_delta;
    cdna_lowq_penalty += cdna_lowq_delta;
  }
  return;
}
#endif


/************************************************************************/

void
Dynprog_init (int maxlookback, int extraquerygap, int maxpeelback,
	      int extramaterial_end, int extramaterial_paired, Mode_T mode) {
  pairdistance_init(mode);
#if 0
  jump_penalty_init(maxlookback,extraquerygap,maxpeelback,
		    extramaterial_end,extramaterial_paired);
#endif
  return;
}

void
Dynprog_term (void) {
  int i;

#if 0
  for (i = 0; i < NJUMPTYPES; i++) {
    FREE(jump_penalty_array[i]);
  }
#endif

  for (i = 0; i < NMISMATCHTYPES; i++) {
    /*
    for (j = 0; j < 128; j++) {
      FREE(pairdistance_array[i][j]);
    }
    */
#ifndef HAVE_SSE4_1
    FREE(pairdistance_array_plus_128[i][0]);
    FREE(pairdistance_array_plus_128[i]);
#endif
    FREE(pairdistance_array[i][0]);
    FREE(pairdistance_array[i]);
  }
  /*
  for (j = 0; j < 128; j++) {
    FREE(consistent_array[j]);
  }
  */
  FREE(consistent_array[0]);
  FREE(consistent_array);

  return;
}

/************************************************************************/

static void
compute_bands (int *lband, int *uband, int rlength, int glength, int extraband, bool widebandp) {
  if (widebandp == false) {
    /* Just go along main diagonal */
    *lband = extraband;
    *uband = extraband;
  } else if (glength >= rlength) {
    /* Widen band to right to reach destination */
    *uband = glength - rlength + extraband;
    *lband = extraband;
  } else {
    /* Widen band to left to reach destination */
    *lband = rlength - glength + extraband;
    *uband = extraband;
  }
  return;
}


#if 0
static void
make_complement_buffered (char *complement, char *sequence, unsigned int length) {
  int i, j;

  /* complement = (char *) CALLOC(length+1,sizeof(char)); */
  for (i = length-1, j = 0; i >= 0; i--, j++) {
    complement[j] = complCode[(int) sequence[i]];
  }
  complement[length] = '\0';
  return;
}
#endif

static void
make_complement_inplace (char *sequence, unsigned int length) {
  char temp;
  unsigned int i, j;

  for (i = 0, j = length-1; i < length/2; i++, j--) {
    temp = complCode[(int) sequence[i]];
    sequence[i] = complCode[(int) sequence[j]];
    sequence[j] = temp;
  }
  if (i == j) {
    sequence[i] = complCode[(int) sequence[i]];
  }

  return;
}

#if 0
static void
get_genomic_seg (Chrpos_T genomicpos, int length, Univcoord_T chroffset,
		 bool watsonp, char *dest) {
  if (watsonp) {
    Genome_fill_buffer_blocks(chroffset + genomicpos,length,dest);

  } else {
    Genome_fill_buffer_blocks(chrhigh - genomicpos,length,dest);
    make_complement_inplace(dest,length);
  }
  return;
}
#endif


#ifdef DEBUG15
/* For debugging of SIMD procedures*/
static void
print_vector_8 (__m128i x, int r, int c, char *label) {
  __m128i a[1];
  Score8_T *s = a;

  _mm_lfence();			/* Needed to print correct values */
  _mm_store_si128(a,x);
  printf("%d,%d %s: %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n",
	 r,c,label,s[0],s[1],s[2],s[3],s[4],s[5],s[6],s[7],s[8],s[9],s[10],s[11],s[12],s[13],s[14],s[15]);
  return;
}

static void
print_vector_16 (__m128i x, int r, int c, char *label) {
  __m128i a[1];
  Score16_T *s = a;

  _mm_lfence();			/* Needed to print correct values */
  _mm_store_si128(a,x);
  printf("%d,%d %s: %d %d %d %d %d %d %d %d\n",r,c,label,s[0],s[1],s[2],s[3],s[4],s[5],s[6],s[7]);
  return;
}
#endif


#ifdef DEBUG14
static void
banded_matrix8_compare (Score8_T **matrix1, Score32_T **matrix2, int rlength, int glength,
			int lband, int uband, char *rsequence, char *gsequence, char *gsequence_alt,
			int goffset, Univcoord_T chroffset, Univcoord_T chrhigh, bool watsonp, bool revp) {
  int r, c, rlo, rhigh;

  for (c = 1; c <= glength; c++) {
    if ((rlo = c - uband) < 1) {
      rlo = 1;
    };

    if ((rhigh = c + lband) > rlength) {
      rhigh = rlength;
    }

    for (r = rlo; r <= rhigh; r++) {
      if (matrix1[c][r] <= NEG_INFINITY_8 + 30 && matrix2[c][r] <= NEG_INFINITY_8 + 30) {
	/* Okay */
      } else if (matrix1[c][r] != matrix2[c][r]) {
	printf("At %d,%d, value %d != value %d\n",r,c,matrix1[c][r],matrix2[c][r]);

	Matrix8_print(matrix1,rlength,glength,rsequence,gsequence,gsequence_alt,
		      goffset,chroffset,chrhigh,watsonp,revp);
	Matrix32_print(matrix2,rlength,glength,rsequence,gsequence,gsequence_alt,
		       goffset,chroffset,chrhigh,watsonp,revp);
	abort();
      }
    }
  }

  return;
}

static void
banded_matrix16_compare (Score16_T **matrix1, Score32_T **matrix2, int rlength, int glength,
			 int lband, int uband, char *rsequence, char *gsequence, char *gsequence_alt,
			 int goffset, Univcoord_T chroffset, Univcoord_T chrhigh, bool watsonp, bool revp) {
  int r, c, rlo, rhigh;

  for (c = 1; c <= glength; c++) {
    if ((rlo = c - uband) < 1) {
      rlo = 1;
    };

    if ((rhigh = c + lband) > rlength) {
      rhigh = rlength;
    }

    for (r = rlo; r <= rhigh; r++) {
      if (matrix1[c][r] <= NEG_INFINITY_16 + 30 && matrix2[c][r] <= NEG_INFINITY_16 + 30) {
	/* Okay */
      } else if (matrix1[c][r] != matrix2[c][r]) {
	printf("At %d,%d, value %d != value %d\n",r,c,matrix1[c][r],matrix2[c][r]);

	Matrix16_print(matrix1,rlength,glength,rsequence,gsequence,gsequence_alt,
		       goffset,chroffset,chrhigh,watsonp,revp);
	Matrix32_print(matrix2,rlength,glength,rsequence,gsequence,gsequence_alt,
		       goffset,chroffset,chrhigh,watsonp,revp);
	abort();
      }
    }
  }

  return;
}
#endif

#ifdef DEBUG14
static void
banded_directions8_compare_nogap (Score8_T **matrix, Direction8_T **directions1, Direction32_T **directions2, int rlength, int glength,
				  int lband, int uband) {
  int r, c, rlo, rhigh;

  for (c = 1; c <= glength; c++) {
    if ((rlo = c - uband) < 1) {
      rlo = 1;
    };

    if ((rhigh = c + lband) > rlength) {
      rhigh = rlength;
    }

    for (r = rlo; r <= rhigh; r++) {
      if (matrix[c][r] < NEG_INFINITY_8 + 30) {
	/* Don't check */

      } else if (directions1[c][r] == 0) {
	if (directions2[c][r] == 0) {
	} else {
	  printf("At %d,%d, nogap dir %d != dir %d\n",r,c,directions1[c][r],directions2[c][r]);
	  abort();
	}

      } else if (directions1[c][r] == 1) {
	if (directions2[c][r] == 1) {
	} else {
	  printf("At %d,%d, nogap dir %d != dir %d\n",r,c,directions1[c][r],directions2[c][r]);
	  abort();
	}

      } else {
	if (directions2[c][r] == 0 || directions2[c][r] == 0) {
	  printf("At %d,%d, nogap dir %d != dir %d\n",r,c,directions1[c][r],directions2[c][r]);
	  abort();
	}
      }
    }
  }

  return;
}


static void
banded_directions16_compare_nogap (Direction16_T **directions1, Direction32_T **directions2, int rlength, int glength,
				   int lband, int uband) {
  int r, c, rlo, rhigh;

  for (c = 1; c <= glength; c++) {
    if ((rlo = c - uband) < 1) {
      rlo = 1;
    };

    if ((rhigh = c + lband) > rlength) {
      rhigh = rlength;
    }

    for (r = rlo; r <= rhigh; r++) {
      if (directions1[c][r] == 0) {
	if (directions2[c][r] == 0) {
	} else {
	  printf("At %d,%d, nogap dir %d != dir %d\n",r,c,directions1[c][r],directions2[c][r]);
	  abort();
	}

      } else if (directions1[c][r] == 1) {
	if (directions2[c][r] == 1) {
	} else {
	  printf("At %d,%d, nogap dir %d != dir %d\n",r,c,directions1[c][r],directions2[c][r]);
	  abort();
	}

      } else {
	if (directions2[c][r] == 0 || directions2[c][r] == 0) {
	  printf("At %d,%d, nogap dir %d != dir %d\n",r,c,directions1[c][r],directions2[c][r]);
	  abort();
	}
      }
    }
  }

  return;
}
#endif

#ifdef DEBUG14
static void
banded_directions8_compare_Egap (Score8_T **matrix1, Direction8_T **directions1, Direction32_T **directions2,
				 int rlength, int glength, int lband, int uband) {
  int r, c, rlo, rhigh, last_check;

  for (c = 1; c <= glength; c++) {
    if ((rlo = c - uband) < 1) {
      rlo = 1;
    };

    if ((rhigh = c + lband) <= rlength) {
      /* Don't check rhigh.  Egap direction derives from a comparison
	 of NEG_INFINITY values, and we should never reach here from
	 directions_nogap anyway. */
      last_check = rhigh - 1;

    } else {
      /* Do check rhigh, which contains instructions for the bottom row */
      rhigh = rlength;
      last_check = rhigh;
    }

    for (r = rlo; r <= last_check; r++) {
      if (matrix1[c][r] < NEG_INFINITY_8 + 30) {
	/* Don't check */

      } else if (directions1[c][r] == 0) {
	if (directions2[c][r] == 0) {
	} else {
	  printf("At %d,%d, Egap dir %d != dir %d\n",r,c,directions1[c][r],directions2[c][r]);
	  abort();
	}

      } else if (directions1[c][r] == 1) {
	if (directions2[c][r] == 1) {
	} else {
	  printf("At %d,%d, Egap dir %d != dir %d\n",r,c,directions1[c][r],directions2[c][r]);
	  abort();
	}

      } else {
	if (directions2[c][r] == 0 || directions2[c][r] == 0) {
	  printf("At %d,%d, Egap dir %d != dir %d\n",r,c,directions1[c][r],directions2[c][r]);
	  abort();
	}
      }
    }
  }

  return;
}

static void
banded_directions16_compare_Egap (Direction16_T **directions1, Direction32_T **directions2, int rlength, int glength,
				  int lband, int uband) {
  int r, c, rlo, rhigh, last_check;

  for (c = 1; c <= glength; c++) {
    if ((rlo = c - uband) < 1) {
      rlo = 1;
    };

    if ((rhigh = c + lband) <= rlength) {
      /* Don't check rhigh.  Egap direction derives from a comparison
	 of NEG_INFINITY values, and we should never reach here from
	 directions_nogap anyway. */
      last_check = rhigh - 1;

    } else {
      /* Do check rhigh, which contains instructions for the bottom row */
      rhigh = rlength;
      last_check = rhigh;
    }

    for (r = rlo; r <= last_check; r++) {
      if (directions1[c][r] == 0) {
	if (directions2[c][r] == 0) {
	} else {
	  printf("At %d,%d, Egap dir %d != dir %d\n",r,c,directions1[c][r],directions2[c][r]);
	  abort();
	}
      } else if (directions1[c][r] == 1) {
	if (directions2[c][r] == 1) {
	} else {
	  printf("At %d,%d, Egap dir %d != dir %d\n",r,c,directions1[c][r],directions2[c][r]);
	  abort();
	}

      } else {
	if (directions2[c][r] == 0 || directions2[c][r] == 0) {
	  printf("At %d,%d, Egap dir %d != dir %d\n",r,c,directions1[c][r],directions2[c][r]);
	  abort();
	}
      }
    }
  }

  return;
}
#endif


#ifdef DEBUG14
static void
banded_directions8_compare_Fgap (Score8_T **matrix1, Direction8_T **directions1, Direction32_T **directions2,
				 int rlength, int glength, int lband, int uband) {
  int r, c, rlo, rhigh, first_check;

  for (c = 1; c <= glength; c++) {
    if ((rlo = c - uband) < 1) {
      first_check = rlo = 1;
    } else {
      first_check = rlo + 1;
    }

    if ((rhigh = c + lband) > rlength) {
      rhigh = rlength;
    }

    for (r = first_check; r <= rhigh; r++) {
      if (matrix1[c][r] < NEG_INFINITY_8 + 30) {
	/* Don't check */

      } else if (directions1[c][r] == 0) {
	if (directions2[c][r] == 0) {
	} else {
	  printf("At %d,%d, Fgap dir %d != dir %d.  Score is %d\n",
		 r,c,directions1[c][r],directions2[c][r],matrix1[c][r]);
	  abort();
	}

      } else if (directions1[c][r] == 1) {
	if (directions2[c][r] == 1) {
	} else {
	  printf("At %d,%d, Fgap dir %d != dir %d.  Score is %d\n",
		 r,c,directions1[c][r],directions2[c][r],matrix1[c][r]);
	  abort();
	}

      } else {
	if (directions2[c][r] == 0 || directions2[c][r] == 0) {
	  printf("At %d,%d, Fgap dir %d != dir %d.  Score is %d\n",
		 r,c,directions1[c][r],directions2[c][r],matrix1[c][r]);
	  abort();
	}
      }
    }
  }

  return;
}

static void
banded_directions16_compare_Fgap (Direction16_T **directions1, Direction32_T **directions2, int rlength, int glength,
				  int lband, int uband) {
  int r, c, rlo, rhigh, first_check;

  for (c = 1; c <= glength; c++) {
    if ((rlo = c - uband) < 1) {
      first_check = rlo = 1;
    } else {
      first_check = rlo + 1;
    }

    if ((rhigh = c + lband) > rlength) {
      rhigh = rlength;
    }

    for (r = first_check; r <= rhigh; r++) {
      if (directions1[c][r] == 0) {
	if (directions2[c][r] == 0) {
	} else {
	  printf("At %d,%d, Fgap dir %d != dir %d\n",r,c,directions1[c][r],directions2[c][r]);
	  abort();
	}
      } else if (directions1[c][r] == 1) {
	if (directions2[c][r] == 1) {
	} else {
	  printf("At %d,%d, Fgap dir %d != dir %d\n",r,c,directions1[c][r],directions2[c][r]);
	  abort();
	}

      } else {
	if (directions2[c][r] == 0 || directions2[c][r] == 0) {
	  printf("At %d,%d, Fgap dir %d != dir %d\n",r,c,directions1[c][r],directions2[c][r]);
	  abort();
	}
      }
    }
  }

  return;
}
#endif


#if defined(HAVE_SSE4_1) || defined(HAVE_SSE2)
/* Makes a matrix of dimensions 0..rlength x 0..glength inclusive */
static Score8_T **
aligned_score8_alloc (int rlength, int glength, void **ptrs, void *space) {
  Score8_T **matrix, *ptr;
  int c;

  matrix = (Score8_T **) ptrs;

  ptr = (Score8_T *) space;
  matrix[0] = &(ptr[SIMD_NCHARS - 1]);	/* Want aligned row to be r = 1, 17, ... */
  for (c = 1; c <= glength; c++) {
    ptr += rlength + SIMD_NCHARS;
    matrix[c] = &(ptr[SIMD_NCHARS - 1]);	/* Want aligned row to be r = 1, 17, ... */
  }
#ifdef DEBUG2
  memset((void *) matrix[0],0,(glength+1)*(rlength+SIMD_NCHARS)*sizeof(Score8_T));
#endif

  return matrix;
}

/* No initialization to DIAG (0), for directions_Egap and directions_nogap */
static Score8_T **
aligned_directions8_alloc (int rlength, int glength, void **ptrs, void *space) {
  Score8_T **matrix, *ptr;
  int c;

  matrix = (Score8_T **) ptrs;

  ptr = (Score8_T *) space;
  matrix[0] = &(ptr[SIMD_NCHARS - 1]);	/* Want aligned row to be r = 1, 17, ... */
  for (c = 1; c <= glength; c++) {
    ptr += rlength + SIMD_NCHARS;
    matrix[c] = &(ptr[SIMD_NCHARS - 1]);	/* Want aligned row to be r = 1, 17, ... */
  }
#ifdef DEBUG2
  memset((void *) matrix[0],/*DIAG*/0,(glength+1)*(rlength+SIMD_NCHARS)*sizeof(Score8_T));
#endif

  return matrix;
}

/* Initialization to DIAG (0), for directions_Fgap */
static Score8_T **
aligned_directions8_calloc (int rlength, int glength, void **ptrs, void *space) {
  Score8_T **matrix, *ptr;
  int c;

  matrix = (Score8_T **) ptrs;

  ptr = (Score8_T *) space;
  matrix[0] = &(ptr[SIMD_NCHARS - 1]);	/* Want aligned row to be r = 1, 17, ... */
  for (c = 1; c <= glength; c++) {
    ptr += rlength + SIMD_NCHARS;
    matrix[c] = &(ptr[SIMD_NCHARS - 1]);	/* Want aligned row to be r = 1, 17, ... */
  }
  memset((void *) matrix[0],/*DIAG*/0,(glength+1)*(rlength+SIMD_NCHARS)*sizeof(Score8_T));

  return matrix;
}



/* Makes a matrix of dimensions 0..rlength x 0..glength inclusive */
static Score16_T **
aligned_score16_alloc (int rlength, int glength, void **ptrs, void *space) {
  Score16_T **matrix, *ptr;
  int c;

  matrix = (Score16_T **) ptrs;

  ptr = (Score16_T *) space;
  matrix[0] = &(ptr[SIMD_NSHORTS - 1]);	/* Want aligned row to be r = 1, 9, 17, ... */
  for (c = 1; c <= glength; c++) {
    ptr += rlength + SIMD_NSHORTS;
    matrix[c] = &(ptr[SIMD_NSHORTS - 1]);	/* Want aligned row to be r = 1, 9, 17, ... */
  }
#ifdef DEBUG2
  memset((void *) matrix[0],0,(glength+1)*(rlength+SIMD_NSHORTS)*sizeof(Score16_T));
#endif
  return matrix;
}

/* No initialization to DIAG (0), for directions_Egap and directions_nogap */
static Score16_T **
aligned_directions16_alloc (int rlength, int glength, void **ptrs, void *space) {
  Score16_T **matrix, *ptr;
  int c;

  matrix = (Score16_T **) ptrs;

  ptr = (Score16_T *) space;
  matrix[0] = &(ptr[SIMD_NSHORTS - 1]);	/* Want aligned row to be r = 1, 9, 17, ... */
  for (c = 1; c <= glength; c++) {
    ptr += rlength + SIMD_NSHORTS;
    matrix[c] = &(ptr[SIMD_NSHORTS - 1]);	/* Want aligned row to be r = 1, 9, 17, ... */
  }
#ifdef DEBUG2
  memset((void *) matrix[0],/*DIAG*/0,(glength+1)*(rlength+SIMD_NSHORTS)*sizeof(Score16_T));
#endif

  return matrix;
}

/* Initialization to DIAG (0), for directions_Fgap */
static Score16_T **
aligned_directions16_calloc (int rlength, int glength, void **ptrs, void *space) {
  Score16_T **matrix, *ptr;
  int c;

  matrix = (Score16_T **) ptrs;

  ptr = (Score16_T *) space;
  matrix[0] = &(ptr[SIMD_NSHORTS - 1]);	/* Want aligned row to be r = 1, 9, 17, ... */
  for (c = 1; c <= glength; c++) {
    ptr += rlength + SIMD_NSHORTS;
    matrix[c] = &(ptr[SIMD_NSHORTS - 1]);	/* Want aligned row to be r = 1, 9, 17, ... */
  }
  memset((void *) matrix[0],/*DIAG*/0,(glength+1)*(rlength+SIMD_NSHORTS)*sizeof(Score16_T));

  return matrix;
}
#endif


#if !defined(HAVE_SSE2) || defined(DEBUG14)
static Score32_T **
compute_scores_standard (Direction32_T ***directions_nogap, Direction32_T ***directions_Egap, Direction32_T ***directions_Fgap,
			 T this, char *rsequence, char *gsequence, char *gsequence_alt,
			 int rlength, int glength,
			 int goffset, Univcoord_T chroffset, Univcoord_T chrhigh, bool watsonp,
			 Mismatchtype_T mismatchtype, Score32_T open, Score32_T extend,
			 int lband, int uband, bool jump_late_p, bool revp) {
#ifdef DEBUG12
  Score32_T bestscore;
  Direction32_T bestdir;
  struct Int3_T **matrix3;
#endif
  Score32_T penalty;
  Score32_T **matrix;
  Score32_T c_gap, *r_gap, *nogap, last_nogap, prev_nogap, first_nogap;
  int r, c, na1, na2;
  char na2_alt;
  Score32_T score, pairscore;
  int rlo, rhigh;
  Pairdistance_T **pairdistance_array_type;

  pairdistance_array_type = pairdistance_array[mismatchtype];

  debug(printf("compute_scores_standard: "));
  debug(printf("Lengths are %d and %d, so bands are %d on left and %d on right\n",rlength,glength,lband,uband));

  matrix = Matrix32_alloc(rlength,glength,this->matrix_ptrs,this->matrix_space);
  *directions_nogap = Directions32_alloc(rlength,glength,this->directions_ptrs_0,this->directions_space_0);
  *directions_Egap = Directions32_alloc(rlength,glength,this->directions_ptrs_1,this->directions_space_1);
  *directions_Fgap = Directions32_alloc(rlength,glength,this->directions_ptrs_2,this->directions_space_2);
  /* (*directions_nogap)[0][0] = STOP; -- Check for r > 0 && c > 0 instead */

  /* Row 0 initialization */
  penalty = open;
  for (c = 1; c <= uband && c <= glength; c++) {
    penalty += extend;
    matrix[c][0] = penalty;
    (*directions_Egap)[c][0] = HORIZ;
    (*directions_nogap)[c][0] = HORIZ;
  }
#if 0
  /* Already initialized to DIAG */
  (*directions_Egap)[1][0] = DIAG; /* previously used STOP */
#endif

  /* Column 0 initialization */
  penalty = open;
  for (r = 1; r <= lband && r <= rlength; r++) {
    penalty += extend;
    matrix[0][r] = penalty;
    (*directions_Fgap)[0][r] = VERT;
    (*directions_nogap)[0][r] = VERT;
  }
#if 0
  /* Already initialized to DIAG */
  (*directions_Fgap)[0][1] = DIAG; /* previously used STOP */
#endif

  r_gap = (Score32_T *) CALLOC(rlength+1,sizeof(Score32_T));
  nogap = (Score32_T *) CALLOC(rlength+1,sizeof(Score32_T));
  nogap[0] = 0;
  penalty = open;
  for (r = 1; r <= lband && r <= rlength; r++) {
    penalty += extend;
    r_gap[r] = NEG_INFINITY_32;
    nogap[r] = penalty;
  }
  for ( ; r <= rlength; r++) {
    r_gap[r] = NEG_INFINITY_32;
    nogap[r] = NEG_INFINITY_32;
  }


#ifdef DEBUG12
  matrix3 = Matrix3_alloc(rlength,glength,this->matrix3_ptrs,this->matrix3_space);
  matrix3[0][0].nogap = 0;
  matrix3[0][0].Egap = matrix3[0][0].Fgap = NEG_INFINITY_32;

  /* Row 0 initialization */
  penalty = open;
  for (c = 1; c <= uband && c <= glength; c++) {
    penalty += extend;
    matrix3[c][0].nogap = penalty;
    matrix3[c][0].Egap = penalty;
    matrix3[c][0].Fgap = NEG_INFINITY_32;
  }

  /* Column 0 initialization */
  penalty = open;
  for (r = 1; r <= lband && r <= rlength; r++) {
    penalty += extend;
    matrix3[0][r].nogap = penalty;
    matrix3[0][r].Egap = NEG_INFINITY_32;
    matrix3[0][r].Fgap = penalty;
  }
#endif


  first_nogap = 0;
  if (jump_late_p) {
    penalty = open + extend;
    for (c = 1; c <= glength; c++) {
      if (gsequence) {
	na2 = revp ? gsequence[1-c] : gsequence[c-1];
	na2_alt = revp ? gsequence_alt[1-c] : gsequence_alt[c-1];
      } else if (revp == false) {
	na2 = get_genomic_nt(&na2_alt,goffset+c-1,chroffset,chrhigh,watsonp);
      } else {
	na2 = get_genomic_nt(&na2_alt,goffset+1-c,chroffset,chrhigh,watsonp);
      }

      c_gap = NEG_INFINITY_32;
      if (c == 1) {
	rlo = 1;
	prev_nogap = 0; /* was nogap[rlo-1] */
	last_nogap = penalty;
      } else if ((rlo = c - uband) < 1) {
	rlo = 1;
	prev_nogap = penalty;
	penalty += extend;
	last_nogap = penalty;
      } else if (rlo == 1) {
	prev_nogap = penalty;
	/* penalty += extend; */
	last_nogap = NEG_INFINITY_32;
#ifdef DEBUG12
	matrix3[c][rlo-1].Fgap = NEG_INFINITY_32;
	matrix3[c][rlo-1].nogap = NEG_INFINITY_32;
#endif
      } else {
	prev_nogap = first_nogap;
	last_nogap = NEG_INFINITY_32;
#ifdef DEBUG12
	matrix3[c][rlo-1].Fgap = NEG_INFINITY_32;
	matrix3[c][rlo-1].nogap = NEG_INFINITY_32;
#endif
      }

      if ((rhigh = c + lband) > rlength) {
	rhigh = rlength;
#ifdef DEBUG12
      } else {
	matrix3[c-1][rhigh].Egap = NEG_INFINITY_32;
	matrix3[c-1][rhigh].nogap = NEG_INFINITY_32;
#endif
      }

      for (r = rlo; r <= rhigh; r++) {
	na1 = revp ? rsequence[1-r] : rsequence[r-1];

#ifdef DEBUG12
	/* FGAP */
	bestscore = matrix3[c][r-1].nogap + open /* + extend */;
	bestdir = DIAG;

	if ((score = matrix3[c][r-1].Fgap /* + extend */) >= bestscore) {  /* Use >= for jump late */
	  bestscore = score;
	  bestdir = VERT;
	}

	matrix3[c][r].Fgap = bestscore + extend;
#endif

	/* FGAP alt */
#ifdef DEBUG12A
	printf("Fgap at r %d, c %d: matrix3[c][r-1].nogap %d vs last_nogap %d\n",r,c,matrix3[c][r-1].nogap,last_nogap);
	printf("Fgap at r %d, c %d: matrix3[c][r-1].Fgap %d vs c_gap %d\n",r,c,matrix3[c][r-1].Fgap,c_gap);
#endif
#ifdef DEBUG12
	assert(matrix3[c][r-1].nogap == last_nogap);
	assert(matrix3[c][r-1].Fgap == c_gap);
#endif
	/* debug2(printf("Fgap at r %d, c %d: c_gap + extend %d vs last_nogap + open + extend %d\n",r,c,c_gap + extend,last_nogap + open + extend)); */
	if (c_gap /* + extend */ >= (score = last_nogap + open /* + extend */)) {  /* Use >= for jump late */
	  c_gap += extend;
	  (*directions_Fgap)[c][r] = VERT;
	} else {
	  c_gap = score + extend;
	  /* bestdir2 = DIAG; -- Already initialized to DIAG */
	}


#ifdef DEBUG12
	/* EGAP */
	bestscore = matrix3[c-1][r].nogap + open /* + extend */;
	bestdir = DIAG;

	if ((score = matrix3[c-1][r].Egap /* + extend */) >= bestscore) {  /* Use >= for jump late */
	  bestscore = score;
	  bestdir = HORIZ;
	}

	matrix3[c][r].Egap = bestscore + extend;
#endif

	/* EGAP alt */
#ifdef DEBUG12A
	printf("Egap at r %d, c %d: matrix3[c-1][r].nogap %d vs nogap[r] %d\n",r,c,matrix3[c-1][r].nogap,nogap[r]);
	printf("Egap at r %d, c %d: matrix3[c-1][r].Egap %d vs r_gap[r] %d\n",r,c,matrix3[c-1][r].Egap,r_gap[r]);
#endif
#ifdef DEBUG12
	assert(matrix3[c-1][r].nogap == nogap[r]);
	assert(matrix3[c-1][r].Egap == r_gap[r]);
#endif
	/* debug2(printf("Egap at r %d, c %d: r_gap[r] %d vs nogap[r] + open %d\n",r,c,r_gap[r],nogap[r]+open)); */
	if (r_gap[r] /* + extend */ >= (score = nogap[r] + open /* + extend */)) {  /* Use >= for jump late */
	  r_gap[r] += extend;
	  (*directions_Egap)[c][r] = HORIZ;
	} else {
	  r_gap[r] = score + extend;
	  /* bestdir2 = DIAG; -- Already initialized to DIAG */
	}


	/* NOGAP */
	pairscore = pairdistance_array_type[na1][na2];
	if ((score = pairdistance_array_type[na1][(int) na2_alt]) > pairscore) {
	  pairscore = score;
	}
#ifdef DEBUG12
	bestscore = matrix3[c-1][r-1].nogap + pairscore;
	bestdir = DIAG;
      
	if ((score = matrix3[c][r].Egap) >= bestscore) {  /* Use >= for jump late */
	  bestscore = score;
	  bestdir = HORIZ;
	}

	if ((score = matrix3[c][r].Fgap) >= bestscore) {  /* Use >= for jump late */
	  bestscore = score;
	  bestdir = VERT;
	}

	matrix3[c][r].nogap = bestscore;
#endif

	/* NOGAP alt */
#ifdef DEBUG12A
	printf("nogap at r %d, c %d: matrix3[c-1][r-1].nogap %d vs prev_nogap %d\n",r,c,matrix3[c-1][r-1].nogap,prev_nogap);
	printf("nogap at r %d, c %d: matrix3[c][r].Fgap %d vs c_gap %d\n",r,c,matrix3[c][r].Fgap,c_gap);
	printf("nogap at r %d, c %d: matrix3[c][r].Egap %d vs r_gap[r] %d\n",r,c,matrix3[c][r].Egap,r_gap[r]);
#endif
#ifdef DEBUG12
	assert(matrix3[c-1][r-1].nogap == prev_nogap);
	assert(matrix3[c][r].Fgap == c_gap);
	assert(matrix3[c][r].Egap == r_gap[r]);
#endif
	last_nogap = prev_nogap + pairscore;
	/* bestdir2 = DIAG; -- Already initialized to DIAG */
	/* debug2(printf("assign nogap at r %d, c %d: H + pairscore %d vs r_horiz + extend %d vs vert + extend %d\n",
	   r,c,last_nogap,r_gap[r],c_gap)); */
	if (r_gap[r] >= last_nogap) {  /* Use >= for jump late */
	  last_nogap = r_gap[r];
	  (*directions_nogap)[c][r] = HORIZ;
	}
	if (c_gap >= last_nogap) {  /* Use >= for jump late */
	  last_nogap = c_gap;
	  (*directions_nogap)[c][r] = VERT;
	}
	/* (*directions_nogap)[c][r] = bestdir2; */

	prev_nogap = nogap[r];	/* Save for next inner loop, before we wipe it out */
	matrix[c][r] = nogap[r] = last_nogap;	/* Save for next outer loop */
	if (r == rlo) {
	  debug12a(printf("At row %d, storing first_nogap to be nogap[r] %d\n",r,nogap[r]));
	  first_nogap = last_nogap;
	}
      }
      debug12a(printf("\n"));
    }

  } else {
    /* Do not jump late */
    penalty = open + extend;
    for (c = 1; c <= glength; c++) {
      if (gsequence) {
	na2 = revp ? gsequence[1-c] : gsequence[c-1];
	na2_alt = revp ? gsequence_alt[1-c] : gsequence_alt[c-1];
      } else if (revp == false) {
	na2 = get_genomic_nt(&na2_alt,goffset+c-1,chroffset,chrhigh,watsonp);
      } else {
	na2 = get_genomic_nt(&na2_alt,goffset+1-c,chroffset,chrhigh,watsonp);
      }

      c_gap = NEG_INFINITY_32;
      if (c == 1) {
	rlo = 1;
	prev_nogap = 0; /* was nogap[rlo-1] */
	last_nogap = penalty;
      } else if ((rlo = c - uband) < 1) {
	rlo = 1;
	prev_nogap = penalty;
	penalty += extend;
	last_nogap = penalty;
      } else if (rlo == 1) {
	prev_nogap = penalty;
	/* penalty += extend; */
	last_nogap = NEG_INFINITY_32;
#ifdef DEBUG12
	matrix3[c][rlo-1].Fgap = NEG_INFINITY_32;
	matrix3[c][rlo-1].nogap = NEG_INFINITY_32;
#endif
      } else {
	prev_nogap = first_nogap;
	last_nogap = NEG_INFINITY_32;
#ifdef DEBUG12
	matrix3[c][rlo-1].Fgap = NEG_INFINITY_32;
	matrix3[c][rlo-1].nogap = NEG_INFINITY_32;
#endif
      }

      if ((rhigh = c + lband) > rlength) {
	rhigh = rlength;
#ifdef DEBUG12
      } else {
	matrix3[c-1][rhigh].Egap = NEG_INFINITY_32;
	matrix3[c-1][rhigh].nogap = NEG_INFINITY_32;
#endif
      }

      for (r = rlo; r <= rhigh; r++) {
	na1 = revp ? rsequence[1-r] : rsequence[r-1];

#ifdef DEBUG12
	/* FGAP */
	bestscore = matrix3[c][r-1].nogap + open /* + extend */;
	bestdir = DIAG;

	if ((score = matrix3[c][r-1].Fgap /* + extend */) > bestscore) {  /* Use > for jump early */
	  bestscore = score;
	  bestdir = VERT;
	}

	matrix3[c][r].Fgap = bestscore + extend;
#endif

	/* FGAP alt */
#ifdef DEBUG12A
	printf("Fgap at r %d, c %d: matrix3[c][r-1].nogap %d vs last_nogap %d\n",r,c,matrix3[c][r-1].nogap,last_nogap);
	printf("Fgap at r %d, c %d: matrix3[c][r-1].Fgap %d vs c_gap %d\n",r,c,matrix3[c][r-1].Fgap,c_gap);
#endif
#ifdef DEBUG12
	assert(matrix3[c][r-1].nogap == last_nogap);
	assert(matrix3[c][r-1].Fgap == c_gap);
#endif
	/* debug2(printf("Fgap at r %d, c %d: c_gap + extend %d vs last_nogap + open + extend %d\n",r,c,c_gap + extend,last_nogap + open + extend)); */
	if (c_gap /* + extend */ > (score = last_nogap + open /* + extend */)) {  /* Use > for jump early */
	  c_gap += extend;
	  (*directions_Fgap)[c][r] = VERT;
	} else {
	  c_gap = score + extend;
	  /* bestdir2 = DIAG; -- Already initialized to DIAG */
	}


#ifdef DEBUG12
	/* EGAP */
	bestscore = matrix3[c-1][r].nogap + open /* + extend */;
	bestdir = DIAG;

	if ((score = matrix3[c-1][r].Egap /* + extend */) > bestscore) {  /* Use > for jump early */
	  bestscore = score;
	  bestdir = HORIZ;
	}

	matrix3[c][r].Egap = bestscore + extend;
#endif

	/* EGAP alt */
#ifdef DEBUG12A
	printf("Egap at r %d, c %d: matrix3[c-1][r].nogap %d vs nogap[r] %d\n",r,c,matrix3[c-1][r].nogap,nogap[r]);
	printf("Egap at r %d, c %d: matrix3[c-1][r].Egap %d vs r_gap[r] %d\n",r,c,matrix3[c-1][r].Egap,r_gap[r]);
#endif
#ifdef DEBUG12
	assert(matrix3[c-1][r].nogap == nogap[r]);
	assert(matrix3[c-1][r].Egap == r_gap[r]);
#endif
	/* debug2(printf("Egap at r %d, c %d: r_gap[r] %d vs nogap[r] + open %d\n",r,c,r_gap[r],nogap[r]+open)); */
	if (r_gap[r] /* + extend */ > (score = nogap[r] + open /* + extend */)) {  /* Use > for jump early */
	  r_gap[r] += extend;
	  (*directions_Egap)[c][r] = HORIZ;
	} else {
	  r_gap[r] = score + extend;
	  /* bestdir2 = DIAG; -- Already initialized to DIAG */
	}


	/* NOGAP */
	pairscore = pairdistance_array_type[na1][na2];
	if ((score = pairdistance_array_type[na1][(int) na2_alt]) > pairscore) {
	  pairscore = score;
	}
#ifdef DEBUG12
	bestscore = matrix3[c-1][r-1].nogap + pairscore;
	bestdir = DIAG;
      
	if ((score = matrix3[c][r].Egap) > bestscore) {  /* Use > for jump early */
	  bestscore = score;
	  bestdir = HORIZ;
	}

	if ((score = matrix3[c][r].Fgap) > bestscore) {  /* Use > for jump early */
	  bestscore = score;
	  bestdir = VERT;
	}

	matrix3[c][r].nogap = bestscore;
#endif

	/* NOGAP alt */
#ifdef DEBUG12A
	printf("nogap at r %d, c %d: matrix3[c-1][r-1].nogap %d vs prev_nogap %d\n",r,c,matrix3[c-1][r-1].nogap,prev_nogap);
	printf("nogap at r %d, c %d: matrix3[c][r].Fgap %d vs c_gap %d\n",r,c,matrix3[c][r].Fgap,c_gap);
	printf("nogap at r %d, c %d: matrix3[c][r].Egap %d vs r_gap[r] %d\n",r,c,matrix3[c][r].Egap,r_gap[r]);
#endif
#ifdef DEBUG12
	assert(matrix3[c-1][r-1].nogap == prev_nogap);
	assert(matrix3[c][r].Fgap == c_gap);
	assert(matrix3[c][r].Egap == r_gap[r]);
#endif
	last_nogap = prev_nogap + pairscore;
	/* bestdir2 = DIAG; -- Already initialized to DIAG */
	/* debug2(printf("assign nogap at r %d, c %d: H + pairscore %d vs r_horiz + extend %d vs vert + extend %d\n",
	   r,c,last_nogap,r_gap[r],c_gap)); */
	if (r_gap[r] > last_nogap) {  /* Use > for jump early */
	  last_nogap = r_gap[r];
	  (*directions_nogap)[c][r] = HORIZ;
	}
	if (c_gap > last_nogap) {  /* Use > for jump early */
	  last_nogap = c_gap;
	  (*directions_nogap)[c][r] = VERT;
	}
	/* (*directions_nogap)[c][r] = bestdir2; */

	prev_nogap = nogap[r];	/* Save for next inner loop, before we wipe it out */
	matrix[c][r] = nogap[r] = last_nogap;	/* Save for next outer loop */
	if (r == rlo) {
	  debug12a(printf("At row %d, storing first_nogap to be nogap[r] %d\n",r,nogap[r]));
	  first_nogap = last_nogap;
	}
      }
      debug12a(printf("\n"));
    }
  }

  debug2(Matrix32_print(matrix,rlength,glength,rsequence,gsequence,gsequence_alt,
			goffset,chroffset,chrhigh,watsonp,revp));
  debug2(Directions32_print(*directions_nogap,*directions_Egap,*directions_Fgap,
			    rlength,glength,rsequence,gsequence,gsequence_alt,
			    goffset,chroffset,chrhigh,watsonp,revp));
  debug12a(Matrix3_print(matrix3,rlength,glength,rsequence,gsequence,gsequence_alt,
			 goffset,chroffset,chrhigh,watsonp,revp));

  FREE(r_gap);
  FREE(nogap);

  return matrix;
}
#endif


#if defined(HAVE_SSE4_1) || defined(HAVE_SSE2)
static Score8_T **
compute_scores_simd_8 (Direction8_T ***directions_nogap, Direction8_T ***directions_Egap, Direction8_T ***directions_Fgap,
		       T this, char *rsequence, char *gsequence, char *gsequence_alt,
		       int rlength, int glength,
#ifdef DEBUG14
		       int goffset, Univcoord_T chroffset, Univcoord_T chrhigh, bool watsonp,
#endif
		       Mismatchtype_T mismatchtype, Score8_T open, Score8_T extend,
		       int lband, int uband, bool jump_late_p, bool revp) {
  int penalty, rpenalty, c_gap, last_nogap, score;		/* Need to have the ability to go past NEG_INFINITY */
  Score8_T **matrix, *score_column, *score_ptr;
  __m128i pairscores_std, pairscores_alt;
#ifndef HAVE_SSE4_1
  __m128i pairscores_best, all_128;
#endif
  __m128i H_nogap_r, X_prev_nogap, E_r_gap, T1, *EE;
  __m128i v_open, v_extend, all_one_bits, end_neg_infinity;
  __m128i dir_horiz;
  __m128i bottom_masks[17], E_mask_bottom;
  int rlength_ceil, r, c;
  int rlo, rlo_floor, rhigh, rhigh_ceil;
  int max_rhigh_filled = 1;
  int na1, na2;
  char na2_alt;
  Score8_T *pairscores[5], *pairscores_std_ptr, *pairscores_alt_ptr;
  Pairdistance_T **pairdistance_array_type;

#ifdef DEBUG14
  Score32_T **matrix_std;
  Direction32_T **directions_nogap_std, **directions_Egap_std, **directions_Fgap_std;
  char na2_single;
#endif


  rlength_ceil = (int) ((rlength + SIMD_NCHARS - 1)/SIMD_NCHARS) * SIMD_NCHARS;

#ifdef HAVE_SSE4_1
  pairdistance_array_type = pairdistance_array[mismatchtype];
#else
  /* Needed to use _mm_max_epu8 and _mm_min_epu8, instead of signed versions */
  pairdistance_array_type = pairdistance_array_plus_128[mismatchtype];
  all_128 = _mm_set1_epi8(128);
#endif
  
  debug(printf("compute_scores_simd_8: "));
  debug(printf("Lengths are %d and %d, so bands are %d on left and %d on right\n",rlength,glength,lband,uband));
  debug(printf("Query length rounded up to %d\n",rlength_ceil));

  matrix = aligned_score8_alloc(rlength_ceil,glength,
				this->aligned_matrix_ptrs,this->aligned_matrix_space);
  *directions_nogap = aligned_directions8_alloc(rlength_ceil,glength,
						this->aligned_directions_ptrs_0,this->aligned_directions_space_0);
  *directions_Egap = aligned_directions8_alloc(rlength_ceil,glength,
					       this->aligned_directions_ptrs_1,this->aligned_directions_space_1);
  *directions_Fgap = aligned_directions8_calloc(rlength_ceil,glength,
						this->aligned_directions_ptrs_2,this->aligned_directions_space_2);

  /* Row 0 initialization */
  /* penalty = open; */
  for (c = 1; c <= uband && c <= glength; c++) {
    /* penalty += extend; */
    (*directions_Egap)[c][0] = HORIZ;
    (*directions_nogap)[c][0] = HORIZ;
  }
#if 1
  /* Already initialized to DIAG.  Actually no longer initializing directions_Egap */
  (*directions_Egap)[1][0] = DIAG; /* previously used STOP */
  (*directions_nogap)[0][0] = DIAG; /* previously used STOP */
#endif

  /* Column 0 initialization */
  /* penalty = open; */
  for (r = 1; r <= lband && r <= rlength; r++) {
    /* penalty += extend; */
    (*directions_Fgap)[0][r] = VERT;
    (*directions_nogap)[0][r] = VERT;
  }
#if 0
  /* Already initialized to DIAG */
  (*directions_Fgap)[0][1] = DIAG; /* previously used STOP */
#endif


  /* Load pairscores.  Store match - mismatch */
  pairscores[0] = (Score8_T *) _mm_malloc(rlength_ceil * sizeof(Score8_T),16);
  pairscores[1] = (Score8_T *) _mm_malloc(rlength_ceil * sizeof(Score8_T),16);
  pairscores[2] = (Score8_T *) _mm_malloc(rlength_ceil * sizeof(Score8_T),16);
  pairscores[3] = (Score8_T *) _mm_malloc(rlength_ceil * sizeof(Score8_T),16);
  pairscores[4] = (Score8_T *) _mm_malloc(rlength_ceil * sizeof(Score8_T),16);

#if 0
  /* Should not be necessary */
  memset((void *) pairscores[0],0,rlength_ceil*sizeof(Score8_T));
  memset((void *) pairscores[1],0,rlength_ceil*sizeof(Score8_T));
  memset((void *) pairscores[2],0,rlength_ceil*sizeof(Score8_T));
  memset((void *) pairscores[3],0,rlength_ceil*sizeof(Score8_T));
  memset((void *) pairscores[4],0,rlength_ceil*sizeof(Score8_T));
#endif

  /* For non-SSE4.1, addition of 128 taken care of by using pairdistance_array_plus_128 above */
  if (revp == false) {
    for (r = 0; r < rlength; r++) {
      na1 = (int) rsequence[r];
      pairscores[0][r] = (Score8_T) pairdistance_array_type[na1][(int) 'A'];
      pairscores[1][r] = (Score8_T) pairdistance_array_type[na1][(int) 'C'];
      pairscores[2][r] = (Score8_T) pairdistance_array_type[na1][(int) 'G'];
      pairscores[3][r] = (Score8_T) pairdistance_array_type[na1][(int) 'T'];
      pairscores[4][r] = (Score8_T) pairdistance_array_type[na1][(int) 'N'];
    }
  } else {
    for (r = 0; r < rlength; r++) {
      na1 = (int) rsequence[-r];
      pairscores[0][r] = (Score8_T) pairdistance_array_type[na1][(int) 'A'];
      pairscores[1][r] = (Score8_T) pairdistance_array_type[na1][(int) 'C'];
      pairscores[2][r] = (Score8_T) pairdistance_array_type[na1][(int) 'G'];
      pairscores[3][r] = (Score8_T) pairdistance_array_type[na1][(int) 'T'];
      pairscores[4][r] = (Score8_T) pairdistance_array_type[na1][(int) 'N'];
    }
  }

  memset((void *) &(pairscores[0][r]),0,(rlength_ceil-r)*sizeof(Score8_T));
  memset((void *) &(pairscores[1][r]),0,(rlength_ceil-r)*sizeof(Score8_T));
  memset((void *) &(pairscores[2][r]),0,(rlength_ceil-r)*sizeof(Score8_T));
  memset((void *) &(pairscores[3][r]),0,(rlength_ceil-r)*sizeof(Score8_T));
  memset((void *) &(pairscores[4][r]),0,(rlength_ceil-r)*sizeof(Score8_T));



  all_one_bits = _mm_set1_epi8(-1);

  end_neg_infinity = _mm_set1_epi8(NEG_INFINITY_8);
  end_neg_infinity = _mm_slli_si128(end_neg_infinity,LAST_CHAR);
  bottom_masks[0] = _mm_set1_epi8(MAX_CHAR);
  bottom_masks[1] = _mm_or_si128(_mm_srli_si128(bottom_masks[0],ONE_CHAR),end_neg_infinity);
  bottom_masks[2] = _mm_or_si128(_mm_srli_si128(bottom_masks[1],ONE_CHAR),end_neg_infinity);
  bottom_masks[3] = _mm_or_si128(_mm_srli_si128(bottom_masks[2],ONE_CHAR),end_neg_infinity);
  bottom_masks[4] = _mm_or_si128(_mm_srli_si128(bottom_masks[3],ONE_CHAR),end_neg_infinity);
  bottom_masks[5] = _mm_or_si128(_mm_srli_si128(bottom_masks[4],ONE_CHAR),end_neg_infinity);
  bottom_masks[6] = _mm_or_si128(_mm_srli_si128(bottom_masks[5],ONE_CHAR),end_neg_infinity);
  bottom_masks[7] = _mm_or_si128(_mm_srli_si128(bottom_masks[6],ONE_CHAR),end_neg_infinity);
  bottom_masks[8] = _mm_or_si128(_mm_srli_si128(bottom_masks[7],ONE_CHAR),end_neg_infinity);
  bottom_masks[9] = _mm_or_si128(_mm_srli_si128(bottom_masks[8],ONE_CHAR),end_neg_infinity);
  bottom_masks[10] = _mm_or_si128(_mm_srli_si128(bottom_masks[9],ONE_CHAR),end_neg_infinity);
  bottom_masks[11] = _mm_or_si128(_mm_srli_si128(bottom_masks[10],ONE_CHAR),end_neg_infinity);
  bottom_masks[12] = _mm_or_si128(_mm_srli_si128(bottom_masks[11],ONE_CHAR),end_neg_infinity);
  bottom_masks[13] = _mm_or_si128(_mm_srli_si128(bottom_masks[12],ONE_CHAR),end_neg_infinity);
  bottom_masks[14] = _mm_or_si128(_mm_srli_si128(bottom_masks[13],ONE_CHAR),end_neg_infinity);
  bottom_masks[15] = _mm_or_si128(_mm_srli_si128(bottom_masks[14],ONE_CHAR),end_neg_infinity);
  bottom_masks[16] = _mm_or_si128(_mm_srli_si128(bottom_masks[15],ONE_CHAR),end_neg_infinity);
#ifndef HAVE_SSE4_1
  for (r = 0; r <= 16; r++) {
    bottom_masks[r] = _mm_add_epi8(bottom_masks[r], all_128);
  }
#endif


  EE = (__m128i *) _mm_malloc(rlength_ceil * sizeof(Score8_T),16);
  for (r = 0; r < rlength_ceil/SIMD_NCHARS; r++) {
    _mm_store_si128(&(EE[r]),_mm_set1_epi8(NEG_INFINITY_8));
  }
    
  v_open = _mm_set1_epi8(open);
  v_extend = _mm_set1_epi8(extend);

  if (jump_late_p) {
    penalty = rpenalty = open + extend;
    for (c = 1; c <= glength; c++) {
      score_column = matrix[c];
      na2 = revp ? gsequence[1-c] : gsequence[c-1];
      na2_alt = revp ? gsequence_alt[1-c] : gsequence_alt[c-1];
#ifdef DEBUG14
      if (revp == false) {
	na2_single = get_genomic_nt(&na2_alt,goffset+c-1,chroffset,chrhigh,watsonp);
      } else {
	na2_single = get_genomic_nt(&na2_alt,goffset+1-c,chroffset,chrhigh,watsonp);
      }
      if (na2 != na2_single) {
	abort();
      }
#endif

      switch (na2) {
      case 'A': pairscores_std_ptr = pairscores[0]; break;
      case 'C': pairscores_std_ptr = pairscores[1]; break;
      case 'G': pairscores_std_ptr = pairscores[2]; break;
      case 'T': pairscores_std_ptr = pairscores[3]; break;
      default: pairscores_std_ptr = pairscores[4];
      }
      switch (na2_alt) {
      case 'A': pairscores_alt_ptr = pairscores[0]; break;
      case 'C': pairscores_alt_ptr = pairscores[1]; break;
      case 'G': pairscores_alt_ptr = pairscores[2]; break;
      case 'T': pairscores_alt_ptr = pairscores[3]; break;
      default: pairscores_alt_ptr = pairscores[4];
      }

      if (c == 1) {
	rlo = 1;
	X_prev_nogap = _mm_set1_epi8(0);
	last_nogap = penalty;

      } else if ((rlo = c - uband) < 1) {
	rlo = 1;
	if (penalty < NEG_INFINITY_8) {
	  X_prev_nogap = _mm_set1_epi8(NEG_INFINITY_8);
	} else {
	  X_prev_nogap = _mm_set1_epi8(penalty);
	}
	X_prev_nogap = _mm_srli_si128(X_prev_nogap,LAST_CHAR);
	penalty += extend;
	last_nogap = penalty;
      
      } else if (rlo == 1) {
	if (penalty < NEG_INFINITY_8) {
	  X_prev_nogap = _mm_set1_epi8(NEG_INFINITY_8);
	} else {
	  X_prev_nogap = _mm_set1_epi8(penalty);
	}
	X_prev_nogap = _mm_srli_si128(X_prev_nogap,LAST_CHAR);
	/* penalty += extend; */
	last_nogap = NEG_INFINITY_8;

      } else if ((rlo_floor = (int)((rlo - 1)/SIMD_NCHARS) * SIMD_NCHARS + 1) == 1) {
	/* first block of 8 */
	X_prev_nogap = _mm_set1_epi8(NEG_INFINITY_8); /* works if we start outside the rlo bounds */
	X_prev_nogap = _mm_srli_si128(X_prev_nogap,LAST_CHAR);
	last_nogap = NEG_INFINITY_8;
	rlo = rlo_floor;

      } else {
	/* second or greater block of 8 */
	X_prev_nogap = _mm_set1_epi8(matrix[c-1][rlo-1]); /* get H from previous block and previous column */
	X_prev_nogap = _mm_srli_si128(X_prev_nogap,LAST_CHAR);
	last_nogap = NEG_INFINITY_8;
	rlo = rlo_floor;
      }

      if ((rhigh = c + lband) > rlength) {
	rhigh = rlength;
	rhigh_ceil = (int) ((rhigh + SIMD_NCHARS - 1)/SIMD_NCHARS) * SIMD_NCHARS;
	E_mask_bottom = bottom_masks[0];
      } else {
	rhigh_ceil = (int) ((rhigh + SIMD_NCHARS - 1)/SIMD_NCHARS) * SIMD_NCHARS;
	E_mask_bottom = bottom_masks[rhigh_ceil-rhigh+1];
      }

      if (rhigh_ceil > max_rhigh_filled) {
	for (r = max_rhigh_filled; r <= rhigh_ceil; r++) {
	  if (r > lband) {
	    matrix[c-1][r] = NEG_INFINITY_8;
	  } else if (rpenalty < NEG_INFINITY_8) {
	    matrix[c-1][r] = NEG_INFINITY_8;
	    rpenalty += extend;
	  } else {
	    matrix[c-1][r] = rpenalty;
	    rpenalty += extend;
	  }
	}
	max_rhigh_filled = r;
      }


      for (r = rlo; r <= rhigh; r += SIMD_NCHARS) {
	/* Load previous E vector at this point in the query sequence */
	/* H vector already loaded before loop or at bottom of loop */
	E_r_gap = _mm_load_si128(&(EE[(r-1)/SIMD_NCHARS]));
	H_nogap_r = _mm_load_si128((__m128i *) &(matrix[c-1][r]));

	/* EGAP */
	T1 = _mm_adds_epi8(H_nogap_r, v_open);
	dir_horiz = _mm_cmplt_epi8(E_r_gap,T1); /* E < H */
	dir_horiz = _mm_andnot_si128(dir_horiz,all_one_bits);	/* E >= H, for jump late */
	_mm_store_si128((__m128i *) &((*directions_Egap)[c][r]),dir_horiz);

#ifdef HAVE_SSE4_1
	E_r_gap = _mm_max_epi8(E_r_gap, T1); /* Compare H + open with vert */
#else
	E_r_gap = _mm_sub_epi8(_mm_max_epu8(_mm_add_epi8(E_r_gap, all_128), _mm_add_epi8(T1, all_128)), all_128);
#endif
	E_r_gap = _mm_adds_epi8(E_r_gap, v_extend); /* Compute scores for Egap (vert + open) */
	/* print_vector_8(E_r_gap,r,c,"E"); */
	if (r + SIMD_NCHARS > rhigh) {
	  /* mask last block (e.g., NM_001193336_71_) */
#ifdef HAVE_SSE4_1
	  E_r_gap = _mm_min_epi8(E_r_gap, E_mask_bottom); /* To handle band limitation from below */
#else
	  /* To handle band limitation from below */
	  /* E_mask_bottom already has 128 added */
	  E_r_gap = _mm_sub_epi8(_mm_min_epu8(_mm_add_epi8(E_r_gap, all_128), E_mask_bottom), all_128);
#endif
	}
	_mm_store_si128(&(EE[(r-1)/SIMD_NCHARS]), E_r_gap);


	/* NOGAP */
	T1 = _mm_srli_si128(H_nogap_r,LAST_CHAR);
	H_nogap_r = _mm_slli_si128(H_nogap_r,ONE_CHAR);
	H_nogap_r = _mm_or_si128(H_nogap_r, X_prev_nogap);
	X_prev_nogap = T1;

	/* Add pairscores, allowing for alternate genomic nt */
#ifdef HAVE_SSE4_1
	pairscores_std = _mm_load_si128((__m128i *) &(pairscores_std_ptr[r-1]));
	pairscores_alt = _mm_load_si128((__m128i *) &(pairscores_alt_ptr[r-1]));
	H_nogap_r = _mm_adds_epi8(H_nogap_r, _mm_max_epi8(pairscores_std,pairscores_alt));
#else
	pairscores_std = _mm_load_si128((__m128i *) &(pairscores_std_ptr[r-1])); /* Has 128 added already */
	pairscores_alt = _mm_load_si128((__m128i *) &(pairscores_alt_ptr[r-1])); /* Has 128 added already */
	pairscores_best = _mm_sub_epi8(_mm_max_epu8(pairscores_std, pairscores_alt), all_128);
	H_nogap_r = _mm_adds_epi8(H_nogap_r, pairscores_best);
#endif
	_mm_clflush(&H_nogap_r); /* Needed for opencc -O3 on AMD */
	/* print_vector_8(H_nogap_r,r,c,"H"); */

	dir_horiz = _mm_cmplt_epi8(E_r_gap,H_nogap_r); /* E < H */
	dir_horiz = _mm_andnot_si128(dir_horiz,all_one_bits);	/* E >= H, for jump late */
	_mm_store_si128((__m128i *) &((*directions_nogap)[c][r]),dir_horiz);

#ifdef HAVE_SSE4_1
	H_nogap_r = _mm_max_epi8(H_nogap_r, E_r_gap); /* Compare H + pairscores with horiz + extend */
#else
	/* Compare H + pairscores with horiz + extend */
	H_nogap_r = _mm_sub_epi8(_mm_max_epu8(_mm_add_epi8(H_nogap_r, all_128), _mm_add_epi8(E_r_gap, all_128)), all_128);
#endif
	_mm_store_si128((__m128i *) &(score_column[r]), H_nogap_r);
      }

      /* Perform F loop on entire column */
      if ((rlo = c - uband) < 1) {
	rlo = 1;
      }
      /* Use same rhigh as computed previously */

      c_gap = NEG_INFINITY_8;
      score_ptr = &(score_column[rlo]);
      for (r = rlo; r <= rhigh; r++) {

	/* FGAP */
	/* debug2(printf("Fgap at r %d, c %d: c_gap + extend %d vs last_nogap + open + extend %d\n",r,c,c_gap + extend,last_nogap + open + extend)); */
	if (c_gap /* + extend */ >= (score = last_nogap + open /* + extend */)) {  /* Use >= for jump late */
	  c_gap += extend;
	  (*directions_Fgap)[c][r] = VERT;
	} else {
	  c_gap = score + extend;
	  /* (*directions_Fgap)[c][r] = DIAG: -- Already initialized to DIAG */
	}

	/* NOGAP */
	last_nogap = *score_ptr;
	/* debug2(printf("assign nogap at r %d, c %d: H/E %d vs vert + extend %d\n",r,c,last_nogap,c_gap)); */
	if (c_gap >= last_nogap) {  /* Use >= for jump late */
	  last_nogap = c_gap;
	  *score_ptr = (c_gap < NEG_INFINITY_8) ? NEG_INFINITY_8 : (Score8_T) c_gap; /* Saturation */
	  (*directions_nogap)[c][r] = VERT;
	}
	score_ptr++;
      }
    }

  } else {
    /* jump early */
    penalty = rpenalty = open + extend;
    for (c = 1; c <= glength; c++) {
      score_column = matrix[c];
      na2 = revp ? gsequence[1-c] : gsequence[c-1];
      na2_alt = revp ? gsequence_alt[1-c] : gsequence_alt[c-1];
#ifdef DEBUG14
      if (revp == false) {
	na2_single = get_genomic_nt(&na2_alt,goffset+c-1,chroffset,chrhigh,watsonp);
      } else {
	na2_single = get_genomic_nt(&na2_alt,goffset+1-c,chroffset,chrhigh,watsonp);
      }
      if (na2 != na2_single) {
	abort();
      }
#endif

      switch (na2) {
      case 'A': pairscores_std_ptr = pairscores[0]; break;
      case 'C': pairscores_std_ptr = pairscores[1]; break;
      case 'G': pairscores_std_ptr = pairscores[2]; break;
      case 'T': pairscores_std_ptr = pairscores[3]; break;
      default: pairscores_std_ptr = pairscores[4];
      }
      switch (na2_alt) {
      case 'A': pairscores_alt_ptr = pairscores[0]; break;
      case 'C': pairscores_alt_ptr = pairscores[1]; break;
      case 'G': pairscores_alt_ptr = pairscores[2]; break;
      case 'T': pairscores_alt_ptr = pairscores[3]; break;
      default: pairscores_alt_ptr = pairscores[4];
      }

      if (c == 1) {
	rlo = 1;
	X_prev_nogap = _mm_set1_epi8(0);
	last_nogap = penalty;

      } else if ((rlo = c - uband) < 1) {
	rlo = 1;
	if (penalty < NEG_INFINITY_8) {
	  X_prev_nogap = _mm_set1_epi8(NEG_INFINITY_8);
	} else {
	  X_prev_nogap = _mm_set1_epi8(penalty);
	}
	X_prev_nogap = _mm_srli_si128(X_prev_nogap,LAST_CHAR);
	penalty += extend;
	last_nogap = penalty;
      
      } else if (rlo == 1) {
	if (penalty < NEG_INFINITY_8) {
	  X_prev_nogap = _mm_set1_epi8(NEG_INFINITY_8);
	} else {
	  X_prev_nogap = _mm_set1_epi8(penalty);
	}
	X_prev_nogap = _mm_srli_si128(X_prev_nogap,LAST_CHAR);
	/* penalty += extend; */
	last_nogap = NEG_INFINITY_8;

      } else if ((rlo_floor = (int)((rlo - 1)/SIMD_NCHARS) * SIMD_NCHARS + 1) == 1) {
	/* first block of 8 */
	X_prev_nogap = _mm_set1_epi8(NEG_INFINITY_8); /* works if we start outside the rlo bounds */
	X_prev_nogap = _mm_srli_si128(X_prev_nogap,LAST_CHAR);
	last_nogap = NEG_INFINITY_8;
	rlo = rlo_floor;

      } else {
	/* second or greater block of 8 */
	X_prev_nogap = _mm_set1_epi8(matrix[c-1][rlo-1]); /* get H from previous block and previous column */
	X_prev_nogap = _mm_srli_si128(X_prev_nogap,LAST_CHAR);
	last_nogap = NEG_INFINITY_8;
	rlo = rlo_floor;
      }

      if ((rhigh = c + lband) > rlength) {
	rhigh = rlength;
	rhigh_ceil = (int) ((rhigh + SIMD_NCHARS - 1)/SIMD_NCHARS) * SIMD_NCHARS;
	E_mask_bottom = bottom_masks[0];
      } else {
	rhigh_ceil = (int) ((rhigh + SIMD_NCHARS - 1)/SIMD_NCHARS) * SIMD_NCHARS;
	E_mask_bottom = bottom_masks[rhigh_ceil-rhigh+1];
      }

      if (rhigh_ceil > max_rhigh_filled) {
	for (r = max_rhigh_filled; r <= rhigh_ceil; r++) {
	  if (r > lband) {
	    matrix[c-1][r] = NEG_INFINITY_8;
	  } else if (rpenalty < NEG_INFINITY_8) {
	    matrix[c-1][r] = NEG_INFINITY_8;
	    rpenalty += extend;
	  } else {
	    matrix[c-1][r] = rpenalty;
	    rpenalty += extend;
	  }
	}
	max_rhigh_filled = r;
      }


      for (r = rlo; r <= rhigh; r += SIMD_NCHARS) {
	/* Load previous E vector at this point in the query sequence */
	/* H vector already loaded before loop or at bottom of loop */
	E_r_gap = _mm_load_si128(&(EE[(r-1)/SIMD_NCHARS]));
	H_nogap_r = _mm_load_si128((__m128i *) &(matrix[c-1][r]));

	/* EGAP */
	T1 = _mm_adds_epi8(H_nogap_r, v_open);
	dir_horiz = _mm_cmpgt_epi8(E_r_gap,T1); /* E > H, for jump early */
	_mm_store_si128((__m128i *) &((*directions_Egap)[c][r]),dir_horiz);

#ifdef HAVE_SSE4_1
	E_r_gap = _mm_max_epi8(E_r_gap, T1); /* Compare H + open with vert */
#else
	/* Compare H + open with vert */
	E_r_gap = _mm_sub_epi8(_mm_max_epu8(_mm_add_epi8(E_r_gap, all_128), _mm_add_epi8(T1, all_128)), all_128);

#endif
	E_r_gap = _mm_adds_epi8(E_r_gap, v_extend); /* Compute scores for Egap (vert + open) */
	/* print_vector_8(E_r_gap,r,c,"E"); */
	if (r + SIMD_NCHARS > rhigh) {
	  /* mask last block (e.g., NM_001193336_71_) */
#ifdef HAVE_SSE4_1
	  E_r_gap = _mm_min_epi8(E_r_gap, E_mask_bottom); /* To handle band limitation from below */
#else
	  /* To handle band limitation from below */
	  /* E_mask_bottom already has 128 added */
	  E_r_gap = _mm_sub_epi8(_mm_min_epu8(_mm_add_epi8(E_r_gap, all_128), E_mask_bottom), all_128);
#endif
	}
	_mm_store_si128(&(EE[(r-1)/SIMD_NCHARS]), E_r_gap);


	/* NOGAP */
	T1 = _mm_srli_si128(H_nogap_r,LAST_CHAR);
	H_nogap_r = _mm_slli_si128(H_nogap_r,ONE_CHAR);
	H_nogap_r = _mm_or_si128(H_nogap_r, X_prev_nogap);
	X_prev_nogap = T1;

	/* Add pairscores, allowing for alternate genomic nt */
#ifdef HAVE_SSE4_1
	pairscores_std = _mm_load_si128((__m128i *) &(pairscores_std_ptr[r-1]));
	pairscores_alt = _mm_load_si128((__m128i *) &(pairscores_alt_ptr[r-1]));
	H_nogap_r = _mm_adds_epi8(H_nogap_r, _mm_max_epi8(pairscores_std,pairscores_alt));
#else
	pairscores_std = _mm_load_si128((__m128i *) &(pairscores_std_ptr[r-1])); /* Has 128 added already */
	pairscores_alt = _mm_load_si128((__m128i *) &(pairscores_alt_ptr[r-1])); /* Has 128 added already */
	pairscores_best = _mm_sub_epi8(_mm_max_epu8(pairscores_std, pairscores_alt), all_128);
	H_nogap_r = _mm_adds_epi8(H_nogap_r, pairscores_best);
#endif
	_mm_clflush(&H_nogap_r); /* Needed for opencc -O3 on AMD */
	/* print_vector_8(H_nogap_r,r,c,"H"); */

	dir_horiz = _mm_cmpgt_epi8(E_r_gap,H_nogap_r); /* E > H, for jump early */
	_mm_store_si128((__m128i *) &((*directions_nogap)[c][r]),dir_horiz);

#ifdef HAVE_SSE4_1
	H_nogap_r = _mm_max_epi8(H_nogap_r, E_r_gap); /* Compare H + pairscores with horiz + extend */
#else
	/* Compare H + pairscores with horiz + extend */
	H_nogap_r = _mm_sub_epi8(_mm_max_epu8(_mm_add_epi8(H_nogap_r, all_128), _mm_add_epi8(E_r_gap, all_128)), all_128);
#endif
	_mm_store_si128((__m128i *) &(score_column[r]), H_nogap_r);
      }

      /* Perform F loop on entire column */
      if ((rlo = c - uband) < 1) {
	rlo = 1;
      }
      /* Use same rhigh as computed previously */

      c_gap = NEG_INFINITY_8;
      score_ptr = &(score_column[rlo]);
      for (r = rlo; r <= rhigh; r++) {

	/* FGAP */
	/* debug2(printf("Fgap at r %d, c %d: c_gap + extend %d vs last_nogap + open + extend %d\n",r,c,c_gap + extend,last_nogap + open + extend)); */
	if (c_gap /* + extend */ > (score = last_nogap + open /* + extend */)) {  /* Use > for jump early */
	  c_gap += extend;
	  (*directions_Fgap)[c][r] = VERT;
	} else {
	  c_gap = score + extend;
	  /* (*directions_Fgap)[c][r] = DIAG: -- Already initialized to DIAG */
	}

	/* NOGAP */
	last_nogap = *score_ptr;
	/* debug2(printf("assign nogap at r %d, c %d: H/E %d vs vert + extend %d\n",r,c,last_nogap,c_gap)); */
	if (c_gap > last_nogap) {  /* Use > for jump early */
	  last_nogap = c_gap;
	  *score_ptr = (c_gap < NEG_INFINITY_8) ? NEG_INFINITY_8 : (Score8_T) c_gap; /* Saturation */
	  (*directions_nogap)[c][r] = VERT;
	}
	score_ptr++;
      }
    }
  }

#ifdef DEBUG2
  printf("SIMD\n");
  Matrix8_print(matrix,rlength,glength,rsequence,gsequence,gsequence_alt,
		goffset,chroffset,chrhigh,watsonp,revp);
  Directions8_print(*directions_nogap,*directions_Egap,*directions_Fgap,
		    rlength,glength,rsequence,gsequence,gsequence_alt,
		    goffset,chroffset,chrhigh,watsonp,revp);
#endif
  
#ifdef DEBUG14
  matrix_std = compute_scores_standard(&directions_nogap_std,&directions_Egap_std,&directions_Fgap_std,
				       this,rsequence,/*gsequence (NULL for debugging)*/NULL,/*gsequence_alt*/NULL,
				       rlength,glength,goffset,chroffset,chrhigh,watsonp,mismatchtype,
				       open,extend,lband,uband,jump_late_p,revp);

#ifdef DEBUG2
  printf("Banded %s\n",revp ? "rev" : "fwd");
  Matrix32_print(matrix_std,rlength,glength,rsequence,gsequence,gsequence_alt,
		 goffset,chroffset,chrhigh,watsonp,revp);
  Directions32_print(directions_nogap_std,directions_Egap_std,directions_Fgap_std,
		     rlength,glength,rsequence,gsequence,gsequence_alt,
		     goffset,chroffset,chrhigh,watsonp,revp);
#endif
  
  banded_matrix8_compare(matrix,matrix_std,rlength,glength,lband,uband,
			 rsequence,gsequence,gsequence_alt,
			 goffset,chroffset,chrhigh,watsonp,revp);

  banded_directions8_compare_nogap(matrix,*directions_nogap,directions_nogap_std,rlength,glength,lband,uband);

  banded_directions8_compare_Egap(matrix,*directions_Egap,directions_Egap_std,rlength,glength,lband,uband);

  banded_directions8_compare_Fgap(matrix,*directions_Fgap,directions_Fgap_std,rlength,glength,lband,uband);
#endif

  _mm_free(EE);
  _mm_free(pairscores[4]);
  _mm_free(pairscores[3]);
  _mm_free(pairscores[2]);
  _mm_free(pairscores[1]);
  _mm_free(pairscores[0]);

  return matrix;
}
#endif


#ifdef HAVE_SSE2
static Score16_T **
compute_scores_simd_16 (Direction16_T ***directions_nogap, Direction16_T ***directions_Egap, Direction16_T ***directions_Fgap,
			T this, char *rsequence, char *gsequence, char *gsequence_alt,
			int rlength, int glength,
#ifdef DEBUG14
			int goffset, Univcoord_T chroffset, Univcoord_T chrhigh, bool watsonp,
#endif
			Mismatchtype_T mismatchtype, Score16_T open, Score16_T extend,
			int lband, int uband, bool jump_late_p, bool revp) {
  int penalty, rpenalty, c_gap, last_nogap, score;		/* Need to have the ability to go past NEG_INFINITY */
  Score16_T **matrix, *score_column, *score_ptr;
  __m128i pairscores_std, pairscores_alt;
  __m128i H_nogap_r, X_prev_nogap, E_r_gap, T1, *EE;
  __m128i v_open, v_extend, all_one_bits, end_neg_infinity;
  __m128i dir_horiz;
  __m128i bottom_masks[9], E_mask_bottom;
  int rlength_ceil, r, c;
  int rlo, rlo_floor, rhigh, rhigh_ceil;
  int max_rhigh_filled = 1;
  int na1, na2;
  char na2_alt;
  Score16_T *pairscores[5], *pairscores_std_ptr, *pairscores_alt_ptr;
  Pairdistance_T **pairdistance_array_type;

#ifdef DEBUG14
  Score32_T **matrix_std;
  Direction32_T **directions_nogap_std, **directions_Egap_std, **directions_Fgap_std;
  char na2_single;
#endif

  rlength_ceil = (int) ((rlength + SIMD_NSHORTS - 1)/SIMD_NSHORTS) * SIMD_NSHORTS;
  pairdistance_array_type = pairdistance_array[mismatchtype];
  
  debug(printf("compute_scores_simd_16: "));
  debug(printf("Lengths are %d and %d, so bands are %d on left and %d on right\n",rlength,glength,lband,uband));
  debug(printf("Query length rounded up to %d\n",rlength_ceil));

  matrix = aligned_score16_alloc(rlength_ceil,glength,
				 this->aligned_matrix_ptrs,this->aligned_matrix_space);
  *directions_nogap = aligned_directions16_alloc(rlength_ceil,glength,
						 this->aligned_directions_ptrs_0,this->aligned_directions_space_0);
  *directions_Egap = aligned_directions16_alloc(rlength_ceil,glength,
						this->aligned_directions_ptrs_1,this->aligned_directions_space_1);
  *directions_Fgap = aligned_directions16_calloc(rlength_ceil,glength,
						 this->aligned_directions_ptrs_2,this->aligned_directions_space_2);

  /* Row 0 initialization */
  /* penalty = open; */
  for (c = 1; c <= uband && c <= glength; c++) {
    /* penalty += extend; */
    (*directions_Egap)[c][0] = HORIZ;
    (*directions_nogap)[c][0] = HORIZ;
  }
#if 1
  /* Already initialized to DIAG.  Actually, no longer initializing directions_Egap */
  (*directions_Egap)[1][0] = DIAG; /* previously used STOP */
  (*directions_nogap)[0][0] = DIAG; /* previously used STOP */
#endif

  /* Column 0 initialization */
  /* penalty = open; */
  for (r = 1; r <= lband && r <= rlength; r++) {
    /* penalty += extend; */
    (*directions_Fgap)[0][r] = VERT;
    (*directions_nogap)[0][r] = VERT;
  }
#if 0
  /* Already initialized to DIAG */
  (*directions_Fgap)[0][1] = DIAG; /* previously used STOP */
#endif


  /* Load pairscores.  Store match - mismatch */
  pairscores[0] = (Score16_T *) _mm_malloc(rlength_ceil * sizeof(Score16_T),16);
  pairscores[1] = (Score16_T *) _mm_malloc(rlength_ceil * sizeof(Score16_T),16);
  pairscores[2] = (Score16_T *) _mm_malloc(rlength_ceil * sizeof(Score16_T),16);
  pairscores[3] = (Score16_T *) _mm_malloc(rlength_ceil * sizeof(Score16_T),16);
  pairscores[4] = (Score16_T *) _mm_malloc(rlength_ceil * sizeof(Score16_T),16);

#if 0
  /* Should not be necessary */
  memset((void *) pairscores[0],0,rlength_ceil*sizeof(Score16_T));
  memset((void *) pairscores[1],0,rlength_ceil*sizeof(Score16_T));
  memset((void *) pairscores[2],0,rlength_ceil*sizeof(Score16_T));
  memset((void *) pairscores[3],0,rlength_ceil*sizeof(Score16_T));
  memset((void *) pairscores[4],0,rlength_ceil*sizeof(Score16_T));
#endif

  if (revp == false) {
    for (r = 0; r < rlength; r++) {
      na1 = (int) rsequence[r];
      pairscores[0][r] = (Score16_T) pairdistance_array_type[na1][(int) 'A'];
      pairscores[1][r] = (Score16_T) pairdistance_array_type[na1][(int) 'C'];
      pairscores[2][r] = (Score16_T) pairdistance_array_type[na1][(int) 'G'];
      pairscores[3][r] = (Score16_T) pairdistance_array_type[na1][(int) 'T'];
      pairscores[4][r] = (Score16_T) pairdistance_array_type[na1][(int) 'N'];
    }
  } else {
    for (r = 0; r < rlength; r++) {
      na1 = (int) rsequence[-r];
      pairscores[0][r] = (Score16_T) pairdistance_array_type[na1][(int) 'A'];
      pairscores[1][r] = (Score16_T) pairdistance_array_type[na1][(int) 'C'];
      pairscores[2][r] = (Score16_T) pairdistance_array_type[na1][(int) 'G'];
      pairscores[3][r] = (Score16_T) pairdistance_array_type[na1][(int) 'T'];
      pairscores[4][r] = (Score16_T) pairdistance_array_type[na1][(int) 'N'];
    }
  }

  memset((void *) &(pairscores[0][r]),0,(rlength_ceil-r)*sizeof(Score16_T));
  memset((void *) &(pairscores[1][r]),0,(rlength_ceil-r)*sizeof(Score16_T));
  memset((void *) &(pairscores[2][r]),0,(rlength_ceil-r)*sizeof(Score16_T));
  memset((void *) &(pairscores[3][r]),0,(rlength_ceil-r)*sizeof(Score16_T));
  memset((void *) &(pairscores[4][r]),0,(rlength_ceil-r)*sizeof(Score16_T));


  all_one_bits = _mm_set1_epi16(-1);

  end_neg_infinity = _mm_set1_epi16(NEG_INFINITY_16);
  end_neg_infinity = _mm_slli_si128(end_neg_infinity,LAST_SHORT);
  bottom_masks[0] = _mm_set1_epi16(MAX_SHORT);
  bottom_masks[1] = _mm_or_si128(_mm_srli_si128(bottom_masks[0],ONE_SHORT),end_neg_infinity);
  bottom_masks[2] = _mm_or_si128(_mm_srli_si128(bottom_masks[1],ONE_SHORT),end_neg_infinity);
  bottom_masks[3] = _mm_or_si128(_mm_srli_si128(bottom_masks[2],ONE_SHORT),end_neg_infinity);
  bottom_masks[4] = _mm_or_si128(_mm_srli_si128(bottom_masks[3],ONE_SHORT),end_neg_infinity);
  bottom_masks[5] = _mm_or_si128(_mm_srli_si128(bottom_masks[4],ONE_SHORT),end_neg_infinity);
  bottom_masks[6] = _mm_or_si128(_mm_srli_si128(bottom_masks[5],ONE_SHORT),end_neg_infinity);
  bottom_masks[7] = _mm_or_si128(_mm_srli_si128(bottom_masks[6],ONE_SHORT),end_neg_infinity);
  bottom_masks[8] = _mm_or_si128(_mm_srli_si128(bottom_masks[7],ONE_SHORT),end_neg_infinity);


  EE = (__m128i *) _mm_malloc(rlength_ceil * sizeof(Score16_T),16);
  for (r = 0; r < rlength_ceil/SIMD_NSHORTS; r++) {
    _mm_store_si128(&(EE[r]),_mm_set1_epi16(NEG_INFINITY_16));
  }
    
  v_open = _mm_set1_epi16(open);
  v_extend = _mm_set1_epi16(extend);

  if (jump_late_p) {
    penalty = rpenalty = open + extend;
    for (c = 1; c <= glength; c++) {
      score_column = matrix[c];
      na2 = revp ? gsequence[1-c] : gsequence[c-1];
      na2_alt = revp ? gsequence_alt[1-c] : gsequence_alt[c-1];
#ifdef DEBUG14
      if (revp == false) {
	na2_single = get_genomic_nt(&na2_alt,goffset+c-1,chroffset,chrhigh,watsonp);
      } else {
	na2_single = get_genomic_nt(&na2_alt,goffset+1-c,chroffset,chrhigh,watsonp);
      }
      if (na2 != na2_single) {
	abort();
      }
#endif

      switch (na2) {
      case 'A': pairscores_std_ptr = pairscores[0]; break;
      case 'C': pairscores_std_ptr = pairscores[1]; break;
      case 'G': pairscores_std_ptr = pairscores[2]; break;
      case 'T': pairscores_std_ptr = pairscores[3]; break;
      default: pairscores_std_ptr = pairscores[4];
      }
      switch (na2_alt) {
      case 'A': pairscores_alt_ptr = pairscores[0]; break;
      case 'C': pairscores_alt_ptr = pairscores[1]; break;
      case 'G': pairscores_alt_ptr = pairscores[2]; break;
      case 'T': pairscores_alt_ptr = pairscores[3]; break;
      default: pairscores_alt_ptr = pairscores[4];
      }

      if (c == 1) {
	rlo = 1;
	X_prev_nogap = _mm_set1_epi16(0);
	last_nogap = penalty;

      } else if ((rlo = c - uband) < 1) {
	rlo = 1;
	if (penalty < NEG_INFINITY_16) {
	  X_prev_nogap = _mm_set1_epi16(NEG_INFINITY_16); /* Unlikely */
	} else {
	  X_prev_nogap = _mm_set1_epi16(penalty);
	}
	X_prev_nogap = _mm_srli_si128(X_prev_nogap,LAST_SHORT);
	penalty += extend;
	last_nogap = penalty;
      
      } else if (rlo == 1) {
	if (penalty < NEG_INFINITY_16) {
	  X_prev_nogap = _mm_set1_epi16(NEG_INFINITY_16);
	} else {
	  X_prev_nogap = _mm_set1_epi16(penalty);
	}
	X_prev_nogap = _mm_srli_si128(X_prev_nogap,LAST_SHORT);
	/* penalty += extend; */
	last_nogap = NEG_INFINITY_16;

      } else if ((rlo_floor = (int)((rlo - 1)/SIMD_NSHORTS) * SIMD_NSHORTS + 1) == 1) {
	/* first block of 8 */
	X_prev_nogap = _mm_set1_epi16(NEG_INFINITY_16); /* works if we start outside the rlo bounds */
	X_prev_nogap = _mm_srli_si128(X_prev_nogap,LAST_SHORT);
	last_nogap = NEG_INFINITY_16;
	rlo = rlo_floor;

      } else {
	/* second or greater block of 8 */
	X_prev_nogap = _mm_set1_epi16(matrix[c-1][rlo-1]); /* get H from previous block and previous column */
	X_prev_nogap = _mm_srli_si128(X_prev_nogap,LAST_SHORT);
	last_nogap = NEG_INFINITY_16;
	rlo = rlo_floor;
      }

      if ((rhigh = c + lband) > rlength) {
	rhigh = rlength;
	rhigh_ceil = (int) ((rhigh + SIMD_NSHORTS - 1)/SIMD_NSHORTS) * SIMD_NSHORTS;
	E_mask_bottom = bottom_masks[0];
      } else {
	rhigh_ceil = (int) ((rhigh + SIMD_NSHORTS - 1)/SIMD_NSHORTS) * SIMD_NSHORTS;
	E_mask_bottom = bottom_masks[rhigh_ceil-rhigh+1];
      }

      if (rhigh_ceil > max_rhigh_filled) {
	for (r = max_rhigh_filled; r <= rhigh_ceil; r++) {
	  if (r > lband) {
	    matrix[c-1][r] = NEG_INFINITY_16;
	  } else if (rpenalty < NEG_INFINITY_16) {
	    /* Unlikely */
	    matrix[c-1][r] = NEG_INFINITY_16;
	    rpenalty += extend;
	  } else {
	    matrix[c-1][r] = rpenalty;
	    rpenalty += extend;
	  }
	}
	max_rhigh_filled = r;
      }


      for (r = rlo; r <= rhigh; r += SIMD_NSHORTS) {
	/* Load previous E vector at this point in the query sequence */
	/* H vector already loaded before loop or at bottom of loop */
	E_r_gap = _mm_load_si128(&(EE[(r-1)/SIMD_NSHORTS]));
	debug15(print_vector_16(E_r_gap,r,c,"E_r_gap"));
	H_nogap_r = _mm_load_si128((__m128i *) &(matrix[c-1][r]));
	debug15(print_vector_16(H_nogap_r,r,c,"H_nogap_r load"));

	/* EGAP */
	T1 = _mm_adds_epi16(H_nogap_r, v_open);
	dir_horiz = _mm_cmplt_epi16(E_r_gap,T1); /* E < H */
	dir_horiz = _mm_andnot_si128(dir_horiz,all_one_bits);	/* E >= H, for jump late */
	_mm_store_si128((__m128i *) &((*directions_Egap)[c][r]),dir_horiz);

	E_r_gap = _mm_max_epi16(E_r_gap, T1); /* Compare H + open with vert */
	E_r_gap = _mm_adds_epi16(E_r_gap, v_extend); /* Compute scores for Egap (vert + open) */
	debug15(print_vector_16(E_r_gap,r,c,"E"));
	if (r + SIMD_NSHORTS > rhigh) {
	  /* mask last block (e.g., NM_001193336_71_) */
	  E_r_gap = _mm_min_epi16(E_r_gap, E_mask_bottom); /* To handle band limitation from below */
	}
	_mm_store_si128(&(EE[(r-1)/SIMD_NSHORTS]), E_r_gap);


	/* NOGAP */
	T1 = _mm_srli_si128(H_nogap_r,LAST_SHORT);
	H_nogap_r = _mm_slli_si128(H_nogap_r,ONE_SHORT);
	H_nogap_r = _mm_or_si128(H_nogap_r, X_prev_nogap);
	X_prev_nogap = T1;

	/* Add pairscores, allowing for alternate genomic nt */
	pairscores_std = _mm_load_si128((__m128i *) &(pairscores_std_ptr[r-1]));
	pairscores_alt = _mm_load_si128((__m128i *) &(pairscores_alt_ptr[r-1]));
	H_nogap_r = _mm_adds_epi16(H_nogap_r, _mm_max_epi16(pairscores_std,pairscores_alt));
	_mm_clflush(&H_nogap_r); /* Needed for opencc -O3 on AMD */
	debug15(print_vector_16(H_nogap_r,r,c,"H"));

	dir_horiz = _mm_cmplt_epi16(E_r_gap,H_nogap_r); /* E < H */
	dir_horiz = _mm_andnot_si128(dir_horiz,all_one_bits);	/* E >= H, for jump late */
	_mm_store_si128((__m128i *) &((*directions_nogap)[c][r]),dir_horiz);

	H_nogap_r = _mm_max_epi16(H_nogap_r, E_r_gap); /* Compare H + pairscores with horiz + extend */
	debug15(print_vector_16(H_nogap_r,r,c,"H_nogap_r store"));
	_mm_store_si128((__m128i *) &(score_column[r]), H_nogap_r);
      }

      /* Perform F loop on entire column */
      if ((rlo = c - uband) < 1) {
	rlo = 1;
      }
      /* Use same rhigh as computed previously */

      c_gap = NEG_INFINITY_16;
      score_ptr = &(score_column[rlo]);
      for (r = rlo; r <= rhigh; r++) {

	/* FGAP */
	/* debug2(printf("Fgap at r %d, c %d: c_gap + extend %d vs last_nogap + open + extend %d\n",r,c,c_gap + extend,last_nogap + open + extend)); */
	if (c_gap /* + extend */ >= (score = last_nogap + open /* + extend */)) {  /* Use >= for jump late */
	  c_gap += extend;
	  (*directions_Fgap)[c][r] = VERT;
	} else {
	  c_gap = score + extend;
	  /* (*directions_Fgap)[c][r] = DIAG: -- Already initialized to DIAG */
	}

	/* NOGAP */
	last_nogap = *score_ptr;
	/* debug2(printf("assign nogap at r %d, c %d: H/E %d vs vert + extend %d\n",r,c,last_nogap,c_gap)); */
	if (c_gap >= last_nogap) {  /* Use >= for jump late */
	  last_nogap = c_gap;
#if 0
	  /* Should not go past NEG_INFINITY_16 for reasonable matrix sizes */
	  *score_ptr = (c_gap < NEG_INFINITY_16) ? NEG_INFINITY_16 : (Score16_T) c_gap; /* Saturation */
#else
	  *score_ptr = c_gap;
#endif
	  (*directions_nogap)[c][r] = VERT;
	}
	score_ptr++;
      }
    }

  } else {
    /* jump early */
    penalty = rpenalty = open + extend;
    for (c = 1; c <= glength; c++) {
      score_column = matrix[c];
      na2 = revp ? gsequence[1-c] : gsequence[c-1];
      na2_alt = revp ? gsequence_alt[1-c] : gsequence_alt[c-1];
#ifdef DEBUG14
      if (revp == false) {
	na2_single = get_genomic_nt(&na2_alt,goffset+c-1,chroffset,chrhigh,watsonp);
      } else {
	na2_single = get_genomic_nt(&na2_alt,goffset+1-c,chroffset,chrhigh,watsonp);
      }
      if (na2 != na2_single) {
	abort();
      }
#endif

      switch (na2) {
      case 'A': pairscores_std_ptr = pairscores[0]; break;
      case 'C': pairscores_std_ptr = pairscores[1]; break;
      case 'G': pairscores_std_ptr = pairscores[2]; break;
      case 'T': pairscores_std_ptr = pairscores[3]; break;
      default: pairscores_std_ptr = pairscores[4];
      }
      switch (na2_alt) {
      case 'A': pairscores_alt_ptr = pairscores[0]; break;
      case 'C': pairscores_alt_ptr = pairscores[1]; break;
      case 'G': pairscores_alt_ptr = pairscores[2]; break;
      case 'T': pairscores_alt_ptr = pairscores[3]; break;
      default: pairscores_alt_ptr = pairscores[4];
      }

      if (c == 1) {
	rlo = 1;
	X_prev_nogap = _mm_set1_epi16(0);
	last_nogap = penalty;

      } else if ((rlo = c - uband) < 1) {
	rlo = 1;
	if (penalty < NEG_INFINITY_16) {
	  X_prev_nogap = _mm_set1_epi16(NEG_INFINITY_16); /* Unlikely */
	} else {
	  X_prev_nogap = _mm_set1_epi16(penalty);
	}
	X_prev_nogap = _mm_srli_si128(X_prev_nogap,LAST_SHORT);
	penalty += extend;
	last_nogap = penalty;
      
      } else if (rlo == 1) {
	if (penalty < NEG_INFINITY_16) {
	  X_prev_nogap = _mm_set1_epi16(NEG_INFINITY_16); /* Unlikely */
	} else {
	  X_prev_nogap = _mm_set1_epi16(penalty);
	}
	X_prev_nogap = _mm_srli_si128(X_prev_nogap,LAST_SHORT);
	/* penalty += extend; */
	last_nogap = NEG_INFINITY_16;

      } else if ((rlo_floor = (int)((rlo - 1)/SIMD_NSHORTS) * SIMD_NSHORTS + 1) == 1) {
	/* first block of 8 */
	X_prev_nogap = _mm_set1_epi16(NEG_INFINITY_16); /* works if we start outside the rlo bounds */
	X_prev_nogap = _mm_srli_si128(X_prev_nogap,LAST_SHORT);
	last_nogap = NEG_INFINITY_16;
	rlo = rlo_floor;

      } else {
	/* second or greater block of 8 */
	X_prev_nogap = _mm_set1_epi16(matrix[c-1][rlo-1]); /* get H from previous block and previous column */
	X_prev_nogap = _mm_srli_si128(X_prev_nogap,LAST_SHORT);
	last_nogap = NEG_INFINITY_16;
	rlo = rlo_floor;
      }

      if ((rhigh = c + lband) > rlength) {
	rhigh = rlength;
	rhigh_ceil = (int) ((rhigh + SIMD_NSHORTS - 1)/SIMD_NSHORTS) * SIMD_NSHORTS;
	E_mask_bottom = bottom_masks[0];
      } else {
	rhigh_ceil = (int) ((rhigh + SIMD_NSHORTS - 1)/SIMD_NSHORTS) * SIMD_NSHORTS;
	E_mask_bottom = bottom_masks[rhigh_ceil-rhigh+1];
      }

      if (rhigh_ceil > max_rhigh_filled) {
	for (r = max_rhigh_filled; r <= rhigh_ceil; r++) {
	  if (r > lband) {
	    matrix[c-1][r] = NEG_INFINITY_16;
	  } else if (rpenalty < NEG_INFINITY_16) {
	    /* Unlikely */
	    matrix[c-1][r] = NEG_INFINITY_16;
	    rpenalty += extend;
	  } else {
	    matrix[c-1][r] = rpenalty;
	    rpenalty += extend;
	  }
	}
	max_rhigh_filled = r;
      }


      for (r = rlo; r <= rhigh; r += SIMD_NSHORTS) {
	/* Load previous E vector at this point in the query sequence */
	/* H vector already loaded before loop or at bottom of loop */
	E_r_gap = _mm_load_si128(&(EE[(r-1)/SIMD_NSHORTS]));
	debug15(print_vector_16(E_r_gap,r,c,"E_r_gap"));
	H_nogap_r = _mm_load_si128((__m128i *) &(matrix[c-1][r]));
	debug15(print_vector_16(H_nogap_r,r,c,"H_nogap_r load"));

	/* EGAP */
	T1 = _mm_adds_epi16(H_nogap_r, v_open);
	dir_horiz = _mm_cmpgt_epi16(E_r_gap,T1); /* E > H, for jump early */
	_mm_store_si128((__m128i *) &((*directions_Egap)[c][r]),dir_horiz);

	E_r_gap = _mm_max_epi16(E_r_gap, T1); /* Compare H + open with vert */
	E_r_gap = _mm_adds_epi16(E_r_gap, v_extend); /* Compute scores for Egap (vert + open) */
	debug15(print_vector_16(E_r_gap,r,c,"E"));
	if (r + SIMD_NSHORTS > rhigh) {
	  /* mask last block (e.g., NM_001193336_71_) */
	  E_r_gap = _mm_min_epi16(E_r_gap, E_mask_bottom); /* To handle band limitation from below */
	}
	_mm_store_si128(&(EE[(r-1)/SIMD_NSHORTS]), E_r_gap);


	/* NOGAP */
	T1 = _mm_srli_si128(H_nogap_r,LAST_SHORT);
	H_nogap_r = _mm_slli_si128(H_nogap_r,ONE_SHORT);
	H_nogap_r = _mm_or_si128(H_nogap_r, X_prev_nogap);
	X_prev_nogap = T1;

	/* Add pairscores, allowing for alternate genomic nt */
	pairscores_std = _mm_load_si128((__m128i *) &(pairscores_std_ptr[r-1]));
	pairscores_alt = _mm_load_si128((__m128i *) &(pairscores_alt_ptr[r-1]));
	H_nogap_r = _mm_adds_epi16(H_nogap_r, _mm_max_epi16(pairscores_std,pairscores_alt));
	_mm_clflush(&H_nogap_r); /* Needed for opencc -O3 on AMD */
	debug15(print_vector_16(H_nogap_r,r,c,"H"));

	dir_horiz = _mm_cmpgt_epi16(E_r_gap,H_nogap_r); /* E > H, for jump early */
	_mm_store_si128((__m128i *) &((*directions_nogap)[c][r]),dir_horiz);

	H_nogap_r = _mm_max_epi16(H_nogap_r, E_r_gap); /* Compare H + pairscores with horiz + extend */
	debug15(print_vector_16(H_nogap_r,r,c,"H_nogap_r store"));
	_mm_store_si128((__m128i *) &(score_column[r]), H_nogap_r);
      }

      /* Perform F loop on entire column */
      if ((rlo = c - uband) < 1) {
	rlo = 1;
      }
      /* Use same rhigh as computed previously */

      c_gap = NEG_INFINITY_16;
      score_ptr = &(score_column[rlo]);
      for (r = rlo; r <= rhigh; r++) {

	/* FGAP */
	/* debug2(printf("Fgap at r %d, c %d: c_gap + extend %d vs last_nogap + open + extend %d\n",r,c,c_gap + extend,last_nogap + open + extend)); */
	if (c_gap /* + extend */ > (score = last_nogap + open /* + extend */)) {  /* Use > for jump early */
	  c_gap += extend;
	  (*directions_Fgap)[c][r] = VERT;
	} else {
	  c_gap = score + extend;
	  /* (*directions_Fgap)[c][r] = DIAG: -- Already initialized to DIAG */
	}

	/* NOGAP */
	last_nogap = *score_ptr;
	/* debug2(printf("assign nogap at r %d, c %d: H/E %d vs vert + extend %d\n",r,c,last_nogap,c_gap)); */
	if (c_gap > last_nogap) {  /* Use > for jump early */
	  last_nogap = c_gap;
#if 0
	  /* Should not go past NEG_INFINITY_16 for reasonable matrix sizes */
	  *score_ptr = (c_gap < NEG_INFINITY_16) ? NEG_INFINITY_16 : (Score16_T) c_gap; /* Saturation */
#else
	  *score_ptr = c_gap;
#endif
	  (*directions_nogap)[c][r] = VERT;
	}
	score_ptr++;
      }
    }
  }

#ifdef DEBUG2
  printf("SIMD\n");
  Matrix16_print(matrix,rlength,glength,rsequence,gsequence,gsequence_alt,
		 goffset,chroffset,chrhigh,watsonp,revp);
  Directions16_print(*directions_nogap,*directions_Egap,*directions_Fgap,
		     rlength,glength,rsequence,gsequence,gsequence_alt,
		     goffset,chroffset,chrhigh,watsonp,revp);
#endif

#ifdef DEBUG14
  matrix_std = compute_scores_standard(&directions_nogap_std,&directions_Egap_std,&directions_Fgap_std,
				       this,rsequence,/*gsequence (NULL for debugging)*/NULL,/*gsequence_alt*/NULL,
				       rlength,glength,
				       goffset,chroffset,chrhigh,watsonp,mismatchtype,
				       open,extend,lband,uband,jump_late_p,revp);

#ifdef DEBUG2
  printf("Banded\n");
  Matrix32_print(matrix_std,rlength,glength,rsequence,gsequence,gsequence_alt,
		 goffset,chroffset,chrhigh,watsonp,revp);
  Directions32_print(directions_nogap_std,directions_Egap_std,directions_Fgap_std,
		     rlength,glength,rsequence,gsequence,gsequence_alt,
		     goffset,chroffset,chrhigh,watsonp,revp);
#endif
  
  banded_matrix16_compare(matrix,matrix_std,rlength,glength,lband,uband,
			  rsequence,gsequence,gsequence_alt,
			  goffset,chroffset,chrhigh,watsonp,revp);

  banded_directions16_compare_nogap(*directions_nogap,directions_nogap_std,rlength,glength,lband,uband);

  banded_directions16_compare_Egap(*directions_Egap,directions_Egap_std,rlength,glength,lband,uband);

  banded_directions16_compare_Fgap(*directions_Fgap,directions_Fgap_std,rlength,glength,lband,uband);
#endif

  _mm_free(EE);
  _mm_free(pairscores[4]);
  _mm_free(pairscores[3]);
  _mm_free(pairscores[2]);
  _mm_free(pairscores[1]);
  _mm_free(pairscores[0]);

  return matrix;
}
#endif


#if 0
static Score16_T **
compute_scores (Direction16_T ***directions, Score_T ***jump, T this, 
		char *rsequence, char *gsequence, bool revp, 
		int rlength, int glength, Mismatchtype_T mismatchtype,
		Score_T open, Score_T extend, bool init_jump_penalty_p,
		int extraband, bool widebandp, bool onesidegapp, bool jump_late_p) {
  Score_T **matrix;
  int r, c, r1, c1, na1, na2, j;
  Score_T bestscore, score, bestjump;
  Direction16_T bestdir;
  int lband, uband, clo, chigh, rlo;
  Pairdistance_T **pairdistance_array_type;

  pairdistance_array_type = pairdistance_array[mismatchtype];

  if (widebandp == false) {
    /* Just go along main diagonal */
    lband = extraband;
    uband = extraband;
  } else if (glength >= rlength) {
    /* Widen band to right to reach destination */
    uband = glength - rlength + extraband;
    lband = extraband;
  } else {
    /* Widen band to left to reach destination */
    lband = rlength - glength + extraband;
    uband = extraband;
  }

  matrix = Matrix16_alloc(rlength,glength,this->matrix_ptrs,this->matrix_space);
  *directions = Directions16_alloc(rlength,glength,this->directions_ptrs,this->directions_space);
  *jump = Matrix16_alloc(rlength,glength,this->jump_ptrs,this->jump_space);

  matrix[0][0] = 0;
  /* (*directions)[0][0] = STOP; -- Check for r > 0 && c > 0 instead */
  (*jump)[0][0] = 0;

  /* Row 0 initialization */
  if (init_jump_penalty_p == true) {
    for (c = 1; c <= uband && c <= glength; c++) {
      matrix[0][c] = jump_penalty(c,open,extend);
      (*directions)[0][c] = HORIZ;
      (*jump)[0][c] = c;
    }
  } else {
    for (c = 1; c <= uband && c <= glength; c++) {
      matrix[0][c] = SINGLE_OPEN_HIGHQ;	/* Needs to be less than 0 to prevent gap and then a single match */
      (*directions)[0][c] = HORIZ;
      (*jump)[0][c] = c;
    }
  }

  /* Column 0 initialization */
  if (init_jump_penalty_p == true) {
    for (r = 1; r <= lband && r <= rlength; r++) {
      matrix[r][0] = jump_penalty(r,open,extend);
      (*directions)[r][0] = VERT;
      (*jump)[r][0] = r;
    }
  } else {
    for (r = 1; r <= lband && r <= rlength; r++) {
      matrix[r][0] = SINGLE_OPEN_HIGHQ;	/* Needs to be less than 0 to prevent gap and then a single match */
      (*directions)[r][0] = VERT;
      (*jump)[r][0] = r;
    }
  }

  for (r = 1; r <= rlength; r++) {
    na1 = revp ? rsequence[1-r] : rsequence[r-1];

    if ((clo = r - lband) < 1) {
      clo = 1;
    }
    if ((chigh = r + uband) > glength) {
      chigh = glength;
    }
    for (c = clo; c <= chigh; c++) {
      na2 = revp ? gsequence[1-c] : gsequence[c-1];

      /* Diagonal case */
      bestscore = matrix[r-1][c-1] + pairdistance_array_type[na1][na2];
      bestdir = DIAG;
      bestjump = 1;
      
      /* Horizontal case */
      if (onesidegapp == true) {
	/* Don't allow horizontal jumps below the main diagonal, and
	   don't cross the main diagonal */
	if (r <= c) {
	  for (c1 = c-1, j = 1; c1 >= r; c1--, j++) {
	    score = matrix[r][c1] + jump_penalty(j,open,extend);
	    if (score > bestscore || (score == bestscore && jump_late_p)) {
	      bestscore = score;
	      bestdir = HORIZ;
	      bestjump = j;
	    }
	  }
	}
      } else {
	for (c1 = c-1, j = 1; c1 >= clo; c1--, j++) {
	  if ((*directions)[r][c1] == DIAG) {
	    score = matrix[r][c1] + jump_penalty(j,open,extend);
	    if (score > bestscore || (score == bestscore && jump_late_p)) {
	      bestscore = score;
	      bestdir = HORIZ;
	      bestjump = j;
	    }
#ifdef RIGHTANGLE
	  } else if ((*directions)[r][c1] == VERT && (*jump)[r][c1] >= 3) {
	    fprintf(stderr,"Error.  Not allowed to check for dir == VERT.\n");
	    abort();
	    score = matrix[r][c1] + jump_penalty(j,open,extend) - open;
	    if (score > bestscore || (score == bestscore && jump_late_p)) {
	      bestscore = score;
	      bestdir = HORIZ;
	      bestjump = j;
	    }
#endif
	  }
	}
      }

      /* Vertical case */
      if (onesidegapp == true) {
	/* Don't allow vertical jumps above the main diagonal, and don't
	   cross the main diagonal */
	if (c <= r) {
	  for (r1 = r-1, j = 1; r1 >= c; r1--, j++) {
	    score = matrix[r1][c] + jump_penalty(j,open,extend);
	    if (score > bestscore || (score == bestscore && jump_late_p)) {
	      bestscore = score;
	      bestdir = VERT;
	      bestjump = j;
	    }
	  }
	}
      } else {
	if ((rlo = c+c-uband-r) < 1) {
	  rlo = 1;
	}
	for (r1 = r-1, j = 1; r1 >= rlo; r1--, j++) {
	  if ((*directions)[r1][c] == DIAG) {
	    score = matrix[r1][c] + jump_penalty(j,open,extend);
	    if (score > bestscore || (score == bestscore && jump_late_p)) {
	      bestscore = score;
	      bestdir = VERT;
	      bestjump = j;
	    }
#ifdef RIGHTANGLE
	  } else if ((*directions)[r1][c] == HORIZ && (*jump)[r1][c] >= 3) {
	    score = matrix[r1][c] + jump_penalty(j,open,extend) - open;
	    if (score > bestscore || (score == bestscore && jump_late_p)) {
	      bestscore = score;
	      bestdir = VERT;
	      bestjump = j;
	    }
#endif
	  }
	}
      }

      /*
	debug(printf("At %d,%d, scoreV = %d, scoreH = %d, scoreD = %d\n",
	r,c,scoreV,scoreH,scoreD));
      */
      
      /* Update */
      matrix[r][c] = bestscore;
      (*directions)[r][c] = bestdir;
      (*jump)[r][c] = bestjump;
    }
  }

  debug2(Matrix16_print(matrix,rlength,glength));
  debug2(Directions16_print(*directions,*jump,rlength,glength,rsequence,gsequence,revp));

  return matrix;
}
#endif


#if defined(HAVE_SSE4_1) || defined(HAVE_SSE2)
static void
find_best_endpoint_8 (int *finalscore, int *bestr, int *bestc, Score8_T **matrix, 
		      int rlength, int glength, int extraband_end_or_paired,
		      bool jump_late_p) {
  Score8_T bestscore = 0;
  int r, c;
  int uband, lband, clo, chigh;
  /* No need for loffset or cmid because they apply only for cdnaend
     == FIVE, which doesn't require searching */

  *bestr = *bestc = 0;

  /* Just go along main diagonal */
  uband = extraband_end_or_paired;
  lband = extraband_end_or_paired;

  if (jump_late_p == false) {
    /* use > bestscore */
    for (r = 1; r <= rlength; r++) {
      if ((clo = r - lband) < 1) {
	clo = 1;
      }
      if ((chigh = r + uband) > glength) {
	chigh = glength;
      }
      for (c = clo; c <= chigh; c++) {
	if (matrix[c][r] > bestscore) {
	  *bestr = r;
	  *bestc = c;
	  bestscore = matrix[c][r];
	}
      }
    }
  } else {
    /* use >= bestscore */
    for (r = 1; r <= rlength; r++) {
      if ((clo = r - lband) < 1) {
	clo = 1;
      }
      if ((chigh = r + uband) > glength) {
	chigh = glength;
      }
      for (c = clo; c <= chigh; c++) {
	if (matrix[c][r] >= bestscore) {
	  *bestr = r;
	  *bestc = c;
	  bestscore = matrix[c][r];
	}
      }
    }
  }


  *finalscore = (int) bestscore;
  return;
}
#endif

static void
find_best_endpoint (int *finalscore, int *bestr, int *bestc,
#ifdef HAVE_SSE2
		    Score16_T **matrix, 
#else
		    Score32_T **matrix,
#endif
		    int rlength, int glength, int extraband_end_or_paired,
		    bool jump_late_p) {
#ifdef HAVE_SSE2
  Score16_T bestscore = 0;
#else
  Score32_T bestscore = 0;
#endif
  int r, c;
  int uband, lband, clo, chigh;
  /* No need for loffset or cmid because they apply only for cdnaend
     == FIVE, which doesn't require searching */

  *bestr = *bestc = 0;

  /* Just go along main diagonal */
  uband = extraband_end_or_paired;
  lband = extraband_end_or_paired;

  if (jump_late_p == false) {
    /* use > bestscore */
    for (r = 1; r <= rlength; r++) {
      if ((clo = r - lband) < 1) {
	clo = 1;
      }
      if ((chigh = r + uband) > glength) {
	chigh = glength;
      }
      for (c = clo; c <= chigh; c++) {
	if (matrix[c][r] > bestscore) {
	  *bestr = r;
	  *bestc = c;
	  bestscore = matrix[c][r];
	}
      }
    }
  } else {
    /* use >= bestscore */
    for (r = 1; r <= rlength; r++) {
      if ((clo = r - lband) < 1) {
	clo = 1;
      }
      if ((chigh = r + uband) > glength) {
	chigh = glength;
      }
      for (c = clo; c <= chigh; c++) {
	if (matrix[c][r] >= bestscore) {
	  *bestr = r;
	  *bestc = c;
	  bestscore = matrix[c][r];
	}
      }
    }
  }


  *finalscore = (int) bestscore;
  return;
}


#if defined(HAVE_SSE4_1) || defined(HAVE_SSE2)
static void
find_best_endpoint_to_queryend_indels_8 (int *finalscore, int *bestr, int *bestc, Score8_T **matrix, 
					 int rlength, int glength, int extraband_end_or_paired,
					 bool jump_late_p) {
  Score8_T bestscore = NEG_INFINITY_8;
  int r, c;
  int uband, lband, clo, chigh;
  /* No need for loffset or cmid because they apply only for cdnaend
     == FIVE, which doesn't require searching */

  if (glength >= rlength) {
    /* Widen band to right to reach destination */
    uband = glength - rlength + extraband_end_or_paired;
    lband = extraband_end_or_paired;
  } else {
    /* Widen band to left to reach destination */
    lband = rlength - glength + extraband_end_or_paired;
    uband = extraband_end_or_paired;
  }

  *bestr = r = rlength;
  *bestc = 0;

  if (jump_late_p == false) {
    if ((clo = r - lband) < 1) {
      clo = 1;
    }
    if ((chigh = r + uband) > glength) {
      chigh = glength;
    }
    if (clo <= chigh) {
      for (c = clo; c <= chigh; c++) {
	/* use > bestscore */
	if (matrix[c][r] > bestscore) {
	  *bestr = r;
	  *bestc = c;
	  bestscore = matrix[c][r];
	}
      }
    }

  } else {
    if ((clo = r - lband) < 1) {
      clo = 1;
    }
    if ((chigh = r + uband) > glength) {
      chigh = glength;
    }
    if (clo <= chigh) {
      for (c = clo; c <= chigh; c++) {
	/* use >= bestscore */
	if (matrix[c][r] >= bestscore) {
	  *bestr = r;
	  *bestc = c;
	  bestscore = matrix[c][r];
	}
      }
    }
  }

  *finalscore = (int) bestscore;
  return;
}
#endif


static void
find_best_endpoint_to_queryend_indels (int *finalscore, int *bestr, int *bestc,
#ifdef HAVE_SSE2
				       Score16_T **matrix, 
#else
				       Score32_T **matrix, 
#endif
				       int rlength, int glength, int extraband_end_or_paired,
				       bool jump_late_p) {
#ifdef HAVE_SSE2
  Score16_T bestscore = NEG_INFINITY_16;
#else
  Score32_T bestscore = NEG_INFINITY_32;
#endif
  int r, c;
  int uband, lband, clo, chigh;
  /* No need for loffset or cmid because they apply only for cdnaend
     == FIVE, which doesn't require searching */

  if (glength >= rlength) {
    /* Widen band to right to reach destination */
    uband = glength - rlength + extraband_end_or_paired;
    lband = extraband_end_or_paired;
  } else {
    /* Widen band to left to reach destination */
    lband = rlength - glength + extraband_end_or_paired;
    uband = extraband_end_or_paired;
  }

  *bestr = r = rlength;
  *bestc = 0;

  if (jump_late_p == false) {
    if ((clo = r - lband) < 1) {
      clo = 1;
    }
    if ((chigh = r + uband) > glength) {
      chigh = glength;
    }
    if (clo <= chigh) {
      for (c = clo; c <= chigh; c++) {
	/* use > bestscore */
	if (matrix[c][r] > bestscore) {
	  *bestr = r;
	  *bestc = c;
	  bestscore = matrix[c][r];
	}
      }
    }

  } else {
    if ((clo = r - lband) < 1) {
      clo = 1;
    }
    if ((chigh = r + uband) > glength) {
      chigh = glength;
    }
    if (clo <= chigh) {
      for (c = clo; c <= chigh; c++) {
	/* use >= bestscore */
	if (matrix[c][r] >= bestscore) {
	  *bestr = r;
	  *bestc = c;
	  bestscore = matrix[c][r];
	}
      }
    }
  }

  *finalscore = (int) bestscore;
  return;
}


static void
find_best_endpoint_to_queryend_nogaps (int *bestr, int *bestc, int rlength, int glength) {
  if (glength < rlength) {
    *bestr = glength;
    *bestc = glength;
  } else {
    *bestr = rlength;
    *bestc = rlength;
  }

  return;
}


static List_T
add_queryskip (List_T pairs, int r, int c, int dist, char *querysequence,
	       int queryoffset, int genomeoffset, Pairpool_T pairpool, bool revp, int dynprogindex) {
  int j;
  char c1;
  int querycoord, genomecoord, step;

  querycoord = r-1;
  genomecoord = c-1;

  if (revp == true) {
    querycoord = -querycoord;
    genomecoord = -genomecoord;
    step = +1;
  } else {
    /* Advance to next genomepos */
    genomecoord++;
    step = -1;
  }

  for (j = 0; j < dist; j++) {
    c1 = querysequence[querycoord];
    debug(printf("Pushing %d,%d [%d,%d] (%c,-), ",r,c,queryoffset+querycoord,genomeoffset+genomecoord,c1));
    pairs = Pairpool_push(pairs,pairpool,queryoffset+querycoord,genomeoffset+genomecoord,
			  c1,INDEL_COMP,/*genome*/' ',/*genomealt*/' ',dynprogindex);
    debug(r--);
    querycoord += step;
  }

  return pairs;
}


static List_T
add_genomeskip (bool *add_dashes_p, List_T pairs, int r, int c, int dist,
		char *genomesequence, char *genomesequenceuc,
		int queryoffset, int genomeoffset, Pairpool_T pairpool, bool revp,
		Univcoord_T chroffset, Univcoord_T chrhigh,
		int cdna_direction, bool watsonp, int dynprogindex, bool use_genomicseg_p) {
  int j;
  char left1, left2, right2, right1, left1_alt, left2_alt, right2_alt, right1_alt, c2, c2_alt;
  int introntype;
  int querycoord, leftgenomecoord, rightgenomecoord, genomecoord, temp, step;

  querycoord = r-1;
  leftgenomecoord = c-dist;
  rightgenomecoord = c-1;

  if (revp == true) {
    querycoord = -querycoord;
    temp = leftgenomecoord;
    leftgenomecoord = -rightgenomecoord;
    rightgenomecoord = -temp;
    step = +1;
  } else {
    /* Advance to next querypos */
    querycoord++;
    step = -1;
  }

#if 0
  if (dist < MICROINTRON_LENGTH) {
    *add_dashes_p = true;
  } else if (cdna_direction == 0) {
    *add_dashes_p = true;
  } else {
    /* Check for intron */
    if (use_genomicseg_p) {
      left1 = left1_alt = genomesequenceuc[leftgenomecoord];
      left2 = left2_alt = genomesequenceuc[leftgenomecoord+1];
      right2 = right2_alt = genomesequenceuc[rightgenomecoord-1];
      right1 = right1_alt = genomesequenceuc[rightgenomecoord];
    } else {
      left1 = get_genomic_nt(&left1_alt,genomeoffset+leftgenomecoord,chroffset,chrhigh,watsonp);
      left2 = get_genomic_nt(&left2_alt,genomeoffset+leftgenomecoord+1,chroffset,chrhigh,watsonp);
      right2 = get_genomic_nt(&right2_alt,genomeoffset+rightgenomecoord-1,chroffset,chrhigh,watsonp);
      right1 = get_genomic_nt(&right1_alt,genomeoffset+rightgenomecoord,chroffset,chrhigh,watsonp);
    }
#ifdef EXTRACT_GENOMICSEG
    assert(left1 == genomesequenceuc[leftgenomecoord]);
    assert(left2 == genomesequenceuc[leftgenomecoord+1]);
    assert(right2 == genomesequenceuc[rightgenomecoord-1]);
    assert(right1 == genomesequenceuc[rightgenomecoord]);
#endif
	
#ifdef PMAP
    introntype = Intron_type(left1,left2,right2,right1,
			     left1_alt,left2_alt,right2_alt,right1_alt,
			     /*cdna_direction*/+1);
#else
    introntype = Intron_type(left1,left2,right2,right1,
			     left1_alt,left2_alt,right2_alt,right1_alt,
			     cdna_direction);
#endif
    if (introntype == NONINTRON) {
      *add_dashes_p = true;
    } else {
      *add_dashes_p = false;
    }
  }
#endif

  *add_dashes_p = true;
  if (*add_dashes_p == true) {
    if (revp == true) {
      genomecoord = leftgenomecoord;
    } else {
      genomecoord = rightgenomecoord;
    }
    for (j = 0; j < dist; j++) {
      if (use_genomicseg_p) {
	c2 = c2_alt = genomesequence[genomecoord];
      } else {
	c2 = get_genomic_nt(&c2_alt,genomeoffset+genomecoord,chroffset,chrhigh,watsonp);
      }
#ifdef EXTRACT_GENOMICSEG
      assert(c2 == genomesequence[genomecoord]);
#endif

      debug(printf("Pushing %d,%d [%d,%d] (-,%c),",r,c,queryoffset+querycoord,genomeoffset+genomecoord,c2));
      pairs = Pairpool_push(pairs,pairpool,queryoffset+querycoord,genomeoffset+genomecoord,
			    ' ',INDEL_COMP,c2,c2_alt,dynprogindex);
      debug(c--);
      genomecoord += step;
    }
  } else {
    debug(printf("Large gap %c%c..%c%c.  Adding gap of type %d.\n",
		 left1,left2,right2,right1,introntype));
#ifndef NOGAPHOLDER
    pairs = Pairpool_push_gapholder(pairs,pairpool,/*queryjump*/0,/*genomejump*/dist,
				    /*leftpair*/NULL,/*rightpair*/NULL,/*knownp*/false);
#endif
  }
  
  return pairs;
}


/* Preferences: Continue in same direction if possible.  This has the
   effect of extending gaps to maximum size. Then, take diagonal if
   possible. Finally, take vertical if possible, because this will use
   up sequence 1, which is the query (cDNA) sequence. */

#if 0
/* revp means both rev1p and rev2p, which must have equal values */
/* Iterative version */
static List_T
traceback_genomicseg (List_T pairs, int *nmatches, int *nmismatches, int *nopens, int *nindels,
		      struct Direction3_T **directions, int r, int c, 
		      char *querysequence, char *querysequenceuc, char *genomesequence, char *genomesequenceuc,
		      int queryoffset, int genomeoffset, Pairpool_T pairpool, bool revp,
		      Univcoord_T chroffset, Univcoord_T chrhigh,
		      int cdna_direction, bool watsonp, int dynprogindex, bool use_genomicseg_p) {
  char c1, c2, c2_alt;
  int dist;
  bool add_dashes_p;
  int querycoord, genomecoord;

  debug(printf("Starting traceback at r=%d,c=%d (roffset=%d, goffset=%d)\n",r,c,queryoffset,genomeoffset));

  while (r > 0 && c > 0) {  /* directions[c][r].nogap != STOP */
    querycoord = r-1;
    genomecoord = c-1;
    if (revp == true) {
      querycoord = -querycoord;
      genomecoord = -genomecoord;
    }

    c1 = querysequence[querycoord];
    if (use_genomicseg_p) {
      c2 = c2_alt = genomesequence[genomecoord];
    } else {
      c2 = get_genomic_nt(&c2_alt,genomeoffset+genomecoord,chroffset,chrhigh,watsonp);
    }
#ifdef EXTRACT_GENOMICSEG
    assert(c2 == genomesequence[genomecoord]);
#endif

    if (c2 == '*') {
      /* Don't push pairs past end of chromosome */

    } else if (/*querysequenceuc[querycoord]*/c1 == c2) {
      debug(printf("D1: Pushing %d,%d [%d,%d] (%c,%c) - match\n",
		   r,c,queryoffset+querycoord,genomeoffset+genomecoord,c1,c2));
      *nmatches += 1;
      pairs = Pairpool_push(pairs,pairpool,queryoffset+querycoord,genomeoffset+genomecoord,
			    c1,DYNPROG_MATCH_COMP,c2,c2_alt,dynprogindex);

    } else if (consistent_array[(int) c1][(int) c2] == true) {
      debug(printf("D1: Pushing %d,%d [%d,%d] (%c,%c) - ambiguous\n",
		   r,c,queryoffset+querycoord,genomeoffset+genomecoord,c1,c2));
      *nmatches += 1;
      pairs = Pairpool_push(pairs,pairpool,queryoffset+querycoord,genomeoffset+genomecoord,
			    c1,AMBIGUOUS_COMP,c2,c2_alt,dynprogindex);

    } else {
      debug(printf("D1: Pushing %d,%d [%d,%d] (%c,%c) - mismatch\n",
		   r,c,queryoffset+querycoord,genomeoffset+genomecoord,c1,c2));
      *nmismatches += 1;
      pairs = Pairpool_push(pairs,pairpool,queryoffset+querycoord,genomeoffset+genomecoord,
			    c1,MISMATCH_COMP,c2,c2_alt,dynprogindex);
    }

    if (directions[c][r].nogap == DIAG) {
      r--; c--;

    } else if (directions[c][r].nogap == HORIZ) {
      dist = 1;
      r--; c--;
      while (c > 1 && directions[c][r].Egap != DIAG) {
	dist++;
	c--;
      }
      c--;

      debug(printf("H%d: ",dist));
      pairs = add_genomeskip(&add_dashes_p,pairs,r,c+dist,dist,genomesequence,genomesequenceuc,
			     queryoffset,genomeoffset,pairpool,revp,chroffset,chrhigh,
			     cdna_direction,watsonp,dynprogindex,use_genomicseg_p);
      if (add_dashes_p == true) {
	*nopens += 1;
	*nindels += dist;
      }
      debug(printf("\n"));

    } else {
      /* Must be VERT */
      dist = 1;
      r--; c--;
      while (r > 1 && directions[c][r].Fgap != DIAG) {
	dist++;
	r--;
      }
      r--;

      debug(printf("V%d: ",dist));
      pairs = add_queryskip(pairs,r+dist,c,dist,querysequence,
			    queryoffset,genomeoffset,pairpool,revp,
			    dynprogindex);
      *nopens += 1;
      *nindels += dist;
      debug(printf("\n"));
    }

  }
  return pairs;
}
#endif


#if 0
static List_T
traceback_prev (List_T pairs, int *nmatches, int *nmismatches, int *nopens, int *nindels,
		struct Direction3_T **directions, int r, int c, 
		char *querysequence, char *querysequenceuc,
		int queryoffset, int genomeoffset, Pairpool_T pairpool, bool revp,
		Univcoord_T chroffset, Univcoord_T chrhigh,
		int cdna_direction, bool watsonp, int dynprogindex) {
  char c1, c2, c2_alt;
  int dist;
  bool add_dashes_p;
  int querycoord, genomecoord;

  debug(printf("Starting traceback at r=%d,c=%d (roffset=%d, goffset=%d)\n",r,c,queryoffset,genomeoffset));

  while (r > 0 && c > 0) {  /* directions[c][r].nogap != STOP */
    querycoord = r-1;
    genomecoord = c-1;
    if (revp == true) {
      querycoord = -querycoord;
      genomecoord = -genomecoord;
    }

    c1 = querysequence[querycoord];
    c2 = get_genomic_nt(&c2_alt,genomeoffset+genomecoord,chroffset,chrhigh,watsonp);
#ifdef EXTRACT_GENOMICSEG
    assert(c2 == genomesequence[genomecoord]);
#endif

    if (c2 == '*') {
      /* Don't push pairs past end of chromosome */
      debug(printf("Don't push pairs past end of chromosome: genomeoffset %u, genomecoord %u, chroffset %u, chrhigh %u, watsonp %d\n",
		   genomeoffset,genomecoord,chroffset,chrhigh,watsonp));

    } else if (/*querysequenceuc[querycoord]*/c1 == c2) {
      debug(printf("D1: Pushing %d,%d [%d,%d] (%c,%c) - match\n",
		   r,c,queryoffset+querycoord,genomeoffset+genomecoord,c1,c2));
      *nmatches += 1;
      pairs = Pairpool_push(pairs,pairpool,queryoffset+querycoord,genomeoffset+genomecoord,
			    c1,DYNPROG_MATCH_COMP,c2,c2_alt,dynprogindex);

    } else if (consistent_array[(int) c1][(int) c2] == true) {
      debug(printf("D1: Pushing %d,%d [%d,%d] (%c,%c) - ambiguous\n",
		   r,c,queryoffset+querycoord,genomeoffset+genomecoord,c1,c2));
      *nmatches += 1;
      pairs = Pairpool_push(pairs,pairpool,queryoffset+querycoord,genomeoffset+genomecoord,
			    c1,AMBIGUOUS_COMP,c2,c2_alt,dynprogindex);

    } else {
      debug(printf("D1: Pushing %d,%d [%d,%d] (%c,%c) - mismatch\n",
		   r,c,queryoffset+querycoord,genomeoffset+genomecoord,c1,c2));
      *nmismatches += 1;
      pairs = Pairpool_push(pairs,pairpool,queryoffset+querycoord,genomeoffset+genomecoord,
			    c1,MISMATCH_COMP,c2,c2_alt,dynprogindex);
    }

    if (directions[c][r].nogap == DIAG) {
      r--; c--;

    } else if (directions[c][r].nogap == HORIZ) {
      dist = 1;
      r--; c--;
      while (c > 1 && directions[c][r].Egap != DIAG) {
	dist++;
	c--;
      }
      c--;

      debug(printf("H%d: ",dist));
      pairs = add_genomeskip(&add_dashes_p,pairs,r,c+dist,dist,
			     /*genomesequence*/NULL,/*genomesequenceuc*/NULL,
			     queryoffset,genomeoffset,pairpool,revp,chroffset,chrhigh,
			     cdna_direction,watsonp,dynprogindex,/*use_genomicseg_p*/false);
      if (add_dashes_p == true) {
	*nopens += 1;
	*nindels += dist;
      }
      debug(printf("\n"));

    } else {
      /* Must be VERT */
      dist = 1;
      r--; c--;
      while (r > 1 && directions[c][r].Fgap != DIAG) {
	dist++;
	r--;
      }
      r--;

      debug(printf("V%d: ",dist));
      pairs = add_queryskip(pairs,r+dist,c,dist,querysequence,
			    queryoffset,genomeoffset,pairpool,revp,
			    dynprogindex);
      *nopens += 1;
      *nindels += dist;
      debug(printf("\n"));
    }

  }
  return pairs;
}
#endif


#if defined(HAVE_SSE4_1) || defined(HAVE_SSE2)
static List_T
traceback_8 (List_T pairs, int *nmatches, int *nmismatches, int *nopens, int *nindels,
	     Direction8_T **directions_nogap, Direction8_T **directions_Egap, Direction8_T **directions_Fgap,
	     int r, int c, char *querysequence, char *querysequenceuc, char *gsequence, char *gsequence_alt,
	     int queryoffset, int genomeoffset, Pairpool_T pairpool, bool revp,
	     Univcoord_T chroffset, Univcoord_T chrhigh,
	     int cdna_direction, bool watsonp, int dynprogindex) {
  char c1, c2, c2_alt;
  int dist;
  bool add_dashes_p;
  int querycoord, genomecoord;
  Direction8_T dir;
#ifdef DEBUG14
  char c2_single;
#endif

  debug(printf("Starting traceback at r=%d,c=%d (roffset=%d, goffset=%d)\n",r,c,queryoffset,genomeoffset));

  /* Handle initial indel */
  if ((dir = directions_nogap[c][r]) == DIAG) {
    /* Not an indel.  Do nothing. */

  } else if (dir == HORIZ) {
    dist = 1;
    while (c > 1 && directions_Egap[c][r] != DIAG) {
      dist++;
      c--;
    }
    c--;
    /* dir = directions_nogap[c][r]; */

    debug(printf("H%d: ",dist));
    pairs = add_genomeskip(&add_dashes_p,pairs,r,c+dist,dist,
			   /*genomesequence*/NULL,/*genomesequenceuc*/NULL,
			   queryoffset,genomeoffset,pairpool,revp,chroffset,chrhigh,
			   cdna_direction,watsonp,dynprogindex,/*use_genomicseg_p*/false);
    if (add_dashes_p == true) {
      *nopens += 1;
      *nindels += dist;
    }
    debug(printf("\n"));

  } else {
    /* Must be VERT */
    dist = 1;
    while (r > 1 && directions_Fgap[c][r] != DIAG) {
      dist++;
      r--;
    }
    r--;
    /* dir = directions_nogap[c][r]; */

    debug(printf("V%d: ",dist));
    pairs = add_queryskip(pairs,r+dist,c,dist,querysequence,
			  queryoffset,genomeoffset,pairpool,revp,
			  dynprogindex);
    *nopens += 1;
    *nindels += dist;
    debug(printf("\n"));
  }

  while (r > 0 && c > 0) {  /* dir != STOP */
    querycoord = r-1;
    genomecoord = c-1;
    if (revp == true) {
      querycoord = -querycoord;
      genomecoord = -genomecoord;
    }

    c1 = querysequence[querycoord];
    c2 = gsequence[genomecoord];
    c2_alt = gsequence_alt[genomecoord];
#ifdef DEBUG14
    c2_single = get_genomic_nt(&c2_alt,genomeoffset+genomecoord,chroffset,chrhigh,watsonp);
    if (c2 != c2_single) {
      abort();
    }
#endif

#ifdef EXTRACT_GENOMICSEG
    assert(c2 == genomesequence[genomecoord]);
#endif

    if (c2 == '*') {
      /* Don't push pairs past end of chromosome */
      debug(printf("Don't push pairs past end of chromosome: genomeoffset %u, genomecoord %u, chroffset %u, chrhigh %u, watsonp %d\n",
		   genomeoffset,genomecoord,chroffset,chrhigh,watsonp));

    } else if (/*querysequenceuc[querycoord]*/c1 == c2) {
      debug(printf("Pushing %d,%d [%d,%d] (%c,%c) - match\n",
		   r,c,queryoffset+querycoord,genomeoffset+genomecoord,c1,c2));
      *nmatches += 1;
      pairs = Pairpool_push(pairs,pairpool,queryoffset+querycoord,genomeoffset+genomecoord,
			    c1,DYNPROG_MATCH_COMP,c2,c2_alt,dynprogindex);

    } else if (consistent_array[(int) c1][(int) c2] == true) {
      debug(printf("Pushing %d,%d [%d,%d] (%c,%c) - ambiguous\n",
		   r,c,queryoffset+querycoord,genomeoffset+genomecoord,c1,c2));
      *nmatches += 1;
      pairs = Pairpool_push(pairs,pairpool,queryoffset+querycoord,genomeoffset+genomecoord,
			    c1,AMBIGUOUS_COMP,c2,c2_alt,dynprogindex);

    } else {
      debug(printf("Pushing %d,%d [%d,%d] (%c,%c) - mismatch\n",
		   r,c,queryoffset+querycoord,genomeoffset+genomecoord,c1,c2));
      *nmismatches += 1;
      pairs = Pairpool_push(pairs,pairpool,queryoffset+querycoord,genomeoffset+genomecoord,
			    c1,MISMATCH_COMP,c2,c2_alt,dynprogindex);
    }

    r--; c--;
    if (r == 0 && c == 0) {
      /* STOP condition.  Do nothing. */

    } else if ((dir = directions_nogap[c][r]) == DIAG) {
      /* Do nothing */

    } else if (dir == HORIZ) {
      dist = 1;
      while (c > 1 && directions_Egap[c][r] != DIAG) {
	dist++;
	c--;
      }
      c--;
      /* dir = directions_nogap[c][r]; */

      debug(printf("H%d: ",dist));
      pairs = add_genomeskip(&add_dashes_p,pairs,r,c+dist,dist,
			     /*genomesequence*/NULL,/*genomesequenceuc*/NULL,
			     queryoffset,genomeoffset,pairpool,revp,chroffset,chrhigh,
			     cdna_direction,watsonp,dynprogindex,/*use_genomicseg_p*/false);
      if (add_dashes_p == true) {
	*nopens += 1;
	*nindels += dist;
      }
      debug(printf("\n"));

    } else {
      /* Must be VERT */
      dist = 1;
      while (r > 1 && directions_Fgap[c][r] != DIAG) {
	dist++;
	r--;
      }
      r--;
      /* dir = directions_nogap[c][r]; */

      debug(printf("V%d: ",dist));
      pairs = add_queryskip(pairs,r+dist,c,dist,querysequence,
			    queryoffset,genomeoffset,pairpool,revp,
			    dynprogindex);
      *nopens += 1;
      *nindels += dist;
      debug(printf("\n"));

    }
  }
  return pairs;
}
#endif


static List_T
traceback (List_T pairs, int *nmatches, int *nmismatches, int *nopens, int *nindels,
#ifdef HAVE_SSE2
	   Direction16_T **directions_nogap, Direction16_T **directions_Egap, Direction16_T **directions_Fgap,
#else
	   Direction32_T **directions_nogap, Direction32_T **directions_Egap, Direction32_T **directions_Fgap,
#endif
	   int r, int c, char *querysequence, char *querysequenceuc, char *gsequence, char *gsequence_alt,
	   int queryoffset, int genomeoffset, Pairpool_T pairpool, bool revp,
	   Univcoord_T chroffset, Univcoord_T chrhigh,
	   int cdna_direction, bool watsonp, int dynprogindex) {
  char c1, c2, c2_alt;
  int dist;
  bool add_dashes_p;
  int querycoord, genomecoord;
#ifdef HAVE_SSE2
  Direction16_T dir;
#else
  Direction32_T dir;
#endif
#ifdef DEBUG14
  char c2_single;
#endif

  debug(printf("Starting traceback at r=%d,c=%d (roffset=%d, goffset=%d)\n",r,c,queryoffset,genomeoffset));

  /* Handle initial indel */
  if ((dir = directions_nogap[c][r]) == DIAG) {
    /* Not an indel.  Do nothing. */

  } else if (dir == HORIZ) {
    dist = 1;
    while (c > 1 && directions_Egap[c][r] != DIAG) {
      dist++;
      c--;
    }
    c--;
    /* dir = directions_nogap[c][r]; */

    debug(printf("H%d: ",dist));
    pairs = add_genomeskip(&add_dashes_p,pairs,r,c+dist,dist,
			   /*genomesequence*/NULL,/*genomesequenceuc*/NULL,
			   queryoffset,genomeoffset,pairpool,revp,chroffset,chrhigh,
			   cdna_direction,watsonp,dynprogindex,/*use_genomicseg_p*/false);
    if (add_dashes_p == true) {
      *nopens += 1;
      *nindels += dist;
    }
    debug(printf("\n"));

  } else {
    /* Must be VERT */
    dist = 1;
    while (r > 1 && directions_Fgap[c][r] != DIAG) {
      dist++;
      r--;
    }
    r--;
    /* dir = directions_nogap[c][r]; */

    debug(printf("V%d: ",dist));
    pairs = add_queryskip(pairs,r+dist,c,dist,querysequence,
			  queryoffset,genomeoffset,pairpool,revp,
			  dynprogindex);
    *nopens += 1;
    *nindels += dist;
    debug(printf("\n"));
  }

  while (r > 0 && c > 0) {  /* dir != STOP */
    querycoord = r-1;
    genomecoord = c-1;
    if (revp == true) {
      querycoord = -querycoord;
      genomecoord = -genomecoord;
    }

    c1 = querysequence[querycoord];
    c2 = gsequence[genomecoord];
    c2_alt = gsequence_alt[genomecoord];
#ifdef DEBUG14
    c2_single = get_genomic_nt(&c2_alt,genomeoffset+genomecoord,chroffset,chrhigh,watsonp);
    if (c2 != c2_single) {
      abort();
    }
#endif

#ifdef EXTRACT_GENOMICSEG
    assert(c2 == genomesequence[genomecoord]);
#endif

    if (c2 == '*') {
      /* Don't push pairs past end of chromosome */
      debug(printf("Don't push pairs past end of chromosome: genomeoffset %u, genomecoord %u, chroffset %u, chrhigh %u, watsonp %d\n",
		   genomeoffset,genomecoord,chroffset,chrhigh,watsonp));

    } else if (/*querysequenceuc[querycoord]*/c1 == c2) {
      debug(printf("Pushing %d,%d [%d,%d] (%c,%c) - match\n",
		   r,c,queryoffset+querycoord,genomeoffset+genomecoord,c1,c2));
      *nmatches += 1;
      pairs = Pairpool_push(pairs,pairpool,queryoffset+querycoord,genomeoffset+genomecoord,
			    c1,DYNPROG_MATCH_COMP,c2,c2_alt,dynprogindex);

    } else if (consistent_array[(int) c1][(int) c2] == true) {
      debug(printf("Pushing %d,%d [%d,%d] (%c,%c) - ambiguous\n",
		   r,c,queryoffset+querycoord,genomeoffset+genomecoord,c1,c2));
      *nmatches += 1;
      pairs = Pairpool_push(pairs,pairpool,queryoffset+querycoord,genomeoffset+genomecoord,
			    c1,AMBIGUOUS_COMP,c2,c2_alt,dynprogindex);

    } else {
      debug(printf("Pushing %d,%d [%d,%d] (%c,%c) - mismatch\n",
		   r,c,queryoffset+querycoord,genomeoffset+genomecoord,c1,c2));
      *nmismatches += 1;
      pairs = Pairpool_push(pairs,pairpool,queryoffset+querycoord,genomeoffset+genomecoord,
			    c1,MISMATCH_COMP,c2,c2_alt,dynprogindex);
    }

    r--; c--;
    if (r == 0 && c == 0) {
      /* STOP condition.  Do nothing. */

    } else if ((dir = directions_nogap[c][r]) == DIAG) {
      /* Do nothing */

    } else if (dir == HORIZ) {
      dist = 1;
      while (c > 1 && directions_Egap[c][r] != DIAG) {
	dist++;
	c--;
      }
      c--;
      /* dir = directions_nogap[c][r]; */

      debug(printf("H%d: ",dist));
      pairs = add_genomeskip(&add_dashes_p,pairs,r,c+dist,dist,
			     /*genomesequence*/NULL,/*genomesequenceuc*/NULL,
			     queryoffset,genomeoffset,pairpool,revp,chroffset,chrhigh,
			     cdna_direction,watsonp,dynprogindex,/*use_genomicseg_p*/false);
      if (add_dashes_p == true) {
	*nopens += 1;
	*nindels += dist;
      }
      debug(printf("\n"));

    } else {
      /* Must be VERT */
      dist = 1;
      while (r > 1 && directions_Fgap[c][r] != DIAG) {
	dist++;
	r--;
      }
      r--;
      /* dir = directions_nogap[c][r]; */

      debug(printf("V%d: ",dist));
      pairs = add_queryskip(pairs,r+dist,c,dist,querysequence,
			    queryoffset,genomeoffset,pairpool,revp,
			    dynprogindex);
      *nopens += 1;
      *nindels += dist;
      debug(printf("\n"));

    }
  }
  return pairs;
}



/* revp means both rev1p and rev2p, which must have equal values */
/* Iterative version */
static List_T
traceback_nogaps (List_T pairs, int *nmatches, int *nmismatches,
		  int r, int c, char *querysequence, char *querysequenceuc,
		  char *gsequence, char *gsequence_alt,
		  int queryoffset, int genomeoffset, Pairpool_T pairpool, 
#ifdef DEBUG14
		  Univcoord_T chroffset, Univcoord_T chrhigh, bool watsonp,
#endif
		  bool revp, int dynprogindex) {
  char c1, c2, c2_alt;
  int querycoord, genomecoord;
#ifdef DEBUG14
  char c2_single;
#endif

  debug11(printf("Starting traceback_nogaps at r=%d,c=%d (roffset=%d, goffset=%d), revp %d\n",
		 r,c,queryoffset,genomeoffset,revp));

  /* printf("genome sequence is %s\n",genomesequence); */
  while (r > 0 && c > 0) {
    querycoord = r-1;
    genomecoord = c-1;
    if (revp == true) {
      querycoord = -querycoord;
      genomecoord = -genomecoord;
    }

    c1 = querysequence[querycoord];
    c2 = gsequence[genomecoord];
    c2_alt = gsequence_alt[genomecoord];

#ifdef DEBUG14
    c2_single = get_genomic_nt(&c2_alt,genomeoffset+genomecoord,chroffset,chrhigh,watsonp);
    if (c2 != c2_single) {
      abort();
    }
#endif

#ifdef EXTRACT_GENOMICSEG
    debug8(printf("genome sequence at %d is %c\n",genomecoord,genomesequence[genomecoord]));
    assert(c2 == genomesequence[genomecoord]);
#endif

    if (c2 == '*') {
      /* Don't push pairs past end of chromosome */

    } else if (/*querysequenceuc[querycoord]*/c1 == c2) {
      debug11(printf("D: Pushing %d,%d [%d,%d] (%c,%c) - match\n",
		   r,c,queryoffset+querycoord,genomeoffset+genomecoord,c1,c2));
      *nmatches += 1;
      pairs = Pairpool_push(pairs,pairpool,queryoffset+querycoord,genomeoffset+genomecoord,
			    c1,DYNPROG_MATCH_COMP,c2,c2_alt,dynprogindex);
      
    } else if (consistent_array[(int) c1][(int) c2] == true) {
      debug11(printf("D: Pushing %d,%d [%d,%d] (%c,%c) - ambiguous\n",
		   r,c,queryoffset+querycoord,genomeoffset+genomecoord,c1,c2));
      *nmatches += 1;
      pairs = Pairpool_push(pairs,pairpool,queryoffset+querycoord,genomeoffset+genomecoord,
			    c1,AMBIGUOUS_COMP,c2,c2_alt,dynprogindex);

    } else {
      debug11(printf("D: Pushing %d,%d [%d,%d] (%c,%c) - mismatch\n",
		   r,c,queryoffset+querycoord,genomeoffset+genomecoord,c1,c2));
      *nmismatches += 1;
      pairs = Pairpool_push(pairs,pairpool,queryoffset+querycoord,genomeoffset+genomecoord,
			    c1,MISMATCH_COMP,c2,c2_alt,dynprogindex);
    }
    r--; c--;
  }

  return pairs;
}


#if defined(HAVE_SSE4_1) || defined(HAVE_SSE2)
static List_T
traceback_local_8 (List_T pairs, int *nmatches, int *nmismatches, int *nopens, int *nindels,
		   Direction8_T **directions_nogap, Direction8_T **directions_Egap, Direction8_T **directions_Fgap,
		   int *r, int *c, int endc, char *querysequence, char *querysequenceuc,
		   char *genomesequence, char *genomesequenceuc, char *genomesequencealt,
		   int queryoffset, int genomeoffset, Pairpool_T pairpool, bool revp,
		   Univcoord_T chroffset, Univcoord_T chrhigh,
		   int cdna_direction, bool watsonp, int dynprogindex) {
  char c1, c2, c2_alt;
  int dist;
  bool add_dashes_p;
  int querycoord, genomecoord;
  Direction8_T dir;

  debug(printf("Starting traceback_local at r=%d,c=%d (roffset=%d, goffset=%d)\n",*r,*c,queryoffset,genomeoffset));

  /* We care only only about genomic coordinate c */

  if (*c <= endc) {
    /* Do nothing */

  } else if ((dir = directions_nogap[*c][*r]) == DIAG) {
    /* Not an indel.  Do nothing. */

  } else if (dir == HORIZ) {
    dist = 1;
    while (*c > 1 && directions_Egap[*c][*r] != DIAG) {
      dist++;
      (*c)--;
    }
    (*c)--;
    /* dir = directions_nogap[*c][*r]; */

    debug(printf("H%d: ",dist));
    pairs = add_genomeskip(&add_dashes_p,pairs,*r,(*c)+dist,dist,genomesequence,genomesequenceuc,
			   queryoffset,genomeoffset,pairpool,revp,chroffset,chrhigh,
			   cdna_direction,watsonp,dynprogindex,/*use_genomicseg_p*/true);
    if (add_dashes_p == true) {
      *nopens += 1;
      *nindels += dist;
    }
    debug(printf("\n"));

  } else {
    /* Must be VERT */
    dist = 1;
    while (*r > 1 && directions_Fgap[*c][*r] != DIAG) {
      dist++;
      (*r)--;
    }
    (*r)--;
    /* dir = directions_nogap[c][r]; */

    debug(printf("V%d: ",dist));
    pairs = add_queryskip(pairs,(*r)+dist,*c,dist,querysequence,
			  queryoffset,genomeoffset,pairpool,revp,
			  dynprogindex);
    *nopens += 1;
    *nindels += dist;
    debug(printf("\n"));
  }

  while (*r > 0 && *c > endc) {
    querycoord = (*r)-1;
    genomecoord = (*c)-1;
    if (revp == true) {
      querycoord = -querycoord;
      genomecoord = -genomecoord;
    }

    c1 = querysequence[querycoord];
    c2 = genomesequence[genomecoord];
    c2_alt = genomesequencealt[genomecoord];

    if (/*querysequenceuc[querycoord]*/c1 == c2) {
      debug(printf("Pushing %d,%d [%d,%d] (%c,%c) - match\n",
		   *r,*c,queryoffset+querycoord,genomeoffset+genomecoord,c1,c2));
      *nmatches += 1;
      pairs = Pairpool_push(pairs,pairpool,queryoffset+querycoord,genomeoffset+genomecoord,
			    c1,DYNPROG_MATCH_COMP,c2,c2_alt,dynprogindex);

    } else if (consistent_array[(int) c1][(int) c2] == true) {
      debug(printf("Pushing %d,%d [%d,%d] (%c,%c) - ambiguous\n",
		   *r,*c,queryoffset+querycoord,genomeoffset+genomecoord,c1,c2));
      *nmatches += 1;
      pairs = Pairpool_push(pairs,pairpool,queryoffset+querycoord,genomeoffset+genomecoord,
			    c1,AMBIGUOUS_COMP,c2,c2_alt,dynprogindex);

    } else {
      debug(printf("Pushing %d,%d [%d,%d] (%c,%c) - mismatch\n",
		   *r,*c,queryoffset+querycoord,genomeoffset+genomecoord,c1,c2));
      *nmismatches += 1;
      pairs = Pairpool_push(pairs,pairpool,queryoffset+querycoord,genomeoffset+genomecoord,
			    c1,MISMATCH_COMP,c2,c2_alt,dynprogindex);
    }

    (*r)--; (*c)--;
    if (*r == 0 && *c == 0) {
      /* STOP condition.  Do nothing. */

    } else if ((dir = directions_nogap[*c][*r]) == DIAG) {
      /* Do nothing */

    } else if (dir == HORIZ) {
      dist = 1;
      while (*c > 1 && directions_Egap[*c][*r] != DIAG) {
	dist++;
	(*c)--;
      }
      (*c)--;
      /* dir = directions_nogap[*c][*r]; */

      debug(printf("H%d: ",dist));
      pairs = add_genomeskip(&add_dashes_p,pairs,*r,(*c)+dist,dist,genomesequence,genomesequenceuc,
			     queryoffset,genomeoffset,pairpool,revp,chroffset,chrhigh,
			     cdna_direction,watsonp,dynprogindex,
			     /*use_genomicseg_p*/true);
      if (add_dashes_p == true) {
	*nopens += 1;
	*nindels += dist;
      }
      debug(printf("\n"));

    } else {
      /* Must be VERT */
      dist = 1;
      while (*r > 1 && directions_Fgap[*c][*r] != DIAG) {
	dist++;
	(*r)--;
      }
      (*r)--;
      /* dir = directions_nogap[*c][*r]; */

      debug(printf("V%d: ",dist));
      pairs = add_queryskip(pairs,(*r)+dist,*c,dist,querysequence,
			    queryoffset,genomeoffset,pairpool,revp,
			    dynprogindex);
      *nopens += 1;
      *nindels += dist;
      debug(printf("\n"));

    }
  }

  return pairs;
}
#endif


static List_T
traceback_local (List_T pairs, int *nmatches, int *nmismatches, int *nopens, int *nindels,
#ifdef HAVE_SSE2
		 Direction16_T **directions_nogap, Direction16_T **directions_Egap, Direction16_T **directions_Fgap,
#else
		 Direction32_T **directions_nogap, Direction32_T **directions_Egap, Direction32_T **directions_Fgap,
#endif
		 int *r, int *c, int endc, char *querysequence, char *querysequenceuc,
		 char *genomesequence, char *genomesequenceuc, char *genomesequencealt,
		 int queryoffset, int genomeoffset, Pairpool_T pairpool, bool revp,
		 Univcoord_T chroffset, Univcoord_T chrhigh,
		 int cdna_direction, bool watsonp, int dynprogindex) {
  char c1, c2, c2_alt;
  int dist;
  bool add_dashes_p;
  int querycoord, genomecoord;
#ifdef HAVE_SSE2
  Direction16_T dir;
#else
  Direction32_T dir;
#endif

  debug(printf("Starting traceback_local at r=%d,c=%d (roffset=%d, goffset=%d)\n",*r,*c,queryoffset,genomeoffset));

  /* We care only only about genomic coordinate c */

  if (*c <= endc) {
    /* Do nothing */

  } else if ((dir = directions_nogap[*c][*r]) == DIAG) {
    /* Not an indel.  Do nothing. */

  } else if (dir == HORIZ) {
    dist = 1;
    while (*c > 1 && directions_Egap[*c][*r] != DIAG) {
      dist++;
      (*c)--;
    }
    (*c)--;
    /* dir = directions_nogap[*c][*r]; */

    debug(printf("H%d: ",dist));
    pairs = add_genomeskip(&add_dashes_p,pairs,*r,(*c)+dist,dist,genomesequence,genomesequenceuc,
			   queryoffset,genomeoffset,pairpool,revp,chroffset,chrhigh,
			   cdna_direction,watsonp,dynprogindex,/*use_genomicseg_p*/true);
    if (add_dashes_p == true) {
      *nopens += 1;
      *nindels += dist;
    }
    debug(printf("\n"));

  } else {
    /* Must be VERT */
    dist = 1;
    while (*r > 1 && directions_Fgap[*c][*r] != DIAG) {
      dist++;
      (*r)--;
    }
    (*r)--;
    /* dir = directions_nogap[c][r]; */

    debug(printf("V%d: ",dist));
    pairs = add_queryskip(pairs,(*r)+dist,*c,dist,querysequence,
			  queryoffset,genomeoffset,pairpool,revp,
			  dynprogindex);
    *nopens += 1;
    *nindels += dist;
    debug(printf("\n"));
  }

  while (*r > 0 && *c > endc) {
    querycoord = (*r)-1;
    genomecoord = (*c)-1;
    if (revp == true) {
      querycoord = -querycoord;
      genomecoord = -genomecoord;
    }

    c1 = querysequence[querycoord];
    c2 = genomesequence[genomecoord];
    c2_alt = genomesequencealt[genomecoord];

    if (/*querysequenceuc[querycoord]*/c1 == c2) {
      debug(printf("Pushing %d,%d [%d,%d] (%c,%c) - match\n",
		   *r,*c,queryoffset+querycoord,genomeoffset+genomecoord,c1,c2));
      *nmatches += 1;
      pairs = Pairpool_push(pairs,pairpool,queryoffset+querycoord,genomeoffset+genomecoord,
			    c1,DYNPROG_MATCH_COMP,c2,c2_alt,dynprogindex);

    } else if (consistent_array[(int) c1][(int) c2] == true) {
      debug(printf("Pushing %d,%d [%d,%d] (%c,%c) - ambiguous\n",
		   *r,*c,queryoffset+querycoord,genomeoffset+genomecoord,c1,c2));
      *nmatches += 1;
      pairs = Pairpool_push(pairs,pairpool,queryoffset+querycoord,genomeoffset+genomecoord,
			    c1,AMBIGUOUS_COMP,c2,c2_alt,dynprogindex);

    } else {
      debug(printf("Pushing %d,%d [%d,%d] (%c,%c) - mismatch\n",
		   *r,*c,queryoffset+querycoord,genomeoffset+genomecoord,c1,c2));
      *nmismatches += 1;
      pairs = Pairpool_push(pairs,pairpool,queryoffset+querycoord,genomeoffset+genomecoord,
			    c1,MISMATCH_COMP,c2,c2_alt,dynprogindex);
    }

    (*r)--; (*c)--;
    if (*r == 0 && *c == 0) {
      /* STOP condition.  Do nothing. */

    } else if ((dir = directions_nogap[*c][*r]) == DIAG) {
      /* Do nothing */

    } else if (dir == HORIZ) {
      dist = 1;
      while (*c > 1 && directions_Egap[*c][*r] != DIAG) {
	dist++;
	(*c)--;
      }
      (*c)--;
      /* dir = directions_nogap[*c][*r]; */

      debug(printf("H%d: ",dist));
      pairs = add_genomeskip(&add_dashes_p,pairs,*r,(*c)+dist,dist,genomesequence,genomesequenceuc,
			     queryoffset,genomeoffset,pairpool,revp,chroffset,chrhigh,
			     cdna_direction,watsonp,dynprogindex,
			     /*use_genomicseg_p*/true);
      if (add_dashes_p == true) {
	*nopens += 1;
	*nindels += dist;
      }
      debug(printf("\n"));

    } else {
      /* Must be VERT */
      dist = 1;
      while (*r > 1 && directions_Fgap[*c][*r] != DIAG) {
	dist++;
	(*r)--;
      }
      (*r)--;
      /* dir = directions_nogap[*c][*r]; */

      debug(printf("V%d: ",dist));
      pairs = add_queryskip(pairs,(*r)+dist,*c,dist,querysequence,
			    queryoffset,genomeoffset,pairpool,revp,
			    dynprogindex);
      *nopens += 1;
      *nindels += dist;
      debug(printf("\n"));

    }
  }

  return pairs;
}




#if 0
/* revp means both rev1p and rev2p, which must have equal values */
/* Iterative version */
static List_T
traceback_local_nogaps (List_T pairs, int *nmatches, int *nmismatches, int r, int c, int endc,
			char *querysequence, char *querysequenceuc, char *genomesequence, char *genomesequenceuc,
			int queryoffset, int genomeoffset, Pairpool_T pairpool, 
			bool revp, int dynprogindex) {
  char c1, c2;
  int querycoord, genomecoord;

  debug11(printf("Starting traceback_local_nogaps at r=%d,c=%d (roffset=%d, goffset=%d), revp %d\n",
		 r,c,queryoffset,genomeoffset,revp));

  /* We care only about the genomic coordinate c */
  while (c > endc) {
    querycoord = r-1;
    genomecoord = c-1;
    if (revp == true) {
      querycoord = -querycoord;
      genomecoord = -genomecoord;
    }

    c1 = querysequence[querycoord];
    c2 = genomesequence[genomecoord];
    c2_alt = genomesequencealt[genomecoord];

    if (querysequenceuc[querycoord] == genomesequenceuc[genomecoord]) {
      debug11(printf("D: Pushing %d,%d [%d,%d] (%c,%c) - match\n",
		   r,c,queryoffset+querycoord,genomeoffset+genomecoord,c1,c2));
      *nmatches += 1;
      pairs = Pairpool_push(pairs,pairpool,queryoffset+querycoord,genomeoffset+genomecoord,
			    c1,DYNPROG_MATCH_COMP,c2,c2_alt,dynprogindex);
      
    } else if (consistent_array[(int) c1][(int) c2] == true) {
      debug11(printf("D: Pushing %d,%d [%d,%d] (%c,%c) - ambiguous\n",
		   r,c,queryoffset+querycoord,genomeoffset+genomecoord,c1,c2));
      *nmatches += 1;
      pairs = Pairpool_push(pairs,pairpool,queryoffset+querycoord,genomeoffset+genomecoord,
			    c1,AMBIGUOUS_COMP,c2,c2_alt,dynprogindex);

    } else {
      debug11(printf("D: Pushing %d,%d [%d,%d] (%c,%c) - mismatch\n",
		   r,c,queryoffset+querycoord,genomeoffset+genomecoord,c1,c2));
      *nmismatches += 1;
      pairs = Pairpool_push(pairs,pairpool,queryoffset+querycoord,genomeoffset+genomecoord,
			    c1,MISMATCH_COMP,c2,c2_alt,dynprogindex);
    }
    r--; c--;
  }

  return pairs;
}
#endif


#if defined(HAVE_SSE4_1) || defined(HAVE_SSE2)
/* Columns are always genomic.  Rows are always query.  Bridging across common columns */
static void
bridge_cdna_gap_8 (int *finalscore, int *bestcL, int *bestcR, int *bestrL, int *bestrR,
		   Score8_T **matrixL, Score8_T **matrixR,
		   int glength, int rlengthL, int rlengthR, int extraband_paired,
		   int open, int extend, int leftoffset, int rightoffset, bool jump_late_p) {
  int bestscore = NEG_INFINITY_16, score, scoreL, scoreR, pen, end_reward = 0;
  int rL, rR, cL, cR;
  int lbandL, ubandL, rloL, rhighL;
  int lbandR, ubandR, rloR, rhighR;

  /* Perform computations */
  lbandL = rlengthL - glength + extraband_paired;
  ubandL = extraband_paired;

  lbandR = rlengthR - glength + extraband_paired;
  ubandR = extraband_paired;

  /* Need a double loop on rows here, in contrast with a single loop
     for introns, because we allow a genomic insertion that doesn't
     match the cDNA.  So we need to add a penalty for a genomic
     insertion */

  if (jump_late_p) {
    for (cL = 1; cL < glength; cL++) {

      /* Note: opening penalty is added at the bottom of the loop */
      for (cR = glength-cL, pen = 0; cR >= 0; cR--, pen += extend) {
	/* debug3(printf("\nAt row %d on left and %d on right\n",cL,cR)); */
	if ((rloL = cL - ubandL) < 1) {
	  rloL = 1;
	}
	if ((rhighL = cL + lbandL) > rlengthL-1) {
	  rhighL = rlengthL-1;
	}

	if ((rloR = cR - ubandR) < 1) {
	  rloR = 1;
	}
	if ((rhighR = cR + lbandR) > rlengthR-1) {
	  rhighR = rlengthR-1;
	}

	for (rL = rloL; rL <= rhighL; rL++) {
	  scoreL = (int) matrixL[cL][rL];
	
	  /* Disallow leftoffset + rL >= rightoffset - rR, or rR >= rightoffset - leftoffset - rL */
	  debug3(printf("  Disallowing rR to be >= %d\n",rightoffset-leftoffset-rL));
	  for (rR = rloR; rR <= rhighR && rR < rightoffset-leftoffset-rL; rR++) {
	    scoreR = (int) matrixR[cR][rR];

	    if ((score = scoreL + scoreR + pen + end_reward) >= bestscore) {  /* Use >= for jump late */
	      debug3(printf("At %d left to %d right, score is (%d)+(%d)+(%d)+(%d) = %d (bestscore)\n",
			    rL,rR,scoreL,scoreR,pen,end_reward,score));

	      bestscore = score;
	      *bestcL = cL;
	      *bestcR = cR;
	      *bestrL = rL;
	      *bestrR = rR;

	    } else {
	      debug3(printf("At %d left to %d right, score is (%d)+(%d)+(%d) = %d\n",
			    rL,rR,scoreL,scoreR,pen,scoreL+scoreR+pen));
	    }
	  }
	}
	pen = open - extend;	/* Subtract extend to compensate for
                                   its addition in the for loop */
      }
    }

  } else {
    /* Do not jump late */
    for (cL = 1; cL < glength; cL++) {

      /* Note: opening penalty is added at the bottom of the loop */
      for (cR = glength-cL, pen = 0; cR >= 0; cR--, pen += extend) {
	/* debug3(printf("\nAt row %d on left and %d on right\n",cL,cR)); */
	if ((rloL = cL - ubandL) < 1) {
	  rloL = 1;
	}
	if ((rhighL = cL + lbandL) > rlengthL-1) {
	  rhighL = rlengthL-1;
	}

	if ((rloR = cR - ubandR) < 1) {
	  rloR = 1;
	}
	if ((rhighR = cR + lbandR) > rlengthR-1) {
	  rhighR = rlengthR-1;
	}

	for (rL = rloL; rL <= rhighL; rL++) {
	  scoreL = (int) matrixL[cL][rL];
	
	  /* Disallow leftoffset + rL >= rightoffset - rR, or rR >= rightoffset - leftoffset - rL */
	  debug3(printf("  Disallowing rR to be >= %d\n",rightoffset-leftoffset-rL));
	  for (rR = rloR; rR <= rhighR && rR < rightoffset-leftoffset-rL; rR++) {
	    scoreR = (int) matrixR[cR][rR];

	    if ((score = scoreL + scoreR + pen + end_reward) > bestscore) {  /* Use > for jump early */
	      debug3(printf("At %d left to %d right, score is (%d)+(%d)+(%d)+(%d) = %d (bestscore)\n",
			    rL,rR,scoreL,scoreR,pen,end_reward,score));

	      bestscore = score;
	      *bestcL = cL;
	      *bestcR = cR;
	      *bestrL = rL;
	      *bestrR = rR;

	    } else {
	      debug3(printf("At %d left to %d right, score is (%d)+(%d)+(%d) = %d\n",
			    rL,rR,scoreL,scoreR,pen,scoreL+scoreR+pen));
	    }
	  }
	}
	pen = open - extend;	/* Subtract extend to compensate for
                                   its addition in the for loop */
      }
    }
  }
      
  *finalscore = (int) bestscore;
  debug3(printf("Returning final score of %d at (%d,%d) left to (%d,%d) right\n",
		*finalscore,*bestcL,*bestrL,*bestcR,*bestrR));

  return;
}
#endif


/* Columns are always genomic.  Rows are always query.  Bridging across common columns */
static void
bridge_cdna_gap (int *finalscore, int *bestcL, int *bestcR, int *bestrL, int *bestrR,
#ifdef HAVE_SSE2
		 Score16_T **matrixL, Score16_T **matrixR,
#else
		 Score32_T **matrixL, Score32_T **matrixR,
#endif
		 int glength, int rlengthL, int rlengthR, int extraband_paired,
		 int open, int extend, int leftoffset, int rightoffset, bool jump_late_p) {
  int bestscore = NEG_INFINITY_32, score, scoreL, scoreR, pen, end_reward = 0;
  int rL, rR, cL, cR;
  int lbandL, ubandL, rloL, rhighL;
  int lbandR, ubandR, rloR, rhighR;

  /* Perform computations */
  lbandL = rlengthL - glength + extraband_paired;
  ubandL = extraband_paired;

  lbandR = rlengthR - glength + extraband_paired;
  ubandR = extraband_paired;

  /* Need a double loop on rows here, in contrast with a single loop
     for introns, because we allow a genomic insertion that doesn't
     match the cDNA.  So we need to add a penalty for a genomic
     insertion */

  if (jump_late_p) {
    for (cL = 1; cL < glength; cL++) {

      /* Note: opening penalty is added at the bottom of the loop */
      for (cR = glength-cL, pen = 0; cR >= 0; cR--, pen += extend) {
	/* debug3(printf("\nAt row %d on left and %d on right\n",cL,cR)); */
	if ((rloL = cL - ubandL) < 1) {
	  rloL = 1;
	}
	if ((rhighL = cL + lbandL) > rlengthL-1) {
	  rhighL = rlengthL-1;
	}

	if ((rloR = cR - ubandR) < 1) {
	  rloR = 1;
	}
	if ((rhighR = cR + lbandR) > rlengthR-1) {
	  rhighR = rlengthR-1;
	}

	for (rL = rloL; rL <= rhighL; rL++) {
	  scoreL = (int) matrixL[cL][rL];
	
	  /* Disallow leftoffset + rL >= rightoffset - rR, or rR >= rightoffset - leftoffset - rL */
	  debug3(printf("  Disallowing rR to be >= %d\n",rightoffset-leftoffset-rL));
	  for (rR = rloR; rR <= rhighR && rR < rightoffset-leftoffset-rL; rR++) {
	    scoreR = (int) matrixR[cR][rR];

	    if ((score = scoreL + scoreR + pen + end_reward) >= bestscore) {  /* Use >= for jump late */
	      debug3(printf("At %d left to %d right, score is (%d)+(%d)+(%d)+(%d) = %d (bestscore)\n",
			    rL,rR,scoreL,scoreR,pen,end_reward,score));

	      bestscore = score;
	      *bestcL = cL;
	      *bestcR = cR;
	      *bestrL = rL;
	      *bestrR = rR;

	    } else {
	      debug3(printf("At %d left to %d right, score is (%d)+(%d)+(%d) = %d\n",
			    rL,rR,scoreL,scoreR,pen,scoreL+scoreR+pen));
	    }
	  }
	}
	pen = open - extend;	/* Subtract extend to compensate for
                                   its addition in the for loop */
      }
    }

  } else {
    /* Do not jump late */
    for (cL = 1; cL < glength; cL++) {

      /* Note: opening penalty is added at the bottom of the loop */
      for (cR = glength-cL, pen = 0; cR >= 0; cR--, pen += extend) {
	/* debug3(printf("\nAt row %d on left and %d on right\n",cL,cR)); */
	if ((rloL = cL - ubandL) < 1) {
	  rloL = 1;
	}
	if ((rhighL = cL + lbandL) > rlengthL-1) {
	  rhighL = rlengthL-1;
	}

	if ((rloR = cR - ubandR) < 1) {
	  rloR = 1;
	}
	if ((rhighR = cR + lbandR) > rlengthR-1) {
	  rhighR = rlengthR-1;
	}

	for (rL = rloL; rL <= rhighL; rL++) {
	  scoreL = (int) matrixL[cL][rL];
	
	  /* Disallow leftoffset + rL >= rightoffset - rR, or rR >= rightoffset - leftoffset - rL */
	  debug3(printf("  Disallowing rR to be >= %d\n",rightoffset-leftoffset-rL));
	  for (rR = rloR; rR <= rhighR && rR < rightoffset-leftoffset-rL; rR++) {
	    scoreR = (int) matrixR[cR][rR];

	    if ((score = scoreL + scoreR + pen + end_reward) > bestscore) {  /* Use > for jump early */
	      debug3(printf("At %d left to %d right, score is (%d)+(%d)+(%d)+(%d) = %d (bestscore)\n",
			    rL,rR,scoreL,scoreR,pen,end_reward,score));

	      bestscore = score;
	      *bestcL = cL;
	      *bestcR = cR;
	      *bestrL = rL;
	      *bestrR = rR;

	    } else {
	      debug3(printf("At %d left to %d right, score is (%d)+(%d)+(%d) = %d\n",
			    rL,rR,scoreL,scoreR,pen,scoreL+scoreR+pen));
	    }
	  }
	}
	pen = open - extend;	/* Subtract extend to compensate for
                                   its addition in the for loop */
      }
    }
  }
      
  *finalscore = (int) bestscore;
  debug3(printf("Returning final score of %d at (%d,%d) left to (%d,%d) right\n",
		*finalscore,*bestcL,*bestrL,*bestcR,*bestrR));

  return;
}


static int
intron_score (int *introntype, int leftdi, int rightdi, int cdna_direction, int canonical_reward, 
	      bool finalp) {
  int scoreI;

#ifdef PMAP
  if ((*introntype = leftdi & rightdi) == NONINTRON) {
    scoreI = 0.0;
  } else {
    switch (*introntype) {
    case GTAG_FWD: scoreI = canonical_reward; break;
    case GCAG_FWD: scoreI = finalp == true ? FINAL_GCAG_INTRON : GCAG_INTRON; break;
    case ATAC_FWD: scoreI = finalp == true ? FINAL_ATAC_INTRON : ATAC_INTRON; break;
    default: *introntype = NONINTRON; scoreI = 0.0;
    }
  }
#else
  if ((*introntype = leftdi & rightdi) == NONINTRON) {
    scoreI = 0.0;
  } else if (cdna_direction > 0) {
    switch (*introntype) {
    case GTAG_FWD: scoreI = canonical_reward; break;
    case GCAG_FWD: scoreI = finalp == true ? FINAL_GCAG_INTRON : GCAG_INTRON; break;
    case ATAC_FWD: scoreI = finalp == true ? FINAL_ATAC_INTRON : ATAC_INTRON; break;
    default: *introntype = NONINTRON; scoreI = 0.0;
    }
  } else if (cdna_direction < 0) {
    switch (*introntype) {
    case GTAG_REV: scoreI = canonical_reward; break;
    case GCAG_REV: scoreI = finalp == true ? FINAL_GCAG_INTRON : GCAG_INTRON; break;
    case ATAC_REV: scoreI = finalp == true ? FINAL_ATAC_INTRON : ATAC_INTRON; break;
    default: *introntype = NONINTRON; scoreI = 0.0;
    }
  } else {
    switch (*introntype) {
    case GTAG_FWD: case GTAG_REV: scoreI = canonical_reward; break;
    case GCAG_FWD: case GCAG_REV: scoreI = finalp == true ? FINAL_GCAG_INTRON : GCAG_INTRON; break;
    case ATAC_FWD: case ATAC_REV: scoreI = finalp == true ? FINAL_ATAC_INTRON : ATAC_INTRON; break;
    default: *introntype = NONINTRON; scoreI = 0.0;
    }
  }
#endif

  return scoreI;
}


static void
get_splicesite_probs (double *left_prob, double *right_prob, int cL, int cR,
		      int *left_known, int *right_known, Univcoord_T leftoffset, Univcoord_T rightoffset,
		      Univcoord_T chroffset, Univcoord_T chrhigh, int cdna_direction, bool watsonp) {
  Univcoord_T splicesitepos;
  
  if (left_known[cL] > 0) {
    debug9(printf("left position is known, so prob is 1.0\n"));
    *left_prob = 1.0;

  } else if (watsonp == true) {
    splicesitepos = chroffset + leftoffset + cL;
    if (cdna_direction > 0) {
      *left_prob = Maxent_hr_donor_prob(splicesitepos,chroffset);
      debug9(printf("1. donor splicesitepos is %u, prob %f, known %d\n",
		    splicesitepos,*left_prob,left_known[cL]));

    } else {
      *left_prob = Maxent_hr_antiacceptor_prob(splicesitepos,chroffset);
      debug9(printf("2. antiacceptor splicesitepos is %u, prob %f, known %d\n",
		    splicesitepos,*left_prob,left_known[cL]));

    }
  } else {
    splicesitepos = chrhigh - leftoffset - cL + 1;
    if (cdna_direction > 0) {
      *left_prob = Maxent_hr_antidonor_prob(splicesitepos,chroffset);
      debug9(printf("3. antidonor splicesitepos is %u, prob %f, known %d\n",
		    splicesitepos,*left_prob,left_known[cL]));

    } else {
      *left_prob = Maxent_hr_acceptor_prob(splicesitepos,chroffset);
      debug9(printf("4. acceptor splicesitepos is %u, prob %f, known %d\n",
		    splicesitepos,*left_prob,left_known[cL]));
    }
  }

  if (right_known[cR] > 0) {
    debug9(printf("right position is known, so prob is 1.0\n"));
    *right_prob = 1.0;

  } else if (watsonp == true) {
    splicesitepos = chroffset + rightoffset - cR + 1;
    if (cdna_direction > 0) {
      *right_prob = Maxent_hr_acceptor_prob(splicesitepos,chroffset);
      debug9(printf("5. acceptor splicesitepos is %u, prob %f, known %d\n",
		    splicesitepos,*right_prob,right_known[cR]));
    } else {
      *right_prob = Maxent_hr_antidonor_prob(splicesitepos,chroffset);
      debug9(printf("6. antidonor splicesitepos is %u, prob %f, known %d\n",
		    splicesitepos,*right_prob,right_known[cR]));

    }
  } else {
    splicesitepos = chrhigh - rightoffset + cR;
    if (cdna_direction > 0) {
      *right_prob = Maxent_hr_antiacceptor_prob(splicesitepos,chroffset);
      debug9(printf("7. antiacceptor splicesitepos is %u, prob %f, known %d\n",
		    splicesitepos,*right_prob,right_known[cR]));

    } else {
      *right_prob = Maxent_hr_donor_prob(splicesitepos,chroffset);
      debug9(printf("8. donor splicesitepos is %u, prob %f, known %d\n",
		    splicesitepos,*right_prob,right_known[cR]));
    }
  }

  return;
}


static void
get_known_splicesites (int *left_known, int *right_known, int glengthL, int glengthR,
		       int leftoffset, int rightoffset, int cdna_direction, bool watsonp,
		       Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh) {

  int *matches, nmatches, i;
  Univcoord_T splicesitepos;

  if (splicing_iit != NULL && donor_typeint >= 0 && acceptor_typeint >= 0) {
    /* Handle known splicing, splice site level */
    if (watsonp == true) {
      if (cdna_direction > 0) {
	/* splicesitepos = leftoffset + cL;  cL = 0 to < glengthL - 1 */
	matches = IIT_get_typed_signed_with_divno(&nmatches,splicing_iit,splicing_divint_crosstable[chrnum],
						  leftoffset+1,leftoffset+glengthL-2,donor_typeint,/*sign*/+1,/*sortp*/false);
	for (i = 0; i < nmatches; i++) {
	  splicesitepos = IIT_interval_low(splicing_iit,matches[i]);
	  debug5(printf("1. Found known donor at %u\n",splicesitepos));
	  left_known[splicesitepos - leftoffset] = KNOWN_SPLICESITE_REWARD;
	}
	FREE(matches);

	/* splicesitepos = rightoffset - cR + 1; cR = 0 to < glengthR - 1 */
	matches = IIT_get_typed_signed_with_divno(&nmatches,splicing_iit,splicing_divint_crosstable[chrnum],
						  rightoffset-glengthR+4,rightoffset+1,acceptor_typeint,/*sign*/+1,/*sortp*/false);
	for (i = 0; i < nmatches; i++) {
	  splicesitepos = IIT_interval_low(splicing_iit,matches[i]);
	  debug5(printf("2. Found known acceptor at %u\n",splicesitepos));
	  right_known[rightoffset - splicesitepos + 1] = KNOWN_SPLICESITE_REWARD;
	}
	FREE(matches);

      } else {
	/* splicesitepos = leftoffset + cL;  cL = 0 to < glengthL - 1 */
	matches = IIT_get_typed_signed_with_divno(&nmatches,splicing_iit,splicing_divint_crosstable[chrnum],
						  leftoffset+1,leftoffset+glengthL-2,acceptor_typeint,/*sign*/-1,/*sortp*/false);
	for (i = 0; i < nmatches; i++) {
	  splicesitepos = IIT_interval_low(splicing_iit,matches[i]);
	  debug5(printf("3. Found known antiacceptor at %u\n",splicesitepos));
	  left_known[splicesitepos - leftoffset] = KNOWN_SPLICESITE_REWARD;
	}
	FREE(matches);

	/* splicesitepos = rightoffset - cR + 1; cR = 0 to < glengthR - 1 */
	matches = IIT_get_typed_signed_with_divno(&nmatches,splicing_iit,splicing_divint_crosstable[chrnum],
						  rightoffset-glengthR+4,rightoffset+1,donor_typeint,/*sign*/-1,/*sortp*/false);
	for (i = 0; i < nmatches; i++) {
	  splicesitepos = IIT_interval_low(splicing_iit,matches[i]);
	  debug5(printf("4. Found known antidonor at %u\n",splicesitepos));
	  right_known[rightoffset - splicesitepos + 1] = KNOWN_SPLICESITE_REWARD;
	}
	FREE(matches);

      }

    } else {
      if (cdna_direction > 0) {
	/* splicesitepos = (chrhigh - chroffset) - leftoffset - cL + 1; cL = 0 to < glengthL - 1 */
	matches = IIT_get_typed_signed_with_divno(&nmatches,splicing_iit,splicing_divint_crosstable[chrnum],
						  (chrhigh - chroffset) - leftoffset - glengthL + 4,
						  (chrhigh - chroffset) - leftoffset + 1,
						  donor_typeint,/*sign*/-1,/*sortp*/false);
	for (i = 0; i < nmatches; i++) {
	  splicesitepos = IIT_interval_low(splicing_iit,matches[i]);
	  debug5(printf("5. Found known antidonor at %u\n",splicesitepos));
	  left_known[(chrhigh - chroffset) - leftoffset - splicesitepos + 1] = KNOWN_SPLICESITE_REWARD;
	}
	FREE(matches);

	/* splicesitepos = (chrhigh - chroffset) - rightoffset + cR; cR = 0 to < glengthR - 1 */
	matches = IIT_get_typed_signed_with_divno(&nmatches,splicing_iit,splicing_divint_crosstable[chrnum],
						  (chrhigh - chroffset) - rightoffset + 1,
						  (chrhigh - chroffset) - rightoffset + glengthR - 2,
						  acceptor_typeint,/*sign*/-1,/*sortp*/false);
	for (i = 0; i < nmatches; i++) {
	  splicesitepos = IIT_interval_low(splicing_iit,matches[i]);
	  debug5(printf("6. Found known antiacceptor at %u\n",splicesitepos));
	  right_known[splicesitepos - (chrhigh - chroffset) + rightoffset] = KNOWN_SPLICESITE_REWARD;
	}
	FREE(matches);

      } else {
	/* splicesitepos = (chrhigh - chroffset) - leftoffset - cL + 1; cL = 0 to < glengthL - 1 */
	matches = IIT_get_typed_signed_with_divno(&nmatches,splicing_iit,splicing_divint_crosstable[chrnum],
						  (chrhigh - chroffset) - leftoffset - glengthL + 4,
						  (chrhigh - chroffset) - leftoffset + 1,
						  acceptor_typeint,/*sign*/+1,/*sortp*/false);
	for (i = 0; i < nmatches; i++) {
	  splicesitepos = IIT_interval_low(splicing_iit,matches[i]);
	  debug5(printf("7. Found known acceptor at %u\n",splicesitepos));
	  left_known[(chrhigh - chroffset) - leftoffset - splicesitepos + 1] = KNOWN_SPLICESITE_REWARD;
	}
	FREE(matches);

	/* splicesitepos = (chrhigh - chroffset) - rightoffset + cR; cR = 0 to < glengthR - 1 */
	matches = IIT_get_typed_signed_with_divno(&nmatches,splicing_iit,splicing_divint_crosstable[chrnum],
						  (chrhigh - chroffset) - rightoffset + 1,
						  (chrhigh - chroffset) - rightoffset + glengthR - 2,
						  donor_typeint,/*sign*/+1,/*sortp*/false);
	for (i = 0; i < nmatches; i++) {
	  splicesitepos = IIT_interval_low(splicing_iit,matches[i]);
	  debug5(printf("8. Found known donor at %u\n",splicesitepos));
	  right_known[splicesitepos - (chrhigh - chroffset) + rightoffset] = KNOWN_SPLICESITE_REWARD;
	}
	FREE(matches);
	  
      }
    }

  } else if (splicing_iit != NULL) {
    /* Handle known splicing, intron level */
    if (watsonp == true) {
      if (cdna_direction > 0) {
	/* splicesitepos = leftoffset + cL; cL = 0 to < glengthL - 1 */
	matches = IIT_get_lows_signed(&nmatches,splicing_iit,splicing_divint_crosstable[chrnum],
				      leftoffset,leftoffset+glengthL-2,/*sign*/+1);
	for (i = 0; i < nmatches; i++) {
	  splicesitepos = IIT_interval_low(splicing_iit,matches[i]);
	  debug5(printf("1. Found known donor at %u\n",splicesitepos));
	  left_known[splicesitepos - leftoffset] = KNOWN_SPLICESITE_REWARD;
	}
	FREE(matches);

	/* splicesitepos+1U = rightoffset - cR + 2; cR = 0 to < glengthR - 1 */
	matches = IIT_get_highs_signed(&nmatches,splicing_iit,splicing_divint_crosstable[chrnum],
				       rightoffset-glengthR+4,rightoffset+2,/*sign*/+1);
	for (i = 0; i < nmatches; i++) {
	  splicesitepos = IIT_interval_high(splicing_iit,matches[i]);
	  debug5(printf("2. Found known acceptor at %u\n",splicesitepos));
	  right_known[rightoffset - splicesitepos + 2] = KNOWN_SPLICESITE_REWARD;
	}
	FREE(matches);

      } else {
	/* splicesitepos = leftoffset + cL; cL = 0 to < glengthL - 1 */
	matches = IIT_get_lows_signed(&nmatches,splicing_iit,splicing_divint_crosstable[chrnum],
				      leftoffset,leftoffset+glengthL-2,/*sign*/-1);
	for (i = 0; i < nmatches; i++) {
	  splicesitepos = IIT_interval_low(splicing_iit,matches[i]);
	  debug5(printf("3. Found known antiacceptor at %u\n",splicesitepos));
	  left_known[splicesitepos - leftoffset] = KNOWN_SPLICESITE_REWARD;
	}
	FREE(matches);

	/* splicesitepos+1U = rightoffset - cR + 2; cR = 0 to < glengthR - 1 */
	matches = IIT_get_highs_signed(&nmatches,splicing_iit,splicing_divint_crosstable[chrnum],
				       rightoffset-glengthR+4,rightoffset+2,/*sign*/-1);
	for (i = 0; i < nmatches; i++) {
	  splicesitepos = IIT_interval_high(splicing_iit,matches[i]);
	  debug5(printf("4. Found known antidonor at %u\n",splicesitepos));
	  right_known[rightoffset - splicesitepos + 2] = KNOWN_SPLICESITE_REWARD;
	}
	FREE(matches);

      }
    } else {
      if (cdna_direction > 0) {
	/* splicesitepos+1U = (chrhigh - chroffset) - leftoffset - cL + 2; cL = 0 to < glengthL - 1 */
	matches = IIT_get_highs_signed(&nmatches,splicing_iit,splicing_divint_crosstable[chrnum],
				       (chrhigh - chroffset) - leftoffset - glengthL + 4,
				       (chrhigh - chroffset) - leftoffset + 2,/*sign*/-1);
	for (i = 0; i < nmatches; i++) {
	  splicesitepos = IIT_interval_high(splicing_iit,matches[i]);
	  debug5(printf("5. Found known antidonor at %u\n",splicesitepos));
	  left_known[(chrhigh - chroffset) - leftoffset - splicesitepos + 2] = KNOWN_SPLICESITE_REWARD;
	}
	FREE(matches);

	/* splicesitepos = (chrhigh - chroffset) - rightoffset + cR; cR = 0 to < glengthR - 1 */
	matches = IIT_get_lows_signed(&nmatches,splicing_iit,splicing_divint_crosstable[chrnum],
				      (chrhigh - chroffset) - rightoffset,
				      (chrhigh - chroffset) - rightoffset + glengthR - 2,/*sign*/-1);
	for (i = 0; i < nmatches; i++) {
	  splicesitepos = IIT_interval_low(splicing_iit,matches[i]);
	  debug5(printf("6. Found known antiacceptor at %u\n",splicesitepos));
	  right_known[splicesitepos - (chrhigh - chroffset) + rightoffset] = KNOWN_SPLICESITE_REWARD;
	}
	FREE(matches);

      } else {
	/* splicesitepos+1U = (chrhigh - chroffset) - leftoffset - cL + 2; cL = 0 to < glengthL - 1 */
	matches = IIT_get_highs_signed(&nmatches,splicing_iit,splicing_divint_crosstable[chrnum],
				       (chrhigh - chroffset) - leftoffset - glengthL + 4,
				       (chrhigh - chroffset) - leftoffset + 2,/*sign*/+1);
	for (i = 0; i < nmatches; i++) {
	  splicesitepos = IIT_interval_high(splicing_iit,matches[i]);
	  debug5(printf("7. Found known acceptor at %u\n",splicesitepos));
	  left_known[(chrhigh - chroffset) - leftoffset - splicesitepos + 2] = KNOWN_SPLICESITE_REWARD;
	}
	FREE(matches);

	/* splicesitepos = (chrhigh - chroffset) - rightoffset + cR; cR = 0 to < glengthR - 1 */
	matches = IIT_get_lows_signed(&nmatches,splicing_iit,splicing_divint_crosstable[chrnum],
				      (chrhigh - chroffset) - rightoffset,
				      (chrhigh - chroffset) - rightoffset + glengthR - 2,/*sign*/+1);
	for (i = 0; i < nmatches; i++) {
	  splicesitepos = IIT_interval_low(splicing_iit,matches[i]);
	  debug5(printf("8. Found known donor at %u\n",splicesitepos));
	  right_known[splicesitepos - (chrhigh - chroffset) + rightoffset] = KNOWN_SPLICESITE_REWARD;
	}
	FREE(matches);
      }
    }
  }

  return;
}


#if defined(HAVE_SSE4_1) || defined(HAVE_SSE2)
static bool
bridge_intron_gap_8 (int *finalscore, int *bestrL, int *bestrR, int *bestcL, int *bestcR,
		     int *best_introntype, double *left_prob, double *right_prob,
		     Score8_T **matrixL, Score8_T **matrixR,
		     Direction8_T **directionsL_nogap, Direction8_T **directionsR_nogap, 
		     char *gsequenceL, char *gsequenceL_alt, char *rev_gsequenceR, char *rev_gsequenceR_alt,
		     int goffsetL, int rev_goffsetR, int rlength, int glengthL, int glengthR,
		     int cdna_direction, bool watsonp, int extraband_paired, double defect_rate, int canonical_reward,
		     int maxhorizjump, int maxvertjump, int leftoffset, int rightoffset,
		     Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
		     bool halfp, bool finalp, bool use_probabilities_p, int score_threshold,
		     bool jump_late_p) {
  bool result;
  int bestscore = NEG_INFINITY_32, score, scoreL, scoreR, scoreI, bestscoreI = NEG_INFINITY_32;
  int bestscore_with_prob = NEG_INFINITY_32;
  int rL, rR, cL, cR;
  int bestrL_with_prob, bestrR_with_prob, bestcL_with_prob, bestcR_with_prob;
  int lbandL, ubandL, cloL, chighL;
  int lbandR, ubandR, cloR, chighR;
  char left1, left2, right2, right1, left1_alt, left2_alt, right2_alt, right1_alt;
  int *leftdi, *rightdi, introntype;
  int *left_known, *right_known;
  double *left_probabilities, *right_probabilities, probL, probR, probL_trunc, probR_trunc, bestprob, bestprob_trunc, prob;
  Univcoord_T splicesitepos, splicesitepos1, splicesitepos2;
  bool bestp;


  debug(printf("Running bridge_intron_gap_8 with use_probabilities_p %d\n",use_probabilities_p));

  if (glengthL+1 <= 0) {
    fprintf(stderr,"Problem with glengthL = %d\n",glengthL);
    abort();
  }

  if (glengthR+1 <= 0) {
    fprintf(stderr,"Problem with glengthR = %d\n",glengthR);
    abort();
  }

  /* Read dinucleotides */
  leftdi = (int *) CALLOC(glengthL+1,sizeof(int));
  rightdi = (int *) CALLOC(glengthR+1,sizeof(int));

  for (cL = 0; cL < glengthL - 1; cL++) {
    left1 = gsequenceL[cL];
    left1_alt = gsequenceL_alt[cL];
    left2 = gsequenceL[cL+1];
    left2_alt = gsequenceL_alt[cL+1];
    assert(left1 == get_genomic_nt(&left1_alt,goffsetL+cL,chroffset,chrhigh,watsonp));
    assert(left2 == get_genomic_nt(&left2_alt,goffsetL+cL+1,chroffset,chrhigh,watsonp));

    if ((left1 == 'G' || left1_alt == 'G') && (left2 == 'T' || left2_alt == 'T')) {
      leftdi[cL] = LEFT_GT;
    } else if ((left1 == 'G' || left1_alt == 'G') && (left2 == 'C' || left2_alt == 'C')) {
      leftdi[cL] = LEFT_GC;
    } else if ((left1 == 'A' || left1_alt == 'A') && (left2 == 'T' || left2_alt == 'T')) {
      leftdi[cL] = LEFT_AT;
#ifndef PMAP
    } else if ((left1 == 'C' || left1_alt == 'C') && (left2 == 'T' || left2_alt == 'T')) {
      leftdi[cL] = LEFT_CT;
#endif
    } else {
      leftdi[cL] = 0x00;
    }
  }
  leftdi[glengthL-1] = leftdi[glengthL] = 0x00;

  for (cR = 0; cR < glengthR - 1; cR++) {
    right2 = rev_gsequenceR[-cR-1];
    right2_alt = rev_gsequenceR_alt[-cR-1];
    right1 = rev_gsequenceR[-cR];
    right1_alt = rev_gsequenceR_alt[-cR];
    assert(right2 == get_genomic_nt(&right2_alt,rev_goffsetR-cR-1,chroffset,chrhigh,watsonp));
    assert(right1 == get_genomic_nt(&right1_alt,rev_goffsetR-cR,chroffset,chrhigh,watsonp));

    if ((right2 == 'A' || right2_alt == 'A') && (right1 == 'G' || right1_alt == 'G')) {
      rightdi[cR] = RIGHT_AG;
    } else if ((right2 == 'A' || right2_alt == 'A') && (right1 == 'C' || right1_alt == 'C')) {
      rightdi[cR] = RIGHT_AC;
#ifndef PMAP
    } else if ((right2 == 'G' || right2_alt == 'G') && (right1 == 'C' || right1_alt == 'C')) {
      rightdi[cR] = RIGHT_GC;
    } else if ((right2 == 'A' || right2_alt == 'A') && (right1 == 'T' || right1_alt == 'T')) {
      rightdi[cR] = RIGHT_AT;
#endif
    } else {
      rightdi[cR] = 0x00;
    }
  }
  rightdi[glengthR-1] = rightdi[glengthR] = 0x00;

  left_known = (int *) CALLOC(glengthL+1,sizeof(int));
  right_known = (int *) CALLOC(glengthR+1,sizeof(int));
  get_known_splicesites(left_known,right_known,glengthL,glengthR,
			/*leftoffset*/goffsetL,/*rightoffset*/rev_goffsetR,
			cdna_direction,watsonp,chrnum,chroffset,chrhigh);

  /* Perform computations */
#if 1
  /* Allows unlimited indel lengths */
  ubandL = glengthL - rlength + extraband_paired;
  lbandL = extraband_paired;

  ubandR = glengthR - rlength + extraband_paired;
  lbandR = extraband_paired;
#else
  /* Limit indels to 3 bp around splice sites.  Doesn't work on PacBio reads. */
  ubandL = 3;
  lbandL = 3;

  ubandR = 3;
  lbandR = 3;
#endif


  if (novelsplicingp == false && splicing_iit != NULL && (donor_typeint < 0 || acceptor_typeint < 0)) {
    /* Constrain to given introns */
    for (rL = 1, rR = rlength-1; rL < rlength; rL++, rR--) {
      debug3(printf("\nGenomic insert: At row %d on left and %d on right\n",rL,rR));
      if ((cloL = rL - lbandL) < 1) {
	cloL = 1;
      }
      if ((chighL = rL + ubandL) > glengthL-1) {
	chighL = glengthL-1;
      }

      if ((cloR = rR - lbandR) < 1) {
	cloR = 1;
      }
      if ((chighR = rR + ubandR) > glengthR-1) {
	chighR = glengthR-1;
      }

      /* Test indels on left and right */
      for (cL = cloL; cL <= chighL; cL++) {
	/* The following check limits genomic inserts (horizontal) and
	   multiple cDNA inserts (vertical). */
	if (left_known[cL] > 0) {
	  scoreL = (int) matrixL[cL][rL];
	  if (directionsL_nogap[cL][rL] != DIAG) {
	    /* Favor gaps away from intron if possible */
	    scoreL -= 1;
	  }

	  /* Disallow leftoffset + cL >= rightoffset - cR, or cR >= rightoffset - leftoffset - cL */
	  for (cR = cloR; cR <= chighR && cR < rightoffset-leftoffset-cL; cR++) {
	    if (right_known[cR] > 0) {
	      scoreR = (int) matrixR[cR][rR];
	      if (directionsR_nogap[cR][rR] != DIAG) {
		/* Favor gaps away from intron if possible */
		scoreR -= 1;
	      }

	      if ((score = scoreL + scoreR) > bestscore ||
		  (score >= bestscore && jump_late_p)) {  /* Use >= for jump late */
		bestp = false;
		if (watsonp == true) {
		  splicesitepos1 = leftoffset + cL;
		  splicesitepos2 = rightoffset - cR + 1;
		  if (IIT_exists_with_divno_signed(splicing_iit,splicing_divint_crosstable[chrnum],
						   splicesitepos1,splicesitepos2+1U,/*sign*/cdna_direction) == true) {
		    bestp = true;
		  }
		} else {
		  splicesitepos1 = (chrhigh - chroffset) - leftoffset - cL + 1;
		  splicesitepos2 = (chrhigh - chroffset) - rightoffset + cR;
		  if (IIT_exists_with_divno_signed(splicing_iit,splicing_divint_crosstable[chrnum],
						   splicesitepos2,splicesitepos1+1U,/*sign*/-cdna_direction) == true) {
		    bestp = true;
		  }
		}
		if (bestp == true) {
		  debug3(printf("At %d left to %d right, score is (%d)+(%d) = %d (bestscore)\n",
				cL,cR,scoreL,scoreR,score));
		  bestscore = score;
		  *bestrL = rL;
		  *bestrR = rR;
		  *bestcL = cL;
		  *bestcR = cR;
		} else {
		  debug3a(printf("At %d left to %d right, score is (%d)+(%d) = %d\n",
				 cL,cR,scoreL,scoreR,score));
		}
	      }
	    }
	  }
	}
      }
    }

    *finalscore = (int) bestscore;
    *best_introntype = NONINTRON;

  } else {
    left_probabilities = (double *) CALLOC(glengthL+1,sizeof(double));
    right_probabilities = (double *) CALLOC(glengthR+1,sizeof(double));

    if (watsonp == true) {
      if (cdna_direction > 0) {
	for (cL = 0; cL < glengthL - 1; cL++) {
	  splicesitepos = chroffset + leftoffset + cL;
	  if (left_known[cL]) {
	    left_probabilities[cL] = 1.0;
	  } else {
	    left_probabilities[cL] = Maxent_hr_donor_prob(splicesitepos,chroffset);
	  }
	}

	for (cR = 0; cR < glengthR - 1; cR++) {
	  splicesitepos = chroffset + rightoffset - cR + 1;
	  if (right_known[cR]) {
	    right_probabilities[cR] = 1.0;
	  } else {
	    right_probabilities[cR] = Maxent_hr_acceptor_prob(splicesitepos,chroffset);
	  }
	}

      } else {
	for (cL = 0; cL < glengthL - 1; cL++) {
	  splicesitepos = chroffset + leftoffset + cL;
	  if (left_known[cL]) {
	    left_probabilities[cL] = 1.0;
	  } else {
	    left_probabilities[cL] = Maxent_hr_antiacceptor_prob(splicesitepos,chroffset);
	  }
	}

	for (cR = 0; cR < glengthR - 1; cR++) {
	  splicesitepos = chroffset + rightoffset - cR + 1;
	  if (right_known[cR]) {
	    right_probabilities[cR] = 1.0;
	  } else {
	    right_probabilities[cR] = Maxent_hr_antidonor_prob(splicesitepos,chroffset);
	  }
	}
      }

    } else {
      if (cdna_direction > 0) {
	for (cL = 0; cL < glengthL - 1; cL++) {
	  splicesitepos = chrhigh - leftoffset - cL + 1;
	  if (left_known[cL]) {
	    left_probabilities[cL] = 1.0;
	  } else {
	    left_probabilities[cL] = Maxent_hr_antidonor_prob(splicesitepos,chroffset);
	  }
	}

	for (cR = 0; cR < glengthR - 1; cR++) {
	  splicesitepos = chrhigh - rightoffset + cR;
	  if (right_known[cR]) {
	    right_probabilities[cR] = 1.0;
	  } else {
	    right_probabilities[cR] = Maxent_hr_antiacceptor_prob(splicesitepos,chroffset);
	  }
	}

      } else {
	for (cL = 0; cL < glengthL - 1; cL++) {
	  splicesitepos = chrhigh - leftoffset - cL + 1;
	  if (left_known[cL]) {
	    left_probabilities[cL] = 1.0;
	  } else {
	    left_probabilities[cL] = Maxent_hr_acceptor_prob(splicesitepos,chroffset);
	  }
	}

	for (cR = 0; cR < glengthR - 1; cR++) {
	  splicesitepos = chrhigh - rightoffset + cR;
	  if (right_known[cR]) {
	    right_probabilities[cR] = 1.0;
	  } else {
	    right_probabilities[cR] = Maxent_hr_donor_prob(splicesitepos,chroffset);
	  }
	}
      }
    }

    /* Search using probs and without simultaneously */
    bestscore = NEG_INFINITY_32;
    bestprob_trunc = 0.0;
    for (rL = 1, rR = rlength-1; rL < rlength; rL++, rR--) {
      debug3(printf("\nAt row %d on left and %d on right\n",rL,rR));
      if ((cloL = rL - lbandL) < 1) {
	cloL = 1;
      }
      if ((chighL = rL + ubandL) > glengthL-1) {
	chighL = glengthL-1;
      }

      if ((cloR = rR - lbandR) < 1) {
	cloR = 1;
      }
      if ((chighR = rR + ubandR) > glengthR-1) {
	chighR = glengthR-1;
      }

      /* Test indels on left and right */
      for (cL = cloL; cL <= chighL; cL++) {
	/* The following check limits genomic inserts (horizontal) and
	   multiple cDNA inserts (vertical). */
	if (1) {
	  probL = left_probabilities[cL];
	  if (probL > PROB_CEILING) {
	    probL_trunc = PROB_CEILING;
	  } else if (probL < PROB_FLOOR) {
	    probL_trunc = PROB_FLOOR;
	  } else {
	    probL_trunc = probL;
	  }
	  scoreL = (int) matrixL[cL][rL];
	  if (directionsL_nogap[cL][rL] != DIAG) {
	    /* Favor gaps away from intron if possible */
	    scoreL -= 1;
	  }

	  /* Disallow leftoffset + cL >= rightoffset - cR, or cR >= rightoffset - leftoffset - cL */
	  for (cR = cloR; cR <= chighR && cR < rightoffset-leftoffset-cL; cR++) {
	    if (1) {
	      probR = right_probabilities[cR];
	      if (probR > PROB_CEILING) {
		probR_trunc = PROB_CEILING;
	      } else if (probR < PROB_FLOOR) {
		probR_trunc = PROB_FLOOR;
	      } else {
		probR_trunc = probR;
	      }
	      scoreR = (int) matrixR[cR][rR];
	      if (directionsR_nogap[cR][rR] != DIAG) {
		/* Favor gaps away from intron if possible */
		scoreR -= 1;
	      }

	      scoreI = intron_score(&introntype,leftdi[cL],rightdi[cR],
				    cdna_direction,canonical_reward,finalp);

	      if ((score = scoreL + scoreI + scoreR) > bestscore) {
		debug3(printf("No prob: At %d left to %d right, score is (%d)+(%d)+(%d) = %d (bestscore)\n",
			      cL,cR,scoreL,scoreI,scoreR,scoreL+scoreI+scoreR));
		bestscore = score;
		*bestrL = rL;
		*bestrR = rR;
		*bestcL = cL;
		*bestcR = cR;
		bestprob = probL + probR;
	      } else if (score == bestscore && probL + probR > bestprob) {
		*bestrL = rL;
		*bestrR = rR;
		*bestcL = cL;
		*bestcR = cR;
		bestprob = probL + probR;
	      }


	      if (probL_trunc + probR_trunc < bestprob_trunc) {
		debug3a(printf("At %d left to %d right, prob is %f + %f = %f\n",
			       cL,cR,probL_trunc,probR_trunc,probL_trunc+probR_trunc));

	      } else if (probL_trunc + probR_trunc == bestprob_trunc) {
		debug3(printf("At %d left to %d right, prob is %f + %f = %f\n",
			      cL,cR,probL_trunc,probR_trunc,probL_trunc+probR_trunc));

		if (scoreL + scoreI + scoreR > bestscore_with_prob) {
		  debug3(printf(" (bestscore %d)\n",scoreL+scoreI+scoreR));
		  bestprob_trunc = probL_trunc + probR_trunc;
		  bestcL_with_prob = cL;
		  bestcR_with_prob = cR;
		  bestrL_with_prob = rL;
		  bestrR_with_prob = rR;
		  bestscore_with_prob = scoreL + scoreI + scoreR;
		}
		  
	      } else {
		/* probL_trunc + probR_trunc > bestprob_trunc */
		debug3(printf("At %d left to %d right, prob is %f + %f = %f\n",
			      cL,cR,probL_trunc,probR_trunc,probL_trunc+probR_trunc));

		debug3(printf(" (bestscore %d)\n",scoreL+scoreI+scoreR));
		bestprob_trunc = probL_trunc + probR_trunc;
		bestcL_with_prob = cL;
		bestcR_with_prob = cR;
		bestrL_with_prob = rL;
		bestrR_with_prob = rR;
		bestscore_with_prob = scoreL + scoreI + scoreR;
	      }
	    }
	  }
	}
      }
    }

    debug(printf("bestscore %d (bestprob %f) vs bestscore_with_prob %d (bestprob_trunc %f)\n",
		 bestscore,bestprob,bestscore_with_prob,bestprob_trunc));
    if (left_probabilities[bestcL_with_prob] < PROB_CEILING || right_probabilities[bestcR_with_prob] < PROB_CEILING) {
      /* Use score without prob */
    } else {
      if (defect_rate < DEFECT_HIGHQ) {
	if (bestscore_with_prob > bestscore - 15) {
	  *bestcL = bestcL_with_prob;
	  *bestcR = bestcR_with_prob;
	  *bestrL = bestrL_with_prob;
	  *bestrR = bestrR_with_prob;
	  bestscore = bestscore_with_prob;
	}

      } else if (defect_rate < DEFECT_MEDQ) {
	if (bestscore_with_prob > bestscore - 25) {
	  *bestcL = bestcL_with_prob;
	  *bestcR = bestcR_with_prob;
	  *bestrL = bestrL_with_prob;
	  *bestrR = bestrR_with_prob;
	  bestscore = bestscore_with_prob;
	}

      } else {
	*bestcL = bestcL_with_prob;
	*bestcR = bestcR_with_prob;
	*bestrL = bestrL_with_prob;
	*bestrR = bestrR_with_prob;
	bestscore = bestscore_with_prob;
      }
    }
    
    scoreI = intron_score(&introntype,leftdi[*bestcL],rightdi[*bestcR],
			  cdna_direction,canonical_reward,finalp);

    if (halfp == true) {
      *finalscore = (int) (bestscore - scoreI/2);
    } else {
      *finalscore = (int) bestscore;
    }

    FREE(left_probabilities);
    FREE(right_probabilities);
  }


  /* Determine if result meets given constraints */
  if (*finalscore < 0) {
    result = false;
  } else if (novelsplicingp == true) {
    result = true;
  } else if (splicing_iit == NULL) {
    result = true;
  } else if (donor_typeint >= 0 && acceptor_typeint >= 0) {
    /* If novelsplicingp is false and using splicing at splice site level, require sites to be known */
    if (left_known[*bestcL] == 0 || right_known[*bestcR] == 0) {
      debug(printf("Novel splicing not allowed, so bridge_intron_gap returning false\n"));
      result = false;
    } else {
      result = true;
    }
  } else {
    /* If novelsplicingp is false and using splicing at splice site level, result was already constrained */
    result = true;
  }


  if (/*finalp == true &&*/ result == true) {
    get_splicesite_probs(&(*left_prob),&(*right_prob),*bestcL,*bestcR,
			 left_known,right_known,leftoffset,rightoffset,chroffset,chrhigh,
			 cdna_direction,watsonp);
  }

  debug3(printf("Returning final score of %d at (%d,%d) left to (%d,%d) right, with probs %f and %f\n",
		*finalscore,*bestrL,*bestcL,*bestrR,*bestcR,*left_prob,*right_prob));
  debug(printf("Returning final score of %d at (%d,%d) left to (%d,%d) right, with probs %f and %f\n",
	       *finalscore,*bestrL,*bestcL,*bestrR,*bestcR,*left_prob,*right_prob));

  FREE(right_known);
  FREE(left_known);

  FREE(rightdi);
  FREE(leftdi);

  return result;
}
#endif


static bool
bridge_intron_gap (int *finalscore, int *bestrL, int *bestrR, int *bestcL, int *bestcR,
		   int *best_introntype, double *left_prob, double *right_prob,
#ifdef HAVE_SSE2
		   Score16_T **matrixL, Score16_T **matrixR,
		   Direction16_T **directionsL_nogap, Direction16_T **directionsR_nogap, 
#else
		   Score32_T **matrixL, Score32_T **matrixR,
		   Direction32_T **directionsL_nogap, Direction32_T **directionsR_nogap, 
#endif
		   char *gsequenceL, char *gsequenceL_alt, char *rev_gsequenceR, char *rev_gsequenceR_alt,
		   int goffsetL, int rev_goffsetR, int rlength, int glengthL, int glengthR,
		   int cdna_direction, bool watsonp, int extraband_paired, double defect_rate, int canonical_reward,
		   int maxhorizjump, int maxvertjump, int leftoffset, int rightoffset,
		   Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
		   bool halfp, bool finalp, bool use_probabilities_p, int score_threshold,
		   bool jump_late_p) {
  bool result;
  int bestscore = NEG_INFINITY_32, score, scoreL, scoreR, scoreI, bestscoreI = NEG_INFINITY_32;
  int bestscore_with_prob = NEG_INFINITY_32;
  int rL, rR, cL, cR;
  int bestrL_with_prob, bestrR_with_prob, bestcL_with_prob, bestcR_with_prob;
  int lbandL, ubandL, cloL, chighL;
  int lbandR, ubandR, cloR, chighR;
  char left1, left2, right2, right1, left1_alt, left2_alt, right2_alt, right1_alt;
  int *leftdi, *rightdi, introntype;
  int *left_known, *right_known;
  double *left_probabilities, *right_probabilities, probL, probR, probL_trunc, probR_trunc, bestprob, bestprob_trunc, prob;
  Univcoord_T splicesitepos, splicesitepos1, splicesitepos2;
  bool bestp;

  debug(printf("Running bridge_intron_gap with use_probabilities_p %d\n",use_probabilities_p));

  if (glengthL+1 <= 0) {
    fprintf(stderr,"Problem with glengthL = %d\n",glengthL);
    abort();
  }

  if (glengthR+1 <= 0) {
    fprintf(stderr,"Problem with glengthR = %d\n",glengthR);
    abort();
  }

  /* Read dinucleotides */
  leftdi = (int *) CALLOC(glengthL+1,sizeof(int));
  rightdi = (int *) CALLOC(glengthR+1,sizeof(int));

  for (cL = 0; cL < glengthL - 1; cL++) {
    left1 = gsequenceL[cL];
    left1_alt = gsequenceL_alt[cL];
    left2 = gsequenceL[cL+1];
    left2_alt = gsequenceL_alt[cL+1];
    assert(left1 == get_genomic_nt(&left1_alt,goffsetL+cL,chroffset,chrhigh,watsonp));
    assert(left2 == get_genomic_nt(&left2_alt,goffsetL+cL+1,chroffset,chrhigh,watsonp));

    if ((left1 == 'G' || left1_alt == 'G') && (left2 == 'T' || left2_alt == 'T')) {
      leftdi[cL] = LEFT_GT;
    } else if ((left1 == 'G' || left1_alt == 'G') && (left2 == 'C' || left2_alt == 'C')) {
      leftdi[cL] = LEFT_GC;
    } else if ((left1 == 'A' || left1_alt == 'A') && (left2 == 'T' || left2_alt == 'T')) {
      leftdi[cL] = LEFT_AT;
#ifndef PMAP
    } else if ((left1 == 'C' || left1_alt == 'C') && (left2 == 'T' || left2_alt == 'T')) {
      leftdi[cL] = LEFT_CT;
#endif
    } else {
      leftdi[cL] = 0x00;
    }
  }
  leftdi[glengthL-1] = leftdi[glengthL] = 0x00;

  for (cR = 0; cR < glengthR - 1; cR++) {
    right2 = rev_gsequenceR[-cR-1];
    right2_alt = rev_gsequenceR_alt[-cR-1];
    right1 = rev_gsequenceR[-cR];
    right1_alt = rev_gsequenceR_alt[-cR];
    assert(right2 == get_genomic_nt(&right2_alt,rev_goffsetR-cR-1,chroffset,chrhigh,watsonp));
    assert(right1 == get_genomic_nt(&right1_alt,rev_goffsetR-cR,chroffset,chrhigh,watsonp));

    if ((right2 == 'A' || right2_alt == 'A') && (right1 == 'G' || right1_alt == 'G')) {
      rightdi[cR] = RIGHT_AG;
    } else if ((right2 == 'A' || right2_alt == 'A') && (right1 == 'C' || right1_alt == 'C')) {
      rightdi[cR] = RIGHT_AC;
#ifndef PMAP
    } else if ((right2 == 'G' || right2_alt == 'G') && (right1 == 'C' || right1_alt == 'C')) {
      rightdi[cR] = RIGHT_GC;
    } else if ((right2 == 'A' || right2_alt == 'A') && (right1 == 'T' || right1_alt == 'T')) {
      rightdi[cR] = RIGHT_AT;
#endif
    } else {
      rightdi[cR] = 0x00;
    }
  }
  rightdi[glengthR-1] = rightdi[glengthR] = 0x00;

  left_known = (int *) CALLOC(glengthL+1,sizeof(int));
  right_known = (int *) CALLOC(glengthR+1,sizeof(int));
  get_known_splicesites(left_known,right_known,glengthL,glengthR,
			/*leftoffset*/goffsetL,/*rightoffset*/rev_goffsetR,
			cdna_direction,watsonp,chrnum,chroffset,chrhigh);

  /* Perform computations */
#if 1
  /* Allows unlimited indel lengths */
  ubandL = glengthL - rlength + extraband_paired;
  lbandL = extraband_paired;

  ubandR = glengthR - rlength + extraband_paired;
  lbandR = extraband_paired;
#else
  /* Limit indels to 3 bp around splice sites.  Doesn't work on PacBio reads. */
  ubandL = 3;
  lbandL = 3;

  ubandR = 3;
  lbandR = 3;
#endif


  if (novelsplicingp == false && splicing_iit != NULL && (donor_typeint < 0 || acceptor_typeint < 0)) {
    /* Constrain to given introns */
    for (rL = 1, rR = rlength-1; rL < rlength; rL++, rR--) {
      debug3(printf("\nGenomic insert: At row %d on left and %d on right\n",rL,rR));
      if ((cloL = rL - lbandL) < 1) {
	cloL = 1;
      }
      if ((chighL = rL + ubandL) > glengthL-1) {
	chighL = glengthL-1;
      }

      if ((cloR = rR - lbandR) < 1) {
	cloR = 1;
      }
      if ((chighR = rR + ubandR) > glengthR-1) {
	chighR = glengthR-1;
      }

      /* Test indels on left and right */
      for (cL = cloL; cL <= chighL; cL++) {
	/* The following check limits genomic inserts (horizontal) and
	   multiple cDNA inserts (vertical). */
	if (left_known[cL] > 0) {
	  scoreL = (int) matrixL[cL][rL];
	  if (directionsL_nogap[cL][rL] != DIAG) {
	    /* Favor gaps away from intron if possible */
	    scoreL -= 1;
	  }

	  /* Disallow leftoffset + cL >= rightoffset - cR, or cR >= rightoffset - leftoffset - cL */
	  for (cR = cloR; cR <= chighR && cR < rightoffset-leftoffset-cL; cR++) {
	    if (right_known[cR] > 0) {
	      scoreR = (int) matrixR[cR][rR];
	      if (directionsR_nogap[cR][rR] != DIAG) {
		/* Favor gaps away from intron if possible */
		scoreR -= 1;
	      }

	      if ((score = scoreL + scoreR) > bestscore ||
		  (score >= bestscore && jump_late_p)) {  /* Use >= for jump late */
		bestp = false;
		if (watsonp == true) {
		  splicesitepos1 = leftoffset + cL;
		  splicesitepos2 = rightoffset - cR + 1;
		  if (IIT_exists_with_divno_signed(splicing_iit,splicing_divint_crosstable[chrnum],
						   splicesitepos1,splicesitepos2+1U,/*sign*/cdna_direction) == true) {
		    bestp = true;
		  }
		} else {
		  splicesitepos1 = (chrhigh - chroffset) - leftoffset - cL + 1;
		  splicesitepos2 = (chrhigh - chroffset) - rightoffset + cR;
		  if (IIT_exists_with_divno_signed(splicing_iit,splicing_divint_crosstable[chrnum],
						   splicesitepos2,splicesitepos1+1U,/*sign*/-cdna_direction) == true) {
		    bestp = true;
		  }
		}
		if (bestp == true) {
		  debug3(printf("At %d left to %d right, score is (%d)+(%d) = %d (bestscore)\n",
				cL,cR,scoreL,scoreR,score));
		  bestscore = score;
		  *bestrL = rL;
		  *bestrR = rR;
		  *bestcL = cL;
		  *bestcR = cR;
		} else {
		  debug3a(printf("At %d left to %d right, score is (%d)+(%d) = %d\n",
				 cL,cR,scoreL,scoreR,score));
		}
	      }
	    }
	  }
	}
      }
    }

    *finalscore = (int) bestscore;
    *best_introntype = NONINTRON;

  } else {
    left_probabilities = (double *) CALLOC(glengthL+1,sizeof(double));
    right_probabilities = (double *) CALLOC(glengthR+1,sizeof(double));

    if (watsonp == true) {
      if (cdna_direction > 0) {
	for (cL = 0; cL < glengthL - 1; cL++) {
	  splicesitepos = chroffset + leftoffset + cL;
	  if (left_known[cL]) {
	    left_probabilities[cL] = 1.0;
	  } else {
	    left_probabilities[cL] = Maxent_hr_donor_prob(splicesitepos,chroffset);
	  }
	}

	for (cR = 0; cR < glengthR - 1; cR++) {
	  splicesitepos = chroffset + rightoffset - cR + 1;
	  if (right_known[cR]) {
	    right_probabilities[cR] = 1.0;
	  } else {
	    right_probabilities[cR] = Maxent_hr_acceptor_prob(splicesitepos,chroffset);
	  }
	}

      } else {
	for (cL = 0; cL < glengthL - 1; cL++) {
	  splicesitepos = chroffset + leftoffset + cL;
	  if (left_known[cL]) {
	    left_probabilities[cL] = 1.0;
	  } else {
	    left_probabilities[cL] = Maxent_hr_antiacceptor_prob(splicesitepos,chroffset);
	  }
	}

	for (cR = 0; cR < glengthR - 1; cR++) {
	  splicesitepos = chroffset + rightoffset - cR + 1;
	  if (right_known[cR]) {
	    right_probabilities[cR] = 1.0;
	  } else {
	    right_probabilities[cR] = Maxent_hr_antidonor_prob(splicesitepos,chroffset);
	  }
	}
      }

    } else {
      if (cdna_direction > 0) {
	for (cL = 0; cL < glengthL - 1; cL++) {
	  splicesitepos = chrhigh - leftoffset - cL + 1;
	  if (left_known[cL]) {
	    left_probabilities[cL] = 1.0;
	  } else {
	    left_probabilities[cL] = Maxent_hr_antidonor_prob(splicesitepos,chroffset);
	  }
	}

	for (cR = 0; cR < glengthR - 1; cR++) {
	  splicesitepos = chrhigh - rightoffset + cR;
	  if (right_known[cR]) {
	    right_probabilities[cR] = 1.0;
	  } else {
	    right_probabilities[cR] = Maxent_hr_antiacceptor_prob(splicesitepos,chroffset);
	  }
	}

      } else {
	for (cL = 0; cL < glengthL - 1; cL++) {
	  splicesitepos = chrhigh - leftoffset - cL + 1;
	  if (left_known[cL]) {
	    left_probabilities[cL] = 1.0;
	  } else {
	    left_probabilities[cL] = Maxent_hr_acceptor_prob(splicesitepos,chroffset);
	  }
	}

	for (cR = 0; cR < glengthR - 1; cR++) {
	  splicesitepos = chrhigh - rightoffset + cR;
	  if (right_known[cR]) {
	    right_probabilities[cR] = 1.0;
	  } else {
	    right_probabilities[cR] = Maxent_hr_donor_prob(splicesitepos,chroffset);
	  }
	}
      }
    }

    /* Search using probs and without simultaneously */
    bestscore = NEG_INFINITY_16;
    bestprob_trunc = 0.0;
    for (rL = 1, rR = rlength-1; rL < rlength; rL++, rR--) {
      debug3(printf("\nAt row %d on left and %d on right\n",rL,rR));
      if ((cloL = rL - lbandL) < 1) {
	cloL = 1;
      }
      if ((chighL = rL + ubandL) > glengthL-1) {
	chighL = glengthL-1;
      }

      if ((cloR = rR - lbandR) < 1) {
	cloR = 1;
      }
      if ((chighR = rR + ubandR) > glengthR-1) {
	chighR = glengthR-1;
      }

      /* Test indels on left and right */
      for (cL = cloL; cL <= chighL; cL++) {
	/* The following check limits genomic inserts (horizontal) and
	   multiple cDNA inserts (vertical). */
	if (1) {
	  probL = left_probabilities[cL];
	  if (probL > PROB_CEILING) {
	    probL_trunc = PROB_CEILING;
	  } else if (probL < PROB_FLOOR) {
	    probL_trunc = PROB_FLOOR;
	  } else {
	    probL_trunc = probL;
	  }
	  scoreL = (int) matrixL[cL][rL];
	  if (directionsL_nogap[cL][rL] != DIAG) {
	    /* Favor gaps away from intron if possible */
	    scoreL -= 1;
	  }

	  /* Disallow leftoffset + cL >= rightoffset - cR, or cR >= rightoffset - leftoffset - cL */
	  for (cR = cloR; cR <= chighR && cR < rightoffset-leftoffset-cL; cR++) {
	    if (1) {
	      probR = right_probabilities[cR];
	      if (probR > PROB_CEILING) {
		probR_trunc = PROB_CEILING;
	      } else if (probR < PROB_FLOOR) {
		probR_trunc = PROB_FLOOR;
	      } else {
		probR_trunc = probR;
	      }
	      scoreR = (int) matrixR[cR][rR];
	      if (directionsR_nogap[cR][rR] != DIAG) {
		/* Favor gaps away from intron if possible */
		scoreR -= 1;
	      }
	      
	      scoreI = intron_score(&introntype,leftdi[cL],rightdi[cR],
				    cdna_direction,canonical_reward,finalp);

	      if ((score = scoreL + scoreI + scoreR) > bestscore) {
		debug3(printf("No prob: At %d left to %d right, score is (%d)+(%d)+(%d) = %d (bestscore)\n",
			      cL,cR,scoreL,scoreI,scoreR,scoreL+scoreI+scoreR));
		bestscore = score;
		*bestrL = rL;
		*bestrR = rR;
		*bestcL = cL;
		*bestcR = cR;
		bestprob = probL + probR;
	      } else if (score == bestscore && probL + probR > bestprob) {
		*bestrL = rL;
		*bestrR = rR;
		*bestcL = cL;
		*bestcR = cR;
		bestprob = probL + probR;
	      }


	      if (probL_trunc + probR_trunc < bestprob_trunc) {
		debug3a(printf("At %d left to %d right, prob is %f + %f = %f\n",
			       cL,cR,probL_trunc,probR_trunc,probL_trunc+probR_trunc));

	      } else if (probL_trunc + probR_trunc == bestprob_trunc) {
		debug3(printf("At %d left to %d right, prob is %f + %f = %f\n",
			      cL,cR,probL_trunc,probR_trunc,probL_trunc+probR_trunc));

		if (scoreL + scoreI + scoreR >= bestscore_with_prob) {
		  debug3(printf(" (bestscore %d)\n",scoreL+scoreI+scoreR));
		  bestprob_trunc = probL_trunc + probR_trunc;
		  bestcL_with_prob = cL;
		  bestcR_with_prob = cR;
		  bestrL_with_prob = rL;
		  bestrR_with_prob = rR;
		  bestscore_with_prob = scoreL + scoreI + scoreR;
		}

	      } else {
		/* probL_trunc + probR_trunc > bestprob_trunc */
		debug3(printf("At %d left to %d right, prob is %f + %f = %f\n",
			      cL,cR,probL_trunc,probR_trunc,probL_trunc+probR_trunc));

		debug3(printf(" (bestscore %d)\n",scoreL+scoreI+scoreR));
		bestprob_trunc = probL_trunc + probR_trunc;
		bestcL_with_prob = cL;
		bestcR_with_prob = cR;
		bestrL_with_prob = rL;
		bestrR_with_prob = rR;
		bestscore_with_prob = scoreL + scoreI + scoreR;
	      }
	    }
	  }
	}
      }
    }

    debug(printf("bestscore %d (bestprob %f) vs bestscore_with_prob %d (bestprob_trunc %f)\n",
		 bestscore,bestprob,bestscore_with_prob,bestprob_trunc));
    if (left_probabilities[bestcL_with_prob] < PROB_CEILING || right_probabilities[bestcR_with_prob] < PROB_CEILING) {
      /* Use score without prob */
    } else {
      if (defect_rate < DEFECT_HIGHQ) {
	if (bestscore_with_prob > bestscore - 15) {
	  *bestcL = bestcL_with_prob;
	  *bestcR = bestcR_with_prob;
	  *bestrL = bestrL_with_prob;
	  *bestrR = bestrR_with_prob;
	  bestscore = bestscore_with_prob;
	}

      } else if (defect_rate < DEFECT_MEDQ) {
	if (bestscore_with_prob > bestscore - 25) {
	  *bestcL = bestcL_with_prob;
	  *bestcR = bestcR_with_prob;
	  *bestrL = bestrL_with_prob;
	  *bestrR = bestrR_with_prob;
	  bestscore = bestscore_with_prob;
	}

      } else {
	*bestcL = bestcL_with_prob;
	*bestcR = bestcR_with_prob;
	*bestrL = bestrL_with_prob;
	*bestrR = bestrR_with_prob;
	bestscore = bestscore_with_prob;
      }
    }

    scoreI = intron_score(&introntype,leftdi[*bestcL],rightdi[*bestcR],
			  cdna_direction,canonical_reward,finalp);

    if (halfp == true) {
      *finalscore = (int) (bestscore - scoreI/2);
    } else {
      *finalscore = (int) bestscore;
    }

    FREE(left_probabilities);
    FREE(right_probabilities);
  }


  /* Determine if result meets given constraints */
  if (*finalscore < 0) {
    result = false;
  } else if (novelsplicingp == true) {
    result = true;
  } else if (splicing_iit == NULL) {
    result = true;
  } else if (donor_typeint >= 0 && acceptor_typeint >= 0) {
    /* If novelsplicingp is false and using splicing at splice site level, require sites to be known */
    if (left_known[*bestcL] == 0 || right_known[*bestcR] == 0) {
      debug(printf("Novel splicing not allowed, so bridge_intron_gap returning false\n"));
      result = false;
    } else {
      result = true;
    }
  } else {
    /* If novelsplicingp is false and using splicing at splice site level, result was already constrained */
    result = true;
  }


  if (/*finalp == true &&*/ result == true) {
    get_splicesite_probs(&(*left_prob),&(*right_prob),*bestcL,*bestcR,
			 left_known,right_known,leftoffset,rightoffset,chroffset,chrhigh,
			 cdna_direction,watsonp);
  }

  debug3(printf("Returning final score of %d at (%d,%d) left to (%d,%d) right, with probs %f and %f\n",
		*finalscore,*bestrL,*bestcL,*bestrR,*bestcR,*left_prob,*right_prob));
  debug(printf("Returning final score of %d at (%d,%d) left to (%d,%d) right, with probs %f and %f\n",
	       *finalscore,*bestrL,*bestcL,*bestrR,*bestcR,*left_prob,*right_prob));

  FREE(right_known);
  FREE(left_known);

  FREE(rightdi);
  FREE(leftdi);

  return result;
}


#if 0
static void
bridge_dual_break_nogap (int *finalscore, int *bestrcL, int *bestrcR,
			 Score16_T **matrixL, Score16_T **matrixR,
			 int diaglength, int goffsetL, int rev_goffsetR,
			 Univcoord_T chroffset, Univcoord_T chrhigh,
			 int cdna_direction, bool watsonp, bool jump_late_p) {
  Score16_T bestscore = NEG_INFINITY_16, score, scoreL, scoreR, scoreI;
  int rcL, rcR;
  char left1, left2, right2, right1, left1_alt, left2_alt, right2_alt, right1_alt;
  int introntype, leftdi, rightdi;

  *bestrcL = *bestrcR = 0;

  for (rcL = 1; rcL <= diaglength; rcL++) {
    rcR = diaglength - rcL;
    left1 = get_genomic_nt(&left1_alt,goffsetL+rcL,chroffset,chrhigh,watsonp);
    left2 = get_genomic_nt(&left2_alt,goffsetL+rcL+1,chroffset,chrhigh,watsonp);
#ifdef EXTRACT_GENOMICSEG
    assert(left1 == gsequence_ucL[rcL]);
    assert(left2 == gsequence_ucL[rcL+1]);
#endif
    if ((left1 == 'G' || left1_alt == 'G') && (left2 == 'T' || left2_alt == 'T')) {
      leftdi = LEFT_GT;
    } else if ((left1 == 'G' || left1_alt == 'G') && (left2 == 'C' || left2_alt == 'C')) {
      leftdi = LEFT_GC;
    } else if ((left1 == 'A' || left1_alt == 'A') && (left2 == 'T' || left2_alt == 'T')) {
      leftdi = LEFT_AT;
    } else if ((left1 == 'C' || left1_alt == 'C') && (left2 == 'T' || left2_alt == 'T')) {
      leftdi = LEFT_CT;
    } else {
      leftdi = 0x00;
    }

    right2 = get_genomic_nt(&right2_alt,rev_goffsetR-rcR-1,chroffset,chrhigh,watsonp);
    right1 = get_genomic_nt(&right1_alt,rev_goffsetR-rcR,chroffset,chrhigh,watsonp);
#ifdef EXTRACT_GENOMICSEG
    assert(right2 == rev_gsequence_ucR[-1-rcR]);
    assert(right1 == rev_gsequence_ucR[-rcR]);
#endif
    if ((right2 == 'A' || right2_alt == 'A') && (right1 == 'G' || right1_alt == 'G')) {
      rightdi = RIGHT_AG;
    } else if ((right2 == 'A' || right2_alt == 'A') && (right1 == 'C' || right1_alt == 'C')) {
      rightdi = RIGHT_AC;
    } else if ((right2 == 'G' || right2_alt == 'G') && (right1 == 'C' || right1_alt == 'C')) {
      rightdi = RIGHT_GC;
    } else if ((right2 == 'A' || right2_alt == 'A') && (right1 == 'T' || right1_alt == 'T')) {
      rightdi = RIGHT_AT;
    } else {
      rightdi = 0x00;
    }
    
    scoreL = matrixL[rcL][rcL];
    scoreI = intron_score(&introntype,leftdi,rightdi,
			  cdna_direction,/*canonical_reward*/FINAL_CANONICAL_INTRON_HIGHQ,/*finalp*/true);
    scoreR = matrixR[rcR][rcR];

    if ((score = scoreL + scoreI + scoreR) > bestscore || (score == bestscore && jump_late_p)) {
      *bestrcL = rcL;
      *bestrcR = rcR;
      bestscore = score;
    }
  }

  *finalscore = (int) bestscore;
  return;
}
#endif



#if 0
static void
bridge_dual_break_fwd (int *finalscore, int *bestrcL, int *bestrcR,
		       Score16_T **matrixL, Score16_T **matrixR,
		       int diaglength, int goffsetL, int rev_goffsetR,
		       Univcoord_T chroffset, Univcoord_T chrhigh, bool watsonp) {
  Score16_T bestscore, scoreL, scoreR, scoreI;
  int rcL, rcR;
  Score16_T bestscoreL_GT, bestscoreL_GC, bestscoreL_AT, bestscoreL_XX;
  Score16_T bestrcL_GT, bestrcL_GC, bestrcL_AT, bestrcL_XX;
  Score16_T bestscoreR_AG, bestscoreR_AC, bestscoreR_XX;
  int bestrcR_AG, bestrcR_AC, bestrcR_XX;
  char left1, left2, right2, right1, left1_alt, left2_alt, right2_alt, right1_alt;
  int introntype;

  *bestrcL = *bestrcR = 0;

  bestrcL_GT = bestrcL_GC = bestrcL_AT = bestrcL_XX = 0;
  bestscoreL_GT = bestscoreL_GC = bestscoreL_AT = bestscoreL_XX = NEG_INFINITY_16;

  for (rcL = 1; rcL <= diaglength; rcL++) {
    if ((scoreL = matrixL[rcL][rcL]) >= bestscoreL_XX) {
      bestscoreL_XX = scoreL;
      bestrcL_XX = rcL;
    }

    left1 = get_genomic_nt(&left1_alt,goffsetL+rcL,chroffset,chrhigh,watsonp);
    left2 = get_genomic_nt(&left2_alt,goffsetL+rcL+1,chroffset,chrhigh,watsonp);
#ifdef EXTRACT_GENOMICSEG
    assert(left1 == gsequence_ucL[rcL]);
    assert(left2 == gsequence_ucL[rcL+1]);
#endif
    if ((left1 == 'G' || left1_alt == 'G') && (left2 == 'T' || left2_alt == 'T')) {
      if (scoreL >= bestscoreL_GT) {
	bestscoreL_GT = scoreL;
	bestrcL_GT = rcL;
      }
    } else if ((left1 == 'G' || left1_alt == 'G') && (left2 == 'C' || left2_alt == 'C')) {
      if (scoreL >= bestscoreL_GC) {
	bestscoreL_GC = scoreL;
	bestrcL_GC = rcL;
      }
    } else if ((left1 == 'A' || left1_alt == 'A') && (left2 == 'T' || left2_alt == 'T')) {
      if (scoreL >= bestscoreL_AT) {
	bestscoreL_AT = scoreL;
	bestrcL_AT = rcL;
      }
    }
  }


  bestrcR_AG = bestrcR_AC = bestrcR_XX = 0;
  bestscoreR_AG = bestscoreR_AC = bestscoreR_XX = NEG_INFINITY_16;

  for (rcR = 1; rcR <= diaglength; rcR++) {
    if ((scoreR = matrixR[rcR][rcR]) >= bestscoreR_XX) {
      bestscoreR_XX = scoreR;
      bestrcR_XX = rcR;
    }
    right2 = get_genomic_nt(&right2_alt,rev_goffsetR-rcR-1,chroffset,chrhigh,watsonp);
    right1 = get_genomic_nt(&right1_alt,rev_goffsetR-rcR,chroffset,chrhigh,watsonp);
#ifdef EXTRACT_GENOMICSEG
    assert(right2 == rev_gsequence_ucR[-1-rcR]);
    assert(right1 == rev_gsequence_ucR[-rcR]);
#endif
    if ((right2 == 'A' || right2_alt == 'A') && (right1 == 'G' || right1_alt == 'G')) {
      if (scoreR >= bestscoreR_AG) {
	bestscoreR_AG = scoreR;
	bestrcR_AG = rcR;
      }
    } else if ((right2 == 'A' || right2_alt == 'A') && (right1 == 'C' || right1_alt == 'C')) {
      if (scoreR >= bestscoreR_AC) {
	bestscoreR_AC = scoreR;
	bestrcR_AC = rcR;
      }
    }
  }

  scoreI = intron_score(&introntype,0x00,0x00,/*cdna_direction*/+1,
			/*canonical_reward*/FINAL_CANONICAL_INTRON_HIGHQ,/*finalp*/true);
  bestscore = bestscoreL_XX + scoreI + bestscoreR_XX;
  *bestrcL = bestrcL_XX;
  *bestrcR = bestrcR_XX;

  scoreI = intron_score(&introntype,LEFT_GT,RIGHT_AG,/*cdna_direction*/+1,
			/*canonical_reward*/FINAL_CANONICAL_INTRON_HIGHQ,/*finalp*/true);
  if (bestscoreL_GT + scoreI + bestscoreR_AG >= bestscore) {
    bestscore = bestscoreL_GT + scoreI + bestscoreR_AG;
    *bestrcL = bestrcL_GT;
    *bestrcR = bestrcR_AG;
  }

  scoreI = intron_score(&introntype,LEFT_GC,RIGHT_AG,/*cdna_direction*/+1,
			/*canonical_reward*/FINAL_CANONICAL_INTRON_HIGHQ,/*finalp*/true);
  if (bestscoreL_GC + scoreI + bestscoreR_AG >= bestscore) {
    bestscore = bestscoreL_GC + scoreI + bestscoreR_AG;
    *bestrcL = bestrcL_GC;
    *bestrcR = bestrcR_AG;
  }

  scoreI = intron_score(&introntype,LEFT_AT,RIGHT_AC,/*cdna_direction*/+1,
			/*canonical_reward*/FINAL_CANONICAL_INTRON_HIGHQ,/*finalp*/true);
  if (bestscoreL_AT + scoreI + bestscoreR_AC >= bestscore) {
    bestscore = bestscoreL_AT + scoreI + bestscoreR_AC;
    *bestrcL = bestrcL_AT;
    *bestrcR = bestrcR_AC;
  }

  *finalscore = (int) bestscore;
  /*
  fprintf(stderr,"%d + %d >=? %d\n",*bestrcL,*bestrcR,diaglength);
  fprintf(stderr,"%d + %d >=? %d\n\n",*bestrcL,*bestrcR,diaglength);
  */
  if (*bestrcL + *bestrcR >= diaglength) {
    bridge_dual_break_nogap(&(*finalscore),&(*bestrcL),&(*bestrcR),matrixL,matrixR,
			    diaglength,goffsetL,rev_goffsetR,chroffset,chrhigh,
			    /*cdna_direction*/+1,watsonp,/*jump_late_p*/false);
  }

  return;
}
#endif


#if 0
static void
bridge_dual_break_rev (int *finalscore, int *bestrcL, int *bestrcR,
		       Score16_T **matrixL, Score16_T **matrixR,
		       int diaglength, int goffsetL, int rev_goffsetR,
		       Univcoord_T chroffset, Univcoord_T chrhigh, bool watsonp) {
  Score16_T bestscore, scoreL, scoreR, scoreI;
  int rcL, rcR;
  Score16_T bestscoreL_CT, bestscoreL_GT, bestscoreL_XX;
  Score16_T bestrcL_CT, bestrcL_GT, bestrcL_XX;
  Score16_T bestscoreR_AC, bestscoreR_GC, bestscoreR_AT, bestscoreR_XX;
  Score16_T bestrcR_AC, bestrcR_GC, bestrcR_AT, bestrcR_XX;
  char left1, left2, right2, right1, left1_alt, left2_alt, right2_alt, right1_alt;
  int introntype;

  *bestrcL = *bestrcR = 0;

  bestrcL_CT = bestrcL_GT = bestrcL_XX = 0;
  bestscoreL_CT = bestscoreL_GT = bestscoreL_XX = NEG_INFINITY_16;

  for (rcL = 1; rcL <= diaglength; rcL++) {
    if ((scoreL = matrixL[rcL][rcL]) >= bestscoreL_XX) {
      bestscoreL_XX = scoreL;
      bestrcL_XX = rcL;
    }

    left1 = get_genomic_nt(&left1_alt,goffsetL+rcL,chroffset,chrhigh,watsonp);
    left2 = get_genomic_nt(&left2_alt,goffsetL+rcL+1,chroffset,chrhigh,watsonp);
#ifdef EXTRACT_GENOMICSEG
    assert(left1 == gsequence_ucL[rcL]);
    assert(left2 == gsequence_ucL[rcL+1]);
#endif
    if ((left1 == 'C' || left1_alt == 'C') && (left2 == 'T' || left2_alt == 'T')) {
      if (scoreL >= bestscoreL_CT) {
	bestscoreL_CT = scoreL;
	bestrcL_CT = rcL;
      }
    } else if ((left1 == 'G' || left1_alt == 'G') && (left2 == 'T' || left2_alt == 'T')) {
      if (scoreL >= bestscoreL_GT) {
	bestscoreL_GT = scoreL;
	bestrcL_GT = rcL;
      }
    }
  }


  bestrcR_AC = bestrcR_GC = bestrcR_AT = bestrcR_XX = 0;
  bestscoreR_AC = bestscoreR_GC = bestscoreR_AT = bestscoreR_XX = NEG_INFINITY_16;

  for (rcR = 1; rcR <= diaglength; rcR++) {
    if ((scoreR = matrixR[rcR][rcR]) >= bestscoreR_XX) {
      bestscoreR_XX = scoreR;
      bestrcR_XX = rcR;
    }

    right2 = get_genomic_nt(&right2_alt,rev_goffsetR-rcR-1,chroffset,chrhigh,watsonp);
    right1 = get_genomic_nt(&right1_alt,rev_goffsetR-rcR,chroffset,chrhigh,watsonp);
#ifdef EXTRACT_GENOMICSEG
    assert(right2 == rev_gsequence_ucR[-1-rcR]);
    assert(right1 == rev_gsequence_ucR[-rcR]);
#endif

    if ((right2 == 'A' || right2_alt == 'A') && (right1 == 'C' || right1_alt == 'C')) {
      if (scoreR >= bestscoreR_AC) {
	bestscoreR_AC = scoreR;
	bestrcR_AC = rcR;
      }
    } else if ((right2 == 'G' || right2_alt == 'G') && (right1 == 'C' || right1_alt == 'C')) {
      if (scoreR >= bestscoreR_GC) {
	bestscoreR_GC = scoreR;
	bestrcR_GC = rcR;
      }
    } else if ((right2 == 'A' || right2_alt == 'A') && (right1 == 'T' || right1_alt == 'T')) {
      if (scoreR >= bestscoreR_AT) {
	bestscoreR_AT = scoreR;
	bestrcR_AT = rcR;
      }
    }
  }

  scoreI = intron_score(&introntype,0x00,0x00,/*cdna_direction*/-1,
			/*canonical_reward*/FINAL_CANONICAL_INTRON_HIGHQ,/*finalp*/true);
  bestscore = bestscoreL_XX + scoreI + bestscoreR_XX;
  *bestrcL = bestrcL_XX;
  *bestrcR = bestrcR_XX;

  scoreI = intron_score(&introntype,LEFT_CT,RIGHT_AC,/*cdna_direction*/-1,
			/*canonical_reward*/FINAL_CANONICAL_INTRON_HIGHQ,/*finalp*/true);
  if (bestscoreL_CT + scoreI + bestscoreR_AC > bestscore) {
    bestscore = bestscoreL_CT + scoreI + bestscoreR_AC;
    *bestrcL = bestrcL_CT;
    *bestrcR = bestrcR_AC;
  }

  scoreI = intron_score(&introntype,LEFT_CT,RIGHT_GC,/*cdna_direction*/-1,
			/*canonical_reward*/FINAL_CANONICAL_INTRON_HIGHQ,/*finalp*/true);
  if (bestscoreL_CT + scoreI + bestscoreR_GC > bestscore) {
    bestscore = bestscoreL_CT + scoreI + bestscoreR_GC;
    *bestrcL = bestrcL_CT;
    *bestrcR = bestrcR_GC;
  }

  scoreI = intron_score(&introntype,LEFT_GT,RIGHT_AT,/*cdna_direction*/-1,
			/*canonical_reward*/FINAL_CANONICAL_INTRON_HIGHQ,/*finalp*/true);
  if (bestscoreL_GT + scoreI + bestscoreR_AT > bestscore) {
    bestscore = bestscoreL_GT + scoreI + bestscoreR_AT;
    *bestrcL = bestrcL_GT;
    *bestrcR = bestrcR_AT;
  }

  *finalscore = (int) bestscore;
  /*
  fprintf(stderr,"%d + %d >=? %d\n",*bestrcL,*bestrcR,diaglength);
  fprintf(stderr,"%d + %d >=? %d\n\n",*bestrcL,*bestrcR,diaglength);
  */
  if (*bestrcL + *bestrcR >= diaglength) {
    bridge_dual_break_nogap(&(*finalscore),&(*bestrcL),&(*bestrcR),matrixL,matrixR,
			    diaglength,goffsetL,rev_goffsetR,chroffset,chrhigh,
			    /*cdna_direction*/-1,watsonp,/*jump_late_p*/true);
  }

  return;
}
#endif



/************************************************************************/

static List_T
single_gap_simple (int *finalscore, int *nmatches, int *nmismatches,
		   char *rsequence, int rlength, char *gsequence, char *gsequence_alt,
		   int roffset, int goffset, Pairpool_T pairpool,
		   int mismatchtype, int dynprogindex) {
  int score;
  List_T pairs = NULL;
  int r;
  int querycoord, genomecoord;
  int c1, c2, c2_alt;
  Pairdistance_T **pairdistance_array_type;

  debug(printf("Starting single_gap_simple\n"));
  pairdistance_array_type = pairdistance_array[mismatchtype];

  *finalscore = 0;
  *nmatches = *nmismatches = 0;

  /* Push from left to right, so we don't need to do List_reverse() later */
  for (r = 1; r <= rlength; r++) {
    querycoord = genomecoord = r-1;

    c1 = rsequence[querycoord];
    c2 = gsequence[genomecoord];
    c2_alt = gsequence_alt[genomecoord];

    if (c2 == '*') {
      /* Don't push pairs past end of chromosome */
      debug(printf("Don't push pairs past end of chromosome: genomeoffset %u, genomecoord %u\n",goffset,genomecoord));

    } else if (c1 == c2 || c1 == c2_alt) {
      debug(printf("Pushing simple %d,%d [%d,%d] (%c,%c) - match\n",
		   r,/*c*/r,roffset+querycoord,goffset+genomecoord,c1,c2));
      score = pairdistance_array_type[c1][c2];
      if (pairdistance_array_type[c1][c2_alt] > score) {
	score = pairdistance_array_type[c1][c2_alt];
      }
      *finalscore += score;
      *nmatches += 1;
      pairs = Pairpool_push(pairs,pairpool,roffset+querycoord,goffset+genomecoord,
			    c1,DYNPROG_MATCH_COMP,c2,c2_alt,dynprogindex);
	
    } else if (consistent_array[(int) c1][(int) c2] == true) {
      debug(printf("Pushing simple %d,%d [%d,%d] (%c,%c) - ambiguous\n",
		   r,/*c*/r,roffset+querycoord,goffset+genomecoord,c1,c2));
      score = pairdistance_array_type[c1][c2];
      if (pairdistance_array_type[c1][c2_alt] > score) {
	score = pairdistance_array_type[c1][c2_alt];
      }
      *finalscore += score;
      *nmatches += 1;
      pairs = Pairpool_push(pairs,pairpool,roffset+querycoord,goffset+genomecoord,
			    c1,AMBIGUOUS_COMP,c2,c2_alt,dynprogindex);
	
    } else {
      debug(printf("Pushing simple %d,%d [%d,%d] (%c,%c) - mismatch\n",
		   r,/*c*/r,roffset+querycoord,goffset+genomecoord,c1,c2));
      score = pairdistance_array_type[c1][c2];
      if (pairdistance_array_type[c1][c2_alt] > score) {
	score = pairdistance_array_type[c1][c2_alt];
      }
      *finalscore += score;
      *nmismatches += 1;
      pairs = Pairpool_push(pairs,pairpool,roffset+querycoord,goffset+genomecoord,
			    c1,MISMATCH_COMP,c2,c2_alt,dynprogindex);
    }
  }

  if (*nmismatches > 1) {
    return (List_T) NULL;
  } else {
    return pairs;
  }
}


static char *
uniq_string (int **nreps, int *uniqlength, char *string, int length) {
  char *uniq, *p, nt, lastnt;
  int i, k, *a;

  *uniqlength = 1;
  lastnt = string[0];
  for (i = 1; i < length; i++) {
    if ((nt = string[i]) != lastnt) {
      (*uniqlength)++;
      lastnt = nt;
    }
  }

  p = uniq = (char *) MALLOC(((*uniqlength) + 1) * sizeof(char));
  a = *nreps = (int *) MALLOC((*uniqlength) * sizeof(int));
  k = 0;

  lastnt = string[0];
  for (i = 1; i < length; i++) {
    if ((nt = string[i]) != lastnt) {
      *p++ = lastnt;
      *a++ = k;
      lastnt = nt;
      k = 0;
    } else {
      k++;
    }
  }
  *p = lastnt;
  *a = k;

#ifdef DEBUG16
  printf("string: %.*s\n",length,string);
  printf("uniq:   %.*s\n",*uniqlength,uniq);
  printf("nreps   ");
  for (i = 0; i < *uniqlength; i++) {
    printf("%d",(*nreps)[i]);
  }
  printf("\n");
#endif

  return uniq;
}


static List_T
augment_pairs (List_T pairs, int *rsequence_nreps, int r_uniqlength, int roffset,
	       int *gsequence_nreps, int g_uniqlength, int goffset,
	       Pairpool_T pairpool, int dynprogindex) {
  List_T augmented = NULL, p;
  Pair_T pair;
  int r, c, r_nreps, c_nreps, i;
  int r_cum = 0, c_cum = 0;
  bool add_dashes_p;

#ifdef DEBUG16
  printf("r_nreps: ");
  for (i = 0; i < r_uniqlength; i++) {
    printf("%d",rsequence_nreps[i]);
  }
  printf("\n");

  printf("g_nreps: ");
  for (i = 0; i < g_uniqlength; i++) {
    printf("%d",gsequence_nreps[i]);
  }
  printf("\n");
#endif

  for (p = pairs; p != NULL; p = p->rest) {
    pair = (Pair_T) p->first;
    r = pair->querypos - roffset;
    c = pair->genomepos - goffset;

    if (pair->comp == INDEL_COMP) {
      augmented = Pairpool_push_existing(augmented,pairpool,pair);
      debug16a(Pair_dump_one(pair,true));
      pair->querypos += r_cum;
      pair->genomepos += c_cum;

      if (pair->cdna == ' ') {
	c_nreps = gsequence_nreps[c];
	debug16a(printf(" genomepos %u, c_cum %d, nreps %d\n",c,c_cum,c_nreps));
	if (c_nreps == 0) {
	  /* Do nothing */
	} else {
	  for (i = 1; i <= c_nreps; i++) {
	    augmented = Pairpool_push(augmented,pairpool,r+r_cum+roffset,c+c_cum+i+goffset,
				      /*cdna*/' ',INDEL_COMP,pair->genome,pair->genomealt,
				      dynprogindex);
	  }
	  c_cum += c_nreps;
	}

      } else if (pair->genome == ' ') {
	r_nreps = rsequence_nreps[r];
	debug16a(printf(" querypos %d, r_cum %d, nreps %d\n",r,r_cum,r_nreps));
	if (r_nreps == 0) {
	  /* Do nothing */
	} else {
	  for (i = 1; i <= r_nreps; i++) {
	    augmented = Pairpool_push(augmented,pairpool,r+r_cum+i+roffset,c+c_cum+goffset,
				      pair->cdna,INDEL_COMP,/*genome*/' ',/*genomealt*/' ',
				      dynprogindex);
	  }
	  r_cum += r_nreps;
	}

      } else {
	fprintf(stderr,"Indel pair is missing both cdna and genome nts\n");
	abort();
      }

    } else {
      r_nreps = rsequence_nreps[r];
      c_nreps = gsequence_nreps[c];
      debug16a(printf(" querypos %d, r_cum %d, nreps %d, genomepos %u, c_cum %d, nreps %d\n",
		      r,r_cum,r_nreps,c,c_cum,c_nreps));
      augmented = Pairpool_push_existing(augmented,pairpool,pair);
      pair->querypos += r_cum;
      pair->genomepos += c_cum;
      if (r_nreps == 0 && c_nreps == 0) {
	/* Do nothing */
      } else if (r_nreps == c_nreps) {
	for (i = 1; i <= r_nreps; i++) {
	  augmented = Pairpool_push_copy(augmented,pairpool,pair);
	  ((Pair_T) augmented->first)->querypos += i;
	  ((Pair_T) augmented->first)->genomepos += i;
	}
	r_cum += r_nreps;
	c_cum += c_nreps;

      } else if (r_nreps < c_nreps) {
	for (i = 1; i <= r_nreps; i++) {
	  augmented = Pairpool_push_copy(augmented,pairpool,pair);
	  ((Pair_T) augmented->first)->querypos += i;
	  ((Pair_T) augmented->first)->genomepos += i;
	}
	r_cum += r_nreps;

	for ( ; i <= c_nreps; i++) {
	  /* Add 1 to r to advance to next coordinate */
	  augmented = Pairpool_push(augmented,pairpool,r+r_cum+roffset + 1,c+c_cum+i+goffset,
				    /*cdna*/' ',INDEL_COMP,pair->genome,pair->genomealt,
				    dynprogindex);
	}
	c_cum += c_nreps;

      } else {
	for (i = 1; i <= c_nreps; i++) {
	  augmented = Pairpool_push_copy(augmented,pairpool,pair);
	  ((Pair_T) augmented->first)->querypos += i;
	  ((Pair_T) augmented->first)->genomepos += i;
	}
	c_cum += c_nreps;

	for ( ; i <= r_nreps; i++) {
	  /* Add 1 to c to advance to next coordinate */
	  augmented = Pairpool_push(augmented,pairpool,r+r_cum+i+roffset,c+c_cum+goffset + 1,
				    pair->cdna,INDEL_COMP,/*genome*/' ',/*genomealt*/' ',
				    dynprogindex);
	}
	r_cum += r_nreps;

      }
    }
  }

  debug16(Pair_dump_list(augmented,true));

  return List_reverse(augmented);
}




List_T
Dynprog_single_gap (int *dynprogindex, int *finalscore, int *nmatches, int *nmismatches, int *nopens, int *nindels,
		    T dynprog, char *rsequence, char *sequenceuc1,
		    int rlength, int glength, int roffset, int goffset,
		    Univcoord_T chroffset, Univcoord_T chrhigh,
#ifdef PMAP
		    char *queryaaseq,
#endif
		    bool watsonp, bool jump_late_p, Pairpool_T pairpool,
		    int extraband_single, double defect_rate, bool widebandp) {
  List_T pairs = NULL;
  char *gsequence, *gsequence_alt;

  char *gsequence_orig, *rsequence_orig;
  int *gsequence_nreps, *rsequence_nreps;
  int glength_orig, rlength_orig;

#ifdef PMAP
  char *inst_rsequence;
#endif
  Mismatchtype_T mismatchtype;
  int lband, uband;
#if defined(HAVE_SSE4_1) || defined(HAVE_SSE2)
  Score8_T **matrix8;
  Direction8_T **directions8_nogap, **directions8_Egap, **directions8_Fgap;

  Score16_T **matrix, open, extend;
  Direction16_T **directions_nogap, **directions_Egap, **directions_Fgap;
#else
  Score32_T **matrix, open, extend;
  Direction32_T **directions_nogap, **directions_Egap, **directions_Fgap;
#endif
  /* bool onesidegapp; */

  if (defect_rate < DEFECT_HIGHQ) {
    mismatchtype = HIGHQ;
    open = SINGLE_OPEN_HIGHQ;
    extend = SINGLE_EXTEND_HIGHQ;
    /* onesidegapp = false; */
  } else if (defect_rate < DEFECT_MEDQ) {
    mismatchtype = MEDQ;
    open = SINGLE_OPEN_MEDQ;
    extend = SINGLE_EXTEND_MEDQ;
    /* onesidegapp = true; */
  } else {
    mismatchtype = LOWQ;
    open = SINGLE_OPEN_LOWQ;
    extend = SINGLE_EXTEND_LOWQ;
    /* onesidegapp = true; */
  }

#if 0
  if (close_indels_mode == +1) {
    /* Allow close indels */
    onesidegapp = false;
  } else if (close_indels_mode == -1) {
    /* Disallow close indels */
    onesidegapp = true;
  } else {
    /* Allow close indels for high quality alignments, as determined above */
  }
#endif    

  /* Rlength: maxlookback+MAXPEELBACK.  Glength +EXTRAMATERIAL */
  debug(
	printf("%c:  ",*dynprogindex > 0 ? (*dynprogindex-1)%26+'a' : (-(*dynprogindex)-1)%26+'A');
	printf("Aligning single gap middle with wideband = %d and extraband %d\n",widebandp,extraband_single);
	printf("At query offset %d-%d, %.*s\n",roffset,roffset+rlength-1,rlength,rsequence));
#ifdef EXTRACT_GENOMICSEG
  debug(printf("At genomic offset %d-%d, %.*s\n",goffset,goffset+glength-1,glength,gsequence));
#endif
  debug(printf("\n"));

#if 0
  /* Can happen in bad genomic regions.  May want to give up, though, if rlength and glength differ greatly. */
  assert(glength < 1000);
#endif

  if (rlength > dynprog->max_rlength || glength > dynprog->max_glength) {
    debug(printf("rlength %d or glength %d is too long.  Returning NULL\n",rlength,glength));
    *finalscore = NEG_INFINITY_32;
    *nmatches = *nmismatches = *nopens = *nindels = 0;
#if 0
    /* Don't push a gapholder for single gap, because gapholder already exists */
    pairs = Pairpool_push_gapholder(NULL,pairpool,rlength,glength,
				    /*leftpair*/NULL,/*rightpair*/NULL,/*knownp*/false);
#endif
    *dynprogindex += (*dynprogindex > 0 ? +1 : -1);
    return (List_T) NULL;
  }

#ifdef PMAP
  inst_rsequence = instantiate_codons(sequenceuc1,queryaaseq,roffset,rlength);
#endif
  
    /* If extraband_single is too large, then gaps may be inserted on
       both sides, like this

       CACCC   AGAGCAGGCACTGCCT
       |||||--- ||| ||---||||| 
       CACCCCAGGGAGGAG   CTGCCC

    */


  if (homopolymerp == true) {
    rsequence_orig = rsequence;
    rlength_orig = rlength;
    rsequence = uniq_string(&rsequence_nreps,&rlength,rsequence_orig,rlength_orig);

    if (watsonp) {
      gsequence_orig = Genome_get_segment_blocks_right(&gsequence_alt,/*left*/chroffset+goffset,glength,chrhigh,/*revcomp*/false);
    } else {
      gsequence_orig = Genome_get_segment_blocks_left(&gsequence_alt,/*left*/chrhigh-goffset+1,glength,chroffset,/*revcomp*/true);
    }
    if (gsequence_orig == NULL) {
      *finalscore = NEG_INFINITY_32;
      *nmatches = *nmismatches = *nopens = *nindels = 0;
      return (List_T) NULL;
    } else {
      glength_orig = glength;
      gsequence = gsequence_alt = uniq_string(&gsequence_nreps,&glength,gsequence_orig,glength_orig);
    }

  } else {
    if (watsonp) {
      gsequence = Genome_get_segment_blocks_right(&gsequence_alt,/*left*/chroffset+goffset,glength,chrhigh,/*revcomp*/false);
    } else {
      gsequence = Genome_get_segment_blocks_left(&gsequence_alt,/*left*/chrhigh-goffset+1,glength,chroffset,/*revcomp*/true);
    }
    if (gsequence == NULL) {
      *finalscore = NEG_INFINITY_32;
      *nmatches = *nmismatches = *nopens = *nindels = 0;
      return (List_T) NULL;
    } else if (glength == rlength &&
	       (pairs = single_gap_simple(&(*finalscore),&(*nmatches),&(*nmismatches),
					  rsequence,rlength,gsequence,gsequence_alt,roffset,goffset,
					  pairpool,mismatchtype,*dynprogindex)) != NULL) {
#ifdef PMAP
      FREE(inst_rsequence);
#endif
      if (gsequence_alt != gsequence) {
	FREE(gsequence_alt);
      }
      FREE(gsequence);

      *nopens = *nindels = 0;
      *dynprogindex += (*dynprogindex > 0 ? +1 : -1);
      return pairs;
    }
  }


  compute_bands(&lband,&uband,rlength,glength,extraband_single,widebandp);
#if defined(HAVE_SSE4_1) || defined(HAVE_SSE2)
  /* Use || because we want the minimum length (which determines the diagonal length) to achieve a score less than 128 */
  if (rlength <= SIMD_MAXLENGTH_EPI8 || glength <= SIMD_MAXLENGTH_EPI8) {
    matrix8 = compute_scores_simd_8(&directions8_nogap,&directions8_Egap,&directions8_Fgap,dynprog,
#ifdef PMAP
				    inst_rsequence,
#else
				    rsequence,
#endif
				    gsequence,gsequence_alt,rlength,glength,
#ifdef DEBUG14
				    goffset,chroffset,chrhigh,watsonp,
#endif
				    mismatchtype,(Score8_T) open,(Score8_T) extend,
				    lband,uband,jump_late_p,/*revp*/false);
    *finalscore = (int) matrix8[glength][rlength];

    *nmatches = *nmismatches = *nopens = *nindels = 0;
    pairs = traceback_8(NULL,&(*nmatches),&(*nmismatches),&(*nopens),&(*nindels),
			directions8_nogap,directions8_Egap,directions8_Fgap,rlength,glength,
#ifdef PMAP
			inst_rsequence,inst_rsequence,
#else
			rsequence,sequenceuc1,
#endif
			gsequence,gsequence_alt,roffset,goffset,pairpool,/*revp*/false,
			chroffset,chrhigh,/*cdna_direction*/0,watsonp,*dynprogindex);

  } else {
    matrix = compute_scores_simd_16(&directions_nogap,&directions_Egap,&directions_Fgap,dynprog,
#ifdef PMAP
				    inst_rsequence,
#else
				    rsequence,
#endif
				    gsequence,gsequence_alt,rlength,glength,
#ifdef DEBUG14
				    goffset,chroffset,chrhigh,watsonp,
#endif
				    mismatchtype,(Score16_T) open, (Score16_T) extend,
				    lband,uband,jump_late_p,/*revp*/false);
    *finalscore = (int) matrix[glength][rlength];

    *nmatches = *nmismatches = *nopens = *nindels = 0;
    pairs = traceback(NULL,&(*nmatches),&(*nmismatches),&(*nopens),&(*nindels),
		      directions_nogap,directions_Egap,directions_Fgap,rlength,glength,
#ifdef PMAP
		      inst_rsequence,inst_rsequence,
#else
		      rsequence,sequenceuc1,
#endif
		      gsequence,gsequence_alt,roffset,goffset,pairpool,/*revp*/false,
		      chroffset,chrhigh,/*cdna_direction*/0,watsonp,*dynprogindex);
  }

#else

  matrix = compute_scores_standard(&directions_nogap,&directions_Egap,&directions_Fgap,dynprog,
#ifdef PMAP
				   inst_rsequence,
#else
				   rsequence,
#endif
				   gsequence,gsequence_alt,rlength,glength,
				   goffset,chroffset,chrhigh,watsonp,
				   mismatchtype,open,extend,
				   lband,uband,jump_late_p,/*revp*/false);
  *finalscore = (int) matrix[glength][rlength];

  *nmatches = *nmismatches = *nopens = *nindels = 0;
  pairs = traceback(NULL,&(*nmatches),&(*nmismatches),&(*nopens),&(*nindels),
		    directions_nogap,directions_Egap,directions_Fgap,rlength,glength,
#ifdef PMAP
		    inst_rsequence,inst_rsequence,
#else
		    rsequence,sequenceuc1,
#endif
		    gsequence,gsequence_alt,roffset,goffset,pairpool,/*revp*/false,
		    chroffset,chrhigh,/*cdna_direction*/0,watsonp,*dynprogindex);
#endif

#ifdef PMAP
  FREE(inst_rsequence);
#endif

  if (homopolymerp == true) {
    pairs = augment_pairs(pairs,rsequence_nreps,rlength,roffset,
			  gsequence_nreps,glength,goffset,pairpool,*dynprogindex);
    FREE(gsequence_nreps);
    FREE(gsequence_orig);
    FREE(gsequence);

    FREE(rsequence_nreps);
    /* Do not free rsequence_orig */
    FREE(rsequence);
    
  } else {
    if (gsequence_alt != gsequence) {
      FREE(gsequence_alt);
    }
    FREE(gsequence);
  }

  /*
  Directions_free(directions);
  Matrix_free(matrix);
  */

  debug(printf("End of dynprog single gap\n"));

  *dynprogindex += (*dynprogindex > 0 ? +1 : -1);
  return List_reverse(pairs);
}


/* Sequences rsequenceL and rsequenceR represent the two ends of the cDNA insertion */
List_T
Dynprog_cdna_gap (int *dynprogindex, int *finalscore, bool *incompletep,
		  T dynprogL, T dynprogR, char *rsequenceL, char *rsequence_ucL, 
		  char *rev_rsequenceR, char *rev_rsequence_ucR,
#if 0
		  char *gsequence, char *gsequence_uc,
#endif
		  int rlengthL, int rlengthR, int glength,
		  int roffsetL, int rev_roffsetR, int goffset,
		  Univcoord_T chroffset, Univcoord_T chrhigh,
#ifdef PMAP
		  char *queryaaseq,
#endif
		  int cdna_direction, bool watsonp, bool jump_late_p, Pairpool_T pairpool,
		  int extraband_paired, double defect_rate) {
  List_T pairs = NULL;
  char *gsequence, *gsequence_alt, *rev_gsequence, *rev_gsequence_alt;
#ifdef PMAP
  char *inst_rsequence, *inst_rsequenceL, *inst_rev_rsequenceR;
#endif
  Mismatchtype_T mismatchtype;
#if defined(HAVE_SSE4_1) || defined(HAVE_SSE2)
  Score8_T **matrix8L, **matrix8R;
  Direction8_T **directions8L_nogap, **directions8L_Egap, **directions8L_Fgap,
    **directions8R_nogap, **directions8R_Egap, **directions8R_Fgap;
  bool use8p;

  Score16_T **matrixL, **matrixR, mismatch, open, extend;
  Direction16_T **directionsL_nogap, **directionsL_Egap, **directionsL_Fgap,
    **directionsR_nogap, **directionsR_Egap, **directionsR_Fgap;
#else
  Score32_T **matrixL, **matrixR, mismatch, open, extend;
  Direction32_T **directionsL_nogap, **directionsL_Egap, **directionsL_Fgap,
    **directionsR_nogap, **directionsR_Egap, **directionsR_Fgap;
#endif
  int rev_goffset, bestrL, bestrR, bestcL, bestcR, lband, uband, k;
  int nmatches, nmismatches, nopens, nindels;
  int queryjump, genomejump;
  char c2, c2_alt;


  if (glength <= 1) {
    return NULL;
  }

  debug(
	printf("\n");
	printf("%c:  ",*dynprogindex > 0 ? (*dynprogindex-1)%26+'a' : (-(*dynprogindex)-1)%26+'A');
	printf("Aligning cdna gap\n");
	printf("At query offset %d-%d, %.*s\n",roffsetL,roffsetL+rlengthL-1,rlengthL,rsequenceL);
	printf("At query offset %d-%d, %.*s\n",rev_roffsetR-rlengthR+1,rev_roffsetR,rlengthR,&(rev_rsequenceR[-rlengthR+1]));
	printf("Whole piece at query offset %d-%d, %.*s\n",roffsetL,rev_roffsetR,rev_roffsetR-roffsetL+1,rsequenceL));
#ifdef EXTRACT_GENOMICSEG
  debug(printf("At genomic offset %d-%d, %.*s\n",goffset,goffset+glength-1,glength,gsequence));
#endif
  debug(printf("\n"));

  /* ?check if offsets are too close.  But this eliminates a segment
     of the cDNA.  Should check in stage 3, and do single gap instead. */
  /*
  if (roffsetL+rlengthL-1 >= rev_roffsetR-rlengthR+1) {
    debug(printf("Bounds don't make sense\n"));
    *finalscore = NEG_INFINITY_16;
    return NULL;
  }
  */

  if (defect_rate < DEFECT_HIGHQ) {
    mismatchtype = HIGHQ;
    mismatch = MISMATCH_HIGHQ;
    open = CDNA_OPEN_HIGHQ;
    extend = CDNA_EXTEND_HIGHQ;
  } else if (defect_rate < DEFECT_MEDQ) {
    mismatchtype = MEDQ;
    mismatch = MISMATCH_MEDQ;
    open = CDNA_OPEN_MEDQ;
    extend = CDNA_EXTEND_MEDQ;
  } else {
    mismatchtype = LOWQ;
    mismatch = MISMATCH_LOWQ;
    open = CDNA_OPEN_LOWQ;
    extend = CDNA_EXTEND_LOWQ;
  }

  if (glength > dynprogR->max_glength || rlengthR > dynprogR->max_rlength) {
    debug(printf("glength %d or rlengthR %d is too long.  Returning NULL\n",glength,rlengthR));
#if 0
    rev_goffset = goffset + glength - 1;
    queryjump = rev_roffsetR - roffsetL + 1;
    genomejump = rev_goffset - goffset + 1;
    pairs = Pairpool_push_gapholder(NULL,pairpool,queryjump,genomejump,
				    /*leftpair*/NULL,/*rightpair*/NULL,/*knownp*/false);
#endif
    *dynprogindex += (*dynprogindex > 0 ? +1 : -1);
    return (List_T) NULL;
  }

  if (glength > dynprogL->max_glength || rlengthL > dynprogL->max_rlength) {
    debug(printf("glength %d or rlengthL %d is too long.  Returning NULL\n",glength,rlengthL));
#if 0
    rev_goffset = goffset + glength - 1;
    queryjump = rev_roffsetR - roffsetL + 1;
    genomejump = rev_goffset - goffset + 1;
    pairs = Pairpool_push_gapholder(NULL,pairpool,queryjump,genomejump,
				    /*leftpair*/NULL,/*rightpair*/NULL,/*knownp*/false);
#endif
    *dynprogindex += (*dynprogindex > 0 ? +1 : -1);
    return (List_T) NULL;
  }

#if 0
  /* Right side looks like 5' end */
  /* Note: sequence 1 and 2 flipped, because 1 has extramaterial */
  rev_gsequence = &(gsequence[glength-1]);
  rev_gsequence_uc = &(gsequence_uc[glength-1]);
#endif
  rev_goffset = goffset+glength-1;

#ifdef PMAP
  inst_rsequence = instantiate_codons(rsequenceL,queryaaseq,roffsetL,rev_roffsetR-roffsetL+1);
  inst_rsequenceL = inst_rsequence;
  inst_rev_rsequenceR = &(inst_rsequence[rev_roffsetR-roffsetL]);
#endif

#if defined(HAVE_SSE4_1) || defined(HAVE_SSE2)
  /* Use || because we want the minimum length (which determines the diagonal length) to achieve a score less than 128 */
  if (glength <= SIMD_MAXLENGTH_EPI8 || (rlengthL <= SIMD_MAXLENGTH_EPI8 && rlengthR <= SIMD_MAXLENGTH_EPI8)) {
    use8p = true;
  } else {
    use8p = false;
  }
#endif


  if (watsonp) {
    rev_gsequence = Genome_get_segment_blocks_left(&rev_gsequence_alt,/*left*/chroffset+rev_goffset+1,glength,chroffset,/*revcomp*/false);
    gsequence = Genome_get_segment_blocks_right(&gsequence_alt,/*left*/chroffset+goffset,glength,chrhigh,/*revcomp*/false);
  } else {
    rev_gsequence = Genome_get_segment_blocks_right(&rev_gsequence_alt,/*left*/chrhigh-rev_goffset,glength,chrhigh,/*revcomp*/true);
    gsequence = Genome_get_segment_blocks_left(&gsequence_alt,/*left*/chrhigh-goffset+1,glength,chroffset,/*revcomp*/true);
  }
  if (gsequence == NULL) {
    FREE(rev_gsequence);
    return (List_T) NULL;
  } else if (rev_gsequence == NULL) {
    FREE(gsequence);
    return (List_T) NULL;
  }

  compute_bands(&lband,&uband,rlengthR,glength,extraband_paired,/*widebandp*/true);
#if defined(HAVE_SSE4_1) || defined(HAVE_SSE2)
  if (use8p == true) {
    matrix8R = compute_scores_simd_8(&directions8R_nogap,&directions8R_Egap,&directions8R_Fgap,dynprogR,
#ifdef PMAP
				     inst_rev_rsequenceR,
#else
				     rev_rsequenceR,
#endif
				     &(rev_gsequence[glength-1]),&(rev_gsequence_alt[glength-1]),
				     rlengthR,glength,
#ifdef DEBUG14
				     rev_goffset,chroffset,chrhigh,watsonp,
#endif
				     mismatchtype,(Score8_T) open,(Score8_T) extend,
				     lband,uband,/*for revp true*/!jump_late_p,/*revp*/true);

  } else {
    matrixR = compute_scores_simd_16(&directionsR_nogap,&directionsR_Egap,&directionsR_Fgap,dynprogR,
#ifdef PMAP
				     inst_rev_rsequenceR,
#else
				     rev_rsequenceR,
#endif
				     &(rev_gsequence[glength-1]),&(rev_gsequence_alt[glength-1]),
				     rlengthR,glength,
#ifdef DEBUG14
				     rev_goffset,chroffset,chrhigh,watsonp,
#endif
				     mismatchtype,(Score16_T) open, (Score16_T) extend,
				     lband,uband,/*for revp true*/!jump_late_p,/*revp*/true);
  }

#else
  matrixR = compute_scores_standard(&directionsR_nogap,&directionsR_Egap,&directionsR_Fgap,dynprogR,
#ifdef PMAP
				    inst_rev_rsequenceR,
#else
				    rev_rsequenceR,
#endif
				    &(rev_gsequence[glength-1]),&(rev_gsequence_alt[glength-1]),
				    rlengthR,glength,
				    rev_goffset,chroffset,chrhigh,watsonp,
				    mismatchtype,open,extend,
				    lband,uband,/*for revp true*/!jump_late_p,/*revp*/true);
#endif


  compute_bands(&lband,&uband,rlengthL,glength,extraband_paired,/*widebandp*/true);
#if defined(HAVE_SSE4_1) || defined(HAVE_SSE2)
  if (use8p == true) {
    matrix8L = compute_scores_simd_8(&directions8L_nogap,&directions8L_Egap,&directions8L_Fgap,dynprogL,
#ifdef PMAP
				     inst_rsequenceL,
#else
				     rsequenceL,
#endif
				     gsequence,gsequence_alt,rlengthL,glength,
#ifdef DEBUG14
				     goffset,chroffset,chrhigh,watsonp,
#endif
				     mismatchtype,(Score8_T) open,(Score8_T) extend,
				     lband,uband,jump_late_p,/*revp*/false);
  } else {
    matrixL = compute_scores_simd_16(&directionsL_nogap,&directionsL_Egap,&directionsL_Fgap,dynprogL,
#ifdef PMAP
				     inst_rsequenceL,
#else
				     rsequenceL,
#endif
				     gsequence,gsequence_alt,rlengthL,glength,
#ifdef DEBUG14
				     goffset,chroffset,chrhigh,watsonp,
#endif
				     mismatchtype,(Score16_T) open, (Score16_T) extend,
				     lband,uband,jump_late_p,/*revp*/false);
  }

#else
  matrixL = compute_scores_standard(&directionsL_nogap,&directionsL_Egap,&directionsL_Fgap,dynprogL,
#ifdef PMAP
				    inst_rsequenceL,
#else
				    rsequenceL,
#endif
				    gsequence,gsequence_alt,rlengthL,glength,
				    goffset,chroffset,chrhigh,watsonp,
				    mismatchtype,open,extend,
				    lband,uband,jump_late_p,/*revp*/false);
#endif

  nmatches = nmismatches = nopens = nindels = 0;
#if defined(HAVE_SSE4_1) || defined(HAVE_SSE2)
  if (use8p == true) {
    bridge_cdna_gap_8(&(*finalscore),&bestcL,&bestcR,&bestrL,&bestrR,matrix8L,matrix8R,
		      glength,rlengthL,rlengthR,extraband_paired,
		      open,extend,roffsetL,rev_roffsetR,jump_late_p);

    pairs = traceback_8(NULL,&nmatches,&nmismatches,&nopens,&nindels,
			directions8R_nogap,directions8R_Egap,directions8R_Fgap,bestrR,bestcR,
#ifdef PMAP
			inst_rev_rsequenceR,inst_rev_rsequenceR,
#else
			rev_rsequenceR,rev_rsequence_ucR,
#endif
			&(rev_gsequence[glength-1]),&(rev_gsequence_alt[glength-1]),
			rev_roffsetR,rev_goffset,pairpool,/*revp*/true,
			chroffset,chrhigh,cdna_direction,watsonp,*dynprogindex);
  } else {
#endif

    bridge_cdna_gap(&(*finalscore),&bestcL,&bestcR,&bestrL,&bestrR,matrixL,matrixR,
		    glength,rlengthL,rlengthR,extraband_paired,
		    open,extend,roffsetL,rev_roffsetR,jump_late_p);

    pairs = traceback(NULL,&nmatches,&nmismatches,&nopens,&nindels,
		      directionsR_nogap,directionsR_Egap,directionsR_Fgap,bestrR,bestcR,
#ifdef PMAP
		      inst_rev_rsequenceR,inst_rev_rsequenceR,
#else
		      rev_rsequenceR,rev_rsequence_ucR,
#endif
		      &(rev_gsequence[glength-1]),&(rev_gsequence_alt[glength-1]),
		      rev_roffsetR,rev_goffset,pairpool,/*revp*/true,
		      chroffset,chrhigh,cdna_direction,watsonp,*dynprogindex);

#if defined(HAVE_SSE4_1) || defined(HAVE_SSE2)
  }
#endif

  pairs = List_reverse(pairs);

  queryjump = (rev_roffsetR-bestrR) - (roffsetL+bestrL) + 1;
  genomejump = (rev_goffset-bestcR) - (goffset+bestcL) + 1;
  /* No need to revise queryjump or genomejump, because the above
     coordinates are internal to the gap. */

  if (queryjump == INSERT_PAIRS && genomejump == INSERT_PAIRS) {
    /* Add cDNA insertion, if any */
    for (k = rev_roffsetR-bestrR; k >= roffsetL+bestrL; k--) {
#ifdef PMAP
      debug(printf("cDNA insertion, Pushing [%d,%d] (%c,-)\n",k,rev_goffset-bestcR+1,inst_rsequenceL[k-roffsetL]));
      pairs = Pairpool_push(pairs,pairpool,k,rev_goffset-bestcR+1,inst_rsequenceL[k-roffsetL],SHORTGAP_COMP,
			    /*genome*/' ',/*genomealt*/' ',*dynprogindex);
#else
      debug(printf("cDNA insertion, Pushing [%d,%d] (%c,-)\n",k,rev_goffset-bestcR+1,rsequenceL[k-roffsetL]));
      pairs = Pairpool_push(pairs,pairpool,k,rev_goffset-bestcR+1,rsequenceL[k-roffsetL],SHORTGAP_COMP,
			    /*genome*/' ',/*genomealt*/' ',*dynprogindex);
#endif
    }
    debug(printf("\n"));


    /* This loop not yet checked for get_genomic_nt giving correct answer */
    for (k = rev_goffset-bestcR; k >= goffset+bestcL; k--) {
      c2 = get_genomic_nt(&c2_alt,k,chroffset,chrhigh,watsonp);
      debug(printf("genome insertion, Pushing [%d,%d] (-,%c)\n",roffsetL+bestrL,k,c2));
#if 0
      assert(c2 == gsequence[k-goffset]);
      pairs = Pairpool_push(pairs,pairpool,roffsetL+bestrL,k,' ',SHORTGAP_COMP,
			    gsequence[k-goffset],/*genomealt*/GENOMEALT_DEFERRED,*dynprogindex);
#else
      pairs = Pairpool_push(pairs,pairpool,roffsetL+bestrL,k,' ',SHORTGAP_COMP,c2,c2_alt,*dynprogindex);
#endif
    }
    debug(printf("\n"));

  } else {

    /* Add gapholder to be solved in the future */
#ifndef NOGAPHOLDER
    debug(printf("Pushing a gap with queryjump = %d, genomejump = %d\n",queryjump,genomejump));
    pairs = Pairpool_push_gapholder(pairs,pairpool,queryjump,genomejump,
				    /*leftpair*/NULL,/*rightpair*/NULL,/*knownp*/false);
#endif
    *incompletep = true;
  }

#if defined(HAVE_SSE4_1) || defined(HAVE_SSE2)
  if (use8p == true) {
    pairs = traceback_8(pairs,&nmatches,&nmismatches,&nopens,&nindels,
			directions8L_nogap,directions8L_Egap,directions8L_Fgap,bestrL,bestcL,
#ifdef PMAP
			inst_rsequenceL,inst_rsequenceL,
#else
			rsequenceL,rsequence_ucL,
#endif
			gsequence,gsequence_alt,roffsetL,goffset,pairpool,/*revp*/false,
			chroffset,chrhigh,cdna_direction,watsonp,*dynprogindex);
  } else {
#endif

    pairs = traceback(pairs,&nmatches,&nmismatches,&nopens,&nindels,
		      directionsL_nogap,directionsL_Egap,directionsL_Fgap,bestrL,bestcL,
#ifdef PMAP
		      inst_rsequenceL,inst_rsequenceL,
#else
		      rsequenceL,rsequence_ucL,
#endif
		      gsequence,gsequence_alt,roffsetL,goffset,pairpool,/*revp*/false,
		      chroffset,chrhigh,cdna_direction,watsonp,*dynprogindex);
#if defined(HAVE_SSE4_1) || defined(HAVE_SSE2)
  }
#endif


#ifdef PMAP
  FREE(inst_rsequence);
#endif

  if (List_length(pairs) == 1) {
    /* Only a gap added */
    pairs = (List_T) NULL;
  }

  if (gsequence_alt != gsequence) {
    FREE(gsequence_alt);
    FREE(rev_gsequence_alt);
  }
  FREE(gsequence);
  FREE(rev_gsequence);

  debug(printf("End of dynprog cDNA gap\n"));

  *dynprogindex += (*dynprogindex > 0 ? +1 : -1);
  return List_reverse(pairs);
}


static List_T
genome_gap_simple (int *finalscore, int *best_introntype, int *new_leftgenomepos, int *new_rightgenomepos,
		   double *left_prob, double *right_prob, int *exonhead, int *nmatches, int *nmismatches,
		   char *rsequence, char *rev_rsequence, int rlength, 
		   char *gsequenceL, char *gsequenceL_alt, char *rev_gsequenceR, char *rev_gsequenceR_alt,
		   int roffset, int rev_roffset, int leftoffset, int rightoffset,
		   Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh, Pairpool_T pairpool,
		   int mismatchtype, int canonical_reward,
		   int cdna_direction, bool watsonp, int dynprogindex, bool halfp) {
  List_T pairs = NULL;
  Pair_T gappair;
  bool result;
  int bestrL, bestrR, rL, rR, r;
  int querycoord, genomecoord;
  int na1, na2, na2_alt, c1, c2, c2_alt;
  char left1, left1_alt, left2, left2_alt, right2, right2_alt, right1, right1_alt;
  int bestscore = 0, bestscoreI = 0, scoreL, scoreI, scoreR, pairscore, score;
  int leftdi, rightdi;
  int *left_known, *right_known;
  int introntype;
  Pairdistance_T **pairdistance_array_type;

  debug(printf("Starting genome_gap_simple with cdna_direction %d and watsonp %d\n",cdna_direction,watsonp));
  pairdistance_array_type = pairdistance_array[mismatchtype];

  left_known = (int *) CALLOC(rlength+1,sizeof(int));
  right_known = (int *) CALLOC(rlength+1,sizeof(int));
  get_known_splicesites(left_known,right_known,/*glengthL*/rlength,/*glengthR*/rlength,
			leftoffset,rightoffset,
			cdna_direction,watsonp,chrnum,chroffset,chrhigh);

  scoreR = 0;
  for (rR = 1; rR < rlength; rR++) {
    na1 = rev_rsequence[1-rR];
    na2 = rev_gsequenceR[1-rR];
    na2_alt = rev_gsequenceR_alt[1-rR];
    pairscore = pairdistance_array_type[na1][na2];
    if ((score = pairdistance_array_type[na1][na2_alt]) > pairscore) {
      pairscore = score;
    }
    scoreR += pairscore;
  }

  scoreL = 0;
  for (rL = 1, rR = rlength-1; rL < rlength; rL++, rR--) {
    na1 = rsequence[rL-1];
    na2 = gsequenceL[rL-1];
    na2_alt = gsequenceL_alt[rL-1];
    pairscore = pairdistance_array_type[na1][na2];
    if ((score = pairdistance_array_type[na1][na2_alt]) > pairscore) {
      pairscore = score;
    }
    scoreL += pairscore;

    /* Read dinucleotides */
    left1 = gsequenceL[rL];
    left1_alt = gsequenceL_alt[rL];
    left2 = gsequenceL[rL+1];
    left2_alt = gsequenceL_alt[rL+1];

    if ((left1 == 'G' || left1_alt == 'G') && (left2 == 'T' || left2_alt == 'T')) {
      leftdi = LEFT_GT;
    } else if ((left1 == 'G' || left1_alt == 'G') && (left2 == 'C' || left2_alt == 'C')) {
      leftdi = LEFT_GC;
    } else if ((left1 == 'A' || left1_alt == 'A') && (left2 == 'T' || left2_alt == 'T')) {
      leftdi = LEFT_AT;
#ifndef PMAP
    } else if ((left1 == 'C' || left1_alt == 'C') && (left2 == 'T' || left2_alt == 'T')) {
      leftdi = LEFT_CT;
#endif
    } else {
      leftdi = 0x00;
    }

    right2 = rev_gsequenceR[-rR-1];
    right2_alt = rev_gsequenceR_alt[-rR-1];
    right1 = rev_gsequenceR[-rR];
    right1_alt = rev_gsequenceR_alt[-rR];
    if ((right2 == 'A' || right2_alt == 'A') && (right1 == 'G' || right1_alt == 'G')) {
      rightdi = RIGHT_AG;
    } else if ((right2 == 'A' || right2_alt == 'A') && (right1 == 'C' || right1_alt == 'C')) {
      rightdi = RIGHT_AC;
#ifndef PMAP
    } else if ((right2 == 'G' || right2_alt == 'G') && (right1 == 'C' || right1_alt == 'C')) {
      rightdi = RIGHT_GC;
    } else if ((right2 == 'A' || right2_alt == 'A') && (right1 == 'T' || right1_alt == 'T')) {
      rightdi = RIGHT_AT;
#endif
    } else {
      rightdi = 0x00;
    }

    scoreI = intron_score(&introntype,leftdi,rightdi,cdna_direction,canonical_reward,/*finalp*/false);
    if ((introntype != NONINTRON || left_known[rL] > 0 || right_known[rR] > 0) &&
	(score = scoreL + left_known[rL] + scoreI + right_known[rR] + scoreR) >= bestscore) {  /* Use >= for jump late */
      debug(printf("At %d left (%c%c) to %d right (%c%c), score is (%d)+(%d)+(%d)+(%d)+(%d) = %d (bestscore)\n",
		   rL,left1,left2,rR,right2,right1,scoreL,left_known[rL],scoreI,right_known[rR],scoreR,
		   scoreL+left_known[rL]+scoreI+right_known[rR]+scoreR));
      bestscore = score;
      bestscoreI = scoreI;
      bestrL = /* *bestcL = */ rL;
      bestrR = /* *bestcR = */ rR;
      *best_introntype = introntype;
    } else {
      debug(printf("At %d left (%c%c) to %d right (%c%c), score is (%d)+(%d)+(%d)+(%d)+(%d) = %d\n",
		   rL,left1,left2,rR,right2,right1,scoreL,left_known[rL],scoreI,right_known[rR],scoreR,
		   scoreL+left_known[rL]+scoreI+right_known[rR]+scoreR));
    }

    /* Subtract pairscore from cumulative scoreR */
    na1 = rev_rsequence[1-rR];
    na2 = rev_gsequenceR[1-rR];
    na2_alt = rev_gsequenceR_alt[1-rR];
    pairscore = pairdistance_array_type[na1][na2];
    if ((score = pairdistance_array_type[na1][na2_alt]) > pairscore) {
      pairscore = score;
    }
    scoreR -= pairscore;
  }

  if (halfp == true) {
    *finalscore = (int) (bestscore - bestscoreI/2);
  } else {
    *finalscore = (int) bestscore;
  }

  if (*finalscore <= 0) {
    result = false;
  } else if (novelsplicingp == true) {
    result = true;
  } else if (splicing_iit == NULL) {
    result = true;
  } else if (donor_typeint >= 0 && acceptor_typeint >= 0) {
    /* If novelsplicingp is false and using splicing at splice site level, require sites to be known */
    if (left_known[bestrL] == 0 || right_known[bestrR] == 0) {
      debug(printf("Novel splicing not allowed, so bridge_intron_gap returning false\n"));
      result = false;
    } else {
      result = true;
    }
  } else {
    /* If novelsplicingp is false and using splicing at splice site level, result was already constrained */
    result = true;
  }

  if (result == true) {
    get_splicesite_probs(&(*left_prob),&(*right_prob),bestrL,bestrR,
			 left_known,right_known,leftoffset,rightoffset,
			 chroffset,chrhigh,cdna_direction,watsonp);
    debug(printf("Probabilities are %f and %f\n",*left_prob,*right_prob));
    if (*left_prob < 0.90 || *right_prob < 0.90) {
      result = false;
    }
  }

  if (result == true) {
    *nmatches = *nmismatches = 0;

    /* Push from left to right, so we don't need to do List_reverse() later */
    for (r = 1; r <= bestrL; r++) {
      querycoord = genomecoord = r-1;

      c1 = rsequence[querycoord];
      c2 = gsequenceL[genomecoord];
      c2_alt = gsequenceL_alt[genomecoord];

      if (c2 == '*') {
	/* Don't push pairs past end of chromosome */
	debug(printf("Don't push pairs past end of chromosome: genomeoffset %u, genomecoord %u, chroffset %u, chrhigh %u, watsonp %d\n",
		     leftoffset,genomecoord,chroffset,chrhigh,watsonp));

      } else if (c1 == c2) {
	debug(printf("Pushing simple %d,%d [%d,%d] (%c,%c) - match\n",
		     r,/*c*/r,roffset+querycoord,leftoffset+genomecoord,c1,c2));
	*nmatches += 1;
	pairs = Pairpool_push(pairs,pairpool,roffset+querycoord,leftoffset+genomecoord,
			      c1,DYNPROG_MATCH_COMP,c2,c2_alt,dynprogindex);
	
      } else if (consistent_array[(int) c1][(int) c2] == true) {
	debug(printf("Pushing simple %d,%d [%d,%d] (%c,%c) - ambiguous\n",
		     r,/*c*/r,roffset+querycoord,leftoffset+genomecoord,c1,c2));
	*nmatches += 1;
	pairs = Pairpool_push(pairs,pairpool,roffset+querycoord,leftoffset+genomecoord,
			      c1,AMBIGUOUS_COMP,c2,c2_alt,dynprogindex);
	
      } else {
	debug(printf("Pushing simple %d,%d [%d,%d] (%c,%c) - mismatch\n",
		     r,/*c*/r,roffset+querycoord,leftoffset+genomecoord,c1,c2));
	*nmismatches += 1;
	pairs = Pairpool_push(pairs,pairpool,roffset+querycoord,leftoffset+genomecoord,
			      c1,MISMATCH_COMP,c2,c2_alt,dynprogindex);
      }
    }

    debug(printf("Pushing a gap\n"));
    *new_leftgenomepos = leftoffset+(bestrL-1);
    *new_rightgenomepos = *exonhead = rightoffset-(bestrR-1);
#ifndef NOGAPHOLDER
    pairs = Pairpool_push_gapholder(pairs,pairpool,/*queryjump*/0,/*genomejump*/(*new_rightgenomepos)-(*new_leftgenomepos)-1,
				    /*leftpair*/NULL,/*rightpair*/NULL,/*knownp*/false);
    gappair = (Pair_T) pairs->first;
    gappair->introntype = introntype;
    gappair->donor_prob = *left_prob;
    gappair->acceptor_prob = *right_prob;
#endif

    for (r = bestrR; r > 0; r--) {
      querycoord = genomecoord = 1-r;

      c1 = rev_rsequence[querycoord];
      c2 = rev_gsequenceR[genomecoord];
      c2_alt = rev_gsequenceR_alt[genomecoord];
      
      if (c2 == '*') {
	/* Don't push pairs past end of chromosome */
	debug(printf("Don't push pairs past end of chromosome: genomeoffset %u, genomecoord %u, chroffset %u, chrhigh %u, watsonp %d\n",
		     rightoffset,genomecoord,chroffset,chrhigh,watsonp));

      } else if (c1 == c2) {
	debug(printf("Pushing simple %d,%d [%d,%d] (%c,%c) - match\n",
		     r,/*c*/r,rev_roffset+querycoord,rightoffset+genomecoord,c1,c2));
	*nmatches += 1;
	pairs = Pairpool_push(pairs,pairpool,rev_roffset+querycoord,rightoffset+genomecoord,
			      c1,DYNPROG_MATCH_COMP,c2,c2_alt,dynprogindex);
	
      } else if (consistent_array[(int) c1][(int) c2] == true) {
	debug(printf("Pushing simple %d,%d [%d,%d] (%c,%c) - ambiguous\n",
		     r,/*c*/r,rev_roffset+querycoord,rightoffset+genomecoord,c1,c2));
	*nmatches += 1;
	pairs = Pairpool_push(pairs,pairpool,rev_roffset+querycoord,rightoffset+genomecoord,
			      c1,AMBIGUOUS_COMP,c2,c2_alt,dynprogindex);
	
      } else {
	debug(printf("Pushing simple %d,%d [%d,%d] (%c,%c) - mismatch\n",
		     r,/*c*/r,rev_roffset+querycoord,rightoffset+genomecoord,c1,c2));
	*nmismatches += 1;
	pairs = Pairpool_push(pairs,pairpool,rev_roffset+querycoord,rightoffset+genomecoord,
			      c1,MISMATCH_COMP,c2,c2_alt,dynprogindex);
      }
    }

  }

  FREE(right_known);
  FREE(left_known);

  return pairs;
}



/* A genome gap is usually an intron.  Sequence 2L and 2R represent
   the two genomic ends of the intron. */
List_T
Dynprog_genome_gap (int *dynprogindex, int *finalscore, int *new_leftgenomepos, int *new_rightgenomepos,
		    double *left_prob, double *right_prob,
		    int *nmatches, int *nmismatches, int *nopens, int *nindels,
		    int *exonhead, int *introntype, T dynprogL, T dynprogR, 
		    char *rsequence, char *sequenceuc1,
		    int rlength, int glengthL, int glengthR, 
		    int roffset, int goffsetL, int rev_goffsetR, 
		    Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
#ifdef PMAP
		    char *queryaaseq,
#endif
		    int cdna_direction, bool watsonp, bool jump_late_p, Pairpool_T pairpool, int extraband_paired,
		    double defect_rate, int maxpeelback, bool halfp, bool finalp, bool use_probabilities_p,
		    int score_threshold, bool splicingp) {
  List_T pairs = NULL;
  Pair_T gappair;
  char *gsequenceL, *gsequenceL_alt, *rev_gsequenceR, *rev_gsequenceR_alt;
#ifdef PMAP
  char *inst_rsequence, *inst_rev_rsequence;
#else
  char *rev_rsequence, *revsequenceuc1;
#endif
  Mismatchtype_T mismatchtype;
  int canonical_reward;
  int rev_roffset, bestrL, bestrR, bestcL, bestcR, lband, uband;
  int maxhorizjump, maxvertjump;
#if defined(HAVE_SSE4_1) || defined(HAVE_SSE2)
  Score8_T **matrix8L, **matrix8R;
  Direction8_T **directions8L_nogap, **directions8L_Egap, **directions8L_Fgap,
    **directions8R_nogap, **directions8R_Egap, **directions8R_Fgap;
  bool use8p;

  Score16_T **matrixL, **matrixR, open, extend;
  Direction16_T **directionsL_nogap, **directionsL_Egap, **directionsL_Fgap,
    **directionsR_nogap, **directionsR_Egap, **directionsR_Fgap;
#else
  Score32_T **matrixL, **matrixR, open, extend;
  Direction32_T **directionsL_nogap, **directionsL_Egap, **directionsL_Fgap,
    **directionsR_nogap, **directionsR_Egap, **directionsR_Fgap;
#endif
  /* int queryjump, genomejump; */

  debug(
	printf("\n");
	printf("%c:  ",*dynprogindex > 0 ? (*dynprogindex-1)%26+'a' : (-(*dynprogindex)-1)%26+'A');
	printf("Aligning genome gap with cdna_direction %d, defect_rate %f\n",cdna_direction,defect_rate);
	printf("At query offset %d-%d, %.*s\n",roffset,roffset+rlength-1,rlength,rsequence));
#ifdef EXTRACT_GENOMICSEG
  debug(printf("At genomic offset %d-%d, %.*s\n",goffsetL,goffsetL+glengthL-1,glengthL,gsequenceL));
  debug(printf("At genomic offset %d-%d, %.*s\n",rev_goffsetR-glengthR+1,rev_goffsetR,glengthR,&(rev_gsequenceR[-glengthR+1])));
#endif
  debug(printf("\n"));

  /* ?check if offsets are too close.  But this eliminates a segment
     of the cDNA.  Should check in stage 3, and do single gap instead. */
  /*
  if (goffsetL+glengthL-1 >= rev_goffsetR-glengthR+1) {
    debug(printf("Bounds don't make sense\n"));
    *finalscore = NEG_INFINITY_16;
    return NULL;
  }
  */

  *nmatches = *nmismatches = *nopens = *nindels = 0;
  *left_prob = *right_prob = 0.0;
  *introntype = NONINTRON;
  if (rlength <= 1) {
    *finalscore = NEG_INFINITY_32;
    return (List_T) NULL;
  }

  if (defect_rate < DEFECT_HIGHQ) {
    mismatchtype = HIGHQ;
    if (rlength > maxpeelback * 4) {
      debug(printf("rlength %d is greater than maxpeelback %d * 4, so using single gap penalties\n",
		   rlength,maxpeelback));
      open = SINGLE_OPEN_HIGHQ;
      extend = SINGLE_EXTEND_HIGHQ;
    } else {
      open = PAIRED_OPEN_HIGHQ;
      extend = PAIRED_EXTEND_HIGHQ;
    }
    if (splicingp == false) {
      canonical_reward = 0;
    } else if (finalp == true) {
      canonical_reward = FINAL_CANONICAL_INTRON_HIGHQ;
    } else {
      canonical_reward = CANONICAL_INTRON_HIGHQ;
    }
    maxhorizjump = MAXHORIZJUMP_HIGHQ;
    maxvertjump = MAXVERTJUMP_HIGHQ;
  } else if (defect_rate < DEFECT_MEDQ) {
    mismatchtype = MEDQ;
    if (rlength > maxpeelback * 4) {
      debug(printf("rlength %d is greater than maxpeelback %d * 4, so using single gap penalties\n",
		   rlength,maxpeelback));
      open = SINGLE_OPEN_MEDQ;
      extend = SINGLE_EXTEND_MEDQ;
    } else {
      open = PAIRED_OPEN_MEDQ;
      extend = PAIRED_EXTEND_MEDQ;
    }
    if (splicingp == false) {
      canonical_reward = 0;
    } else if (finalp == true) {
      canonical_reward = FINAL_CANONICAL_INTRON_MEDQ;
    } else {
      canonical_reward = CANONICAL_INTRON_MEDQ;
    }
    maxhorizjump = MAXHORIZJUMP_MEDQ;
    maxvertjump = MAXVERTJUMP_MEDQ;
  } else {
    mismatchtype = LOWQ;
    if (rlength > maxpeelback * 4) {
      debug(printf("rlength %d is greater than maxpeelback %d * 4, so using single gap penalties\n",
		   rlength,maxpeelback));
      open = SINGLE_OPEN_LOWQ;
      extend = SINGLE_EXTEND_LOWQ;
    } else {
      open = PAIRED_OPEN_LOWQ;
      extend = PAIRED_EXTEND_LOWQ;
    }
    if (splicingp == false) {
      canonical_reward = 0;
    } else if (finalp == true) {
      canonical_reward = FINAL_CANONICAL_INTRON_LOWQ;
    } else {
      canonical_reward = CANONICAL_INTRON_LOWQ;
    }
    maxhorizjump = MAXHORIZJUMP_LOWQ;
    maxvertjump = MAXVERTJUMP_LOWQ;
  }

  if (rlength > dynprogL->max_rlength || glengthL > dynprogL->max_glength) {
    debug(printf("rlength %d or glengthL %d is too long.  Returning NULL\n",rlength,glengthL));
    *new_leftgenomepos = goffsetL-1;
    *new_rightgenomepos = rev_goffsetR+1;
    *exonhead = rev_roffset = roffset+rlength-1;
#ifndef NOGAPHOLDER
    /*
    queryjump = rev_roffset - roffset + 1;
    genomejump = rev_goffsetR - goffsetL + 1;
    pairs = Pairpool_push_gapholder(NULL,pairpool,queryjump,genomejump,false);
    */
#endif
    *dynprogindex += (*dynprogindex > 0 ? +1 : -1);
    *finalscore = NEG_INFINITY_32;
    *introntype = NONINTRON;
    return (List_T) NULL;
  }

  if (rlength > dynprogR->max_rlength || glengthR > dynprogR->max_glength) {
    debug(printf("rlength %d or glengthR %d is too long.  Returning NULL\n",rlength,glengthR));
    *new_leftgenomepos = goffsetL-1;
    *new_rightgenomepos = rev_goffsetR+1;
    *exonhead = rev_roffset = roffset+rlength-1;
#ifndef NOGAPHOLDER
    /*
    queryjump = rev_roffset - roffset + 1;
    genomejump = rev_goffsetR - goffsetL + 1;
    pairs = Pairpool_push_gapholder(NULL,pairpool,queryjump,genomejump,false);
    */
#endif
    *dynprogindex += (*dynprogindex > 0 ? +1 : -1);
    *finalscore = NEG_INFINITY_32;
    *introntype = NONINTRON;
    return (List_T) NULL;
  }

#ifdef PMAP
  inst_rsequence = instantiate_codons(rsequence,queryaaseq,roffset,rlength);
  inst_rev_rsequence = &(inst_rsequence[rlength-1]);
#else
  rev_rsequence = &(rsequence[rlength-1]);
  revsequenceuc1 = &(sequenceuc1[rlength-1]);
#endif
  rev_roffset = roffset+rlength-1;

  if (watsonp) {
    gsequenceL = Genome_get_segment_blocks_right(&gsequenceL_alt,/*left*/chroffset+goffsetL,glengthL,chrhigh,/*revcomp*/false);
    rev_gsequenceR = Genome_get_segment_blocks_left(&rev_gsequenceR_alt,/*left*/chroffset+rev_goffsetR+1,glengthR,chroffset,/*revcomp*/false);
  } else {
    gsequenceL = Genome_get_segment_blocks_left(&gsequenceL_alt,/*left*/chrhigh-goffsetL+1,glengthL,chroffset,/*revcomp*/true);
    rev_gsequenceR = Genome_get_segment_blocks_right(&rev_gsequenceR_alt,/*left*/chrhigh-rev_goffsetR,glengthR,chrhigh,/*revcomp*/true);
  }
  if (gsequenceL == NULL) {
    FREE(rev_gsequenceR);
    *finalscore = NEG_INFINITY_32;
    return (List_T) NULL;
  } else if (rev_gsequenceR == NULL) {
    FREE(gsequenceL);
    *finalscore = NEG_INFINITY_32;
    return (List_T) NULL;
  }


  debug(printf("At genomic offset %d-%d, %.*s\n",goffsetL,goffsetL+glengthL-1,glengthL,gsequenceL));
  debug(printf("At genomic offset %d-%d, %.*s\n",rev_goffsetR-glengthR+1,rev_goffsetR,glengthR,rev_gsequenceR));

  /* In low-identity alignments, the simple procedure can
     lead to multiple mismatches, which will invalidate the intron
     because of its neighborhood */
  if (defect_rate < DEFECT_MEDQ &&
      (pairs = genome_gap_simple(&(*finalscore),&(*introntype),&(*new_leftgenomepos),&(*new_rightgenomepos),
				 &(*left_prob),&(*right_prob),&(*exonhead),&(*nmatches),&(*nmismatches),
				 rsequence,rev_rsequence,rlength,
				 gsequenceL,gsequenceL_alt,&(rev_gsequenceR[glengthR-1]),&(rev_gsequenceR_alt[glengthR-1]),
				 roffset,rev_roffset,goffsetL,rev_goffsetR,
				 chrnum,chroffset,chrhigh,pairpool,mismatchtype,canonical_reward,
				 cdna_direction,watsonp,*dynprogindex,halfp)) != NULL) {
    debug(printf("simple procedure worked\n"));
#ifdef PMAP
    FREE(inst_rsequence);
#endif

    if (gsequenceL_alt != gsequenceL) {
      FREE(gsequenceL_alt);
      FREE(rev_gsequenceR_alt);
    }
    FREE(gsequenceL);
    FREE(rev_gsequenceR);

    *dynprogindex += (*dynprogindex > 0 ? +1 : -1);
    return pairs;		/* Already reversed */
  }


#if defined(HAVE_SSE4_1) || defined(HAVE_SSE2)
  /* Use || because we want the minimum length (which determines the diagonal length) to achieve a score less than 128 */
  if (rlength <= SIMD_MAXLENGTH_EPI8 || (glengthL <= SIMD_MAXLENGTH_EPI8 && glengthR <= SIMD_MAXLENGTH_EPI8)) {
    use8p = true;
  } else {
    use8p = false;
  }
#endif

  compute_bands(&lband,&uband,rlength,glengthL,extraband_paired,/*widebandp*/true);
#if defined(HAVE_SSE4_1) || defined(HAVE_SSE2)
  if (use8p == true) {
    matrix8L = compute_scores_simd_8(&directions8L_nogap,&directions8L_Egap,&directions8L_Fgap,dynprogL,
#ifdef PMAP
				     inst_rsequence,
#else
				     rsequence,
#endif
				     gsequenceL,gsequenceL_alt,rlength,glengthL,
#ifdef DEBUG14
				     goffsetL,chroffset,chrhigh,watsonp,
#endif
				     mismatchtype,(Score8_T) open,(Score8_T) extend,
				     lband,uband,jump_late_p,/*revp*/false);
  } else {
    matrixL = compute_scores_simd_16(&directionsL_nogap,&directionsL_Egap,&directionsL_Fgap,dynprogL,
#ifdef PMAP
				     inst_rsequence,
#else
				     rsequence,
#endif
				     gsequenceL,gsequenceL_alt,rlength,glengthL,
#ifdef DEBUG14
				     goffsetL,chroffset,chrhigh,watsonp,
#endif
				     mismatchtype,(Score16_T) open, (Score16_T) extend,
				     lband,uband,jump_late_p,/*revp*/false);
  }

#else
  matrixL = compute_scores_standard(&directionsL_nogap,&directionsL_Egap,&directionsL_Fgap,dynprogL,
#ifdef PMAP
				    inst_rsequence,
#else
				    rsequence,
#endif
				    gsequenceL,gsequenceL_alt,rlength,glengthL,
				    goffsetL,chroffset,chrhigh,watsonp,
				    mismatchtype,open,extend,lband,uband,
				    jump_late_p,/*revp*/false);
#endif
  
  compute_bands(&lband,&uband,rlength,glengthR,extraband_paired,/*widebandp*/true);
#if defined(HAVE_SSE4_1) || defined(HAVE_SSE2)
  if (use8p == true) {
    matrix8R = compute_scores_simd_8(&directions8R_nogap,&directions8R_Egap,&directions8R_Fgap,dynprogR,
#ifdef PMAP
				     inst_rev_rsequence,
#else
				     rev_rsequence,
#endif
				     &(rev_gsequenceR[glengthR-1]),&(rev_gsequenceR_alt[glengthR-1]),
				     rlength,glengthR,
#ifdef DEBUG14
				     rev_goffsetR,chroffset,chrhigh,watsonp,
#endif
				     mismatchtype,(Score8_T) open,(Score8_T) extend,
				     lband,uband,/*for revp true*/!jump_late_p,/*revp*/true);
  } else {
    matrixR = compute_scores_simd_16(&directionsR_nogap,&directionsR_Egap,&directionsR_Fgap,dynprogR,
#ifdef PMAP
				     inst_rev_rsequence,
#else
				     rev_rsequence,
#endif
				     &(rev_gsequenceR[glengthR-1]),&(rev_gsequenceR_alt[glengthR-1]),
				     rlength,glengthR,
#ifdef DEBUG14
				     rev_goffsetR,chroffset,chrhigh,watsonp,
#endif
				     mismatchtype,(Score16_T) open, (Score16_T) extend,
				     lband,uband,/*for revp true*/!jump_late_p,/*revp*/true);
  }

#else
  matrixR = compute_scores_standard(&directionsR_nogap,&directionsR_Egap,&directionsR_Fgap,dynprogR,
#ifdef PMAP
				    inst_rev_rsequence,
#else
				    rev_rsequence,
#endif
				    &(rev_gsequenceR[glengthR-1]),&(rev_gsequenceR_alt[glengthR-1]),
				    rlength,glengthR,
				    rev_goffsetR,chroffset,chrhigh,watsonp,
				    mismatchtype,open,extend,
				    lband,uband,/*for revp true*/!jump_late_p,/*revp*/true);
#endif

#if defined(HAVE_SSE4_1) || defined(HAVE_SSE2)
  if (use8p == true) {
    if (bridge_intron_gap_8(&(*finalscore),&bestrL,&bestrR,&bestcL,&bestcR,
			    &(*introntype),&(*left_prob),&(*right_prob),
			    matrix8L,matrix8R,directions8L_nogap,directions8R_nogap,
			    gsequenceL,gsequenceL_alt,&(rev_gsequenceR[glengthR-1]),&(rev_gsequenceR_alt[glengthR-1]),
			    goffsetL,rev_goffsetR,rlength,glengthL,glengthR,
			    cdna_direction,watsonp,extraband_paired,defect_rate,
			    canonical_reward,maxhorizjump,maxvertjump,goffsetL,rev_goffsetR,
			    chrnum,chroffset,chrhigh,halfp,finalp,use_probabilities_p,
			    score_threshold,jump_late_p) == false) {
      if (gsequenceL_alt != gsequenceL) {
	FREE(gsequenceL_alt);
	FREE(rev_gsequenceR_alt);
      }
      FREE(gsequenceL);
      FREE(rev_gsequenceR);

      return (List_T) NULL;
    } else {
      *new_leftgenomepos = goffsetL+(bestcL-1);
      *new_rightgenomepos = rev_goffsetR-(bestcR-1);
      debug(printf("New leftgenomepos = %d, New rightgenomepos = %d\n",*new_leftgenomepos,*new_rightgenomepos));

      *exonhead = rev_roffset-(bestrR-1);

      pairs = traceback_8(NULL,&(*nmatches),&(*nmismatches),&(*nopens),&(*nindels),
			  directions8R_nogap,directions8R_Egap,directions8R_Fgap,bestrR,bestcR,
#ifdef PMAP
			  inst_rev_rsequence,inst_rev_rsequence,
#else
			  rev_rsequence,revsequenceuc1,
#endif
			  &(rev_gsequenceR[glengthR-1]),&(rev_gsequenceR_alt[glengthR-1]),
			  rev_roffset,rev_goffsetR,pairpool,/*revp*/true,
			  chroffset,chrhigh,cdna_direction,watsonp,*dynprogindex);
      pairs = List_reverse(pairs);

      /* queryjump = (rev_roffset-bestrR) - (roffset+bestrL) + 1; */
      /* genomejump = (rev_goffsetR-bestcR) - (goffsetL+bestcL) + 1; */
      /* No need to revise queryjump or genomejump, because the above
	 coordinates are internal to the gap. */

      debug(printf("Pushing a gap with genomejump %d, introntype %s, prob %f and %f\n",
		   (*new_rightgenomepos)-(*new_leftgenomepos)-1,
		   Intron_type_string(*introntype),*left_prob,*right_prob));
#ifndef NOGAPHOLDER
      pairs = Pairpool_push_gapholder(pairs,pairpool,/*queryjump*/(rev_roffset-bestrR) - (roffset+bestrL) + 1,
				      /*genomejump*/(*new_rightgenomepos)-(*new_leftgenomepos)-1,
				      /*leftpair*/NULL,/*rightpair*/NULL,/*knownp*/false);
      gappair = (Pair_T) pairs->first;
      gappair->introntype = *introntype;
      gappair->donor_prob = *left_prob;
      gappair->acceptor_prob = *right_prob;
#endif

      pairs = traceback_8(pairs,&(*nmatches),&(*nmismatches),&(*nopens),&(*nindels),
			  directions8L_nogap,directions8L_Egap,directions8L_Fgap,bestrL,bestcL,
#ifdef PMAP
			  inst_rsequence,inst_rsequence,
#else
			  rsequence,sequenceuc1,
#endif
			  gsequenceL,gsequenceL_alt,
			  roffset,goffsetL,pairpool,/*revp*/false,
			  chroffset,chrhigh,cdna_direction,watsonp,*dynprogindex);

#ifdef PMAP
      FREE(inst_rsequence);
#endif

      if (List_length(pairs) == 1) {
	/* Only a gap inserted */
	pairs = (List_T) NULL;
      }

      if (gsequenceL_alt != gsequenceL) {
	FREE(gsequenceL_alt);
	FREE(rev_gsequenceR_alt);
      }
      FREE(gsequenceL);
      FREE(rev_gsequenceR);

      debug(printf("End of dynprog genome gap\n"));

      *dynprogindex += (*dynprogindex > 0 ? +1 : -1);
      return List_reverse(pairs);
    }

  } else {
#endif  /* HAVE_SSE4_1 || HAVE_SSE2 */

    if (bridge_intron_gap(&(*finalscore),&bestrL,&bestrR,&bestcL,&bestcR,
			  &(*introntype),&(*left_prob),&(*right_prob),
			  matrixL,matrixR,directionsL_nogap,directionsR_nogap,
			  gsequenceL,gsequenceL_alt,&(rev_gsequenceR[glengthR-1]),&(rev_gsequenceR_alt[glengthR-1]),
			  goffsetL,rev_goffsetR,rlength,glengthL,glengthR,
			  cdna_direction,watsonp,extraband_paired,defect_rate,
			  canonical_reward,maxhorizjump,maxvertjump,goffsetL,rev_goffsetR,
			  chrnum,chroffset,chrhigh,halfp,finalp,use_probabilities_p,
			  score_threshold,jump_late_p) == false) {
      if (gsequenceL_alt != gsequenceL) {
	FREE(gsequenceL_alt);
	FREE(rev_gsequenceR_alt);
      }
      FREE(gsequenceL);
      FREE(rev_gsequenceR);

      return (List_T) NULL;

    } else {
      *new_leftgenomepos = goffsetL+(bestcL-1);
      *new_rightgenomepos = rev_goffsetR-(bestcR-1);
      debug(printf("New leftgenomepos = %d, New rightgenomepos = %d\n",*new_leftgenomepos,*new_rightgenomepos));

      *exonhead = rev_roffset-(bestrR-1);

      pairs = traceback(NULL,&(*nmatches),&(*nmismatches),&(*nopens),&(*nindels),
			directionsR_nogap,directionsR_Egap,directionsR_Fgap,bestrR,bestcR,
#ifdef PMAP
			inst_rev_rsequence,inst_rev_rsequence,
#else
			rev_rsequence,revsequenceuc1,
#endif
			&(rev_gsequenceR[glengthR-1]),&(rev_gsequenceR_alt[glengthR-1]),
			rev_roffset,rev_goffsetR,pairpool,/*revp*/true,
			chroffset,chrhigh,cdna_direction,watsonp,*dynprogindex);
      pairs = List_reverse(pairs);

      /* queryjump = (rev_roffset-bestrR) - (roffset+bestrL) + 1; */
      /* genomejump = (rev_goffsetR-bestcR) - (goffsetL+bestcL) + 1; */
      /* No need to revise queryjump or genomejump, because the above
	 coordinates are internal to the gap. */

      debug(printf("Pushing a gap with genomejump %d, introntype %s, prob %f and %f\n",
		   (*new_rightgenomepos)-(*new_leftgenomepos)-1,
		   Intron_type_string(*introntype),*left_prob,*right_prob));
#ifndef NOGAPHOLDER
      pairs = Pairpool_push_gapholder(pairs,pairpool,/*queryjump*/(rev_roffset-bestrR) - (roffset+bestrL) + 1,
				      /*genomejump*/(*new_rightgenomepos)-(*new_leftgenomepos)-1,
				      /*leftpair*/NULL,/*rightpair*/NULL,/*knownp*/false);
      gappair = (Pair_T) pairs->first;
      gappair->introntype = *introntype;
      gappair->donor_prob = *left_prob;
      gappair->acceptor_prob = *right_prob;
#endif

      pairs = traceback(pairs,&(*nmatches),&(*nmismatches),&(*nopens),&(*nindels),
			directionsL_nogap,directionsL_Egap,directionsL_Fgap,bestrL,bestcL,
#ifdef PMAP
			inst_rsequence,inst_rsequence,
#else
			rsequence,sequenceuc1,
#endif
			gsequenceL,gsequenceL_alt,roffset,goffsetL,pairpool,/*revp*/false,
			chroffset,chrhigh,cdna_direction,watsonp,*dynprogindex);
#ifdef PMAP
      FREE(inst_rsequence);
#endif

      if (List_length(pairs) == 1) {
	/* Only a gap inserted */
	pairs = (List_T) NULL;
      }

      if (gsequenceL_alt != gsequenceL) {
	FREE(gsequenceL_alt);
	FREE(rev_gsequenceR_alt);
      }
      FREE(gsequenceL);
      FREE(rev_gsequenceR);

      debug(printf("End of dynprog genome gap\n"));

      *dynprogindex += (*dynprogindex > 0 ? +1 : -1);
      return List_reverse(pairs);
    }

#if defined(HAVE_SSE4_1) || defined(HAVE_SSE2)
  }
#endif

}





static int
binary_search (int lowi, int highi, Univcoord_T *positions, Univcoord_T goal) {
  int middlei;

  debug10(printf("entered binary search with lowi=%d, highi=%d, goal=%u\n",lowi,highi,goal));

  while (lowi < highi) {
    middlei = (lowi+highi)/2;
    debug10(printf("  binary: %d:%u %d:%u %d:%u   vs. %u\n",
		  lowi,positions[lowi],middlei,positions[middlei],
		  highi,positions[highi],goal));
    if (goal < positions[middlei]) {
      highi = middlei;
    } else if (goal > positions[middlei]) {
      lowi = middlei + 1;
    } else {
      debug10(printf("binary search returns %d\n",middlei));
      return middlei;
    }
  }

  debug10(printf("binary search returns %d\n",highi));
  return highi;
}



List_T
Dynprog_end5_gap (int *dynprogindex, int *finalscore, int *nmatches, int *nmismatches, 
		  int *nopens, int *nindels, T dynprog, 
		  char *rev_rsequence, char *revsequenceuc1,
		  int rlength, int glength, int rev_roffset, int rev_goffset, 
		  Univcoord_T chroffset, Univcoord_T chrhigh,
#ifdef PMAP
		  char *queryaaseq,
#endif
		  int cdna_direction, bool watsonp, bool jump_late_p, Pairpool_T pairpool,
		  int extraband_end, double defect_rate, Endalign_T endalign) {
  List_T pairs = NULL;
  char *rev_gsequence, *rev_gsequence_alt;
#ifdef PMAP
  char *inst_rsequence, *inst_rev_rsequence;
#endif
  Pair_T pair;
  Mismatchtype_T mismatchtype; 
  int bestr, bestc, lband, uband;
#if defined(HAVE_SSE4_1) || defined(HAVE_SSE2)
  Direction8_T **directions8_nogap, **directions8_Egap, **directions8_Fgap;
  Score8_T **matrix8;
  bool use8p = false;

  Score16_T **matrix, open, extend;
  Direction16_T **directions_nogap, **directions_Egap, **directions_Fgap;
#else
  Score32_T **matrix, open, extend;
  Direction32_T **directions_nogap, **directions_Egap, **directions_Fgap;
#endif
#ifdef PMAP
  int initpos, initmod;
#endif

  debug6(
	printf("%c:  ",*dynprogindex > 0 ? (*dynprogindex-1)%26+'a' : (-(*dynprogindex)-1)%26+'A');
	printf("Aligning 5' end gap with endalign %d\n",endalign);
	);

  mismatchtype = ENDQ;
  if (defect_rate < DEFECT_HIGHQ) {
    open = END_OPEN_HIGHQ;
    extend = END_EXTEND_HIGHQ;
  } else if (defect_rate < DEFECT_MEDQ) {
    open = END_OPEN_MEDQ;
    extend = END_EXTEND_MEDQ;
  } else {
    open = END_OPEN_LOWQ;
    extend = END_EXTEND_LOWQ;
  }

  /* We can just chop lengths to work, since we're not constrained on 5' end */
  if (rlength <= 0) {
    /* Needed to avoid abort by Matrix16_alloc */
    debug6(printf("rlength %d <= 0, so returning NULL\n",rlength));
    *nmatches = *nmismatches = *nopens = *nindels = 0;
    *finalscore = 0;
    return (List_T) NULL;
  } else if (endalign == QUERYEND_NOGAPS) {
    /* Don't shorten rlength */
  } else if (rlength > dynprog->max_rlength) {
    debug6(printf("rlength %d is too long.  Chopping to %d\n",rlength,dynprog->max_rlength));
    rlength = dynprog->max_rlength;
  }
  if (glength <= 0) {
    /* Needed to avoid abort by Matrix16_alloc */
    debug6(printf("glength %d <= 0, so returning NULL\n",glength));
    *nmatches = *nmismatches = *nopens = *nindels = 0;
    *finalscore = 0;
    return (List_T) NULL;
  } else if (endalign == QUERYEND_NOGAPS) {
    /* Don't shorten glength */
  } else if (glength > dynprog->max_glength) {
    debug6(printf("glength %d is too long.  Chopping to %d\n",glength,dynprog->max_glength));
    glength = dynprog->max_glength;
  }

#ifdef PMAP
  inst_rsequence = instantiate_codons(&(rev_rsequence[-rlength+1]),queryaaseq,rev_roffset-rlength+1,rlength);
  inst_rev_rsequence = &(inst_rsequence[rlength-1]);
  debug6(printf("At query offset %d-%d, %.*s\n",rev_roffset-rlength+1,rev_roffset,rlength,&(inst_rev_rsequence[-rlength+1])));
#else
  debug6(printf("At query offset %d-%d, %.*s\n",rev_roffset-rlength+1,rev_roffset,rlength,&(rev_rsequence[-rlength+1])));
#endif

#ifdef EXTRACT_GENOMICSEG
  debug6(printf("At genomic offset %d-%d, %.*s\n",
		rev_goffset-glength+1,rev_goffset,glength,&(rev_gsequence[-glength+1])));
#endif


  if (watsonp) {
    rev_gsequence = Genome_get_segment_blocks_left(&rev_gsequence_alt,/*left*/chroffset+rev_goffset+1,glength,chroffset,/*revcomp*/false);
  } else {
    rev_gsequence = Genome_get_segment_blocks_right(&rev_gsequence_alt,/*left*/chrhigh-rev_goffset,glength,chrhigh,/*revcomp*/true);
  }
  if (rev_gsequence == NULL) {
    *nmatches = *nmismatches = *nopens = *nindels = 0;
    *finalscore = 0;
    return (List_T) NULL;
  }

  if (endalign == QUERYEND_GAP || endalign == BEST_LOCAL) {
    compute_bands(&lband,&uband,rlength,glength,extraband_end,/*widebandp*/true);
#if defined(HAVE_SSE4_1) || defined(HAVE_SSE2)
    /* Use || because we want the minimum length (which determines the diagonal length) to achieve a score less than 128 */
    if (rlength <= SIMD_MAXLENGTH_EPI8 || glength <= SIMD_MAXLENGTH_EPI8) {
      use8p = true;
      matrix8 = compute_scores_simd_8(&directions8_nogap,&directions8_Egap,&directions8_Fgap,dynprog,
#ifdef PMAP
				      inst_rev_rsequence,
#else
				      rev_rsequence,
#endif
				      &(rev_gsequence[glength-1]),&(rev_gsequence_alt[glength-1]),
				      rlength,glength,
#ifdef DEBUG14
				      rev_goffset,chroffset,chrhigh,watsonp,
#endif
				      mismatchtype,(Score8_T) open,(Score8_T) extend,
				      lband,uband,/*for revp true*/!jump_late_p,/*revp*/true);
      find_best_endpoint_8(&(*finalscore),&bestr,&bestc,matrix8,rlength,glength,extraband_end,
			   !jump_late_p);

    } else {
      matrix = compute_scores_simd_16(&directions_nogap,&directions_Egap,&directions_Fgap,dynprog,
#ifdef PMAP
				      inst_rev_rsequence,
#else
				      rev_rsequence,
#endif
				      &(rev_gsequence[glength-1]),&(rev_gsequence_alt[glength-1]),
				      rlength,glength,
#ifdef DEBUG14
				      rev_goffset,chroffset,chrhigh,watsonp,
#endif
				      mismatchtype,(Score16_T) open, (Score16_T) extend,
				      lband,uband,/*for revp true*/!jump_late_p,/*revp*/true);
      find_best_endpoint(&(*finalscore),&bestr,&bestc,matrix,rlength,glength,extraband_end,
			 !jump_late_p);
    }

#else

    matrix = compute_scores_standard(&directions_nogap,&directions_Egap,&directions_Fgap,dynprog,
#ifdef PMAP
				     inst_rev_rsequence,
#else
				     rev_rsequence,
#endif
				     &(rev_gsequence[glength-1]),&(rev_gsequence_alt[glength-1]),
				     rlength,glength,
				     rev_goffset,chroffset,chrhigh,watsonp,
				     mismatchtype,open,extend,
				     lband,uband,/*for revp true*/!jump_late_p,/*revp*/true);
    find_best_endpoint(&(*finalscore),&bestr,&bestc,matrix,rlength,glength,extraband_end,
		       !jump_late_p);
#endif

  } else if (endalign == QUERYEND_INDELS) {
    compute_bands(&lband,&uband,rlength,glength,extraband_end,/*widebandp*/true);
#if defined(HAVE_SSE4_1) || defined(HAVE_SSE2)
    /* Use || because we want the minimum length (which determines the diagonal length) to achive a score less than 128 */
    if (rlength <= SIMD_MAXLENGTH_EPI8 || glength <= SIMD_MAXLENGTH_EPI8) {
      use8p = true;
      matrix8 = compute_scores_simd_8(&directions8_nogap,&directions8_Egap,&directions8_Fgap,dynprog,
#ifdef PMAP
				      inst_rev_rsequence,
#else
				      rev_rsequence,
#endif
				      &(rev_gsequence[glength-1]),&(rev_gsequence_alt[glength-1]),
				      rlength,glength,
#ifdef DEBUG14
				      rev_goffset,chroffset,chrhigh,watsonp,
#endif
				      mismatchtype,(Score8_T) open,(Score8_T) extend,
				      lband,uband,/*for revp true*/!jump_late_p,/*revp*/true);
      find_best_endpoint_to_queryend_indels_8(&(*finalscore),&bestr,&bestc,matrix8,rlength,glength,extraband_end,
					      !jump_late_p);
      /* *finalscore = 0 -- Splicetrie procedures need to know finalscore */

    } else {
      matrix = compute_scores_simd_16(&directions_nogap,&directions_Egap,&directions_Fgap,dynprog,
#ifdef PMAP
				      inst_rev_rsequence,
#else
				      rev_rsequence,
#endif
				      &(rev_gsequence[glength-1]),&(rev_gsequence_alt[glength-1]),
				      rlength,glength,
#ifdef DEBUG14
				      rev_goffset,chroffset,chrhigh,watsonp,
#endif
				      mismatchtype,(Score16_T) open, (Score16_T) extend,
				      lband,uband,/*for revp true*/!jump_late_p,/*revp*/true);
      find_best_endpoint_to_queryend_indels(&(*finalscore),&bestr,&bestc,matrix,rlength,glength,extraband_end,
					    !jump_late_p);
      /* *finalscore = 0 -- Splicetrie procedures need to know finalscore */
    }

#else
    matrix = compute_scores_standard(&directions_nogap,&directions_Egap,&directions_Fgap,dynprog,
#ifdef PMAP
				     inst_rev_rsequence,
#else
				     rev_rsequence,
#endif
				     &(rev_gsequence[glength-1]),&(rev_gsequence_alt[glength-1]),
				     rlength,glength,
				     rev_goffset,chroffset,chrhigh,watsonp,
				     mismatchtype,open,extend,lband,uband,
				     /*for revp true*/!jump_late_p,/*revp*/true);
    find_best_endpoint_to_queryend_indels(&(*finalscore),&bestr,&bestc,matrix,rlength,glength,extraband_end,
					  !jump_late_p);
    /* *finalscore = 0 -- Splicetrie procedures need to know finalscore */

#endif


  } else if (endalign == QUERYEND_NOGAPS) {
    find_best_endpoint_to_queryend_nogaps(&bestr,&bestc,rlength,glength);
    /* *finalscore = 0;	-- Splicetrie procedures need to know finalscore */

  } else {
    fprintf(stderr,"Unexpected endalign value %d\n",endalign);
    abort();
  }


#ifdef PMAP
  initpos = rev_roffset-(bestc-1);
  debug6(printf("Initial query pos is %d\n",initpos));
  if ((initmod = initpos % 3) > 0) {
    if (bestr + initmod < rlength && bestc + initmod < glength) {
      debug6(printf("Rounding down by %d\n",initmod));
      bestr += initmod;
      bestc += initmod;
    }
  }
#endif

  *nmatches = *nmismatches = *nopens = *nindels = 0;
  if (endalign == QUERYEND_NOGAPS) {
    pairs = traceback_nogaps(NULL,&(*nmatches),&(*nmismatches),bestr,bestc,
#ifdef PMAP
			     inst_rev_rsequence,inst_rev_rsequence,
#else
			     rev_rsequence,revsequenceuc1,
#endif
			     &(rev_gsequence[glength-1]),&(rev_gsequence_alt[glength-1]),
			     rev_roffset,rev_goffset,pairpool,
#ifdef DEBUG14
			     chroffset,chrhigh,watsonp,
#endif
			     /*revp*/true,*dynprogindex);
    *finalscore = (*nmatches)*FULLMATCH + (*nmismatches)*MISMATCH_ENDQ;

#if defined(HAVE_SSE4_1) || defined(HAVE_SSE2)
  } else if (use8p == true) {
    pairs = traceback_8(NULL,&(*nmatches),&(*nmismatches),&(*nopens),&(*nindels),
			directions8_nogap,directions8_Egap,directions8_Fgap,bestr,bestc,
#ifdef PMAP
			inst_rev_rsequence,inst_rev_rsequence,
#else
			rev_rsequence,revsequenceuc1,
#endif
			&(rev_gsequence[glength-1]),&(rev_gsequence_alt[glength-1]),
			rev_roffset,rev_goffset,pairpool,/*revp*/true,
			chroffset,chrhigh,cdna_direction,watsonp,*dynprogindex);
#endif

  } else {
    pairs = traceback(NULL,&(*nmatches),&(*nmismatches),&(*nopens),&(*nindels),
		      directions_nogap,directions_Egap,directions_Fgap,bestr,bestc,
#ifdef PMAP
		      inst_rev_rsequence,inst_rev_rsequence,
#else
		      rev_rsequence,revsequenceuc1,
#endif
		      &(rev_gsequence[glength-1]),&(rev_gsequence_alt[glength-1]),
		      rev_roffset,rev_goffset,pairpool,/*revp*/true,
		      chroffset,chrhigh,cdna_direction,watsonp,*dynprogindex);
  }


  if ((endalign == QUERYEND_GAP || endalign == BEST_LOCAL) && (*nmatches + 1) < *nmismatches) {
    *finalscore = 0;
    /* No need to free pairs */
    pairs = NULL;
  } else {
    /* Add 1 to count the match already in the alignment */
    pairs = List_reverse(pairs); /* Look at 5' end to remove excess gaps */
    while (pairs != NULL && (pair = List_head(pairs)) && pair->comp == INDEL_COMP) {
      pairs = List_next(pairs);
    }
  }

  /*
    Directions_free(directions);
    Matrix_free(matrix);
  */
  
#ifdef PMAP
  FREE(inst_rsequence);
#endif

  if (rev_gsequence_alt != rev_gsequence) {
    FREE(rev_gsequence_alt);
  }
  FREE(rev_gsequence);

  debug6(printf("End of dynprog end5 gap\n\n"));

  *dynprogindex += (*dynprogindex > 0 ? +1 : -1);
  return List_reverse(pairs);
}



/* rev_gsequence is the splicejunction */
List_T
Dynprog_end5_splicejunction (int *dynprogindex, int *finalscore, int *missscore,
			     int *nmatches, int *nmismatches, int *nopens, int *nindels, T dynprog, 
			     char *rev_rsequence, char *revsequenceuc1,
			     char *rev_gsequence, char *rev_gsequence_uc, char *rev_gsequence_alt,
			     int rlength, int glength, int rev_roffset, int rev_goffset_anchor, int rev_goffset_far,
			     Univcoord_T chroffset, Univcoord_T chrhigh,
#ifdef PMAP
			     char *queryaaseq,
#endif
			     int cdna_direction, bool watsonp, bool jump_late_p, Pairpool_T pairpool,
			     int extraband_end, double defect_rate, int contlength) {
  List_T pairs = NULL;
#ifdef PMAP
  char *inst_rsequence, *inst_rev_rsequence;
#endif
  Pair_T pair;
  Mismatchtype_T mismatchtype;
  int bestr, bestc, lband, uband;
#if defined(HAVE_SSE4_1) || defined(HAVE_SSE2)
  Score8_T **matrix8;
  Direction8_T **directions8_nogap, **directions8_Egap, **directions8_Fgap;
  bool use8p = false;

  Score16_T **matrix, open, extend;
  Direction16_T **directions_nogap, **directions_Egap, **directions_Fgap;
#else
  Score32_T **matrix, open, extend;
  Direction32_T **directions_nogap, **directions_Egap, **directions_Fgap;
#endif
#ifdef PMAP
  int initpos, initmod;
#endif

  debug6(
	printf("%c:  ",*dynprogindex > 0 ? (*dynprogindex-1)%26+'a' : (-(*dynprogindex)-1)%26+'A');
	printf("Aligning 5' end gap with endalign QUERYEND_NOGAPS\n");
	);

  mismatchtype = ENDQ;
  if (defect_rate < DEFECT_HIGHQ) {
    open = END_OPEN_HIGHQ;
    extend = END_EXTEND_HIGHQ;
  } else if (defect_rate < DEFECT_MEDQ) {
    open = END_OPEN_MEDQ;
    extend = END_EXTEND_MEDQ;
  } else {
    open = END_OPEN_LOWQ;
    extend = END_EXTEND_LOWQ;
  }

  /* We can just chop lengths to work, since we're not constrained on 5' end */
  if (rlength <= 0 || rlength > dynprog->max_rlength) {
    /* Needed to avoid abort by Matrix16_alloc */
    debug6(printf("rlength %d <= 0, so returning NULL\n",rlength));
    *nmatches = *nmismatches = *nopens = *nindels = 0;
    *finalscore = 0;
    *missscore = -100;
    return (List_T) NULL;
  }
  if (glength <= 0 || glength > dynprog->max_glength) {
    /* Needed to avoid abort by Matrix16_alloc */
    debug6(printf("glength %d <= 0, so returning NULL\n",glength));
    *nmatches = *nmismatches = *nopens = *nindels = 0;
    *finalscore = 0;
    *missscore = -100;
    return (List_T) NULL;
  }

#ifdef PMAP
  inst_rsequence = instantiate_codons(&(rev_rsequence[-rlength+1]),queryaaseq,rev_roffset-rlength+1,rlength);
  inst_rev_rsequence = &(inst_rsequence[rlength-1]);
  debug6(printf("At query offset %d-%d, %.*s\n",rev_roffset-rlength+1,rev_roffset,rlength,&(inst_rev_rsequence[-rlength+1])));
#else
  debug6(printf("At query offset %d-%d, %.*s\n",rev_roffset-rlength+1,rev_roffset,rlength,&(rev_rsequence[-rlength+1])));
#endif

  debug6(printf("At genomic offset %d-%d, %.*s\n",
		rev_goffset_anchor-glength+1,rev_goffset_anchor,glength,&(rev_gsequence[-glength+1])));
  
#ifdef PMAP
  initpos = rev_roffset-(bestc-1);
  debug6(printf("Initial query pos is %d\n",initpos));
  if ((initmod = initpos % 3) > 0) {
    if (bestr + initmod < rlength && bestc + initmod < glength) {
      debug6(printf("Rounding down by %d\n",initmod));
      bestr += initmod;
      bestc += initmod;
    }
  }
#endif

  compute_bands(&lband,&uband,rlength,glength,extraband_end,/*widebandp*/true);
#if defined(HAVE_SSE4_1) || defined(HAVE_SSE2)
  /* Use || because we want the minimum length (which determines the diagonal length) to achieve a score less than 128 */
  if (rlength <= SIMD_MAXLENGTH_EPI8 || glength <= SIMD_MAXLENGTH_EPI8) {
    use8p = true;
    matrix8 = compute_scores_simd_8(&directions8_nogap,&directions8_Egap,&directions8_Fgap,dynprog,
#ifdef PMAP
				    inst_rev_rsequence,rev_gsequence,rev_gsequence_alt,
#else
				    rev_rsequence,rev_gsequence_uc,rev_gsequence_alt,
#endif
				    rlength,glength,
#ifdef DEBUG14
				    /*goffset*/0,chroffset,chrhigh,watsonp,
#endif
				    mismatchtype,(Score8_T) open,(Score8_T) extend,
				    lband,uband,/*for revp true*/!jump_late_p,/*revp*/true);
    find_best_endpoint_to_queryend_indels_8(&(*finalscore),&bestr,&bestc,matrix8,rlength,glength,extraband_end,
					    !jump_late_p);

  } else {
    matrix = compute_scores_simd_16(&directions_nogap,&directions_Egap,&directions_Fgap,dynprog,
#ifdef PMAP
				    inst_rev_rsequence,rev_gsequence,rev_gsequence_alt,
#else
				    rev_rsequence,rev_gsequence_uc,rev_gsequence_alt,
#endif
				    rlength,glength,
#ifdef DEBUG14
				    /*goffset*/0,chroffset,chrhigh,watsonp,
#endif
				    mismatchtype,(Score16_T) open, (Score16_T) extend,
				    lband,uband,/*for revp true*/!jump_late_p,/*revp*/true);
    find_best_endpoint_to_queryend_indels(&(*finalscore),&bestr,&bestc,matrix,rlength,glength,extraband_end,
					  !jump_late_p);
  }

#else
  matrix = compute_scores_standard(&directions_nogap,&directions_Egap,&directions_Fgap,dynprog,
#ifdef PMAP
				   inst_rev_rsequence,rev_gsequence,rev_gsequence_alt,
#else
				   rev_rsequence,rev_gsequence,rev_gsequence_alt,
#endif
				   rlength,glength,
				   /*goffset*/0,chroffset,chrhigh,watsonp,
				   mismatchtype,open,extend,
				   lband,uband,/*for revp true*/!jump_late_p,/*revp*/true);
  find_best_endpoint_to_queryend_indels(&(*finalscore),&bestr,&bestc,matrix,rlength,glength,extraband_end,
					!jump_late_p);
#endif

  *nmatches = *nmismatches = *nopens = *nindels = 0;
#if defined(HAVE_SSE4_1) || defined(HAVE_SSE2)
  if (use8p == true) {
    pairs = traceback_local_8(NULL,&(*nmatches),&(*nmismatches),&(*nopens),&(*nindels),
			      directions8_nogap,directions8_Egap,directions8_Fgap,&bestr,&bestc,/*endc*/contlength,
#ifdef PMAP
			      inst_rev_rsequence,inst_rev_rsequence,
#else
			      rev_rsequence,revsequenceuc1,
#endif
			      rev_gsequence,rev_gsequence_uc,rev_gsequence_alt,
			      rev_roffset,rev_goffset_far,pairpool,/*revp*/true,
			      chroffset,chrhigh,cdna_direction,watsonp,*dynprogindex);

    pairs = Pairpool_push_gapholder(pairs,pairpool,/*queryjump*/0,/*genomejump*/rev_goffset_anchor - rev_goffset_far,
				    /*leftpair*/NULL,/*rightpair*/NULL,/*knownp*/true);

    pairs = traceback_local_8(pairs,&(*nmatches),&(*nmismatches),&(*nopens),&(*nindels),
			      directions8_nogap,directions8_Egap,directions8_Fgap,&bestr,&bestc,/*endc*/0,
#ifdef PMAP
			      inst_rev_rsequence,inst_rev_rsequence,
#else
			      rev_rsequence,revsequenceuc1,
#endif
			      rev_gsequence,rev_gsequence_uc,rev_gsequence_alt,
			      rev_roffset,rev_goffset_anchor,pairpool,/*revp*/true,
			      chroffset,chrhigh,cdna_direction,watsonp,*dynprogindex);

  } else {
#endif

    pairs = traceback_local(NULL,&(*nmatches),&(*nmismatches),&(*nopens),&(*nindels),
			    directions_nogap,directions_Egap,directions_Fgap,&bestr,&bestc,/*endc*/contlength,
#ifdef PMAP
			    inst_rev_rsequence,inst_rev_rsequence,
#else
			    rev_rsequence,revsequenceuc1,
#endif
			    rev_gsequence,rev_gsequence_uc,rev_gsequence_alt,
			    rev_roffset,rev_goffset_far,pairpool,/*revp*/true,
			    chroffset,chrhigh,cdna_direction,watsonp,*dynprogindex);

    pairs = Pairpool_push_gapholder(pairs,pairpool,/*queryjump*/0,/*genomejump*/rev_goffset_anchor - rev_goffset_far,
				    /*leftpair*/NULL,/*rightpair*/NULL,/*knownp*/true);

    pairs = traceback_local(pairs,&(*nmatches),&(*nmismatches),&(*nopens),&(*nindels),
			    directions_nogap,directions_Egap,directions_Fgap,&bestr,&bestc,/*endc*/0,
#ifdef PMAP
			    inst_rev_rsequence,inst_rev_rsequence,
#else
			    rev_rsequence,revsequenceuc1,
#endif
			    rev_gsequence,rev_gsequence_uc,rev_gsequence_alt,
			    rev_roffset,rev_goffset_anchor,pairpool,/*revp*/true,
			    chroffset,chrhigh,cdna_direction,watsonp,*dynprogindex);
#if defined(HAVE_SSE4_1) || defined(HAVE_SSE2)
  }
#endif

  /* Score compared with perfect score, so heavy weight on mismatches may not be necessary */
  *finalscore = (*nmatches)*FULLMATCH + (*nmismatches)*MISMATCH_ENDQ + (*nopens)*open + (*nindels)*extend;
  *missscore = (*finalscore) - rlength*FULLMATCH;
  debug6(printf("finalscore %d = %d*%d matches + %d*%d mismatches + %d*%d opens + %d*%d extends\n",
		*finalscore,FULLMATCH,*nmatches,MISMATCH_ENDQ,*nmismatches,open,*nopens,extend,*nindels));
  debug6(printf("missscore = %d\n",*missscore));

  /* Add 1 to count the match already in the alignment */
  pairs = List_reverse(pairs); /* Look at 5' end to remove excess gaps */
  while (pairs != NULL && (pair = List_head(pairs)) && pair->comp == INDEL_COMP) {
    pairs = List_next(pairs);
  }

#ifdef PMAP
  FREE(inst_rsequence);
#endif

  debug6(Pair_dump_list(pairs,true));
  debug6(printf("End of dynprog end5 gap splicejunction\n\n"));

  *dynprogindex += (*dynprogindex > 0 ? +1 : -1);
  return List_reverse(pairs);
}



List_T
Dynprog_end3_gap (int *dynprogindex, int *finalscore, int *nmatches, int *nmismatches, 
		  int *nopens, int *nindels, T dynprog, 
		  char *rsequence, char *sequenceuc1,
		  int rlength, int glength, int roffset, int goffset, 
		  Univcoord_T chroffset, Univcoord_T chrhigh,
#ifdef PMAP
		  char *queryaaseq,
#endif
		  int cdna_direction, bool watsonp, bool jump_late_p, Pairpool_T pairpool,
		  int extraband_end, double defect_rate, Endalign_T endalign) {
  List_T pairs = NULL;
  char *gsequence, *gsequence_alt;
#ifdef PMAP
  char *inst_rsequence;
#endif
  Pair_T pair;
  Mismatchtype_T mismatchtype;
  int bestr, bestc, lband, uband;
#if defined(HAVE_SSE4_1) || defined(HAVE_SSE2)
  Score8_T **matrix8;
  Direction8_T **directions8_nogap, **directions8_Egap, **directions8_Fgap;
  bool use8p = false;

  Score16_T **matrix, open, extend;
  Direction16_T **directions_nogap, **directions_Egap, **directions_Fgap;
#else
  Score32_T **matrix, open, extend;
  Direction32_T **directions_nogap, **directions_Egap, **directions_Fgap;
#endif
#ifdef PMAP
  int termpos, termmod;
#endif

  debug6(
	printf("%c:  ",*dynprogindex > 0 ? (*dynprogindex-1)%26+'a' : (-(*dynprogindex)-1)%26+'A');
	printf("Aligning 3' end gap with endalign = %d\n",endalign);
	);

  mismatchtype = ENDQ;
  if (defect_rate < DEFECT_HIGHQ) {
    open = END_OPEN_HIGHQ;
    extend = END_EXTEND_HIGHQ;
  } else if (defect_rate < DEFECT_MEDQ) {
    open = END_OPEN_MEDQ;
    extend = END_EXTEND_MEDQ;
  } else {
    open = END_OPEN_LOWQ;
    extend = END_EXTEND_LOWQ;
  }

  /* We can just chop lengths to work, since we're not constrained on 3' end */
  if (rlength <= 0) {
    /* Needed to avoid abort by Matrix16_alloc */
    *nmatches = *nmismatches = *nopens = *nindels = 0;
    *finalscore = 0;
    return (List_T) NULL;
  } else if (endalign == QUERYEND_NOGAPS) {
    /* Don't shorten rlength */
  } else if (rlength > dynprog->max_rlength) {
    debug6(printf("rlength %d is too long.  Chopping to %d\n",rlength,dynprog->max_rlength));
    rlength = dynprog->max_rlength;
  }
  if (glength <= 0) {
    /* Needed to avoid abort by Matrix16_alloc */
    *nmatches = *nmismatches = *nopens = *nindels = 0;
    *finalscore = 0;
    return (List_T) NULL;
  } else if (endalign == QUERYEND_NOGAPS) {
    /* Don't shorten glength */
  } else if (glength > dynprog->max_glength) {
    debug6(printf("glength %d is too long.  Chopping to %d\n",glength,dynprog->max_glength));
    glength = dynprog->max_glength;
  }

#ifdef PMAP
  inst_rsequence = instantiate_codons(rsequence,queryaaseq,roffset,rlength);
  debug6(printf("At query offset %d-%d, %.*s\n",roffset,roffset+rlength-1,rlength,inst_rsequence));
#else
  debug6(printf("At query offset %d-%d, %.*s\n",roffset,roffset+rlength-1,rlength,rsequence));
#endif
	
#ifdef EXTRACT_GENOMICSEG
  debug6(printf("At genomic offset %d-%d, %.*s\n",goffset,goffset+glength-1,glength,gsequence));
#endif


  if (watsonp) {
    gsequence = Genome_get_segment_blocks_right(&gsequence_alt,/*left*/chroffset+goffset,glength,chrhigh,/*revcomp*/false);
  } else {
    gsequence = Genome_get_segment_blocks_left(&gsequence_alt,/*left*/chrhigh-goffset+1,glength,chroffset,/*revcomp*/true);
  }
  if (gsequence == NULL) {
    *nmatches = *nmismatches = *nopens = *nindels = 0;
    *finalscore = 0;
    return (List_T) NULL;
  }

  if (endalign == QUERYEND_GAP || endalign == BEST_LOCAL) {
    compute_bands(&lband,&uband,rlength,glength,extraband_end,/*widebandp*/true);
#if defined(HAVE_SSE4_1) || defined(HAVE_SSE2)
    /* Use || because we want the minimum length (which determines the diagonal length) to achieve a score less than 128 */
    if (rlength <= SIMD_MAXLENGTH_EPI8 || glength <= SIMD_MAXLENGTH_EPI8) {
      use8p = true;
      matrix8 = compute_scores_simd_8(&directions8_nogap,&directions8_Egap,&directions8_Fgap,dynprog,
#ifdef PMAP
				      inst_rsequence,
#else
				      rsequence,
#endif
				      gsequence,gsequence_alt,rlength,glength,
#ifdef DEBUG14
				      goffset,chroffset,chrhigh,watsonp,
#endif
				      mismatchtype,(Score8_T) open,(Score8_T) extend,
				      lband,uband,jump_late_p,/*revp*/false);
      find_best_endpoint_8(&(*finalscore),&bestr,&bestc,matrix8,rlength,glength,extraband_end,
			   jump_late_p);
    } else {
      matrix = compute_scores_simd_16(&directions_nogap,&directions_Egap,&directions_Fgap,dynprog,
#ifdef PMAP
				      inst_rsequence,
#else
				      rsequence,
#endif
				      gsequence,gsequence_alt,rlength,glength,
#ifdef DEBUG14
				      goffset,chroffset,chrhigh,watsonp,
#endif
				      mismatchtype,(Score16_T) open, (Score16_T) extend,
				      lband,uband,jump_late_p,/*revp*/false);
      find_best_endpoint(&(*finalscore),&bestr,&bestc,matrix,rlength,glength,extraband_end,
			 jump_late_p);
    }

#else

    matrix = compute_scores_standard(&directions_nogap,&directions_Egap,&directions_Fgap,dynprog,
#ifdef PMAP
				     inst_rsequence,
#else
				     rsequence,
#endif
				     gsequence,gsequence_alt,rlength,glength,
				     goffset,chroffset,chrhigh,watsonp,
				     mismatchtype,open,extend,
				     lband,uband,jump_late_p,/*revp*/false);
    find_best_endpoint(&(*finalscore),&bestr,&bestc,matrix,rlength,glength,extraband_end,
		       jump_late_p);
#endif

  } else if (endalign == QUERYEND_INDELS) {
    compute_bands(&lband,&uband,rlength,glength,extraband_end,/*widebandp*/true);
#if defined(HAVE_SSE4_1) || defined(HAVE_SSE2)
    /* Use || because we want the minimum length (which determines the diagonal length) to achieve a score less than 128 */
    if (rlength <= SIMD_MAXLENGTH_EPI8 || glength <= SIMD_MAXLENGTH_EPI8) {
      use8p = true;
      matrix8 = compute_scores_simd_8(&directions8_nogap,&directions8_Egap,&directions8_Fgap,dynprog,
#ifdef PMAP
				      inst_rsequence,
#else
				      rsequence,
#endif
				      gsequence,gsequence_alt,rlength,glength,
#ifdef DEBUG14
				      goffset,chroffset,chrhigh,watsonp,
#endif
				      mismatchtype,(Score8_T) open,(Score8_T) extend,
				      lband,uband,jump_late_p,/*revp*/false);
      find_best_endpoint_to_queryend_indels_8(&(*finalscore),&bestr,&bestc,matrix8,rlength,glength,extraband_end,
					      jump_late_p);
      /* *finalscore = 0; -- Splicetrie procedures need to know finalscore */

    } else {
      matrix = compute_scores_simd_16(&directions_nogap,&directions_Egap,&directions_Fgap,dynprog,
#ifdef PMAP
				      inst_rsequence,
#else
				      rsequence,
#endif
				      gsequence,gsequence_alt,rlength,glength,
#ifdef DEBUG14
				      goffset,chroffset,chrhigh,watsonp,
#endif
				      mismatchtype,(Score16_T) open, (Score16_T) extend,
				      lband,uband,jump_late_p,/*revp*/false);
      find_best_endpoint_to_queryend_indels(&(*finalscore),&bestr,&bestc,matrix,rlength,glength,extraband_end,
					    jump_late_p);
      /* *finalscore = 0; -- Splicetrie procedures need to know finalscore */
    }

#else

    matrix = compute_scores_standard(&directions_nogap,&directions_Egap,&directions_Fgap,dynprog,
#ifdef PMAP
				     inst_rsequence,
#else
				     rsequence,
#endif
				     gsequence,gsequence_alt,rlength,glength,
				     goffset,chroffset,chrhigh,watsonp,
				     mismatchtype,open,extend,
				     lband,uband,jump_late_p,/*revp*/false);
    find_best_endpoint_to_queryend_indels(&(*finalscore),&bestr,&bestc,matrix,rlength,glength,extraband_end,
					  jump_late_p);
    /* *finalscore = 0; -- Splicetrie procedures need to know finalscore */
#endif

  } else if (endalign == QUERYEND_NOGAPS) {
    find_best_endpoint_to_queryend_nogaps(&bestr,&bestc,rlength,glength);
    /* *finalscore = 0; -- Splicetrie procedures need to know finalscore */

  } else {
    fprintf(stderr,"Unexpected endalign value %d\n",endalign);
    abort();
  }

#ifdef PMAP
  termpos = roffset+(bestc-1);
  debug6(printf("Final query pos is %d\n",termpos));
  if ((termmod = termpos % 3) < 2) {
    if (bestr + (2 - termmod) < rlength && bestc + (2 - termmod) < glength) {
      debug6(printf("Rounding up by %d\n",2 - termmod));
      bestr += 2 - termmod;
      bestc += 2 - termmod;
    }
  }
#endif

  *nmatches = *nmismatches = *nopens = *nindels = 0;
  if (endalign == QUERYEND_NOGAPS) {
    pairs = traceback_nogaps(NULL,&(*nmatches),&(*nmismatches),bestr,bestc,
#ifdef PMAP
			     inst_rsequence,inst_rsequence,
#else
			     rsequence,sequenceuc1,
#endif
			     gsequence,gsequence_alt,roffset,goffset,pairpool,
#ifdef DEBUG14
			     chroffset,chrhigh,watsonp,
#endif
			     /*revp*/false,*dynprogindex);
    *finalscore = (*nmatches)*FULLMATCH + (*nmismatches)*MISMATCH_ENDQ;

#if defined(HAVE_SSE4_1) || defined(HAVE_SSE2)
  } else if (use8p == true) {
    pairs = traceback_8(NULL,&(*nmatches),&(*nmismatches),&(*nopens),&(*nindels),
			directions8_nogap,directions8_Egap,directions8_Fgap,bestr,bestc,
#ifdef PMAP
			inst_rsequence,inst_rsequence,
#else
			rsequence,sequenceuc1,
#endif
			gsequence,gsequence_alt,roffset,goffset,pairpool,/*revp*/false,
			chroffset,chrhigh,cdna_direction,watsonp,*dynprogindex);
#endif

  } else {
    pairs = traceback(NULL,&(*nmatches),&(*nmismatches),&(*nopens),&(*nindels),
		      directions_nogap,directions_Egap,directions_Fgap,bestr,bestc,
#ifdef PMAP
		      inst_rsequence,inst_rsequence,
#else
		      rsequence,sequenceuc1,
#endif
		      gsequence,gsequence_alt,roffset,goffset,pairpool,/*revp*/false,
		      chroffset,chrhigh,cdna_direction,watsonp,*dynprogindex);
  }

  if ((endalign == QUERYEND_GAP || endalign == BEST_LOCAL) && (*nmatches + 1) < *nmismatches) {
    *finalscore = 0;
    /* No need to free pairs */
    pairs = NULL;

  } else {
    /* Add 1 to count the match already in the alignment */
    pairs = List_reverse(pairs); /* Look at 3' end to remove excess gaps */
    while (pairs != NULL && (pair = List_head(pairs)) && pair->comp == INDEL_COMP) {
      pairs = List_next(pairs);
    }
  }

  /*
    Directions_free(directions);
    Matrix_free(matrix);
  */

#ifdef PMAP
  FREE(inst_rsequence);
#endif

  if (gsequence_alt != gsequence) {
    FREE(gsequence_alt);
  }
  FREE(gsequence);

  debug6(printf("End of dynprog end3 gap\n\n"));

  *dynprogindex += (*dynprogindex > 0 ? +1 : -1);
  return pairs;			/* not List_reverse(pairs) */
}


/* gsequence is the splicejunction */
List_T
Dynprog_end3_splicejunction (int *dynprogindex, int *finalscore, int *missscore,
			     int *nmatches, int *nmismatches, int *nopens, int *nindels, T dynprog, 
			     char *rsequence, char *sequenceuc1,
			     char *gsequence, char *gsequence_uc, char *gsequence_alt,
			     int rlength, int glength, int roffset, int goffset_anchor, int goffset_far,
			     Univcoord_T chroffset, Univcoord_T chrhigh,
#ifdef PMAP
			     char *queryaaseq,
#endif
			     int cdna_direction, bool watsonp, bool jump_late_p, Pairpool_T pairpool,
			     int extraband_end, double defect_rate, int contlength) {
  List_T pairs = NULL;
#ifdef PMAP
  char *inst_rsequence;
#endif
  Pair_T pair;
  Mismatchtype_T mismatchtype;
  int bestr, bestc, lband, uband;
#if defined(HAVE_SSE4_1) || defined(HAVE_SSE2)
  Score8_T **matrix8;
  Direction8_T **directions8_nogap, **directions8_Egap, **directions8_Fgap;
  bool use8p = false;

  Score16_T **matrix, open, extend;
  Direction16_T **directions_nogap, **directions_Egap, **directions_Fgap;
#else
  Score32_T **matrix, open, extend;
  Direction32_T **directions_nogap, **directions_Egap, **directions_Fgap;
#endif
#ifdef PMAP
  int termpos, termmod;
#endif

  debug6(
	printf("%c:  ",*dynprogindex > 0 ? (*dynprogindex-1)%26+'a' : (-(*dynprogindex)-1)%26+'A');
	printf("Aligning 3' end gap splicejunction\n");
	);

  mismatchtype = ENDQ;
  if (defect_rate < DEFECT_HIGHQ) {
    open = END_OPEN_HIGHQ;
    extend = END_EXTEND_HIGHQ;
  } else if (defect_rate < DEFECT_MEDQ) {
    open = END_OPEN_MEDQ;
    extend = END_EXTEND_MEDQ;
  } else {
    open = END_OPEN_LOWQ;
    extend = END_EXTEND_LOWQ;
  }


  /* We can just chop lengths to work, since we're not constrained on 3' end */
  if (rlength <= 0 || rlength > dynprog->max_rlength) {
    /* Needed to avoid abort by Matrix16_alloc */
    *nmatches = *nmismatches = *nopens = *nindels = 0;
    *finalscore = 0;
    *missscore = -100;
    return (List_T) NULL;
  }
  if (glength <= 0 || glength > dynprog->max_glength) {
    /* Needed to avoid abort by Matrix16_alloc */
    *nmatches = *nmismatches = *nopens = *nindels = 0;
    *finalscore = 0;
    *missscore = -100;
    return (List_T) NULL;
  }

#ifdef PMAP
  inst_rsequence = instantiate_codons(rsequence,queryaaseq,roffset,rlength);
  debug6(printf("At query offset %d-%d, %.*s\n",roffset,roffset+rlength-1,rlength,inst_rsequence));
#else
  debug6(printf("At query offset %d-%d, %.*s\n",roffset,roffset+rlength-1,rlength,rsequence));
#endif
	
  debug6(printf("At genomic offset %d-%d, %.*s\n",
		goffset_anchor,goffset_anchor+glength-1,glength,gsequence));


  /* find_best_endpoint_to_queryend_nogaps(&bestr,&bestc,rlength,glength); */
  /* bestr = bestc = rlength; */
  /* *finalscore = 0; -- Splicetrie procedures need to know finalscore */

#ifdef PMAP
  termpos = roffset+(bestc-1);
  debug6(printf("Final query pos is %d\n",termpos));
  if ((termmod = termpos % 3) < 2) {
    if (bestr + (2 - termmod) < rlength && bestc + (2 - termmod) < glength) {
      debug6(printf("Rounding up by %d\n",2 - termmod));
      bestr += 2 - termmod;
      bestc += 2 - termmod;
    }
  }
#endif

  compute_bands(&lband,&uband,rlength,glength,extraband_end,/*widebandp*/true);
#if defined(HAVE_SSE4_1) || defined(HAVE_SSE2)
  /* Use || because we want the minimum length (which determines the diagonal length) to achieve a score less than 128 */
  if (rlength <= SIMD_MAXLENGTH_EPI8 || glength <= SIMD_MAXLENGTH_EPI8) {
    use8p = true;
    matrix8 = compute_scores_simd_8(&directions8_nogap,&directions8_Egap,&directions8_Fgap,dynprog,
#ifdef PMAP
				    inst_rsequence,gsequence,gsequence_alt,
#else
				    rsequence,gsequence_uc,gsequence_alt,
#endif
				    rlength,glength,
#ifdef DEBUG14
				    /*goffset*/0,chroffset,chrhigh,watsonp,
#endif
				    mismatchtype,(Score8_T) open,(Score8_T) extend,
				    lband,uband,jump_late_p,/*revp*/false);
    find_best_endpoint_to_queryend_indels_8(&(*finalscore),&bestr,&bestc,matrix8,rlength,glength,extraband_end,
					    jump_late_p);

  } else {
    matrix = compute_scores_simd_16(&directions_nogap,&directions_Egap,&directions_Fgap,dynprog,
#ifdef PMAP
				    inst_rsequence,gsequence,gsequence_alt,
#else
				    rsequence,gsequence_uc,gsequence_alt,
#endif
				    rlength,glength,
#ifdef DEBUG14
				    /*goffset*/0,chroffset,chrhigh,watsonp,
#endif
				    mismatchtype,(Score16_T) open, (Score16_T) extend,
				    lband,uband,jump_late_p,/*revp*/false);
    find_best_endpoint_to_queryend_indels(&(*finalscore),&bestr,&bestc,matrix,rlength,glength,extraband_end,
					  jump_late_p);
  }

#else

  matrix = compute_scores_standard(&directions_nogap,&directions_Egap,&directions_Fgap,dynprog,
#ifdef PMAP
				   inst_rsequence,gsequence,gsequence_alt,
#else
				   rsequence,gsequence,gsequence_alt,
#endif
				   rlength,glength,
				   /*goffset*/0,chroffset,chrhigh,watsonp,
				   mismatchtype,open,extend,
				   lband,uband,jump_late_p,/*revp*/false);
  find_best_endpoint_to_queryend_indels(&(*finalscore),&bestr,&bestc,matrix,rlength,glength,extraband_end,
					jump_late_p);
#endif

  *nmatches = *nmismatches = *nopens = *nindels = 0;
#if defined(HAVE_SSE4_1) || defined(HAVE_SSE2)
  if (use8p == true) {
    pairs = traceback_local_8(NULL,&(*nmatches),&(*nmismatches),&(*nopens),&(*nindels),
			      directions8_nogap,directions8_Egap,directions8_Fgap,&bestr,&bestc,/*endc*/contlength,
#ifdef PMAP
			      inst_rsequence,inst_rsequence,
#else
			      rsequence,sequenceuc1,
#endif
			      gsequence,gsequence_uc,gsequence_alt,
			      roffset,goffset_far,pairpool,/*revp*/false,
			      chroffset,chrhigh,cdna_direction,watsonp,*dynprogindex);

    pairs = Pairpool_push_gapholder(pairs,pairpool,/*queryjump*/0,/*genomejump*/goffset_far - goffset_anchor,
				    /*leftpair*/NULL,/*rightpair*/NULL,/*knownp*/true);

    pairs = traceback_local_8(pairs,&(*nmatches),&(*nmismatches),&(*nopens),&(*nindels),
			      directions8_nogap,directions8_Egap,directions8_Fgap,&bestr,&bestc,/*endc*/0,
#ifdef PMAP
			      inst_rsequence,inst_rsequence,
#else
			      rsequence,sequenceuc1,
#endif
			      gsequence,gsequence_uc,gsequence_alt,
			      roffset,goffset_anchor,pairpool,/*revp*/false,
			      chroffset,chrhigh,cdna_direction,watsonp,*dynprogindex);

  } else {
#endif

    pairs = traceback_local(NULL,&(*nmatches),&(*nmismatches),&(*nopens),&(*nindels),
			    directions_nogap,directions_Egap,directions_Fgap,&bestr,&bestc,/*endc*/contlength,
#ifdef PMAP
			    inst_rsequence,inst_rsequence,
#else
			    rsequence,sequenceuc1,
#endif
			    gsequence,gsequence_uc,gsequence_alt,
			    roffset,goffset_far,pairpool,/*revp*/false,
			    chroffset,chrhigh,cdna_direction,watsonp,*dynprogindex);

    pairs = Pairpool_push_gapholder(pairs,pairpool,/*queryjump*/0,/*genomejump*/goffset_far - goffset_anchor,
				    /*leftpair*/NULL,/*rightpair*/NULL,/*knownp*/true);

    pairs = traceback_local(pairs,&(*nmatches),&(*nmismatches),&(*nopens),&(*nindels),
			    directions_nogap,directions_Egap,directions_Fgap,&bestr,&bestc,/*endc*/0,
#ifdef PMAP
			    inst_rsequence,inst_rsequence,
#else
			    rsequence,sequenceuc1,
#endif
			    gsequence,gsequence_uc,gsequence_alt,
			    roffset,goffset_anchor,pairpool,/*revp*/false,
			    chroffset,chrhigh,cdna_direction,watsonp,*dynprogindex);
#if defined(HAVE_SSE4_1) || defined(HAVE_SSE2)
  }
#endif

  /* Score compared with perfect score, so heavy weight on mismatches may not be necessary */
  *finalscore = (*nmatches)*FULLMATCH + (*nmismatches)*MISMATCH_ENDQ + (*nopens)*open + (*nindels)*extend;
  *missscore = (*finalscore) - rlength*FULLMATCH;
  debug6(printf("finalscore %d = %d*%d matches + %d*%d mismatches + %d*%d opens + %d*%d extends\n",
		*finalscore,FULLMATCH,*nmatches,MISMATCH_ENDQ,*nmismatches,open,*nopens,extend,*nindels));
  debug6(printf("missscore = %d\n",*missscore));

  /* Add 1 to count the match already in the alignment */
  pairs = List_reverse(pairs); /* Look at 3' end to remove excess gaps */
  while (pairs != NULL && (pair = List_head(pairs)) && pair->comp == INDEL_COMP) {
    pairs = List_next(pairs);
  }

#ifdef PMAP
  FREE(inst_rsequence);
#endif

  debug6(Pair_dump_list(pairs,true));
  debug6(printf("End of dynprog end3 gap splicejunction\n\n"));

  *dynprogindex += (*dynprogindex > 0 ? +1 : -1);
  return pairs;			/* not List_reverse(pairs) */
}


static void
make_contjunction_5 (char *splicejunction, char *splicejunction_alt, Univcoord_T splicecoord,
		     int splicelength, int contlength, Splicetype_T anchor_splicetype,
		     bool watsonp) {
  char *proximal, *proximal_alt;

  debug7(printf("make_contjunction_5 at %u, splice (%s), contlength %d, splicelength %d:",
		splicecoord,Splicetype_string(anchor_splicetype),contlength, splicelength));

  proximal = &(splicejunction[splicelength]);
  proximal_alt = &(splicejunction_alt[splicelength]);

  if (anchor_splicetype == ACCEPTOR) {
    Genome_fill_buffer_blocks_noterm(splicecoord,contlength,proximal,proximal_alt);

  } else if (anchor_splicetype == ANTIDONOR) {
    Genome_fill_buffer_blocks_noterm(splicecoord,contlength,proximal,proximal_alt);

  } else if (anchor_splicetype == ANTIACCEPTOR) {
    Genome_fill_buffer_blocks_noterm(splicecoord-contlength,contlength,proximal,proximal_alt);

  } else if (anchor_splicetype == DONOR) {
    Genome_fill_buffer_blocks_noterm(splicecoord-contlength,contlength,proximal,proximal_alt);
    
  } else {
    fprintf(stderr,"Unexpected anchor_splicetype value %d\n",anchor_splicetype);
    abort();
  }

  if (watsonp == false) {
    make_complement_inplace(proximal,contlength);
    make_complement_inplace(proximal_alt,contlength);
  }

#ifdef DEBUG7
  if (watsonp == true) {
    printf(" (fwd)  contjunction    : %.*s\n",contlength,proximal);
    printf(" (fwd)  contjunction_alt: %.*s\n",contlength,proximal_alt);
  } else {
    printf(" (rev)  contjunction    : %.*s\n",contlength,proximal);
    printf(" (rev)  contjunction_alt: %.*s\n",contlength,proximal_alt);
  }
#endif

  return;
}



/* Fills in just the distal part, keeping the proximal part same for contlength */
bool
Dynprog_make_splicejunction_5 (char *splicejunction, char *splicejunction_alt, Univcoord_T splicecoord,
			       int splicelength, int contlength, Splicetype_T far_splicetype,
			       bool watsonp) {
  char *distal, *distal_alt;

  debug7(printf("make_splicejunction_5 at %u, splice (%s), contlength %d, splicelength %d:\n",
		splicecoord,Splicetype_string(far_splicetype),contlength, splicelength));

  distal = &(splicejunction[0]);
  distal_alt = &(splicejunction_alt[0]);

  if (far_splicetype == ACCEPTOR) {
    Genome_fill_buffer_blocks_noterm(splicecoord,splicelength,distal,distal_alt);

  } else if (far_splicetype == ANTIDONOR) {
    Genome_fill_buffer_blocks_noterm(splicecoord,splicelength,distal,distal_alt);

  } else if (splicecoord <= splicelength) {
    return false;

  } else if (far_splicetype == ANTIACCEPTOR) {
    Genome_fill_buffer_blocks_noterm(splicecoord-splicelength,splicelength,distal,distal_alt);

  } else if (far_splicetype == DONOR) {
    Genome_fill_buffer_blocks_noterm(splicecoord-splicelength,splicelength,distal,distal_alt);
    
  } else {
    fprintf(stderr,"Unexpected far_splicetype value %d\n",far_splicetype);
    abort();
  }

  if (watsonp == false) {
    make_complement_inplace(distal,splicelength);
    make_complement_inplace(distal_alt,splicelength);
  }

#ifdef DEBUG7
  if (watsonp == true) {
    printf(" (fwd)  splicejunction    : %s\n",splicejunction);
    printf(" (fwd)  splicejunction_alt: %s\n",splicejunction_alt);
  } else {
    printf(" (rev)  splicejunction    : %s\n",splicejunction);
    printf(" (rev)  splicejunction_alt: %s\n",splicejunction_alt);
  }
#endif

  return true;
}


static void
make_contjunction_3 (char *splicejunction, char *splicejunction_alt, Univcoord_T splicecoord,
		     int splicelength, int contlength, Splicetype_T anchor_splicetype,
		     bool watsonp) {
  char *proximal, *proximal_alt;

  debug7(printf("make_contjunction_3 at %u, splice (%s), contlength %d, splicelength %d:\n",
		splicecoord,Splicetype_string(anchor_splicetype),contlength,splicelength));

  proximal = &(splicejunction[0]);
  proximal_alt = &(splicejunction_alt[0]);

  if (anchor_splicetype == DONOR) {
    Genome_fill_buffer_blocks_noterm(splicecoord-contlength,contlength,proximal,proximal_alt);

  } else if (anchor_splicetype == ANTIACCEPTOR) {
    Genome_fill_buffer_blocks_noterm(splicecoord-contlength,contlength,proximal,proximal_alt);

  } else if (anchor_splicetype == ANTIDONOR) {
    Genome_fill_buffer_blocks_noterm(splicecoord,contlength,proximal,proximal_alt);

  } else if (anchor_splicetype == ACCEPTOR) {
    Genome_fill_buffer_blocks_noterm(splicecoord,contlength,proximal,proximal_alt);
    
  } else {
    fprintf(stderr,"Unexpected anchor_splicetype value %d\n",anchor_splicetype);
    abort();
  }

  if (watsonp == false) {
    make_complement_inplace(proximal,contlength);
    make_complement_inplace(proximal_alt,contlength);
  }

#ifdef DEBUG7
  if (watsonp == true) {
    printf(" (fwd)  contjunction    : %.*s\n",contlength,proximal);
    printf(" (fwd)  contjunction_alt: %.*s\n",contlength,proximal_alt);
  } else {
    printf(" (rev)  contjunction    : %.*s\n",contlength,proximal);
    printf(" (rev)  contjunction_alt: %.*s\n",contlength,proximal_alt);
  }
#endif

  return;
}



/* Fills in just the distal part, keeping the proximal part same for contlength */
bool
Dynprog_make_splicejunction_3 (char *splicejunction, char *splicejunction_alt, Univcoord_T splicecoord,
			       int splicelength, int contlength, Splicetype_T far_splicetype,
			       bool watsonp) {
  char *distal, *distal_alt;

  debug7(printf("make_splicejunction_3 at %u, splice (%s), contlength %d, splicelength %d:\n",
		splicecoord,Splicetype_string(far_splicetype),contlength,splicelength));

  distal = &(splicejunction[contlength]);
  distal_alt = &(splicejunction_alt[contlength]);

  if (far_splicetype == ANTIDONOR) {
    Genome_fill_buffer_blocks_noterm(splicecoord,splicelength,distal,distal_alt);

  } else if (far_splicetype == ACCEPTOR) {
    Genome_fill_buffer_blocks_noterm(splicecoord,splicelength,distal,distal_alt);
    
  } else if (splicecoord <= splicelength) {
    return false;

  } else if (far_splicetype == DONOR) {
    Genome_fill_buffer_blocks_noterm(splicecoord-splicelength,splicelength,distal,distal_alt);

  } else if (far_splicetype == ANTIACCEPTOR) {
    Genome_fill_buffer_blocks_noterm(splicecoord-splicelength,splicelength,distal,distal_alt);


  } else {
    fprintf(stderr,"Unexpected far_splicetype value %d\n",far_splicetype);
    abort();
  }

  if (watsonp == false) {
    make_complement_inplace(distal,splicelength);
    make_complement_inplace(distal_alt,splicelength);
  }

#ifdef DEBUG7
  if (watsonp == true) {
    printf(" (fwd)  splicejunction    : %s\n",splicejunction);
    printf(" (fwd)  splicejunction_alt: %s\n",splicejunction_alt);
  } else {
    printf(" (rev)  splicejunction    : %s\n",splicejunction);
    printf(" (rev)  splicejunction_alt: %s\n",splicejunction_alt);
  }
#endif

  return true;
}


List_T
Dynprog_end5_known (bool *knownsplicep, int *dynprogindex, int *finalscore,
		    int *ambig_end_length, Splicetype_T *ambig_splicetype,
		    int *nmatches, int *nmismatches, int *nopens, int *nindels, T dynprog, 
		    char *rev_rsequence, char *revsequenceuc1,
		    int rlength, int glength, int rev_roffset, int rev_goffset, 
		    Univcoord_T chroffset, Univcoord_T chrhigh,
		    Univcoord_T knownsplice_limit_low, Univcoord_T knownsplice_limit_high,
#ifdef PMAP
		    char *queryaaseq,
#endif
		    int cdna_direction, bool watsonp, bool jump_late_p, Pairpool_T pairpool,
		    int extraband_end, double defect_rate) {
  List_T best_pairs = NULL, orig_pairs;
  Pair_T pair;
  Univcoord_T low, high, far_limit_low, far_limit_high;
  Splicetype_T anchor_splicetype, far_splicetype;
  int contlength, splicelength, endlength;
  char *splicejunction, *splicejunction_alt;
#ifdef EXTRACT_GENOMICSEG
  char *splicejunction_test;
#endif
  int jstart, j;

  int orig_score, threshold_miss_score, perfect_score;
  int obsmax_penalty;


  assert(glength >= rlength);

  debug7(
	printf("%c:  ",*dynprogindex > 0 ? (*dynprogindex-1)%26+'a' : (-(*dynprogindex)-1)%26+'A');
	printf("Aligning 5' end gap, known\n")
	);

  *ambig_end_length = 0;

  /* We can just chop lengths to work, since we're not constrained on 5' end */
  if (rlength <= 0) {
    /* Needed to avoid abort by Matrix16_alloc */
    *finalscore = 0;
    *knownsplicep = false;
    return (List_T) NULL;
  }
  if (glength <= 0) {
    /* Needed to avoid abort by Matrix16_alloc */
    *finalscore = 0;
    *knownsplicep = false;
    return (List_T) NULL;
  }

  debug7(printf("At query offset %d-%d, %.*s\n",rev_roffset-rlength+1,rev_roffset,rlength,&(rev_rsequence[-rlength+1])));
#ifdef EXTRACT_GENOMICSEG
  debug7(printf("At genomic offset %d-%d, %.*s\n",rev_goffset-glength+1,rev_goffset,glength,&(rev_gsequence[-glength+1])));
#endif

  perfect_score = rlength*FULLMATCH;

  /* Try without splicing, all the way to query end */
  best_pairs = Dynprog_end5_gap(&(*dynprogindex),&(*finalscore),&(*nmatches),&(*nmismatches),
				&(*nopens),&(*nindels),dynprog,rev_rsequence,revsequenceuc1,
				rlength,glength,rev_roffset,rev_goffset,chroffset,chrhigh,
#ifdef PMAP
				queryaaseq,
#endif
				cdna_direction,watsonp,jump_late_p,pairpool,
				extraband_end,defect_rate,/*endalign*/QUERYEND_NOGAPS);
  if (*finalscore < 0) {
    orig_score = 0;
    orig_pairs = best_pairs = (List_T) NULL;
  } else {
    orig_score = *finalscore;
    orig_pairs = best_pairs;
  }
  threshold_miss_score = orig_score - perfect_score;
  debug7(printf("score %d - perfect score %d = threshold %d",
		orig_score,perfect_score,threshold_miss_score));
  if (threshold_miss_score < -2*FULLMATCH) {
    /* Don't allow more than 2 mismatches in a distant splice */
    threshold_miss_score = -2*FULLMATCH;
    debug7(printf(", but revising to %d\n",threshold_miss_score));
  }
  debug7(printf("\n"));
  *knownsplicep = false;


  if (threshold_miss_score < 0 && glength > 0) {
    /* Try known splicing */
    splicejunction = (char *) CALLOC(glength+1,sizeof(char));
    splicejunction_alt = (char *) CALLOC(glength+1,sizeof(char));
#ifdef EXTRACT_GENOMICSEG
    splicejunction_test = (char *) CALLOC(glength+1,sizeof(char));
#endif

    endlength = rlength;
    if (watsonp == true) {
      low = chroffset + rev_goffset-endlength + 2;
      high = chroffset + rev_goffset + 1;
      debug7(printf("5' watson\n"));
      debug7(printf("Calculating low %u (%u) = %u + %d-%d + 2\n",
		    low,low-chroffset,chroffset,rev_goffset,endlength));
      debug7(printf("Calculating high %u (%u) = %u + %d + 1\n",
		    high,high-chroffset,chroffset,rev_goffset));
      if (cdna_direction > 0) {
	anchor_splicetype = ACCEPTOR;
	far_splicetype = DONOR;
      } else {
	anchor_splicetype = ANTIDONOR;
	far_splicetype = ANTIACCEPTOR;
      }
    } else {
      low = chrhigh - rev_goffset;
      high = chrhigh - (rev_goffset-endlength) - 1;
      debug7(printf("5' crick\n"));
      debug7(printf("Calculating low %u (%u) = %u - %d\n",
		    low,low-chroffset,chrhigh,rev_goffset));
      debug7(printf("Calculating high %u (%u) = %u - (%d-%d) - 1\n",
		    high,high-chroffset,chrhigh,rev_goffset,endlength));
      if (cdna_direction > 0) {
	anchor_splicetype = ANTIACCEPTOR;
	far_splicetype = ANTIDONOR;
      } else {
	anchor_splicetype = DONOR;
	far_splicetype = ACCEPTOR;
      }
    }

    far_limit_low = knownsplice_limit_low;
    far_limit_high = knownsplice_limit_high;
    debug7(printf("Genomic positions: %u..%u (%u..%u), looking for anchor splicetype %s\n",
		  low,high,low-chroffset,high-chroffset,Splicetype_string(anchor_splicetype)));
    j = jstart = binary_search(0,nsplicesites,splicesites,low);
    while (j < nsplicesites && splicesites[j] <= high) {
      if (splicetypes[j] == anchor_splicetype) {
	debug7(printf("Found one at %u (%u)\n",splicesites[j],splicesites[j]-chroffset));
	if (watsonp == true) {
	  contlength = high - splicesites[j];
	} else {
	  contlength = splicesites[j] - low;
	}
	debug7(printf("contlength %d, splicelength %d, rlength %d, glength %d\n",
		      contlength,glength-contlength,rlength,glength));
	assert(contlength >= 0 && contlength < rlength);

#ifdef EXTRACT_GENOMICSEG
	debug7(printf("cont: %.*s\n",contlength,&(rev_gsequence[-contlength+1])));
#endif
	splicelength = glength - contlength;
	assert(splicelength > 0);
	debug7(printf("  Saw %u (%u) of type %s (cont length %d, splice length %d)\n",
		      splicesites[j],splicesites[j]-chroffset,Splicetype_string(splicetypes[j]),contlength,splicelength));

	make_contjunction_5(splicejunction,splicejunction_alt,splicesites[j],splicelength,contlength,anchor_splicetype,watsonp);
#ifdef EXTRACT_GENOMICSEG
	strncpy(&(splicejunction_test[splicelength]),&(rev_gsequence[-contlength+1]),contlength);
	debug7(printf("contjunction_gen:  %s\n",&(splicejunction[splicelength])));
	debug7(printf("contjunction_test: %s\n",&(splicejunction_test[splicelength])));
	assert(!strncmp(&(splicejunction[splicelength]),&(splicejunction_test[splicelength]),contlength));
#endif

	if (watsonp) {
	  far_limit_high = splicesites[j];
	} else {
	  far_limit_low = splicesites[j];
	}

	obsmax_penalty = 0;
	if (trieoffsets_obs != NULL) {
	  debug7(printf("  Running Splicetrie_solve_end5 on observed splice sites with rev_goffset %d\n",rev_goffset));
	  best_pairs = Splicetrie_solve_end5(best_pairs,triecontents_obs,trieoffsets_obs,j,
					     far_limit_low,far_limit_high,
					     &(*finalscore),&(*nmatches),&(*nmismatches),
					     &(*nopens),&(*nindels),&(*knownsplicep),&(*ambig_end_length),
					     &threshold_miss_score,/*obsmax_penalty*/0,perfect_score,
					     /*anchor_splicesite*/splicesites[j],splicejunction,splicejunction_alt,
					     splicelength,contlength,far_splicetype,
					     chroffset,chrhigh,&(*dynprogindex),dynprog,
					     rev_rsequence,revsequenceuc1,rlength,glength,rev_roffset,rev_goffset,
#ifdef PMAP
					     queryaaseq,
#endif
					     cdna_direction,watsonp,jump_late_p,pairpool,extraband_end,defect_rate);
	  debug7(printf("  Result on obs with ambig_end_length_5 %d\n",*ambig_end_length));
	  debug7(Pair_dump_list(best_pairs,/*zerobasedp*/true));
	  obsmax_penalty += FULLMATCH;
	}

	if (threshold_miss_score + obsmax_penalty < 0 && trieoffsets_max != NULL) {
	  debug7(printf("  Running Splicetrie_solve_end5 on maxdistance splice sites with rev_goffset %d\n",rev_goffset));
	  best_pairs = Splicetrie_solve_end5(best_pairs,triecontents_max,trieoffsets_max,j,
					     far_limit_low,far_limit_high,
					     &(*finalscore),&(*nmatches),&(*nmismatches),
					     &(*nopens),&(*nindels),&(*knownsplicep),&(*ambig_end_length),
					     &threshold_miss_score,obsmax_penalty,perfect_score,
					     /*anchor_splicesite*/splicesites[j],splicejunction,splicejunction_alt,
					     splicelength,contlength,far_splicetype,
					     chroffset,chrhigh,&(*dynprogindex),dynprog,
					     rev_rsequence,revsequenceuc1,rlength,glength,rev_roffset,rev_goffset,
#ifdef PMAP
					     queryaaseq,
#endif
					     cdna_direction,watsonp,jump_late_p,pairpool,extraband_end,defect_rate);
	  debug7(printf("  Result on max with ambig_end_length_5 %d\n",*ambig_end_length));
	  debug7(Pair_dump_list(best_pairs,/*zerobasedp*/true));
	}
      }
      j++;
    }

#ifdef EXTRACT_GENOMICSEG
    FREE(splicejunction_test);
#endif
    FREE(splicejunction_alt);
    FREE(splicejunction);
  }


  if (best_pairs == NULL) {
    if (*ambig_end_length == 0) {
      /* Don't go to query end this time */
      if (rlength > dynprog->max_rlength) {
	debug7(printf("rlength %d is too long.  Chopping to %d\n",rlength,dynprog->max_rlength));
	rlength = dynprog->max_rlength;
      }
      if (glength > dynprog->max_glength) {
	debug7(printf("glength %d is too long.  Chopping to %d\n",glength,dynprog->max_glength));
	glength = dynprog->max_glength;
      }
      orig_pairs = Dynprog_end5_gap(&(*dynprogindex),&(*finalscore),&(*nmatches),&(*nmismatches),
				    &(*nopens),&(*nindels),dynprog,rev_rsequence,revsequenceuc1,
				    rlength,glength,rev_roffset,rev_goffset,chroffset,chrhigh,
#ifdef PMAP
				    queryaaseq,
#endif
				    cdna_direction,watsonp,jump_late_p,pairpool,
				    extraband_end,defect_rate,/*endalign*/BEST_LOCAL);
      debug7(Pair_dump_list(orig_pairs,/*zerobasedp*/true));
      debug7(printf("End of dynprog end5 known\n"));
      *knownsplicep = false;
      return orig_pairs;

    } else {
      *ambig_splicetype = anchor_splicetype;
      debug7(printf("Final result: best_pairs is NULL.  ambig_end_length is %d.  ambig_splicetype is %s.  Result after truncate:\n",
		    *ambig_end_length,Splicetype_string(*ambig_splicetype)));
      /* Truncate ambiguous part.  querypos is decreasing. */
      orig_pairs = List_reverse(orig_pairs);
      while (orig_pairs != NULL && ((Pair_T) orig_pairs->first)->querypos < *ambig_end_length) {
	orig_pairs = Pairpool_pop(orig_pairs,&pair);
      }
      orig_pairs = List_reverse(orig_pairs);
      debug7(Pair_dump_list(orig_pairs,/*zerobasedp*/true));
      *knownsplicep = false;
      *finalscore = orig_score;
      debug7(printf("End of dynprog end5 known\n"));
      return orig_pairs;
    }

  } else {
    debug7(printf("Found a best splice\n"));
    *ambig_end_length = 0;
    debug7(printf("End of dynprog end5 known\n"));
    if (*knownsplicep == true) {
      return Pair_protect_end5(best_pairs,pairpool);
    } else {
      return best_pairs;
    }
  }
}


List_T
Dynprog_end3_known (bool *knownsplicep, int *dynprogindex, int *finalscore,
		    int *ambig_end_length, Splicetype_T *ambig_splicetype,
		    int *nmatches, int *nmismatches, int *nopens, int *nindels, T dynprog, 
		    char *rsequence, char *sequenceuc1,
		    int rlength, int glength, int roffset, int goffset, int querylength,
		    Univcoord_T chroffset, Univcoord_T chrhigh,
		    Univcoord_T knownsplice_limit_low, Univcoord_T knownsplice_limit_high,
#ifdef PMAP
		    char *queryaaseq,
#endif
		    int cdna_direction, bool watsonp, bool jump_late_p, Pairpool_T pairpool,
		    int extraband_end, double defect_rate) {
  List_T best_pairs = NULL, orig_pairs;
  Pair_T pair;
  Univcoord_T low, high, far_limit_low, far_limit_high;
  Splicetype_T anchor_splicetype, far_splicetype;
  int contlength, splicelength, endlength;
  char *splicejunction, *splicejunction_alt;
#ifdef EXTRACT_GENOMICSEG
  char *splicejunction_test;
#endif
  int jstart, j;

  int orig_score, threshold_miss_score, perfect_score;
  int obsmax_penalty;


  assert(glength >= rlength);

  debug7(
	printf("%c:  ",*dynprogindex > 0 ? (*dynprogindex-1)%26+'a' : (-(*dynprogindex)-1)%26+'A');
	printf("Aligning 3' end gap, known\n")
	);

  *ambig_end_length = 0;

  /* We can just chop lengths to work, since we're not constrained on 3' end */
  if (rlength <= 0) {
    /* Needed to avoid abort by Matrix16_alloc */
    *finalscore = 0;
    *knownsplicep = false;
    return (List_T) NULL;
  }
  if (glength <= 0) {
    /* Needed to avoid abort by Matrix16_alloc */
    *finalscore = 0;
    *knownsplicep = false;
    return (List_T) NULL;
  }


  debug7(printf("At query offset %d-%d, %.*s\n",roffset,roffset+rlength-1,rlength,rsequence));
#ifdef EXTRACT_GENOMICSEG
  debug7(printf("At genomic offset %d-%d, %.*s\n",goffset,goffset+glength-1,glength,gsequence));
#endif

  perfect_score = rlength*FULLMATCH;

  /* Try without splicing, all the way to query end */
  best_pairs = Dynprog_end3_gap(&(*dynprogindex),&(*finalscore),&(*nmatches),&(*nmismatches),
				&(*nopens),&(*nindels),dynprog,rsequence,sequenceuc1,
				rlength,glength,roffset,goffset,chroffset,chrhigh,
#ifdef PMAP
				queryaaseq,
#endif
				cdna_direction,watsonp,jump_late_p,pairpool,
				extraband_end,defect_rate,/*endalign*/QUERYEND_NOGAPS);
  if (*finalscore < 0) {
    orig_score = 0;
    orig_pairs = best_pairs = (List_T) NULL;
  } else {
    orig_score = *finalscore;
    orig_pairs = best_pairs;
  }
  threshold_miss_score = orig_score - perfect_score;
  debug7(printf("score %d - perfect score %d = threshold %d",
		orig_score,perfect_score,threshold_miss_score));
  if (threshold_miss_score < -2*FULLMATCH) {
    /* Don't allow more than 2 mismatches in a distant splice */
    threshold_miss_score = -2*FULLMATCH;
    debug7(printf(", but revising to %d\n",threshold_miss_score));
  }
  debug7(printf("\n"));
  *knownsplicep = false;


  if (threshold_miss_score < 0 && glength > 0) {
    /* Try known splicing */
    splicejunction = (char *) CALLOC(glength+1,sizeof(char));
    splicejunction_alt = (char *) CALLOC(glength+1,sizeof(char));
#ifdef EXTRACT_GENOMICSEG
    splicejunction_test = (char *) CALLOC(glength+1,sizeof(char));
#endif

    endlength = rlength;
    if (watsonp == true) {
      low = chroffset + goffset;
      high = chroffset + goffset+endlength - 1;
      debug7(printf("3' watson\n"));
      debug7(printf("Calculating low %u (%u) = %u + %d\n",
		    low,low-chroffset,chroffset,goffset));
      debug7(printf("Calculating high %u (%u) = %u + %d+%d - 1\n",
		    high,high-chroffset,chroffset,goffset,endlength));
      if (cdna_direction > 0) {
	anchor_splicetype = DONOR;
	far_splicetype = ACCEPTOR;
      } else {
	anchor_splicetype = ANTIACCEPTOR;
	far_splicetype = ANTIDONOR;
      }
    } else {
      low = chrhigh - (goffset+endlength) + 2;
      high = chrhigh - goffset + 1;
      debug7(printf("3' crick\n"));
      debug7(printf("Calculating low %u (%u) = %u - (%d+%d)\n",
		    low,low-chroffset,chrhigh,goffset,endlength));
      debug7(printf("Calculating high %u (%u) = %u - %d + 1\n",
		    high,high-chroffset,chrhigh,goffset));
      if (cdna_direction > 0) {
	anchor_splicetype = ANTIDONOR;
	far_splicetype = ANTIACCEPTOR;
      } else {
	anchor_splicetype = ACCEPTOR;
	far_splicetype = DONOR;
      }
    }

    far_limit_low = knownsplice_limit_low;
    far_limit_high = knownsplice_limit_high;
    debug7(printf("Genomic positions: %u..%u (%u..%u), looking for anchor splicetype %s\n",
		  low,high,low-chroffset,high-chroffset,Splicetype_string(anchor_splicetype)));
    j = jstart = binary_search(0,nsplicesites,splicesites,low);
    while (j < nsplicesites && splicesites[j] <= high) {
      if (splicetypes[j] == anchor_splicetype) {
	debug7(printf("Found one at %u (%u)\n",splicesites[j],splicesites[j]-chroffset));
	if (watsonp == true) {
	  contlength = splicesites[j] - low;
	} else {
	  contlength = high - splicesites[j];
	}
	debug7(printf("contlength %d, splicelength %d, rlength %d, glength %d\n",
		      contlength,glength-contlength,rlength,glength));
	assert(contlength >= 0 && contlength < rlength);

#ifdef EXTRACT_GENOMICSEG
	debug7(printf("cont: %.*s\n",contlength,gsequence));
#endif
	splicelength = glength - contlength;
	assert(splicelength > 0);
	debug7(printf("  Saw %u (%u) of type %s (cont length %d, splice length %d)\n",
		      splicesites[j],splicesites[j]-chroffset,Splicetype_string(splicetypes[j]),contlength,splicelength));

	make_contjunction_3(splicejunction,splicejunction_alt,splicesites[j],splicelength,contlength,anchor_splicetype,watsonp);
#ifdef EXTRACT_GENOMICSEG
	strncpy(splicejunction_test,gsequence,contlength);
	debug7(printf("contjunction_gen:  %s\n",splicejunction));
	debug7(printf("contjunction_test: %s\n",splicejunction_test));
	assert(!strncmp(splicejunction,splicejunction_test,contlength));
#endif

	if (watsonp) {
	  far_limit_low = splicesites[j];
	} else {
	  far_limit_high = splicesites[j];
	}

	obsmax_penalty = 0;
	if (trieoffsets_obs != NULL) {
	  debug7(printf("  Running Splicetrie_solve_end3 on observed splice sites with goffset %d\n",goffset));
	  best_pairs = Splicetrie_solve_end3(best_pairs,triecontents_obs,trieoffsets_obs,j,
					     far_limit_low,far_limit_high,
					     &(*finalscore),&(*nmatches),&(*nmismatches),
					     &(*nopens),&(*nindels),&(*knownsplicep),&(*ambig_end_length),
					     &threshold_miss_score,/*obsmax_penalty*/0,perfect_score,
					     /*anchor_splicesite*/splicesites[j],splicejunction,splicejunction_alt,
					     splicelength,contlength,far_splicetype,
					     chroffset,chrhigh,&(*dynprogindex),dynprog,
					     rsequence,sequenceuc1,rlength,glength,roffset,goffset,
#ifdef PMAP
					     queryaaseq,
#endif
					     cdna_direction,watsonp,jump_late_p,pairpool,extraband_end,defect_rate);
	  debug7(printf("  Result on obs with ambig_end_length_3 %d\n",*ambig_end_length));
	  debug7(Pair_dump_list(best_pairs,/*zerobasedp*/true));
	  obsmax_penalty += FULLMATCH;
	}

	if (threshold_miss_score + obsmax_penalty < 0 && trieoffsets_max != NULL) {
	  debug7(printf("  Running Splicetrie_solve_end3 on maxdistance splice sites with goffset %d\n",goffset));
	  best_pairs = Splicetrie_solve_end3(best_pairs,triecontents_max,trieoffsets_max,j,
					     far_limit_low,far_limit_high,
					     &(*finalscore),&(*nmatches),&(*nmismatches),
					     &(*nopens),&(*nindels),&(*knownsplicep),&(*ambig_end_length),
					     &threshold_miss_score,obsmax_penalty,perfect_score,
					     /*anchor_splicesite*/splicesites[j],splicejunction,splicejunction_alt,
					     splicelength,contlength,far_splicetype,
					     chroffset,chrhigh,&(*dynprogindex),dynprog,
					     rsequence,sequenceuc1,rlength,glength,roffset,goffset,
#ifdef PMAP
					     queryaaseq,
#endif
					     cdna_direction,watsonp,jump_late_p,pairpool,extraband_end,defect_rate);
	  debug7(printf("  Result on max with ambig_end_length_3 %d\n",*ambig_end_length));
	  debug7(Pair_dump_list(best_pairs,/*zerobasedp*/true));
	}
      }
      j++;
    }

#ifdef EXTRACT_GENOMICSEG
    FREE(splicejunction_test);
#endif
    FREE(splicejunction_alt);
    FREE(splicejunction);
  }


  if (best_pairs == NULL) {
    if (*ambig_end_length == 0) {
      /* Don't go to query end this time */
      if (rlength > dynprog->max_rlength) {
	debug7(printf("rlength %d is too long.  Chopping to %d\n",rlength,dynprog->max_rlength));
	rlength = dynprog->max_rlength;
      }
      if (glength > dynprog->max_glength) {
	debug7(printf("glength %d is too long.  Chopping to %d\n",glength,dynprog->max_glength));
	glength = dynprog->max_glength;
      }
      orig_pairs = Dynprog_end3_gap(&(*dynprogindex),&(*finalscore),&(*nmatches),&(*nmismatches),
				    &(*nopens),&(*nindels),dynprog,rsequence,sequenceuc1,
				    rlength,glength,roffset,goffset,chroffset,chrhigh,
#ifdef PMAP
				    queryaaseq,
#endif
				    cdna_direction,watsonp,jump_late_p,pairpool,
				    extraband_end,defect_rate,/*endalign*/BEST_LOCAL);
      debug7(Pair_dump_list(orig_pairs,/*zerobasedp*/true));
      *knownsplicep = false;
      debug7(printf("End of dynprog end5 known\n"));
      return orig_pairs;

    } else {
      *ambig_splicetype = anchor_splicetype;
      debug7(printf("Final result: best_pairs is NULL.  ambig_end_length is %d.  ambig_splicetype is %s.  Result after truncate:\n",
		    *ambig_end_length,Splicetype_string(*ambig_splicetype)));
      /* Truncate ambiguous part.  querypos is decreasing */
      while (orig_pairs != NULL && ((Pair_T) orig_pairs->first)->querypos >= querylength - *ambig_end_length) {
	orig_pairs = Pairpool_pop(orig_pairs,&pair);
      }
      debug7(Pair_dump_list(orig_pairs,/*zerobasedp*/true));
      *knownsplicep = false;
      *finalscore = orig_score;
      debug7(printf("End of dynprog end5 known\n"));
      return orig_pairs;
    }

  } else {
    debug7(printf("Found a best splice\n"));
    *ambig_end_length = 0;
    debug7(printf("End of dynprog end3 known\n"));
    if (*knownsplicep == true) {
      return Pair_protect_end3(best_pairs,pairpool);
    } else {
      return best_pairs;
    }
  }
}




static const Except_T microexon_error = {"Microexon error"};

static List_T
make_microexon_pairs_double (int roffsetL, int roffsetM, int roffsetR,
			     int goffsetL, int goffsetM, int goffsetR,
			     int lengthL, int lengthM, int lengthR,
			     char *queryseq, char *queryuc,
			     Univcoord_T chroffset, Univcoord_T chrhigh, bool watsonp,
			     Pairpool_T pairpool, char gapchar, int dynprogindex) {
  List_T pairs = NULL;
  Pair_T gappair;
  char c1, c2, c2_alt;
  int i;

  /* Left segment */
  for (i = 0; i < lengthL; i++) {
    c1 = queryseq[roffsetL+i];

    c2 = get_genomic_nt(&c2_alt,goffsetL+i,chroffset,chrhigh,watsonp);
#ifdef EXTRACT_GENOMICSEG
    assert(c2 == genomicseg[goffsetL+i]);
#endif

    if (queryuc[roffsetL+i] == c2) {
      pairs = Pairpool_push(pairs,pairpool,roffsetL+i,goffsetL+i,c1,DYNPROG_MATCH_COMP,c2,c2_alt,
			    dynprogindex);
    } else if (consistent_array[(int) c1][(int) c2] == true) {
      pairs = Pairpool_push(pairs,pairpool,roffsetL+i,goffsetL+i,c1,AMBIGUOUS_COMP,c2,c2_alt,
			    dynprogindex);
    } else {
      pairs = Pairpool_push(pairs,pairpool,roffsetL+i,goffsetL+i,c1,MISMATCH_COMP,c2,c2_alt,
			    dynprogindex);
    }
  }

  /* First gap */
  /* Don't have to adjust querypos/genomepos, since no cdna/genome skips allowed */
  pairs = Pairpool_push_gapholder(pairs,pairpool,/*queryjump*/0,/*genomejump*/goffsetM-(goffsetL+lengthL),
				  /*leftpair*/NULL,/*rightpair*/NULL,/*knownp*/false);
  
  /* Assign pair->comp, because might occur after assign_gap_types */
  gappair = (Pair_T) List_head(pairs);
  gappair->comp = gapchar;


  /* Microexon */
  for (i = 0; i < lengthM; i++) {
    c1 = queryseq[roffsetM+i];

    c2 = get_genomic_nt(&c2_alt,goffsetM+i,chroffset,chrhigh,watsonp);
#ifdef EXTRACT_GENOMICSEG
    assert(c2 == genomicseg[goffsetM+i]);
#endif

    if (queryuc[roffsetM+i] == c2) {
      pairs = Pairpool_push(pairs,pairpool,roffsetM+i,goffsetM+i,c1,DYNPROG_MATCH_COMP,c2,c2_alt,
			    dynprogindex);
    } else if (consistent_array[(int) c1][(int) c2] == true) {
      pairs = Pairpool_push(pairs,pairpool,roffsetM+i,goffsetM+i,c1,AMBIGUOUS_COMP,c2,c2_alt,
			    dynprogindex);
    } else {
      pairs = Pairpool_push(pairs,pairpool,roffsetM+i,goffsetM+i,c1,MISMATCH_COMP,c2,c2_alt,
			    dynprogindex);
    }
  }

  /* Second gap */
  /* Don't have to adjust querypos/genomepos, since no cdna/genome skips allowed */
  if (lengthR == 0) {
    /* If lengthR is zero, then we will have a gap after a gap */
    Except_raise(&microexon_error,__FILE__,__LINE__);
  } else {
    pairs = Pairpool_push_gapholder(pairs,pairpool,/*queryjump*/0,/*genomejump*/goffsetR-(goffsetM+lengthM),
				    /*leftpair*/NULL,/*rightpair*/NULL,/*knownp*/false);
  }

  /* Assign pair->comp, because might occur after assign_gap_types */
  gappair = (Pair_T) List_head(pairs);
  gappair->comp = gapchar;

  
  /* Right segment */
  for (i = 0; i < lengthR; i++) {
    c1 = queryseq[roffsetR+i];

    c2 = get_genomic_nt(&c2_alt,goffsetR+i,chroffset,chrhigh,watsonp);
#ifdef EXTRACT_GENOMICSEG
    assert(c2 == genomicseg[goffsetR+i]);
#endif

    if (queryuc[roffsetR+i] == c2) {
      pairs = Pairpool_push(pairs,pairpool,roffsetR+i,goffsetR+i,c1,DYNPROG_MATCH_COMP,c2,c2_alt,
			    dynprogindex);
    } else if (consistent_array[(int) c1][(int) c2] == true) {
      pairs = Pairpool_push(pairs,pairpool,roffsetR+i,goffsetR+i,c1,AMBIGUOUS_COMP,c2,c2_alt,
			    dynprogindex);
    } else {
      pairs = Pairpool_push(pairs,pairpool,roffsetR+i,goffsetR+i,c1,MISMATCH_COMP,c2,c2_alt,
			    dynprogindex);
    }
  }

  return pairs;
}


#if 0
static List_T
make_microexon_pairs_single (int roffsetL, int roffsetR,
			     int goffsetL, int goffsetR,
			     int lengthL, int lengthR, char *queryseq, char *queryuc,
			     Univcoord_T chroffset, Univcoord_T chrhigh, bool watsonp,
			     Pairpool_T pairpool, char gapchar, int dynprogindex) {
  List_T pairs = NULL;
  Pair_T gappair;
  char c1, c2, c2_alt;
  int i;

  /* Microexon */
  for (i = 0; i < lengthL; i++) {
    c1 = queryseq[roffsetL+i];

    c2 = get_genomic_nt(&c2_alt,goffsetL+i,chroffset,chrhigh,watsonp);
#ifdef EXTRACT_GENOMICSEG
    assert(c2 == genomicseg[goffsetL+i]);
#endif

    if (queryuc[roffsetL+i] == c2) {
      pairs = Pairpool_push(pairs,pairpool,roffsetL+i,goffsetL+i,c1,DYNPROG_MATCH_COMP,c2,c2_alt,
			    dynprogindex);
    } else if (consistent_array[(int) c1][(int) c2] == true) {
      pairs = Pairpool_push(pairs,pairpool,roffsetL+i,goffsetL+i,c1,AMBIGUOUS_COMP,c2,c2_alt,
			    dynprogindex);
    } else {
      pairs = Pairpool_push(pairs,pairpool,roffsetL+i,goffsetL+i,c1,MISMATCH_COMP,c2,c2_alt,
			    dynprogindex);
    }
  }

  /* Gap */
  /* Don't have to adjust querypos/genomepos, since no cdna/genome skips allowed */
  pairs = Pairpool_push_gapholder(pairs,pairpool,/*queryjump*/0,/*genomejump*/goffsetR-(goffsetL+lengthL),
				  /*leftpair*/NULL,/*rightpair*/NULL,/*knownp*/false);

  /* Assign pair->comp, because might occur after assign_gap_types */
  gappair = (Pair_T) List_head(pairs);
  gappair->comp = gapchar;
  

  /* Right segment */
  for (i = 0; i < lengthR; i++) {
    c1 = queryseq[roffsetR+i];

    c2 = get_genomic_nt(&c2_alt,goffsetR+i,chroffset,chrhigh,watsonp);
#ifdef EXTRACT_GENOMICSEG
    assert(c2 == genomicseg[goffsetR+i]);
#endif

    if (queryuc[roffsetR+i] == c2) {
      pairs = Pairpool_push(pairs,pairpool,roffsetR+i,goffsetR+i,c1,DYNPROG_MATCH_COMP,c2,c2_alt,
			    dynprogindex);
    } else if (consistent_array[(int) c1][(int) c2] == true) {
      pairs = Pairpool_push(pairs,pairpool,roffsetR+i,goffsetR+i,c1,AMBIGUOUS_COMP,c2,c2_alt,
			    dynprogindex);
    } else {
      pairs = Pairpool_push(pairs,pairpool,roffsetR+i,goffsetR+i,c1,MISMATCH_COMP,c2,c2_alt,
			    dynprogindex);
    }
  }

  return pairs;
}
#endif


List_T
Dynprog_microexon_int (double *bestprob2, double *bestprob3, int *dynprogindex, int *microintrontype,
		       char *rsequence, char *sequenceuc1,
		       int rlength, int glengthL, int glengthR,
		       int roffset, int goffsetL, int rev_goffsetR, int cdna_direction,
#ifdef PMAP
		       char *queryaaseq, char *genomicuc,
#endif
		       char *queryseq, char *queryuc,
		       Univcoord_T chroffset, Univcoord_T chrhigh, bool watsonp,
		       Pairpool_T pairpool, double defect_rate) {
  List_T pairs = NULL;
#ifdef PMAP
  char *inst_rsequence;
#endif
  Intlist_T hits = NULL, p;
#ifdef EXTRACT_GENOMICSEG
  Intlist_T hits_old;
#endif
  int bestcL = -1, bestcR = -1, best_middlelength;
  int middlelength, cL, cR, mincR, maxcR, leftbound, rightbound, textleft, textright, candidate, i;
  int min_microexon_length, span, nmismatches;
  char left1, left2, right2, right1, left1_alt, left2_alt, right2_alt, right1_alt;
  char c, c_alt;
  char c1_alt, c2_alt, c3_alt, c4_alt;
  char intron1, intron2, intron3, intron4, gapchar;
  float pvalue, bestprob = 0.0, prob2, prob3;
  Univcoord_T splicesitepos;


  *bestprob2 = *bestprob3 = 0.0;

  if (defect_rate < DEFECT_HIGHQ) {
    pvalue = MICROEXON_PVALUE_HIGHQ;
  } else if (defect_rate < DEFECT_MEDQ) {
    pvalue = MICROEXON_PVALUE_MEDQ;
  } else {
    pvalue = MICROEXON_PVALUE_LOWQ;
  }

#ifdef PMAP
  intron1 = 'G';
  intron2 = 'T';
  intron3 = 'A';
  intron4 = 'G';
  gapchar = FWD_CANONICAL_INTRON_COMP;
  *microintrontype = GTAG_FWD;
#else
  if (cdna_direction > 0) {
    intron1 = 'G';
    intron2 = 'T';
    intron3 = 'A';
    intron4 = 'G';
    gapchar = FWD_CANONICAL_INTRON_COMP;
    *microintrontype = GTAG_FWD;
  } else if (cdna_direction < 0) {
    intron1 = 'C';
    intron2 = 'T';
    intron3 = 'A';
    intron4 = 'C';
    gapchar = REV_CANONICAL_INTRON_COMP;
    *microintrontype = GTAG_REV;
  } else {
    /* Can occur when called by Stage3_merge_local_splice */
    /* fprintf(stderr,"cdna_direction is 0 in Dynprog_microexon_int\n"); */
    *microintrontype = NONINTRON;
    return NULL;
  }
#endif

#ifdef EXTRACT_GENOMICSEG
  debug(printf("Begin microexon search for %.*s and %.*s\n",
	       glengthL,gsequenceL,glengthR,&(rev_gsequenceR[-glengthR+1])));
#else
  debug(printf("Begin microexon search\n"));
#endif

#ifdef PMAP
  inst_rsequence = instantiate_codons(rsequence,queryaaseq,roffset,rlength);
  debug(printf("  Query sequence is %.*s\n",rlength,inst_rsequence));
#else
  debug(printf("  Query sequence is %.*s\n",rlength,rsequence));
#endif
  span = rev_goffsetR-goffsetL;
  debug(printf("  Genomic span is of length %d\n",span));

  if (span <= 0) {
    fprintf(stderr,"Bug in Dynprog_microexon_int.  span %d <= 0.  Please report to twu@gene.com\n",span);
    abort();
  } else {
    min_microexon_length = ceilf(-fasterlog(1.0-powf(1.0-pvalue,1.0/(float) span)) / /*log(4)*/1.386294);
  }
  min_microexon_length -= 8;	/* Two donor-acceptor pairs */
  debug(printf("  Min microexon length is %d\n",min_microexon_length));
  if (min_microexon_length > MAX_MICROEXON_LENGTH) {
#ifdef PMAP
    FREE(inst_rsequence);
#endif
    *microintrontype = NONINTRON;
    return NULL;
  } else if (min_microexon_length < MIN_MICROEXON_LENGTH) {
    min_microexon_length = MIN_MICROEXON_LENGTH;
  }

  debug(printf("\nFinding starting boundary on left\n"));
  leftbound = 0;
  nmismatches = 0;
  while (leftbound < rlength - 1 && nmismatches <= 1) {
    debug(printf("  leftbound = %d, nmismatches = %d.",leftbound,nmismatches));
#ifdef PMAP
    c = get_genomic_nt(&c_alt,goffsetL+leftbound,chroffset,chrhigh,watsonp);
    debug(printf("  Comparing %c with %c\n",inst_rsequence[leftbound],c));
    if (matchtable[inst_rsequence[leftbound]-'A'][c-'A'] == false) {
      nmismatches++;
    }
#else
    c = get_genomic_nt(&c_alt,goffsetL+leftbound,chroffset,chrhigh,watsonp);
#ifdef EXTRACT_GENOMICSEG
    assert(c == gsequence_ucL[leftbound]);
#endif
    debug(printf("  Comparing %c with %c\n",sequenceuc1[leftbound],c));
    if (sequenceuc1[leftbound] != c) {
      nmismatches++;
    }
#endif
    leftbound++;
  }
  leftbound--;			/* This is where the leftmost mismatch occurred */

  debug(printf("\nFinding starting boundary on right\n"));
  rightbound = 0;
  i = rlength-1;
  nmismatches = 0;
  while (i >= 0 && nmismatches <= 1) {
    debug(printf("  rightbound = %d, nmismatches = %d.",rightbound,nmismatches));
#ifdef PMAP
    c = get_genomic_nt(&c_alt,rev_goffsetR-rightbound,chroffset,chrhigh,watsonp);
    debug(printf("  Comparing %c with %c\n",inst_rsequence[i],c));
    if (matchtable[inst_rsequence[i]-'A'][c-'A'] == false) {
      nmismatches++;
    }
#else
    c = get_genomic_nt(&c_alt,rev_goffsetR-rightbound,chroffset,chrhigh,watsonp);
#ifdef EXTRACT_GENOMICSEG
    assert(c == rev_gsequence_ucR[-rightbound]);
#endif
    debug(printf("  Comparing %c with %c\n",sequenceuc1[i],c));
    if (sequenceuc1[i] != c) {
      nmismatches++;
    }
#endif
    rightbound++;
    i--;
  }
  rightbound--;			/* This is where the rightmost mismatch occurred */

  debug(printf("  Left must start before %d from left end of query.  Right must start after %d from right end of query\n",
	       leftbound,rightbound));

  /* We require that cL >= 1 and cR >= 1 so that lengthL and lengthR are >= 1 */
  for (cL = 1; cL <= leftbound; cL++) {
    left1 = get_genomic_nt(&left1_alt,goffsetL+cL,chroffset,chrhigh,watsonp);
    left2 = get_genomic_nt(&left2_alt,goffsetL+cL+1,chroffset,chrhigh,watsonp);
#ifdef EXTRACT_GENOMICSEG
    assert(left1 == gsequence_ucL[cL]);
    assert(left2 == gsequence_ucL[cL+1]);
#endif

    debug(printf("  %d: %c%c\n",cL,left1,left2));
    if (left1 == intron1 && left2 == intron2) {
      mincR = rlength - MAX_MICROEXON_LENGTH - cL;
      if (mincR < 1) {
	mincR = 1;
      }
      maxcR = rlength - min_microexon_length - cL;
      if (maxcR > rightbound) {
	maxcR = rightbound;
	}
      debug(printf("  Found left GT at %d.  Scanning from %d - cL - (1-7), or %d to %d\n",
		   cL,rlength,mincR,maxcR));
      for (cR = mincR; cR <= maxcR; cR++) {
	right2 = get_genomic_nt(&right2_alt,rev_goffsetR-cR-1,chroffset,chrhigh,watsonp);
	right1 = get_genomic_nt(&right1_alt,rev_goffsetR-cR,chroffset,chrhigh,watsonp);
#ifdef EXTRACT_GENOMICSEG
	assert(right2 == rev_gsequence_ucR[-cR-1]);
	assert(right1 == rev_gsequence_ucR[-cR]);
#endif
	debug(printf("   Checking %d: %c%c\n",cR,right2,right1));
	if (right2 == intron3 && right1 == intron4) {
	  middlelength = rlength - cL - cR;
#ifdef PMAP
	  debug(printf("  Found pair at %d to %d, length %d.  Middle sequence is %.*s\n",
		       cL,cR,middlelength,middlelength,&(inst_rsequence[cL])));
#else
	  debug(printf("  Found pair at %d to %d, length %d.  Middle sequence is %.*s\n",
		       cL,cR,middlelength,middlelength,&(rsequence[cL])));
#endif
	  
	  textleft = goffsetL + cL + MICROINTRON_LENGTH;
	  textright = rev_goffsetR - cR - MICROINTRON_LENGTH;

	  if (textright >= textleft + middlelength) {
#ifdef PMAP
	    hits = BoyerMoore(&(inst_rsequence[cL]),middlelength,&(genomicuc[textleft]),textright-textleft);
#else
	    hits = BoyerMoore_nt(&(sequenceuc1[cL]),middlelength,textleft,textright-textleft,
				 chroffset,chrhigh,watsonp);
#ifdef EXTRACT_GENOMICSEG
	    hits_old = BoyerMoore(&(sequenceuc1[cL]),middlelength,&(genomicuc[textleft]),textright-textleft);
	    assert(Intlist_equal(hits,hits_old));
	    Intlist_free(&hits_old);
#endif
#endif
	    for (p = hits; p != NULL; p = Intlist_next(p)) {
	      candidate = textleft + Intlist_head(p);
#ifdef EXTRACT_GENOMICSEG
	      assert(get_genomic_nt(candidate-2,chroffset,chrhigh,watsonp) == genomicuc[candidate - 2]);
	      assert(get_genomic_nt(candidate-1,chroffset,chrhigh,watsonp) == genomicuc[candidate - 1]);
	      assert(get_genomic_nt(candidate+middlelength,chroffset,chrhigh,watsonp) == genomicuc[candidate + middlelength]);
	      assert(get_genomic_nt(candidate+middlelength+1,chroffset,chrhigh,watsonp) == genomicuc[candidate + middlelength+1]);
#endif
	      if (/*genomicuc[candidate - 2]*/ get_genomic_nt(&c3_alt,candidate-2,chroffset,chrhigh,watsonp) == intron3 &&
		  /*genomicuc[candidate - 1]*/ get_genomic_nt(&c4_alt,candidate-1,chroffset,chrhigh,watsonp)  == intron4 &&
		  /*genomicuc[candidate + middlelength]*/ get_genomic_nt(&c1_alt,candidate+middlelength,chroffset,chrhigh,watsonp) == intron1 &&
		  /*genomicuc[candidate + middlelength + 1]*/ get_genomic_nt(&c2_alt,candidate+middlelength+1,chroffset,chrhigh,watsonp) == intron2) {
		debug(printf("  Successful microexon at %d >>> %d..%d >>> %d\n",goffsetL+cL,candidate,candidate+middlelength,rev_goffsetR-cR));

		/* Not handling known splice sites yet */
		if (watsonp == true) {
		  if (cdna_direction > 0) {
		    splicesitepos = chroffset + (candidate-1) + 1;
		    prob2 = Maxent_hr_acceptor_prob(splicesitepos,chroffset);
		    splicesitepos = chroffset + candidate+middlelength;
		    prob3 = Maxent_hr_donor_prob(splicesitepos,chroffset);
		  } else {
		    splicesitepos = chroffset + (candidate-1) + 1;
		    prob2 = Maxent_hr_antidonor_prob(splicesitepos,chroffset);
		    splicesitepos = chroffset + candidate+middlelength;
		    prob3 = Maxent_hr_antiacceptor_prob(splicesitepos,chroffset);
		  }
		} else {
		  if (cdna_direction > 0) {
		    splicesitepos = chrhigh - (candidate-1);
		    prob2 = Maxent_hr_antiacceptor_prob(splicesitepos,chroffset);
		    splicesitepos = chrhigh - (candidate+middlelength) + 1;
		    prob3 = Maxent_hr_antidonor_prob(splicesitepos,chroffset);
		  } else {
		    splicesitepos = chrhigh - (candidate-1);
		    prob2 = Maxent_hr_donor_prob(splicesitepos,chroffset);
		    splicesitepos = chrhigh - (candidate+middlelength) + 1;
		    prob3 = Maxent_hr_acceptor_prob(splicesitepos,chroffset);
		  }
		}
	      
		debug(printf("microexon probabilities: prob2 = %f, prob3 = %f\n",prob2,prob3));
		if (prob2 + prob3 > bestprob) {
		  bestcL = cL;
		  bestcR = cR;
		  best_middlelength = middlelength;
		  *bestprob2 = prob2;
		  *bestprob3 = prob3;
		  bestprob = prob2 + prob3;
		}
	      }
	    }
	    Intlist_free(&hits);
	  }

	}
      }
    }
  }

  if (bestcL < 0 || bestcR < 0) {
    debug(printf("End of dynprog microexon int\n"));

#ifdef PMAP
    FREE(inst_rsequence);
#endif

    *microintrontype = NONINTRON;
    return NULL;

  } else {
    pairs = make_microexon_pairs_double(roffset,roffset+bestcL,roffset+bestcL+best_middlelength,
					goffsetL,candidate,rev_goffsetR-bestcR+1,
					/*lengthL*/bestcL,/*lengthM*/best_middlelength,/*lengthR*/bestcR,
#ifdef PMAP
					&(inst_rsequence[-roffset]),&(inst_rsequence[-roffset]),
#else
					queryseq,queryuc,
#endif
					chroffset,chrhigh,watsonp,pairpool,gapchar,*dynprogindex);
#ifdef PMAP
    FREE(inst_rsequence);
#endif
    *dynprogindex += (*dynprogindex > 0 ? +1 : -1);
    return pairs;
  }
}


#if 0
/* Based on probability of seeing a pattern of length n in L is
   1-(1-p1)^L, where p1 is 4^n.  We determine L so chance probability
   is less than ENDSEQUENCE_PVALUE */
static int
search_length (int endlength, int maxlength, bool end_microexons_p) {
  double p1;
  int effective_maxlength, extrant, result;

  if (end_microexons_p == true) {
    extrant = 4;		/* Count the four nucleotides at the intron bounds */
    effective_maxlength = maxlength;
  } else {
    extrant = 0;		/* Don't count the four nucleotides */
    effective_maxlength = 5000;
    if (maxlength < effective_maxlength) {
      effective_maxlength = maxlength;
    }
  }

  if (endlength + extrant > 12) {
    debug(printf("  Search length for endlength of %d is maxlength %d\n",endlength,effective_maxlength));
    return effective_maxlength;
  } else {
    p1 = 1.0/pow(4.0,(double) (endlength + extrant));
    result = (int) (fasterlog(1.0-ENDSEQUENCE_PVALUE)/fasterlog(1-p1));
    debug(printf("  Search length for endlength of %d plus extra nt of %d is %d\n",endlength,extrant,result));
    if (result > effective_maxlength) {
      return effective_maxlength;
    } else {
      return result;
    }
  }
}
#endif


#if 0
/* Not currently used */
List_T
Dynprog_microexon_5 (int *dynprogindex, int *microintrontype, int *microexonlength,
		     char *rev_rsequence, char *revsequenceuc1, char *rev_gsequence, char *rev_gsequence_uc,
		     int rlength, int glength, int rev_roffset, int rev_goffset, int cdna_direction,
#ifdef PMAP
		     char *queryaaseq,
#endif
		     char *queryseq, char *queryuc, char *genomicseg, char *genomicuc,
		     Pairpool_T pairpool, bool end_microexons_p) {
  List_T pairs = NULL;
#ifdef PMAP
  char *inst_rsequence, *inst_rev_rsequence;
#endif
  Intlist_T hits = NULL, p;
  int endlength, maxc, c, textleft, textright, candidate, nmismatches = 0;
  char right2, right1;
  char intron1, intron2, intron3, intron4, gapchar;

#ifdef PMAP
  intron1 = 'G';
  intron2 = 'T';
  intron3 = 'A';
  intron4 = 'G';
  gapchar = FWD_CANONICAL_INTRON_COMP;
  *microintrontype = GTAG_FWD;
#else
  if (cdna_direction > 0) {
    intron1 = 'G';
    intron2 = 'T';
    intron3 = 'A';
    intron4 = 'G';
    gapchar = FWD_CANONICAL_INTRON_COMP;
    *microintrontype = GTAG_FWD;
  } else if (cdna_direction < 0) {
    intron1 = 'C';
    intron2 = 'T';
    intron3 = 'A';
    intron4 = 'C';
    gapchar = REV_CANONICAL_INTRON_COMP;
    *microintrontype = GTAG_REV;
  } else {
    *microintrontype = NONINTRON;
    return (List_T) NULL;
    abort();
  }
#endif

#ifdef EXTRACT_GENOMICSEG
  debug(printf("Begin microexon search at 5' for %.*s\n",
	       glength,&(rev_gsequence[-glength+1])));
#else
  debug(printf("Begin microexon search at 5'\n"));
#endif

#ifdef PMAP
  inst_rsequence = instantiate_codons(&(rev_rsequence[-rlength+1]),queryaaseq,rev_roffset-rlength+1,rlength);
  inst_rev_rsequence = &(inst_rsequence[rlength-1]);
  debug(printf("  Query sequence is %.*s\n",rlength,&(inst_rev_rsequence[-rlength+1])));
#else
  debug(printf("  Query sequence is %.*s\n",rlength,&(rev_rsequence[-rlength+1])));
#endif

  *microexonlength = 0;
  if (glength < rlength) {
    maxc = glength - MIN_MICROEXON_LENGTH;
  } else {
    maxc = rlength - MIN_MICROEXON_LENGTH;
  }
  for (c = 0; c < maxc; c++) {
    right2 = rev_gsequence_uc[-c-1];
    right1 = rev_gsequence_uc[-c];
    debug(printf("   Checking %c%c\n",right2,right1));
#ifdef PMAP
    if (c > 0 && matchtable[inst_rev_rsequence[-c+1]-'A'][rev_gsequence_uc[-c+1]-'A'] == false) {
      nmismatches++;
    }
#else
    if (c > 0 && revsequenceuc1[-c+1] != rev_gsequence_uc[-c+1]) {
      nmismatches++;
    }
#endif
    if (nmismatches > 1) {
#ifdef PMAP
      debug(printf("   Aborting at %c !~ %c\n",inst_rev_rsequence[-c+1],rev_gsequence[-c+1]));
      FREE(inst_rsequence);
#else
      debug(printf("   Aborting at %c != %c\n",rev_rsequence[-c+1],rev_gsequence[-c+1]));
#endif
      *microintrontype = NONINTRON;
      return NULL;
    }
    if (right2 == intron3 && right1 == intron4) {
      endlength = rlength - c;
#ifdef PMAP
      debug(printf("  Found acceptor at %d, length %d.  End sequence is %.*s\n",
		       c,endlength,endlength,&(inst_rev_rsequence[-endlength+1])));
#else
      debug(printf("  Found acceptor at %d, length %d.  End sequence is %.*s\n",
		       c,endlength,endlength,&(rev_rsequence[-endlength+1])));
#endif

      textright = rev_goffset - c - MICROINTRON_LENGTH;
      textleft = textright - search_length(endlength,textright,end_microexons_p) + MICROINTRON_LENGTH;

      if (textright >= textleft + endlength) {
#ifdef PMAP
	hits = BoyerMoore(&(inst_rev_rsequence[-c-endlength+1]),endlength,&(genomicuc[textleft]),textright-textleft);
#else
	hits = BoyerMoore(&(revsequenceuc1[-c-endlength+1]),endlength,&(genomicuc[textleft]),textright-textleft);
#endif
	for (p = hits; p != NULL; p = Intlist_next(p)) {
	  candidate = textleft + Intlist_head(p);
	  if (genomicseg[candidate + endlength] == intron1 &&
	      genomicseg[candidate + endlength + 1] == intron2) {
	    debug(printf("  Successful microexon at %d\n",candidate));

	    Intlist_free(&hits);
	    *microexonlength = endlength;
	    pairs = make_microexon_pairs_single(rev_roffset-c-endlength+1,rev_roffset-c+1,
						candidate,rev_goffset-c+1,endlength,c,
#ifdef PMAP
						&(inst_rev_rsequence[-rev_roffset]),&(inst_rev_rsequence[-rev_roffset]),
#else
						queryseq,queryuc,
#endif
						chroffset,watsonp,pairpool,gapchar,*dynprogindex);
#ifdef PMAP
	    FREE(inst_rsequence);
#endif
	    *dynprogindex += (*dynprogindex > 0 ? +1 : -1);
	    return pairs;
	  }
	}
	Intlist_free(&hits);
      }

    }
  }

#ifdef PMAP
  FREE(inst_rsequence);
#endif

  debug(printf("End of dynprog microexon 5\n"));

  *microintrontype = NONINTRON;
  return NULL;
}
#endif


#if 0
/* Not currently used */
List_T
Dynprog_microexon_3 (int *dynprogindex, int *microintrontype, int *microexonlength, 
		     char *rsequence, char *sequenceuc1, char *gsequence, char *gsequence_uc,
		     int rlength, int glength, int roffset, int goffset, int cdna_direction,
#ifdef PMAP
		     char *queryaaseq,
#endif
		     char *queryseq, char *queryuc, char *genomicseg, char *genomicuc,
		     int genomiclength, Pairpool_T pairpool, bool end_microexons_p) {
  List_T pairs = NULL;
#ifdef PMAP
  char *inst_rsequence;
#endif
  Intlist_T hits = NULL, p;
  int endlength, maxc, c, textleft, textright, candidate, nmismatches = 0;
  char left1, left2;
  char intron1, intron2, intron3, intron4, gapchar;

#ifdef PMAP
  intron1 = 'G';
  intron2 = 'T';
  intron3 = 'A';
  intron4 = 'G';
  gapchar = FWD_CANONICAL_INTRON_COMP;
  *microintrontype = GTAG_FWD;
#else
  if (cdna_direction > 0) {
    intron1 = 'G';
    intron2 = 'T';
    intron3 = 'A';
    intron4 = 'G';
    gapchar = FWD_CANONICAL_INTRON_COMP;
    *microintrontype = GTAG_FWD;
  } else if (cdna_direction < 0) {
    intron1 = 'C';
    intron2 = 'T';
    intron3 = 'A';
    intron4 = 'C';
    gapchar = REV_CANONICAL_INTRON_COMP;
    *microintrontype = GTAG_REV;
  } else {
    *microintrontype = NONINTRON;
    return (List_T) NULL;
    abort();
  }
#endif

#ifdef EXTRACT_GENOMICSEG
  debug(printf("Begin microexon search at 3' for %.*s\n",glength,gsequence));
#else
  debug(printf("Begin microexon search at 3'\n"));
#endif

#ifdef PMAP
  inst_rsequence = instantiate_codons(rsequence,queryaaseq,roffset,rlength);
  debug(printf("  Query sequence is %.*s\n",rlength,inst_rsequence));
#else
  debug(printf("  Query sequence is %.*s\n",rlength,rsequence));
#endif

  *microexonlength = 0;
  if (glength < rlength) {
    maxc = glength - MIN_MICROEXON_LENGTH;
  } else {
    maxc = rlength - MIN_MICROEXON_LENGTH;
  }
  for (c = 0; c < maxc; c++) {
    left1 = gsequence_uc[c];
    left2 = gsequence_uc[c+1];
    debug(printf("   Checking %c%c\n",left1,left2));
#ifdef PMAP
    if (c > 0 && matchtable[inst_rsequence[c-1]-'A'][gsequence_uc[c-1]-'A'] == false) {
      nmismatches++;
    }
#else
    if (c > 0 && sequenceuc1[c-1] != gsequence_uc[c-1]) {
      nmismatches++;
    }
#endif
    if (nmismatches > 1) {
#ifdef PMAP
      debug(printf("   Aborting at %c !~ %c\n",inst_rsequence[c-1],gsequence[c-1]));
      FREE(inst_rsequence);
#else
      debug(printf("   Aborting at %c != %c\n",rsequence[c-1],gsequence[c-1]));
#endif
      *microintrontype = NONINTRON;
      return NULL;
    }
    if (left1 == intron1 && left2 == intron2) {
      endlength = rlength - c;
#ifdef PMAP
      debug(printf("  Found donor at %d, length %d.  End sequence is %.*s\n",
		   c,endlength,endlength,&(inst_rsequence[c])));
#else
      debug(printf("  Found donor at %d, length %d.  End sequence is %.*s\n",
		   c,endlength,endlength,&(rsequence[c])));
#endif

      textleft = goffset + c;
      textright = textleft + search_length(endlength,genomiclength-textleft,end_microexons_p);
      
      if (textright >= textleft + endlength) {
#ifdef PMAP
	hits = BoyerMoore(&(inst_rsequence[c]),endlength,&(genomicuc[textleft]),textright-textleft);
#else
	hits = BoyerMoore(&(sequenceuc1[c]),endlength,&(genomicuc[textleft]),textright-textleft);
#endif
	for (p = hits; p != NULL; p = Intlist_next(p)) {
	  candidate = textleft + Intlist_head(p);
	  if (genomicseg[candidate - 2] == intron3 &&
	      genomicseg[candidate - 1] == intron4) {
	    debug(printf("  Successful microexon at %d\n",candidate));

	    Intlist_free(&hits);
	    *microexonlength = endlength;
	    pairs = make_microexon_pairs_single(roffset,roffset+c,
						goffset,candidate,c,endlength,
#ifdef PMAP
						&(inst_rsequence[-roffset]),&(inst_rsequence[-roffset]),
#else
						queryseq,queryuc,
#endif
						genomicseg,genomicuc,
						pairpool,gapchar,*dynprogindex);
#ifdef PMAP
	    FREE(inst_rsequence);
#endif
	    *dynprogindex += (*dynprogindex > 0 ? +1 : -1);
	    return pairs;
	  }
	}
	Intlist_free(&hits);
      }

    }
  }

#ifdef PMAP
  FREE(inst_rsequence);
#endif

  debug(printf("End of dynprog microexon 3\n"));

  *microintrontype = NONINTRON;
  return NULL;
}
#endif

