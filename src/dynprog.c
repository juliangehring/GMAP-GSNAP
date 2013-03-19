static char rcsid[] = "$Id: dynprog.c 87907 2013-03-05 02:23:47Z twu $";

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "dynprog.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>		/* For ceil, log, pow */
#include <ctype.h>		/* For tolower */
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






#ifdef DEBUG2
#define NEG_INFINITY -99
#else
#define NEG_INFINITY -1000000
#endif

#define ONESIDEGAP 1

/*
#define RIGHTANGLE 1
*/

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
#define MISMATCH_MEDQ -2
#define MISMATCH_LOWQ -1

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

#define PAIRED_OPEN_HIGHQ -18
#define PAIRED_OPEN_MEDQ -18
#define PAIRED_OPEN_LOWQ -18

#define PAIRED_EXTEND_HIGHQ -3
#define PAIRED_EXTEND_MEDQ -3
#define PAIRED_EXTEND_LOWQ -3

#define SINGLE_OPEN_HIGHQ -10
#define SINGLE_OPEN_MEDQ -10
#define SINGLE_OPEN_LOWQ -10

#define SINGLE_EXTEND_HIGHQ -3
#define SINGLE_EXTEND_MEDQ -3
#define SINGLE_EXTEND_LOWQ -3


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
#define END_OPEN_HIGHQ -12
#define END_OPEN_MEDQ -12
#define END_OPEN_LOWQ -12

#define END_EXTEND_HIGHQ -1
#define END_EXTEND_MEDQ -1
#define END_EXTEND_LOWQ -1


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


typedef char Direction_T;
#define VERT 4
#define HORIZ 2
#define DIAG 1
#define STOP 0


static IIT_T splicing_iit;
static int *splicing_divint_crosstable;
static int donor_typeint;
static int acceptor_typeint;

static Genomicpos_T *splicesites;
static Splicetype_T *splicetypes;
static Genomicpos_T *splicedists;
static int nsplicesites;
static unsigned int *trieoffsets_obs;
static unsigned int *triecontents_obs;
static unsigned int *trieoffsets_max;
static unsigned int *triecontents_max;

static bool novelsplicingp;


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
	       Genomicpos_T *splicesites_in, Splicetype_T *splicetypes_in,
	       Genomicpos_T *splicedists_in, int nsplicesites_in,
	       unsigned int *trieoffsets_obs_in, unsigned int *triecontents_obs_in,
	       unsigned int *trieoffsets_max_in, unsigned int *triecontents_max_in) {
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
get_genomic_nt (char *g_alt, int genomicpos, Genomicpos_T chroffset, Genomicpos_T chrhigh,
		bool watsonp) {
  char c2, c2_alt;
  Genomicpos_T pos;

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
      return '*';

    } else if (pos >= chrhigh) {
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

/* Makes a matrix of dimensions 0..length1 x 0..length2 inclusive */
static int **
Matrix_alloc (int length1, int length2, int **ptrs, int *space) {
  int **matrix, i;

  if (length1 <= 0 || length2 <= 0) {
    fprintf(stderr,"dynprog: lengths are negative: %d %d\n",length1,length2);
    abort();
  }

  matrix = ptrs;
  matrix[0] = space;
  for (i = 1; i <= length1; i++) {
    matrix[i] = &(matrix[i-1][length2 + 1]);
  }

  /* Clear memory only around the band, but doesn't work with Gotoh P1
     and Q1 matrices */
  /*
  for (r = 0; r <= length1; r++) {
    if ((clo = r - lband - 1) >= 0) {
      matrix[r][clo] = 0;
    }
    if ((chigh = r + rband + 1) <= length2) {
      matrix[r][chigh] = 0;
    }
  }
  */
  memset((void *) space,0,(length1+1)*(length2+1)*sizeof(int));

  return matrix;
}


struct Int3_T {
  int gap1;
  int gap2;
  int nogap;
};


/* Makes a matrix of dimensions 0..length1 x 0..length2 inclusive */
static struct Int3_T **
Matrix3_alloc (int length1, int length2, struct Int3_T **ptrs, struct Int3_T *space) {
  struct Int3_T **matrix;
  int i;

  if (length1 <= 0 || length2 <= 0) {
    fprintf(stderr,"dynprog: lengths are negative: %d %d\n",length1,length2);
    abort();
  }

  matrix = ptrs;
  matrix[0] = space;
  for (i = 1; i <= length1; i++) {
    matrix[i] = &(matrix[i-1][length2 + 1]);
  }

  /* Clear memory only around the band, but doesn't work with Gotoh P1
     and Q1 matrices */
  /*
  for (r = 0; r <= length1; r++) {
    if ((clo = r - lband - 1) >= 0) {
      matrix[r][clo] = 0;
    }
    if ((chigh = r + rband + 1) <= length2) {
      matrix[r][chigh] = 0;
    }
  }
  */
  memset((void *) space,0,(length1+1)*(length2+1)*sizeof(struct Int3_T));

  return matrix;
}


static void
Matrix3_print (struct Int3_T **matrix, int length1, int length2, char *sequence1,
	       int offset2, Genomicpos_T chroffset, Genomicpos_T chrhigh,
	       bool watsonp, bool revp) {
  int i, j;
  char g_alt;

  printf("G1");
  for (j = 0; j <= length2; ++j) {
    if (j == 0) {
      printf("    ");
    } else {
#ifdef EXTRACT_GENOMICSEG
      printf("  %c ",revp ? sequence2[-j+1] : sequence2[j-1]);
#else
      if (revp == false) {
	printf("  %c ",get_genomic_nt(&g_alt,offset2+j-1,chroffset,chrhigh,watsonp));
      } else {
	printf("  %c ",get_genomic_nt(&g_alt,offset2+1-j,chroffset,chrhigh,watsonp));
      }
#endif
    }
  }
  printf("\n");

  for (i = 0; i <= length1; ++i) {
    if (i == 0) {
      printf("  ");
    } else {
      printf("%c ",revp ? sequence1[-i+1] : sequence1[i-1]);
    }
    for (j = 0; j <= length2; ++j) {
      if (matrix[i][j].gap1 < NEG_INFINITY) {
	printf("%3d ",NEG_INFINITY);
      } else {
	printf("%3d ",matrix[i][j].gap1);
      }
    }
    printf("\n");
  }
  printf("\n");


  printf("NG");
  for (j = 0; j <= length2; ++j) {
    if (j == 0) {
      printf("    ");
    } else {
#ifdef EXTRACT_GENOMICSEG
      printf("  %c ",revp ? sequence2[-j+1] : sequence2[j-1]);
#else
      if (revp == false) {
	printf("  %c ",get_genomic_nt(&g_alt,offset2+j-1,chroffset,chrhigh,watsonp));
      } else {
	printf("  %c ",get_genomic_nt(&g_alt,offset2+1-j,chroffset,chrhigh,watsonp));
      }
#endif
    }
  }
  printf("\n");

  for (i = 0; i <= length1; ++i) {
    if (i == 0) {
      printf("  ");
    } else {
      printf("%c ",revp ? sequence1[-i+1] : sequence1[i-1]);
    }
    for (j = 0; j <= length2; ++j) {
      if (matrix[i][j].nogap < NEG_INFINITY) {
	printf("%3d ",NEG_INFINITY);
      } else {
	printf("%3d ",matrix[i][j].nogap);
      }
    }
    printf("\n");
  }
  printf("\n");


  printf("G2");
  for (j = 0; j <= length2; ++j) {
    if (j == 0) {
      printf("    ");
    } else {
#ifdef EXTRACT_GENOMICSEG
      printf("  %c ",revp ? sequence2[-j+1] : sequence2[j-1]);
#else
      if (revp == false) {
	printf("  %c ",get_genomic_nt(&g_alt,offset2+j-1,chroffset,chrhigh,watsonp));
      } else {
	printf("  %c ",get_genomic_nt(&g_alt,offset2+1-j,chroffset,chrhigh,watsonp));
      }
#endif
    }
  }
  printf("\n");

  for (i = 0; i <= length1; ++i) {
    if (i == 0) {
      printf("  ");
    } else {
      printf("%c ",revp ? sequence1[-i+1] : sequence1[i-1]);
    }
    for (j = 0; j <= length2; ++j) {
      if (matrix[i][j].gap2 < NEG_INFINITY) {
	printf("%3d ",NEG_INFINITY);
      } else {
	printf("%3d ",matrix[i][j].gap2);
      }
    }
    printf("\n");
  }
  printf("\n");

  return;
}


/************************************************************************/
/*  Directions  */
/************************************************************************/

/* Makes a matrix of dimensions 0..length1 x 0..length2 inclusive */
static Direction_T **
Directions_alloc (int length1, int length2, Direction_T **ptrs, Direction_T *space) {
  Direction_T **directions;
  int i;

  directions = ptrs;
  directions[0] = space;
  for (i = 1; i <= length1; i++) {
    directions[i] = &(directions[i-1][length2 + 1]);
  }

  /* Clear memory only around the band, but may not work with Gotoh
     method */
  /*
    for (r = 0; r <= length1; r++) {
    if ((clo = r - lband - 1) >= 0) {
    directions[r][clo] = STOP;
    }
    if ((chigh = r + rband + 1) <= length2) {
    directions[r][chigh] = STOP;
    }
    }
  */
  memset((void *) space,0,(length1+1)*(length2+1)*sizeof(Direction_T));

  return directions;
}

static void
Directions_print (Direction_T **directions, int **jump, int length1, int length2, 
		  char *sequence1, bool revp) {
  int i, j;
  char buffer[4];

#ifdef EXTRACT_GENOMICSEG
  printf("  ");
  for (j = 0; j <= length2; ++j) {
    if (j == 0) {
      printf("    ");
    } else {
      printf("  %c ",revp ? sequence2[-j+1] : sequence2[j-1]);
    }
  }
  printf("\n");
#endif

  for (i = 0; i <= length1; ++i) {
    if (i == 0) {
      printf("  ");
    } else {
      printf("%c ",revp ? sequence1[-i+1] : sequence1[i-1]);
    }
    for (j = 0; j <= length2; ++j) {
      if (directions[i][j] == DIAG) {
	sprintf(buffer,"D%d",jump[i][j]);
      } else if (directions[i][j] == HORIZ) {
	sprintf(buffer,"H%d",jump[i][j]);
      } else if (directions[i][j] == VERT) {
	sprintf(buffer,"V%d",jump[i][j]);
      } else {
	sprintf(buffer,"S%d",0);
      }
      printf("%3s ",buffer);
    }
    printf("\n");
  }
  printf("\n");
  return;
}


struct Direction3_T {
  Direction_T gap1;
  Direction_T gap2;
  Direction_T nogap;
};


/* Makes a matrix of dimensions 0..length1 x 0..length2 inclusive */
static struct Direction3_T **
Directions3_alloc (int length1, int length2, struct Direction3_T **ptrs, struct Direction3_T *space) {
  struct Direction3_T **directions;
  int i;

  directions = ptrs;
  directions[0] = space;
  for (i = 1; i <= length1; i++) {
    directions[i] = &(directions[i-1][length2 + 1]);
  }

  /* Clear memory only around the band, but may not work with Gotoh
     method */
  /*
    for (r = 0; r <= length1; r++) {
    if ((clo = r - lband - 1) >= 0) {
    directions[r][clo] = STOP;
    }
    if ((chigh = r + rband + 1) <= length2) {
    directions[r][chigh] = STOP;
    }
    }
  */
  memset((void *) space,0,(length1+1)*(length2+1)*sizeof(struct Direction3_T));

  return directions;
}


static void
Directions3_print (struct Direction3_T **directions, int length1, int length2, 
		   char *sequence1, bool revp) {
  int i, j;

#ifdef EXTRACT_GENOMICSEG
  printf("  ");
  for (j = 0; j <= length2; ++j) {
    if (j == 0) {
      printf("    ");
    } else {
      printf("  %c   ",revp ? sequence2[-j+1] : sequence2[j-1]);
    }
  }
  printf("\n");
#endif

  for (i = 0; i <= length1; ++i) {
    if (i == 0) {
      printf("  ");
    } else {
      printf("%c ",revp ? sequence1[-i+1] : sequence1[i-1]);
    }
    for (j = 0; j <= length2; ++j) {
      if (directions[i][j].gap1 == DIAG) {
	printf("D");
      } else if (directions[i][j].gap1 == HORIZ) {
	printf("H");
      } else if (directions[i][j].gap1 == VERT) {
	printf("V");
      } else {
	printf("S");
      }
      printf("|");
      if (directions[i][j].nogap == DIAG) {
	printf("D");
      } else if (directions[i][j].nogap == HORIZ) {
	printf("H");
      } else if (directions[i][j].nogap == VERT) {
	printf("V");
      } else {
	printf("S");
      }
      printf("|");
      if (directions[i][j].gap2 == DIAG) {
	printf("D");
      } else if (directions[i][j].gap2 == HORIZ) {
	printf("H");
      } else if (directions[i][j].gap2 == VERT) {
	printf("V");
      } else {
	printf("S");
      }
      printf(" ");
    }
    printf("\n");
  }
  printf("\n");
  return;
}




#define QUERY_MAXLENGTH 500
#define GENOMIC_MAXLENGTH 2000


#define T Dynprog_T
struct T {
  int maxlength1;
  int maxlength2;

  struct Int3_T **matrix_ptrs, *matrix_space;
  struct Direction3_T **directions_ptrs, *directions_space;
};

static void
compute_maxlengths (int *maxlength1, int *maxlength2,
		    int maxlookback, int extraquerygap, int maxpeelback,
		    int extramaterial_end, int extramaterial_paired) {
  *maxlength1 = maxlookback + maxpeelback;
  if (*maxlength1 < QUERY_MAXLENGTH) {
    *maxlength1 = QUERY_MAXLENGTH;
  }

  *maxlength2 = *maxlength1 + extraquerygap;
  if (extramaterial_end > extramaterial_paired) {
    *maxlength2 += extramaterial_end;
  } else {
    *maxlength2 += extramaterial_paired;
  }

  if (*maxlength2 < GENOMIC_MAXLENGTH) {
    *maxlength2 = GENOMIC_MAXLENGTH;
  }

  return;
}


T
Dynprog_new (int maxlookback, int extraquerygap, int maxpeelback,
	     int extramaterial_end, int extramaterial_paired) {
  T new = (T) MALLOC(sizeof(*new));
  int maxlength1, maxlength2;

  compute_maxlengths(&maxlength1,&maxlength2,
		     maxlookback,extraquerygap,maxpeelback,
		     extramaterial_end,extramaterial_paired);
  new->maxlength1 = maxlength1;
  new->maxlength2 = maxlength2;

  new->matrix_ptrs = (struct Int3_T **) CALLOC(maxlength1+1,sizeof(struct Int3_T *));
  new->matrix_space = (struct Int3_T *) CALLOC((maxlength1+1)*(maxlength2+1),sizeof(struct Int3_T));
  new->directions_ptrs = (struct Direction3_T **) CALLOC(maxlength1+1,sizeof(struct Direction3_T *));
  new->directions_space = (struct Direction3_T *) CALLOC((maxlength1+1)*(maxlength2+1),sizeof(struct Direction3_T));

  return new;
}


void
Dynprog_free (T *old) {
  if (*old) {
    FREE((*old)->matrix_ptrs);
    FREE((*old)->matrix_space);
    FREE((*old)->directions_ptrs);
    FREE((*old)->directions_space);

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

static int **pairdistance_array[NMISMATCHTYPES];
static bool **consistent_array;

int
Dynprog_pairdistance (int c1, int c2) {
  return pairdistance_array[HIGHQ][c1][c2];
}

static void
permute_cases (int NA1, int NA2, int score) {
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
permute_cases_oneway (int NA1, int NA2, int score) {
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
    pairdistance_array[i] = (int **) CALLOC(128,sizeof(int *));
    pairdistance_array[i][0] = (int *) CALLOC(128*128,sizeof(int));
    ptr = 0;
    for (j = 1; j < 128; j++) {
      ptr += 128;
      pairdistance_array[i][j] = &(pairdistance_array[i][0][ptr]);
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

  return;
}


/************************************************************************/

/* static int max_jump_penalty_lookup; */
static int *jump_penalty_array[NJUMPTYPES];


#if 0
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


static void
jump_penalty_init (int maxlookback, int extraquerygap, int maxpeelback,
		   int extramaterial_end, int extramaterial_paired) {
  int max_lookup, maxlength1, maxlength2, i;
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

  compute_maxlengths(&maxlength1,&maxlength2,
		     maxlookback,extraquerygap,maxpeelback,
		     extramaterial_end,extramaterial_paired);

  if (maxlength1 > maxlength2) {
    max_lookup = maxlength1;
  } else {
    max_lookup = maxlength2;
  }
  if ((remainder = max_lookup % 3) != 0) {
    max_lookup += (3 - remainder);
  }

  /* Set global */
  /* max_jump_penalty_lookup = max_lookup; */

  for (i = 0; i < NJUMPTYPES; i++) {
    jump_penalty_array[i] = (int *) CALLOC(max_lookup+1,sizeof(int));
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

/************************************************************************/

void
Dynprog_init (int maxlookback, int extraquerygap, int maxpeelback,
	      int extramaterial_end, int extramaterial_paired, Mode_T mode) {
  pairdistance_init(mode);
  jump_penalty_init(maxlookback,extraquerygap,maxpeelback,
		    extramaterial_end,extramaterial_paired);
  return;
}

void
Dynprog_term (void) {
  int i;

  for (i = 0; i < NJUMPTYPES; i++) {
    FREE(jump_penalty_array[i]);
  }

  for (i = 0; i < NMISMATCHTYPES; i++) {
    /*
    for (j = 0; j < 128; j++) {
      FREE(pairdistance_array[i][j]);
    }
    */
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
get_genomic_seg (Genomicpos_T genomicpos, int length, Genomicpos_T chroffset,
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


static struct Int3_T **
compute_scores_lookup_fwd (struct Direction3_T ***directions, T this, 
			   char *sequence1, int offset2, int length1, int length2, 
			   Genomicpos_T chroffset, Genomicpos_T chrhigh, bool watsonp,
			   Mismatchtype_T mismatchtype, int open, int extend,
			   int extraband, bool widebandp, bool jump_late_p) {
  struct Int3_T **matrix;
  int r, c, na1, na2;
  char na2_alt;
  int bestscore, score, score_alt;
  Direction_T bestdir;
  int lband, rband, rlo, rhigh;
  int **pairdistance_array_type;
  int penalty;

  pairdistance_array_type = pairdistance_array[mismatchtype];

  if (widebandp == false) {
    /* Just go along main diagonal */
    lband = extraband;
    rband = extraband;
  } else if (length2 >= length1) {
    /* Widen band to right to reach destination */
    rband = length2 - length1 + extraband;
    lband = extraband;
  } else {
    /* Widen band to left to reach destination */
    lband = length1 - length2 + extraband;
    rband = extraband;
  }
  debug(printf("compute_scores_lookup_fwd: "));
  debug(printf("Lengths are %d and %d, so bands are %d on left and %d on right\n",length1,length2,lband,rband));

  matrix = Matrix3_alloc(length1,length2,this->matrix_ptrs,this->matrix_space);
  *directions = Directions3_alloc(length1,length2,this->directions_ptrs,this->directions_space);

  matrix[0][0].nogap = 0;
  (*directions)[0][0].nogap = STOP;
  matrix[0][0].gap1 = matrix[0][0].gap2 = NEG_INFINITY;

  /* Row 0 initialization */
  penalty = open;
  for (c = 1; c <= rband && c <= length2; c++) {
    penalty += extend;
    matrix[0][c].nogap = NEG_INFINITY;

    matrix[0][c].gap1 = penalty;
    (*directions)[0][c].gap1 = HORIZ;

    matrix[0][c].gap2 = NEG_INFINITY;
  }
  (*directions)[0][1].gap1 = STOP;

  /* Column 0 initialization */
  penalty = open;
  for (r = 1; r <= lband && r <= length1; r++) {
    penalty += extend;
    matrix[r][0].nogap = NEG_INFINITY;

    matrix[r][0].gap1 = NEG_INFINITY;

    matrix[r][0].gap2 = penalty;
    (*directions)[r][0].gap2 = VERT;
  }
  (*directions)[1][0].gap2 = STOP;


  for (c = 1; c <= length2; c++) {
    na2 = get_genomic_nt(&na2_alt,offset2+c-1,chroffset,chrhigh,watsonp);
#ifdef EXTRACT_GENOMICSEG
    assert(na2 == sequence2[c-1]); /* na2 = revp ? sequence2[1-c] : sequence2[c-1]; */
#endif

    if ((rlo = c - rband) < 1) {
      rlo = 1;
    } else {
      matrix[rlo-1][c].gap2 = NEG_INFINITY;
      matrix[rlo-1][c].nogap = NEG_INFINITY;
    }

    if ((rhigh = c + lband) > length1) {
      rhigh = length1;
    } else {
      matrix[rhigh][c-1].gap1 = NEG_INFINITY;
      matrix[rhigh][c-1].nogap = NEG_INFINITY;
    }

    
    for (r = rlo; r <= rhigh; r++) {
      na1 = sequence1[r-1]; /* na1 = revp ? sequence1[1-r] : sequence1[r-1]; */

      /* GAP1 */
      bestscore = matrix[r][c-1].nogap + open /* + extend */;
      bestdir = DIAG;

      if ((score = matrix[r][c-1].gap1 /* + extend */) > bestscore || (score == bestscore && jump_late_p)) {
	bestscore = score;
	bestdir = HORIZ;
      }

      matrix[r][c].gap1 = bestscore + extend;
      (*directions)[r][c].gap1 = bestdir;


      /* GAP2 */
      bestscore = matrix[r-1][c].nogap + open /* + extend */;
      bestdir = DIAG;

      if ((score = matrix[r-1][c].gap2 /* + extend */) > bestscore || (score == bestscore && jump_late_p)) {
	bestscore = score;
	bestdir = VERT;
      }

      matrix[r][c].gap2 = bestscore + extend;
      (*directions)[r][c].gap2 = bestdir;


      /* NOGAP */
      /* pairscore = pairdistance_array_type[na1][na2]; */
      bestscore = matrix[r-1][c-1].nogap /* + pairscore */;
      bestdir = DIAG;
      
      if ((score = matrix[r-1][c-1].gap1 /* + pairscore */) > bestscore || (score == bestscore && jump_late_p)) {
	bestscore = score;
	bestdir = HORIZ;
      }

      if ((score = matrix[r-1][c-1].gap2 /* + pairscore */) > bestscore || (score == bestscore && jump_late_p)) {
	bestscore = score;
	bestdir = VERT;
      }

      score = pairdistance_array_type[na1][na2];
      score_alt = pairdistance_array_type[na1][(int) na2_alt];
      if (score >= score_alt) {
	matrix[r][c].nogap = bestscore + score;
      } else {
	matrix[r][c].nogap = bestscore + score_alt;
      }
      (*directions)[r][c].nogap = bestdir;
    }
  }

  /*
  debug2(Matrix_print(P0,length1,length2));
  debug2(Matrix_print(P1,length1,length2));
  debug2(Matrix_print(Q0,length1,length2));
  debug2(Matrix_print(Q1,length1,length2));
  */
  debug2(Matrix3_print(matrix,length1,length2,sequence1,
		       offset2,chroffset,chrhigh,watsonp,
		       /*revp*/false));
  debug2(Directions3_print(*directions,length1,length2,sequence1,/*revp*/false));

  return matrix;
}


static struct Int3_T **
compute_scores_genomicseg_fwd (struct Direction3_T ***directions, T this, 
			       char *sequence1, char *sequence2, char *sequencealt2,
			       int length1, int length2, 
			       Mismatchtype_T mismatchtype, int open, int extend,
			       int extraband, bool widebandp, bool jump_late_p) {
  struct Int3_T **matrix;
  int r, c, na1, na2;
  char na2_alt;
  int bestscore, score, score_alt;
  Direction_T bestdir;
  int lband, rband, rlo, rhigh;
  int **pairdistance_array_type;
  int penalty;

  pairdistance_array_type = pairdistance_array[mismatchtype];

  if (widebandp == false) {
    /* Just go along main diagonal */
    lband = extraband;
    rband = extraband;
  } else if (length2 >= length1) {
    /* Widen band to right to reach destination */
    rband = length2 - length1 + extraband;
    lband = extraband;
  } else {
    /* Widen band to left to reach destination */
    lband = length1 - length2 + extraband;
    rband = extraband;
  }
  debug(printf("compute_scores_genomicseg_fwd: "));
  debug(printf("Lengths are %d and %d, so bands are %d on left and %d on right\n",length1,length2,lband,rband));

  matrix = Matrix3_alloc(length1,length2,this->matrix_ptrs,this->matrix_space);
  *directions = Directions3_alloc(length1,length2,this->directions_ptrs,this->directions_space);

  matrix[0][0].nogap = 0;
  (*directions)[0][0].nogap = STOP;
  matrix[0][0].gap1 = matrix[0][0].gap2 = NEG_INFINITY;

  /* Row 0 initialization */
  penalty = open;
  for (c = 1; c <= rband && c <= length2; c++) {
    penalty += extend;
    matrix[0][c].nogap = NEG_INFINITY;

    matrix[0][c].gap1 = penalty;
    (*directions)[0][c].gap1 = HORIZ;

    matrix[0][c].gap2 = NEG_INFINITY;
  }
  (*directions)[0][1].gap1 = STOP;

  /* Column 0 initialization */
  penalty = open;
  for (r = 1; r <= lband && r <= length1; r++) {
    penalty += extend;
    matrix[r][0].nogap = NEG_INFINITY;

    matrix[r][0].gap1 = NEG_INFINITY;

    matrix[r][0].gap2 = penalty;
    (*directions)[r][0].gap2 = VERT;
  }
  (*directions)[1][0].gap2 = STOP;


  for (c = 1; c <= length2; c++) {
    na2 = sequence2[c-1];
    na2_alt = sequencealt2[c-1];
#ifdef EXTRACT_GENOMICSEG
    assert(na2 == sequence2[c-1]); /* na2 = revp ? sequence2[1-c] : sequence2[c-1]; */
#endif

    if ((rlo = c - rband) < 1) {
      rlo = 1;
    } else {
      matrix[rlo-1][c].gap2 = NEG_INFINITY;
      matrix[rlo-1][c].nogap = NEG_INFINITY;
    }

    if ((rhigh = c + lband) > length1) {
      rhigh = length1;
    } else {
      matrix[rhigh][c-1].gap1 = NEG_INFINITY;
      matrix[rhigh][c-1].nogap = NEG_INFINITY;
    }

    
    for (r = rlo; r <= rhigh; r++) {
      na1 = sequence1[r-1]; /* na1 = revp ? sequence1[1-r] : sequence1[r-1]; */

      /* GAP1 */
      bestscore = matrix[r][c-1].nogap + open /* + extend */;
      bestdir = DIAG;

      if ((score = matrix[r][c-1].gap1 /* + extend */) > bestscore || (score == bestscore && jump_late_p)) {
	bestscore = score;
	bestdir = HORIZ;
      }

      matrix[r][c].gap1 = bestscore + extend;
      (*directions)[r][c].gap1 = bestdir;


      /* GAP2 */
      bestscore = matrix[r-1][c].nogap + open /* + extend */;
      bestdir = DIAG;

      if ((score = matrix[r-1][c].gap2 /* + extend */) > bestscore || (score == bestscore && jump_late_p)) {
	bestscore = score;
	bestdir = VERT;
      }

      matrix[r][c].gap2 = bestscore + extend;
      (*directions)[r][c].gap2 = bestdir;


      /* NOGAP */
      /* pairscore = pairdistance_array_type[na1][na2]; */
      bestscore = matrix[r-1][c-1].nogap /* + pairscore */;
      bestdir = DIAG;
      
      if ((score = matrix[r-1][c-1].gap1 /* + pairscore */) > bestscore || (score == bestscore && jump_late_p)) {
	bestscore = score;
	bestdir = HORIZ;
      }

      if ((score = matrix[r-1][c-1].gap2 /* + pairscore */) > bestscore || (score == bestscore && jump_late_p)) {
	bestscore = score;
	bestdir = VERT;
      }

      score = pairdistance_array_type[na1][na2];
      score_alt = pairdistance_array_type[na1][(int) na2_alt];
      if (score >= score_alt) {
	matrix[r][c].nogap = bestscore + score;
      } else {
	matrix[r][c].nogap = bestscore + score_alt;
      }
      (*directions)[r][c].nogap = bestdir;
    }
  }

  /*
  debug2(Matrix_print(P0,length1,length2));
  debug2(Matrix_print(P1,length1,length2));
  debug2(Matrix_print(Q0,length1,length2));
  debug2(Matrix_print(Q1,length1,length2));
  */
#if 0
  debug2(Matrix3_print(matrix,length1,length2,sequence1,
		       offset2,chroffset,chrhigh,watsonp,
		       /*revp*/false));
  debug2(Directions3_print(*directions,length1,length2,sequence1,/*revp*/false));
#endif

  return matrix;
}


static struct Int3_T **
compute_scores_lookup_rev (struct Direction3_T ***directions, T this, 
			   char *sequence1, int offset2, int length1, int length2, 
			   Genomicpos_T chroffset, Genomicpos_T chrhigh, bool watsonp,
			   Mismatchtype_T mismatchtype, int open, int extend,
			   int extraband, bool widebandp, bool jump_late_p) {
  struct Int3_T **matrix;
  int r, c, na1, na2;
  char na2_alt;
  int bestscore, score, score_alt;
  Direction_T bestdir;
  int lband, rband, rlo, rhigh;
  int **pairdistance_array_type;
  int penalty;

  pairdistance_array_type = pairdistance_array[mismatchtype];

  if (widebandp == false) {
    /* Just go along main diagonal */
    lband = extraband;
    rband = extraband;
  } else if (length2 >= length1) {
    /* Widen band to right to reach destination */
    rband = length2 - length1 + extraband;
    lband = extraband;
  } else {
    /* Widen band to left to reach destination */
    lband = length1 - length2 + extraband;
    rband = extraband;
  }
  debug(printf("compute_scores_lookup_rev: "));
  debug(printf("Lengths are %d and %d, so bands are %d on left and %d on right\n",length1,length2,lband,rband));

  matrix = Matrix3_alloc(length1,length2,this->matrix_ptrs,this->matrix_space);
  *directions = Directions3_alloc(length1,length2,this->directions_ptrs,this->directions_space);

  matrix[0][0].nogap = 0;
  (*directions)[0][0].nogap = STOP;
  matrix[0][0].gap1 = matrix[0][0].gap2 = NEG_INFINITY;

  /* Row 0 initialization */
  penalty = open;
  for (c = 1; c <= rband && c <= length2; c++) {
    penalty += extend;
    matrix[0][c].nogap = NEG_INFINITY;

    matrix[0][c].gap1 = penalty;
    (*directions)[0][c].gap1 = HORIZ;

    matrix[0][c].gap2 = NEG_INFINITY;
  }
  (*directions)[0][1].gap1 = STOP;

  /* Column 0 initialization */
  penalty = open;
  for (r = 1; r <= lband && r <= length1; r++) {
    penalty += extend;
    matrix[r][0].nogap = NEG_INFINITY;

    matrix[r][0].gap1 = NEG_INFINITY;

    matrix[r][0].gap2 = penalty;
    (*directions)[r][0].gap2 = VERT;
  }
  (*directions)[1][0].gap2 = STOP;


  for (c = 1; c <= length2; c++) {
    na2 = get_genomic_nt(&na2_alt,offset2+1-c,chroffset,chrhigh,watsonp);
#ifdef EXTRACT_GENOMICSEG
    assert(na2 == sequence2[1-c]);  /* na2 = revp ? sequence2[1-c] : sequence2[c-1]; */
#endif

    if ((rlo = c - rband) < 1) {
      rlo = 1;
    } else {
      matrix[rlo-1][c].gap2 = NEG_INFINITY;
      matrix[rlo-1][c].nogap = NEG_INFINITY;
    }

    if ((rhigh = c + lband) > length1) {
      rhigh = length1;
    } else {
      matrix[rhigh][c-1].gap1 = NEG_INFINITY;
      matrix[rhigh][c-1].nogap = NEG_INFINITY;
    }

    
    for (r = rlo; r <= rhigh; r++) {
      na1 = sequence1[1-r]; /* na1 = revp ? sequence1[1-r] : sequence1[r-1]; */

      /* GAP1 */
      bestscore = matrix[r][c-1].nogap + open /* + extend */;
      bestdir = DIAG;

      if ((score = matrix[r][c-1].gap1 /* + extend */) > bestscore || (score == bestscore && jump_late_p)) {
	bestscore = score;
	bestdir = HORIZ;
      }

      matrix[r][c].gap1 = bestscore + extend;
      (*directions)[r][c].gap1 = bestdir;


      /* GAP2 */
      bestscore = matrix[r-1][c].nogap + open /* + extend */;
      bestdir = DIAG;

      if ((score = matrix[r-1][c].gap2 /* + extend */) > bestscore || (score == bestscore && jump_late_p)) {
	bestscore = score;
	bestdir = VERT;
      }

      matrix[r][c].gap2 = bestscore + extend;
      (*directions)[r][c].gap2 = bestdir;


      /* NOGAP */
      /* pairscore = pairdistance_array_type[na1][na2]; */
      bestscore = matrix[r-1][c-1].nogap /* + pairscore */;
      bestdir = DIAG;
      
      if ((score = matrix[r-1][c-1].gap1 /* + pairscore */) > bestscore || (score == bestscore && jump_late_p)) {
	bestscore = score;
	bestdir = HORIZ;
      }

      if ((score = matrix[r-1][c-1].gap2 /* + pairscore */) > bestscore || (score == bestscore && jump_late_p)) {
	bestscore = score;
	bestdir = VERT;
      }

      score = pairdistance_array_type[na1][na2];
      score_alt = pairdistance_array_type[na1][(int) na2_alt];
      if (score >= score_alt) {
	matrix[r][c].nogap = bestscore + score;
      } else {
	matrix[r][c].nogap = bestscore + score_alt;
      }
      (*directions)[r][c].nogap = bestdir;
    }
  }


  /*
  debug2(Matrix_print(P0,length1,length2));
  debug2(Matrix_print(P1,length1,length2));
  debug2(Matrix_print(Q0,length1,length2));
  debug2(Matrix_print(Q1,length1,length2));
  */
  debug2(Matrix3_print(matrix,length1,length2,sequence1,
		       offset2,chroffset,chrhigh,watsonp,
		       /*revp*/true));
  debug2(Directions3_print(*directions,length1,length2,sequence1,/*revp*/true));

  return matrix;
}


static struct Int3_T **
compute_scores_genomicseg_rev (struct Direction3_T ***directions, T this, 
			       char *sequence1, char *sequence2, char *sequencealt2,
			       int length1, int length2, 
			       Mismatchtype_T mismatchtype, int open, int extend,
			       int extraband, bool widebandp, bool jump_late_p) {
  struct Int3_T **matrix;
  int r, c, na1, na2, na2_alt;
  int bestscore, score, score_alt;
  Direction_T bestdir;
  int lband, rband, rlo, rhigh;
  int **pairdistance_array_type;
  int penalty;

  pairdistance_array_type = pairdistance_array[mismatchtype];

  if (widebandp == false) {
    /* Just go along main diagonal */
    lband = extraband;
    rband = extraband;
  } else if (length2 >= length1) {
    /* Widen band to right to reach destination */
    rband = length2 - length1 + extraband;
    lband = extraband;
  } else {
    /* Widen band to left to reach destination */
    lband = length1 - length2 + extraband;
    rband = extraband;
  }
  debug(printf("compute_scores_genomicseg_rev: "));
  debug(printf("Lengths are %d and %d, so bands are %d on left and %d on right\n",length1,length2,lband,rband));

  matrix = Matrix3_alloc(length1,length2,this->matrix_ptrs,this->matrix_space);
  *directions = Directions3_alloc(length1,length2,this->directions_ptrs,this->directions_space);

  matrix[0][0].nogap = 0;
  (*directions)[0][0].nogap = STOP;
  matrix[0][0].gap1 = matrix[0][0].gap2 = NEG_INFINITY;

  /* Row 0 initialization */
  penalty = open;
  for (c = 1; c <= rband && c <= length2; c++) {
    penalty += extend;
    matrix[0][c].nogap = NEG_INFINITY;

    matrix[0][c].gap1 = penalty;
    (*directions)[0][c].gap1 = HORIZ;

    matrix[0][c].gap2 = NEG_INFINITY;
  }
  (*directions)[0][1].gap1 = STOP;

  /* Column 0 initialization */
  penalty = open;
  for (r = 1; r <= lband && r <= length1; r++) {
    penalty += extend;
    matrix[r][0].nogap = NEG_INFINITY;

    matrix[r][0].gap1 = NEG_INFINITY;

    matrix[r][0].gap2 = penalty;
    (*directions)[r][0].gap2 = VERT;
  }
  (*directions)[1][0].gap2 = STOP;


  for (c = 1; c <= length2; c++) {
    na2 = sequence2[1-c];
    na2_alt = sequencealt2[1-c];
#ifdef EXTRACT_GENOMICSEG
    assert(na2 == sequence2[1-c]);  /* na2 = revp ? sequence2[1-c] : sequence2[c-1]; */
#endif

    if ((rlo = c - rband) < 1) {
      rlo = 1;
    } else {
      matrix[rlo-1][c].gap2 = NEG_INFINITY;
      matrix[rlo-1][c].nogap = NEG_INFINITY;
    }

    if ((rhigh = c + lband) > length1) {
      rhigh = length1;
    } else {
      matrix[rhigh][c-1].gap1 = NEG_INFINITY;
      matrix[rhigh][c-1].nogap = NEG_INFINITY;
    }

    
    for (r = rlo; r <= rhigh; r++) {
      na1 = sequence1[1-r]; /* na1 = revp ? sequence1[1-r] : sequence1[r-1]; */

      /* GAP1 */
      bestscore = matrix[r][c-1].nogap + open /* + extend */;
      bestdir = DIAG;

      if ((score = matrix[r][c-1].gap1 /* + extend */) > bestscore || (score == bestscore && jump_late_p)) {
	bestscore = score;
	bestdir = HORIZ;
      }

      matrix[r][c].gap1 = bestscore + extend;
      (*directions)[r][c].gap1 = bestdir;


      /* GAP2 */
      bestscore = matrix[r-1][c].nogap + open /* + extend */;
      bestdir = DIAG;

      if ((score = matrix[r-1][c].gap2 /* + extend */) > bestscore || (score == bestscore && jump_late_p)) {
	bestscore = score;
	bestdir = VERT;
      }

      matrix[r][c].gap2 = bestscore + extend;
      (*directions)[r][c].gap2 = bestdir;


      /* NOGAP */
      /* pairscore = pairdistance_array_type[na1][na2]; */
      bestscore = matrix[r-1][c-1].nogap /* + pairscore */;
      bestdir = DIAG;
      
      if ((score = matrix[r-1][c-1].gap1 /* + pairscore */) > bestscore || (score == bestscore && jump_late_p)) {
	bestscore = score;
	bestdir = HORIZ;
      }

      if ((score = matrix[r-1][c-1].gap2 /* + pairscore */) > bestscore || (score == bestscore && jump_late_p)) {
	bestscore = score;
	bestdir = VERT;
      }

      score = pairdistance_array_type[na1][na2];
      score_alt = pairdistance_array_type[na1][(int) na2_alt];
      if (score >= score_alt) {
	matrix[r][c].nogap = bestscore + score;
      } else {
	matrix[r][c].nogap = bestscore + score_alt;
      }
      (*directions)[r][c].nogap = bestdir;
    }
  }


  /*
  debug2(Matrix_print(P0,length1,length2));
  debug2(Matrix_print(P1,length1,length2));
  debug2(Matrix_print(Q0,length1,length2));
  debug2(Matrix_print(Q1,length1,length2));
  */
#if 0
  debug2(Matrix3_print(matrix,length1,length2,sequence1,
		       offset2,chroffset,chrhigh,watsonp,
		       /*revp*/true));
  debug2(Directions3_print(*directions,length1,length2,sequence1,/*revp*/true));
#endif

  return matrix;
}


/* Puts sequence1 (genomic) in outer loop for cDNA gaps */
static struct Int3_T **
compute_scores_lookup_fwd_12 (struct Direction3_T ***directions, T this, 
			      char *sequence2, int offset1, int length1, int length2, 
			      Genomicpos_T chroffset, Genomicpos_T chrhigh, bool watsonp,
			      Mismatchtype_T mismatchtype, int open, int extend,
			      int extraband, bool widebandp, bool jump_late_p) {
  struct Int3_T **matrix;
  int r, c, na1, na2;
  char na1_alt;
  int bestscore, score, score_alt;
  Direction_T bestdir;
  int lband, rband, clo, chigh;
  int **pairdistance_array_type;
  int penalty;

  pairdistance_array_type = pairdistance_array[mismatchtype];

  if (widebandp == false) {
    /* Just go along main diagonal */
    lband = extraband;
    rband = extraband;
  } else if (length2 >= length1) {
    /* Widen band to right to reach destination */
    rband = length2 - length1 + extraband;
    lband = extraband;
  } else {
    /* Widen band to left to reach destination */
    lband = length1 - length2 + extraband;
    rband = extraband;
  }
  debug(printf("compute_scores_lookup_fwd_12: "));
  debug(printf("Lengths are %d and %d, so bands are %d on left and %d on right\n",length1,length2,lband,rband));

  matrix = Matrix3_alloc(length1,length2,this->matrix_ptrs,this->matrix_space);
  *directions = Directions3_alloc(length1,length2,this->directions_ptrs,this->directions_space);

  matrix[0][0].nogap = 0;
  (*directions)[0][0].nogap = STOP;
  matrix[0][0].gap1 = matrix[0][0].gap2 = NEG_INFINITY;

  /* Row 0 initialization */
  penalty = open;
  for (c = 1; c <= rband && c <= length2; c++) {
    penalty += extend;
    matrix[0][c].nogap = NEG_INFINITY;

    matrix[0][c].gap1 = penalty;
    (*directions)[0][c].gap1 = HORIZ;

    matrix[0][c].gap2 = NEG_INFINITY;
  }
  (*directions)[0][1].gap1 = STOP;

  /* Column 0 initialization */
  penalty = open;
  for (r = 1; r <= lband && r <= length1; r++) {
    penalty += extend;
    matrix[r][0].nogap = NEG_INFINITY;

    matrix[r][0].gap1 = NEG_INFINITY;

    matrix[r][0].gap2 = penalty;
    (*directions)[r][0].gap2 = VERT;
  }
  (*directions)[1][0].gap2 = STOP;


  for (r = 1; r <= length1; r++) {
    na1 = get_genomic_nt(&na1_alt,offset1+r-1,chroffset,chrhigh,watsonp);
#ifdef EXTRACT_GENOMICSEG
    assert(na1 == sequence1[r-1]);  /* na1 = revp ? sequence1[1-r] : sequence1[r-1]; */
#endif
    
    if ((clo = r - lband) < 1) {
      clo = 1;
    } else {
      matrix[r][clo-1].gap1 = NEG_INFINITY;
      matrix[r][clo-1].nogap = NEG_INFINITY;
    }

    if ((chigh = r + rband) > length2) {
      chigh = length2;
    } else {
      matrix[r-1][chigh].gap2 = NEG_INFINITY;
      matrix[r-1][chigh].nogap = NEG_INFINITY;
    }

    
    for (c = clo; c <= chigh; c++) {
      na2 = sequence2[c-1]; /* na2 = revp ? sequence2[1-c] : sequence2[c-1]; */

      /* GAP1 */
      bestscore = matrix[r][c-1].nogap + open /* + extend */;
      bestdir = DIAG;

      if ((score = matrix[r][c-1].gap1 /* + extend */) > bestscore || (score == bestscore && jump_late_p)) {
	bestscore = score;
	bestdir = HORIZ;
      }

      matrix[r][c].gap1 = bestscore + extend;
      (*directions)[r][c].gap1 = bestdir;


      /* GAP2 */
      bestscore = matrix[r-1][c].nogap + open /* + extend */;
      bestdir = DIAG;

      if ((score = matrix[r-1][c].gap2 /* + extend */) > bestscore || (score == bestscore && jump_late_p)) {
	bestscore = score;
	bestdir = VERT;
      }

      matrix[r][c].gap2 = bestscore + extend;
      (*directions)[r][c].gap2 = bestdir;


      /* NOGAP */
      /* pairscore = pairdistance_array_type[na2][na1]; -- Note swap */
      bestscore = matrix[r-1][c-1].nogap /* + pairscore */;
      bestdir = DIAG;
      
      if ((score = matrix[r-1][c-1].gap1 /* + pairscore */) > bestscore || (score == bestscore && jump_late_p)) {
	bestscore = score;
	bestdir = HORIZ;
      }

      if ((score = matrix[r-1][c-1].gap2 /* + pairscore */) > bestscore || (score == bestscore && jump_late_p)) {
	bestscore = score;
	bestdir = VERT;
      }

      score = pairdistance_array_type[na2][na1]; /* Note swap */
      score_alt = pairdistance_array_type[na2][(int) na1_alt]; /* Note swap */
      if (score >= score_alt) {
	matrix[r][c].nogap = bestscore + score;
      } else {
	matrix[r][c].nogap = bestscore + score_alt;
      }
      (*directions)[r][c].nogap = bestdir;
    }
  }

  /*
  debug2(Matrix_print(P0,length1,length2));
  debug2(Matrix_print(P1,length1,length2));
  debug2(Matrix_print(Q0,length1,length2));
  debug2(Matrix_print(Q1,length1,length2));
  */
  debug2(Matrix3_print(matrix,length2,length1,sequence2,
		       offset1,chroffset,chrhigh,watsonp,
		       /*revp*/false));
  debug2(Directions3_print(*directions,length2,length1,sequence2,/*revp*/false));

  return matrix;
}


/* Puts sequence1 (genomic) in outer loop for cDNA gaps */
static struct Int3_T **
compute_scores_lookup_rev_12 (struct Direction3_T ***directions, T this, 
			      char *sequence2, int offset1, int length1, int length2, 
			      Genomicpos_T chroffset, Genomicpos_T chrhigh, bool watsonp,
			      Mismatchtype_T mismatchtype, int open, int extend,
			      int extraband, bool widebandp, bool jump_late_p) {
  struct Int3_T **matrix;
  int r, c, na1, na2;
  char na1_alt;
  int bestscore, score, score_alt;
  Direction_T bestdir;
  int lband, rband, clo, chigh;
  int **pairdistance_array_type;
  int penalty;

  pairdistance_array_type = pairdistance_array[mismatchtype];

  if (widebandp == false) {
    /* Just go along main diagonal */
    lband = extraband;
    rband = extraband;
  } else if (length2 >= length1) {
    /* Widen band to right to reach destination */
    rband = length2 - length1 + extraband;
    lband = extraband;
  } else {
    /* Widen band to left to reach destination */
    lband = length1 - length2 + extraband;
    rband = extraband;
  }
  debug(printf("compute_scores_lookup_rev_12: "));
  debug(printf("Lengths are %d and %d, so bands are %d on left and %d on right\n",length1,length2,lband,rband));

  matrix = Matrix3_alloc(length1,length2,this->matrix_ptrs,this->matrix_space);
  *directions = Directions3_alloc(length1,length2,this->directions_ptrs,this->directions_space);

  matrix[0][0].nogap = 0;
  (*directions)[0][0].nogap = STOP;
  matrix[0][0].gap1 = matrix[0][0].gap2 = NEG_INFINITY;

  /* Row 0 initialization */
  penalty = open;
  for (c = 1; c <= rband && c <= length2; c++) {
    penalty += extend;
    matrix[0][c].nogap = NEG_INFINITY;

    matrix[0][c].gap1 = penalty;
    (*directions)[0][c].gap1 = HORIZ;

    matrix[0][c].gap2 = NEG_INFINITY;
  }
  (*directions)[0][1].gap1 = STOP;

  /* Column 0 initialization */
  penalty = open;
  for (r = 1; r <= lband && r <= length1; r++) {
    penalty += extend;
    matrix[r][0].nogap = NEG_INFINITY;

    matrix[r][0].gap1 = NEG_INFINITY;

    matrix[r][0].gap2 = penalty;
    (*directions)[r][0].gap2 = VERT;
  }
  (*directions)[1][0].gap2 = STOP;


  for (r = 1; r <= length1; r++) {
    na1 = get_genomic_nt(&na1_alt,offset1+1-r,chroffset,chrhigh,watsonp);
#ifdef EXTRACT_GENOMICSEG
    assert(na1 == sequence1[1-r]); /* na1 = revp ? sequence1[1-r] : sequence1[r-1]; */
#endif

    if ((clo = r - lband) < 1) {
      clo = 1;
    } else {
      matrix[r][clo-1].gap1 = NEG_INFINITY;
      matrix[r][clo-1].nogap = NEG_INFINITY;
    }

    if ((chigh = r + rband) > length2) {
      chigh = length2;
    } else {
      matrix[r-1][chigh].gap2 = NEG_INFINITY;
      matrix[r-1][chigh].nogap = NEG_INFINITY;
    }

    
    for (c = clo; c <= chigh; c++) {
      na2 = sequence2[1-c]; /* na2 = revp ? sequence2[1-c] : sequence2[c-1]; */

      /* GAP1 */
      bestscore = matrix[r][c-1].nogap + open /* + extend */;
      bestdir = DIAG;

      if ((score = matrix[r][c-1].gap1 /* + extend */) > bestscore || (score == bestscore && jump_late_p)) {
	bestscore = score;
	bestdir = HORIZ;
      }

      matrix[r][c].gap1 = bestscore + extend;
      (*directions)[r][c].gap1 = bestdir;


      /* GAP2 */
      bestscore = matrix[r-1][c].nogap + open /* + extend */;
      bestdir = DIAG;

      if ((score = matrix[r-1][c].gap2 /* + extend */) > bestscore || (score == bestscore && jump_late_p)) {
	bestscore = score;
	bestdir = VERT;
      }

      matrix[r][c].gap2 = bestscore + extend;
      (*directions)[r][c].gap2 = bestdir;


      /* NOGAP */
      /* pairscore = pairdistance_array_type[na2][na1]; -- Note swap */
      bestscore = matrix[r-1][c-1].nogap /* + pairscore */;
      bestdir = DIAG;
      
      if ((score = matrix[r-1][c-1].gap1 /* + pairscore */) > bestscore || (score == bestscore && jump_late_p)) {
	bestscore = score;
	bestdir = HORIZ;
      }

      if ((score = matrix[r-1][c-1].gap2 /* + pairscore */) > bestscore || (score == bestscore && jump_late_p)) {
	bestscore = score;
	bestdir = VERT;
      }

      score = pairdistance_array_type[na2][na1]; /* Note swap */
      score_alt = pairdistance_array_type[na2][(int) na1_alt]; /* Note swap */
      if (score >= score_alt) {
	matrix[r][c].nogap = bestscore + score;
      } else {
	matrix[r][c].nogap = bestscore + score_alt;
      }
      (*directions)[r][c].nogap = bestdir;
    }
  }


  /*
  debug2(Matrix_print(P0,length1,length2));
  debug2(Matrix_print(P1,length1,length2));
  debug2(Matrix_print(Q0,length1,length2));
  debug2(Matrix_print(Q1,length1,length2));
  */
  debug2(Matrix3_print(matrix,length2,length1,sequence2,
		       offset1,chroffset,chrhigh,watsonp,
		       /*revp*/true));
  debug2(Directions3_print(*directions,length2,length1,sequence2,/*revp*/true));

  return matrix;
}



#if 0
static int **
compute_scores (Direction_T ***directions, int ***jump, T this, 
		char *sequence1, char *sequence2, bool revp, 
		int length1, int length2, Mismatchtype_T mismatchtype,
		int open, int extend, bool init_jump_penalty_p,
		int extraband, bool widebandp, bool onesidegapp, bool jump_late_p) {
  int **matrix;
  int r, c, r1, c1, na1, na2;
  int bestscore, score, bestjump, j;
  Direction_T bestdir;
  int lband, rband, clo, chigh, rlo;
  int **pairdistance_array_type;

  pairdistance_array_type = pairdistance_array[mismatchtype];

  if (widebandp == false) {
    /* Just go along main diagonal */
    lband = extraband;
    rband = extraband;
  } else if (length2 >= length1) {
    /* Widen band to right to reach destination */
    rband = length2 - length1 + extraband;
    lband = extraband;
  } else {
    /* Widen band to left to reach destination */
    lband = length1 - length2 + extraband;
    rband = extraband;
  }

  matrix = Matrix_alloc(length1,length2,this->matrix_ptrs,this->matrix_space);
  *directions = Directions_alloc(length1,length2,this->directions_ptrs,this->directions_space);
  *jump = Matrix_alloc(length1,length2,this->jump_ptrs,this->jump_space);

  matrix[0][0] = 0;
  (*directions)[0][0] = STOP;
  (*jump)[0][0] = 0;

  /* Row 0 initialization */
  if (init_jump_penalty_p == true) {
    for (c = 1; c <= rband && c <= length2; c++) {
      matrix[0][c] = jump_penalty(c,open,extend);
      (*directions)[0][c] = HORIZ;
      (*jump)[0][c] = c;
    }
  } else {
    for (c = 1; c <= rband && c <= length2; c++) {
      matrix[0][c] = SINGLE_OPEN_HIGHQ;	/* Needs to be less than 0 to prevent gap and then a single match */
      (*directions)[0][c] = HORIZ;
      (*jump)[0][c] = c;
    }
  }

  /* Column 0 initialization */
  if (init_jump_penalty_p == true) {
    for (r = 1; r <= lband && r <= length1; r++) {
      matrix[r][0] = jump_penalty(r,open,extend);
      (*directions)[r][0] = VERT;
      (*jump)[r][0] = r;
    }
  } else {
    for (r = 1; r <= lband && r <= length1; r++) {
      matrix[r][0] = SINGLE_OPEN_HIGHQ;	/* Needs to be less than 0 to prevent gap and then a single match */
      (*directions)[r][0] = VERT;
      (*jump)[r][0] = r;
    }
  }

  for (r = 1; r <= length1; r++) {
    na1 = revp ? sequence1[1-r] : sequence1[r-1];

    if ((clo = r - lband) < 1) {
      clo = 1;
    }
    if ((chigh = r + rband) > length2) {
      chigh = length2;
    }
    for (c = clo; c <= chigh; c++) {
      na2 = revp ? sequence2[1-c] : sequence2[c-1];

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
	if ((rlo = c+c-rband-r) < 1) {
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

  /*
  debug2(Matrix_print(P0,length1,length2));
  debug2(Matrix_print(P1,length1,length2));
  debug2(Matrix_print(Q0,length1,length2));
  debug2(Matrix_print(Q1,length1,length2));
  */
  debug2(Matrix_print(matrix,length1,length2,sequence1,sequence2,revp));
  debug2(Directions_print(*directions,*jump,length1,length2,
			  sequence1,sequence2,revp));

  return matrix;
}
#endif


static void
find_best_endpoint (int *finalscore, int *bestr, int *bestc, struct Int3_T **matrix, 
		    int length1, int length2, int extraband_end_or_paired,
		    bool jump_late_p) {
  int bestscore = 0;
  int r, c;
  int rband, lband, clo, chigh;
  /* No need for loffset or cmid because they apply only for cdnaend
     == FIVE, which doesn't require searching */

  *bestr = *bestc = 0;

  /* Just go along main diagonal */
  rband = extraband_end_or_paired;
  lband = extraband_end_or_paired;

  if (jump_late_p == false) {
    /* use > bestscore */
    for (r = 1; r <= length1; r++) {
      if ((clo = r - lband) < 1) {
	clo = 1;
      }
      if ((chigh = r + rband) > length2) {
	chigh = length2;
      }
      for (c = clo; c <= chigh; c++) {
	if (matrix[r][c].nogap > bestscore) {
	  *bestr = r;
	  *bestc = c;
	  bestscore = matrix[r][c].nogap;
	}
      }
    }
  } else {
    /* use >= bestscore */
    for (r = 1; r <= length1; r++) {
      if ((clo = r - lband) < 1) {
	clo = 1;
      }
      if ((chigh = r + rband) > length2) {
	chigh = length2;
      }
      for (c = clo; c <= chigh; c++) {
	if (matrix[r][c].nogap >= bestscore) {
	  *bestr = r;
	  *bestc = c;
	  bestscore = matrix[r][c].nogap;
	}
      }
    }
  }


  *finalscore = bestscore;
  return;
}


static void
find_best_endpoint_to_queryend_indels (int *finalscore, int *bestr, int *bestc, struct Int3_T **matrix, 
				       int length1, int length2, int extraband_end_or_paired,
				       bool jump_late_p) {
  int bestscore = NEG_INFINITY;
  int r, c;
  int rband, lband, clo, chigh;
  /* No need for loffset or cmid because they apply only for cdnaend
     == FIVE, which doesn't require searching */

  if (length2 >= length1) {
    /* Widen band to right to reach destination */
    rband = length2 - length1 + extraband_end_or_paired;
    lband = extraband_end_or_paired;
  } else {
    /* Widen band to left to reach destination */
    lband = length1 - length2 + extraband_end_or_paired;
    rband = extraband_end_or_paired;
  }

  *bestr = r = length1;
  *bestc = 0;

  if (jump_late_p == false) {
    if ((clo = r - lband) < 1) {
      clo = 1;
    }
    if ((chigh = r + rband) > length2) {
      chigh = length2;
    }
    if (clo <= chigh) {
      for (c = clo; c <= chigh; c++) {
	/* use > bestscore */
	if (matrix[r][c].nogap > bestscore) {
	  *bestr = r;
	  *bestc = c;
	  bestscore = matrix[r][c].nogap;
	}
      }
    }

  } else {
    if ((clo = r - lband) < 1) {
      clo = 1;
    }
    if ((chigh = r + rband) > length2) {
      chigh = length2;
    }
    if (clo <= chigh) {
      for (c = clo; c <= chigh; c++) {
	/* use >= bestscore */
	if (matrix[r][c].nogap >= bestscore) {
	  *bestr = r;
	  *bestc = c;
	  bestscore = matrix[r][c].nogap;
	}
      }
    }
  }

  *finalscore = bestscore;
  return;
}


static void
find_best_endpoint_to_queryend_nogaps (int *bestr, int *bestc, int length1, int length2) {
  if (length2 < length1) {
    *bestr = length2;
    *bestc = length2;
  } else {
    *bestr = length1;
    *bestc = length1;
  }

  return;
}


static List_T
add_queryskip (List_T pairs, int r, int c, int dist, char *querysequence,
	       int queryoffset, int genomeoffset, Pairpool_T pairpool, bool revp, bool cdna_gap_p, int dynprogindex) {
  int j;
  char c1;
  int querycoord, genomecoord, step;

  if (cdna_gap_p == false) {
    querycoord = r-1;
    genomecoord = c-1;
  } else {
    querycoord = c-1;
    genomecoord = r-1;
  }

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
    debug(
	  if (cdna_gap_p == false) {
	    r--;
	  } else {
	    c--;
	  }
	  );
    querycoord += step;
  }

  return pairs;
}


static List_T
add_genomeskip (bool *add_dashes_p, List_T pairs, int r, int c, int dist,
		char *genomesequence, char *genomesequenceuc,
		int queryoffset, int genomeoffset, Pairpool_T pairpool, bool revp,
		Genomicpos_T chroffset, Genomicpos_T chrhigh,
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

  if (dist < MICROINTRON_LENGTH) {
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
#if 0
    pairs = Pairpool_push_gapholder(pairs,pairpool,/*queryjump*/0,/*genomejump*/dist,/*knownp*/false);
#endif
    pairs = Pairpool_push_gapholder(pairs,pairpool,/*queryjump*/UNKNOWNJUMP,/*genomejump*/UNKNOWNJUMP,/*knownp*/false);
#endif
  }
  
  return pairs;
}


static List_T
add_genomeskip_cdna (bool *add_dashes_p, List_T pairs, int r, int c, int dist,
		     int queryoffset, int genomeoffset, Pairpool_T pairpool, bool revp,
		     Genomicpos_T chroffset, Genomicpos_T chrhigh,
		     int cdna_direction, bool watsonp, int dynprogindex) {
  int j;
  char left1, left2, right2, right1, left1_alt, left2_alt, right2_alt, right1_alt, c2, c2_alt;
  int introntype;
  int querycoord, leftgenomecoord, rightgenomecoord, genomecoord, temp, step;

  querycoord = c-1;
  leftgenomecoord = r-dist;
  rightgenomecoord = r-1;

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

  if (dist < MICROINTRON_LENGTH) {
    *add_dashes_p = true;
  } else {
    /* Check for intron */
    left1 = get_genomic_nt(&left1_alt,genomeoffset+leftgenomecoord,chroffset,chrhigh,watsonp);
    left2 = get_genomic_nt(&left2_alt,genomeoffset+leftgenomecoord+1,chroffset,chrhigh,watsonp);
    right2 = get_genomic_nt(&right2_alt,genomeoffset+rightgenomecoord-1,chroffset,chrhigh,watsonp);
    right1 = get_genomic_nt(&right1_alt,genomeoffset+rightgenomecoord,chroffset,chrhigh,watsonp);
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

  if (*add_dashes_p == true) {
    if (revp == true) {
      genomecoord = leftgenomecoord;
    } else {
      genomecoord = rightgenomecoord;
    }
    for (j = 0; j < dist; j++) {
      c2 = get_genomic_nt(&c2_alt,genomeoffset+genomecoord,chroffset,chrhigh,watsonp);
#ifdef EXTRACT_GENOMICSEG
      assert(c2 == genomesequence[genomecoord]);
#endif

      debug(printf("Pushing %d,%d [%d,%d] (-,%c),",r,c,queryoffset+querycoord,genomeoffset+genomecoord,c2));
      pairs = Pairpool_push(pairs,pairpool,queryoffset+querycoord,genomeoffset+genomecoord,
			    ' ',INDEL_COMP,c2,c2_alt,dynprogindex);
      debug(r--);
      genomecoord += step;
    }
  } else {
    debug(printf("Large gap %c%c..%c%c.  Adding gap of type %d.\n",
		 left1,left2,right2,right1,introntype));
#ifndef NOGAPHOLDER
#if 0
    pairs = Pairpool_push_gapholder(pairs,pairpool,/*queryjump*/0,/*genomejump*/dist,/*knownp*/false);
#endif
    pairs = Pairpool_push_gapholder(pairs,pairpool,/*queryjump*/UNKNOWNJUMP,/*genomejump*/UNKNOWNJUMP,/*knownp*/false);
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
		      Genomicpos_T chroffset, Genomicpos_T chrhigh,
		      int cdna_direction, bool watsonp, int dynprogindex, bool use_genomicseg_p) {
  char c1, c2, c2_alt;
  int dist;
  bool add_dashes_p;
  int querycoord, genomecoord;

  debug(printf("Starting traceback at r=%d,c=%d (offset1=%d, offset2=%d)\n",r,c,queryoffset,genomeoffset));

  while (directions[r][c].nogap != STOP) {
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

    } else if (querysequenceuc[querycoord] == c2) {
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

    if (directions[r][c].nogap == DIAG) {
      r--; c--;

    } else if (directions[r][c].nogap == HORIZ) {
      dist = 1;
      r--; c--;
      while (directions[r][c].gap1 == HORIZ) {
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

    } else /* if (directions[r][c].nogap == VERT) */ {
      dist = 1;
      r--; c--;
      while (directions[r][c].gap2 == VERT) {
	dist++;
	r--;
      }
      r--;

      debug(printf("V%d: ",dist));
      pairs = add_queryskip(pairs,r+dist,c,dist,querysequence,
			    queryoffset,genomeoffset,pairpool,revp,/*cdna_gap_p*/false,
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
	   struct Direction3_T **directions, int r, int c, 
	   char *querysequence, char *querysequenceuc,
	   int queryoffset, int genomeoffset, Pairpool_T pairpool, bool revp,
	   Genomicpos_T chroffset, Genomicpos_T chrhigh,
	   int cdna_direction, bool watsonp, int dynprogindex) {
  char c1, c2, c2_alt;
  int dist;
  bool add_dashes_p;
  int querycoord, genomecoord;

  debug(printf("Starting traceback at r=%d,c=%d (offset1=%d, offset2=%d)\n",r,c,queryoffset,genomeoffset));

  while (directions[r][c].nogap != STOP) {
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

    } else if (querysequenceuc[querycoord] == c2) {
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

    if (directions[r][c].nogap == DIAG) {
      r--; c--;

    } else if (directions[r][c].nogap == HORIZ) {
      dist = 1;
      r--; c--;
      while (directions[r][c].gap1 == HORIZ) {
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

    } else /* if (directions[r][c].nogap == VERT) */ {
      dist = 1;
      r--; c--;
      while (directions[r][c].gap2 == VERT) {
	dist++;
	r--;
      }
      r--;

      debug(printf("V%d: ",dist));
      pairs = add_queryskip(pairs,r+dist,c,dist,querysequence,
			    queryoffset,genomeoffset,pairpool,revp,/*cdna_gap_p*/false,
			    dynprogindex);
      *nopens += 1;
      *nindels += dist;
      debug(printf("\n"));
    }

  }
  return pairs;
}


static List_T
traceback_cdna (List_T pairs, int *nmatches, int *nmismatches, int *nopens, int *nindels,
		struct Direction3_T **directions, int r, int c, 
		char *querysequence, char *querysequenceuc,
		int queryoffset, int genomeoffset, Pairpool_T pairpool, bool revp,
		Genomicpos_T chroffset, Genomicpos_T chrhigh,
		int cdna_direction, bool watsonp, int dynprogindex) {
  char c1, c2, c2_alt;
  int dist;
  bool add_dashes_p;
  int querycoord, genomecoord;

  debug(printf("Starting traceback_cdna at r=%d,c=%d (offset1=%d, offset2=%d)\n",r,c,queryoffset,genomeoffset));

  while (directions[r][c].nogap != STOP) {
    querycoord = c-1;
    genomecoord = r-1;
    if (revp == true) {
      querycoord = -querycoord;
      genomecoord = -genomecoord;
    }

    c1 = querysequence[querycoord];
    c2 = get_genomic_nt(&c2_alt,genomeoffset+genomecoord,chroffset,chrhigh,watsonp);
#ifdef EXTRACT_GENOMICSEG
    assert(c2 == genomesequence[genomecoord]);
#endif

    /* Note swap on consistent_array */
    if (querysequenceuc[querycoord] == c2) {
      debug(printf("D1: Pushing %d,%d [%d,%d] (%c,%c) - match\n",
		   r,c,queryoffset+querycoord,genomeoffset+genomecoord,c1,c2));
      *nmatches += 1;
      pairs = Pairpool_push(pairs,pairpool,queryoffset+querycoord,genomeoffset+genomecoord,
			    c1,DYNPROG_MATCH_COMP,c2,c2_alt,dynprogindex);

    } else if (consistent_array[(int) c2][(int) c1] == true) {
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

    if (directions[r][c].nogap == DIAG) {
      r--; c--;

    } else if (directions[r][c].nogap == HORIZ) {
      dist = 1;
      r--; c--;
      while (directions[r][c].gap1 == HORIZ) {
	dist++;
	c--;
      }
      c--;

      debug(printf("H%d: ",dist));
      pairs = add_queryskip(pairs,r,c+dist,dist,querysequence,
			    queryoffset,genomeoffset,pairpool,revp,/*cdna_gap_p*/true,
			    dynprogindex);
      *nopens += 1;
      *nindels += dist;
      debug(printf("\n"));

    } else /* if (directions[r][c].nogap == VERT) */ {
      dist = 1;
      r--; c--;
      while (directions[r][c].gap2 == VERT) {
	dist++;
	r--;
      }
      r--;

      debug(printf("V%d: ",dist));
      pairs = add_genomeskip_cdna(&add_dashes_p,pairs,r+dist,c,dist,
				  queryoffset,genomeoffset,pairpool,revp,
				  chroffset,chrhigh,
				  cdna_direction,watsonp,dynprogindex);
      if (add_dashes_p == true) {
	*nopens += 1;
	*nindels += dist;
      }
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
		  int queryoffset, int genomeoffset, Pairpool_T pairpool, 
		  Genomicpos_T chroffset, Genomicpos_T chrhigh,
		  bool revp, bool watsonp, int dynprogindex) {
  char c1, c2, c2_alt;
  int querycoord, genomecoord;

  debug11(printf("Starting traceback_nogaps at r=%d,c=%d (offset1=%d, offset2=%d), revp %d\n",
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
    c2 = get_genomic_nt(&c2_alt,genomeoffset+genomecoord,chroffset,chrhigh,watsonp);
#ifdef EXTRACT_GENOMICSEG
    debug8(printf("genome sequence at %d is %c\n",genomecoord,genomesequence[genomecoord]));
    assert(c2 == genomesequence[genomecoord]);
#endif

    if (c2 == '*') {
      /* Don't push pairs past end of chromosome */

    } else if (querysequenceuc[querycoord] == c2) {
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


static List_T
traceback_local (List_T pairs, int *nmatches, int *nmismatches, int *nopens, int *nindels,
		 struct Direction3_T **directions, int *r, int *c, int endc,
		 char *querysequence, char *querysequenceuc,
		 char *genomesequence, char *genomesequenceuc, char *genomesequencealt,
		 int queryoffset, int genomeoffset, Pairpool_T pairpool, bool revp,
		 Genomicpos_T chroffset, Genomicpos_T chrhigh,
		 int cdna_direction, bool watsonp, int dynprogindex) {
  char c1, c2, c2_alt;
  int dist;
  bool add_dashes_p;
  int querycoord, genomecoord;

  debug(printf("Starting traceback_local at r=%d,c=%d (offset1=%d, offset2=%d)\n",*r,*c,queryoffset,genomeoffset));

  /* We care only only about genomic coordinate c */
  while (*c > endc) {
    querycoord = (*r)-1;
    genomecoord = (*c)-1;
    if (revp == true) {
      querycoord = -querycoord;
      genomecoord = -genomecoord;
    }

    c1 = querysequence[querycoord];
    c2 = genomesequence[genomecoord];
    c2_alt = genomesequencealt[genomecoord];

    if (querysequenceuc[querycoord] == c2) {
      debug(printf("D1: Pushing %d,%d [%d,%d] (%c,%c) - match\n",
		   *r,*c,queryoffset+querycoord,genomeoffset+genomecoord,c1,c2));
      *nmatches += 1;
      pairs = Pairpool_push(pairs,pairpool,queryoffset+querycoord,genomeoffset+genomecoord,
			    c1,DYNPROG_MATCH_COMP,c2,c2_alt,dynprogindex);

    } else if (consistent_array[(int) c1][(int) c2] == true) {
      debug(printf("D1: Pushing %d,%d [%d,%d] (%c,%c) - ambiguous\n",
		   *r,*c,queryoffset+querycoord,genomeoffset+genomecoord,c1,c2));
      *nmatches += 1;
      pairs = Pairpool_push(pairs,pairpool,queryoffset+querycoord,genomeoffset+genomecoord,
			    c1,AMBIGUOUS_COMP,c2,c2_alt,dynprogindex);

    } else {
      debug(printf("D1: Pushing %d,%d [%d,%d] (%c,%c) - mismatch\n",
		   *r,*c,queryoffset+querycoord,genomeoffset+genomecoord,c1,c2));
      *nmismatches += 1;
      pairs = Pairpool_push(pairs,pairpool,queryoffset+querycoord,genomeoffset+genomecoord,
			    c1,MISMATCH_COMP,c2,c2_alt,dynprogindex);
    }

    if (directions[*r][*c].nogap == DIAG) {
      (*r)--; (*c)--;

    } else if (directions[*r][*c].nogap == HORIZ) {
      dist = 1;
      (*r)--; (*c)--;
      while (directions[*r][*c].gap1 == HORIZ) {
	dist++;
	(*c)--;
      }
      (*c)--;

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

    } else /* if (directions[*r][*c].nogap == VERT) */ {
      dist = 1;
      (*r)--; (*c)--;
      while (directions[*r][*c].gap2 == VERT) {
	dist++;
	(*r)--;
      }
      (*r)--;

      debug(printf("V%d: ",dist));
      pairs = add_queryskip(pairs,(*r)+dist,*c,dist,querysequence,
			    queryoffset,genomeoffset,pairpool,revp,/*cdna_gap_p*/false,
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

  debug11(printf("Starting traceback_local_nogaps at r=%d,c=%d (offset1=%d, offset2=%d), revp %d\n",
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


/* We have switched length2 for columns, and length 1L and 1R for rows */
static void
bridge_cdna_gap (int *finalscore, int *bestrL, int *bestrR, int *bestcL, int *bestcR,
		 struct Int3_T **matrixL, struct Int3_T **matrixR,
		 int length2, int length1L, int length1R, int extraband_paired,
		 int open, int extend, int leftoffset, int rightoffset, bool jump_late_p) {
  int bestscore = -100000, score, scoreL, scoreR, pen, end_reward = 0;
  int rL, rR, cL, cR;
  int lbandL, rbandL, cloL, chighL;
  int lbandR, rbandR, cloR, chighR;

  /* Perform computations */
  rbandL = length1L - length2 + extraband_paired;
  lbandL = extraband_paired;

  rbandR = length1R - length2 + extraband_paired;
  lbandR = extraband_paired;

  /* Need a double loop on rows here, in contrast with a single loop
     for introns, because we allow a genomic insertion that doesn't
     match the cDNA.  So we need to add a penalty for a genomic
     insertion */

  for (rL = 1; rL < length2; rL++) {

    /* Note: opening penalty is added at the bottom of the loop */
    for (rR = length2-rL, pen = 0; rR >= 0; rR--, pen += extend) {
      /* debug3(printf("\nAt row %d on left and %d on right\n",rL,rR)); */
      if ((cloL = rL - lbandL) < 1) {
	cloL = 1;
      }
      if ((chighL = rL + rbandL) > length1L-1) {
	chighL = length1L-1;
      }

      if ((cloR = rR - lbandR) < 1) {
	cloR = 1;
      }
      if ((chighR = rR + rbandR) > length1R-1) {
	chighR = length1R-1;
      }

      for (cL = cloL; cL <= chighL; cL++) {
	scoreL = matrixL[rL][cL].nogap;
	
	/* Disallow leftoffset + cL >= rightoffset - cR, or cR >= rightoffset - leftoffset - cL */
	/* debug3(printf("  Disallowing cR to be >= %d\n",rightoffset-leftoffset-cL)); */
	for (cR = cloR; cR <= chighR && cR < rightoffset-leftoffset-cL; cR++) {
	  scoreR = matrixR[rR][cR].nogap;

	  if ((score = scoreL + scoreR + pen + end_reward) > bestscore || (score == bestscore && jump_late_p)) {
	    /*
	    debug3(printf("At %d left to %d right, score is (%d)+(%d)+(%d)+(%d) = %d (bestscore)\n",
			  cL,cR,scoreL,scoreR,pen,end_reward,score));
	    */

	    bestscore = score;
	    *bestrL = rL;
	    *bestrR = rR;
	    *bestcL = cL;
	    *bestcR = cR;

	  } else {
	    /*
	      debug3(printf("At %d left to %d right, score is (%d)+(%d)+(%d) = %d\n",
	      cL,cR,scoreL,scoreR,pen,scoreL+scoreR+pen));
	    */
	  }
	}
      }
      pen = open - extend;	/* Subtract extend to compensate for
                                   its addition in the for loop */
    }
  }
      
  *finalscore = bestscore;
  debug3(printf("Returning final score of %d at (%d,%d) left to (%d,%d) right\n",
		*finalscore,*bestrL,*bestcL,*bestrR,*bestcR));

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
		      int *left_known, int *right_known, Genomicpos_T leftoffset, Genomicpos_T rightoffset,
		      Genomicpos_T chroffset, Genomicpos_T chrhigh, int cdna_direction, bool watsonp) {
  Genomicpos_T splicesitepos;
  
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


static bool
bridge_intron_gap (int *finalscore, int *bestrL, int *bestrR, int *bestcL, int *bestcR,
		   int *best_introntype, double *left_prob, double *right_prob,
		   struct Int3_T **matrixL, struct Int3_T **matrixR,
		   struct Direction3_T **directionsL, struct Direction3_T **directionsR, 
		   int offset2L, int revoffset2R, int length1, int length2L, int length2R,
		   int cdna_direction, bool watsonp, int extraband_paired, int canonical_reward,
		   int maxhorizjump, int maxvertjump, int leftoffset, int rightoffset,
		   Chrnum_T chrnum, Genomicpos_T chroffset, Genomicpos_T chrhigh,
		   bool halfp, bool finalp, bool use_probabilities_p, int score_threshold,
		   bool jump_late_p) {
  bool result;
  int bestscore = -100000, score, scoreL, scoreR, scoreI, bestscoreI = -100000;
  int rL, rR, cL, cR;
  int lbandL, rbandL, cloL, chighL;
  int lbandR, rbandR, cloR, chighR;
  char left1, left2, right2, right1, left1_alt, left2_alt, right2_alt, right1_alt;
  int *leftdi, *rightdi, introntype;
  int *left_known, *right_known;
  double *left_probabilities, *right_probabilities, probL, probR, bestprob;
  Genomicpos_T splicesitepos, splicesitepos1, splicesitepos2;
  bool bestp;


  debug(printf("Running bridge_intron_gap with use_probabilities_p %d\n",use_probabilities_p));

  if (length2L+1 <= 0) {
    fprintf(stderr,"Problem with length2L = %d\n",length2L);
    abort();
  }

  if (length2R+1 <= 0) {
    fprintf(stderr,"Problem with length2R = %d\n",length2R);
    abort();
  }

  /* Read dinucleotides */
  leftdi = (int *) CALLOC(length2L+1,sizeof(int));
  rightdi = (int *) CALLOC(length2R+1,sizeof(int));

  for (cL = 0; cL < length2L - 1; cL++) {
    left1 = get_genomic_nt(&left1_alt,offset2L+cL,chroffset,chrhigh,watsonp);
    left2 = get_genomic_nt(&left2_alt,offset2L+cL+1,chroffset,chrhigh,watsonp);
#ifdef EXTRACT_GENOMICSEG
    assert(left1 == sequenceuc2L[cL]);
    assert(left2 == sequenceuc2L[cL+1]);
#endif

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
  leftdi[length2L-1] = leftdi[length2L] = 0x00;

  for (cR = 0; cR < length2R - 1; cR++) {
    right2 = get_genomic_nt(&right2_alt,revoffset2R-cR-1,chroffset,chrhigh,watsonp);
    right1 = get_genomic_nt(&right1_alt,revoffset2R-cR,chroffset,chrhigh,watsonp);
#ifdef EXTRACT_GENOMICSEG
    assert(right2 == revsequenceuc2R[-cR-1]);
    assert(right1 == revsequenceuc2R[-cR]);
#endif

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
  rightdi[length2R-1] = rightdi[length2R] = 0x00;

  left_known = (int *) CALLOC(length2L+1,sizeof(int));
  right_known = (int *) CALLOC(length2R+1,sizeof(int));
  if (splicing_iit != NULL && donor_typeint >= 0 && acceptor_typeint >= 0) {
    /* Handle known splicing, splice site level */
    if (watsonp == true) {
      if (cdna_direction > 0) {
	for (cL = 0; cL < length2L - 1; cL++) {
	  splicesitepos = leftoffset + cL;
	  if (IIT_exists_with_divno_typed_signed(splicing_iit,splicing_divint_crosstable[chrnum],
						 splicesitepos,splicesitepos+1U,donor_typeint,/*sign*/+1) == true) {
	    debug5(printf("Found known donor at %u\n",splicesitepos));
	    left_known[cL] = KNOWN_SPLICESITE_REWARD;
	  }
	}

	for (cR = 0; cR < length2R - 1; cR++) {
	  splicesitepos = rightoffset - cR + 1;
	  if (IIT_exists_with_divno_typed_signed(splicing_iit,splicing_divint_crosstable[chrnum],
						 splicesitepos,splicesitepos+1U,acceptor_typeint,/*sign*/+1) == true) {
	    debug5(printf("Found known acceptor at %u\n",splicesitepos));
	    right_known[cR] = KNOWN_SPLICESITE_REWARD;
	  }
	}

      } else {
	for (cL = 0; cL < length2L - 1; cL++) {
	  splicesitepos = leftoffset + cL;
	  if (IIT_exists_with_divno_typed_signed(splicing_iit,splicing_divint_crosstable[chrnum],
						 splicesitepos,splicesitepos+1U,acceptor_typeint,/*sign*/-1) == true) {
	    debug5(printf("Found known antiacceptor at %u\n",splicesitepos));
	    left_known[cL] = KNOWN_SPLICESITE_REWARD;
	  }
	}

	for (cR = 0; cR < length2R - 1; cR++) {
	  splicesitepos = rightoffset - cR + 1;
	  if (IIT_exists_with_divno_typed_signed(splicing_iit,splicing_divint_crosstable[chrnum],
						 splicesitepos,splicesitepos+1U,donor_typeint,/*sign*/-1) == true) {
	    debug5(printf("Found known antidonor at %u\n",splicesitepos));
	    right_known[cR] = KNOWN_SPLICESITE_REWARD;
	  }
	}
      }

    } else {
      if (cdna_direction > 0) {
	for (cL = 0; cL < length2L - 1; cL++) {
	  splicesitepos = (chrhigh - chroffset) - leftoffset - cL + 1;
	  if (IIT_exists_with_divno_typed_signed(splicing_iit,splicing_divint_crosstable[chrnum],
						 splicesitepos,splicesitepos+1U,donor_typeint,/*sign*/-1) == true) {
	    debug5(printf("Found known antidonor at %u\n",splicesitepos));
	    left_known[cL] = KNOWN_SPLICESITE_REWARD;
	  }
	}

	for (cR = 0; cR < length2R - 1; cR++) {
	  splicesitepos = (chrhigh - chroffset) - rightoffset + cR;
	  if (IIT_exists_with_divno_typed_signed(splicing_iit,splicing_divint_crosstable[chrnum],
						 splicesitepos,splicesitepos+1U,acceptor_typeint,/*sign*/-1) == true) {
	    debug5(printf("Found known antiacceptor at %u\n",splicesitepos));
	    right_known[cR] = KNOWN_SPLICESITE_REWARD;
	  }
	}

      } else {
	for (cL = 0; cL < length2L - 1; cL++) {
	  splicesitepos = (chrhigh - chroffset) - leftoffset - cL + 1;
	  if (IIT_exists_with_divno_typed_signed(splicing_iit,splicing_divint_crosstable[chrnum],
						 splicesitepos,splicesitepos+1U,acceptor_typeint,/*sign*/+1) == true) {
	    debug5(printf("Found known acceptor at %u\n",splicesitepos));
	    left_known[cL] = KNOWN_SPLICESITE_REWARD;
	  }
	}

	for (cR = 0; cR < length2R - 1; cR++) {
	  splicesitepos = (chrhigh - chroffset) - rightoffset + cR;
	  if (IIT_exists_with_divno_typed_signed(splicing_iit,splicing_divint_crosstable[chrnum],
						 splicesitepos,splicesitepos+1U,donor_typeint,/*sign*/+1) == true) {
	    debug5(printf("Found known donor at %u\n",splicesitepos));
	    right_known[cR] = KNOWN_SPLICESITE_REWARD;
	  }
	}
      }
    }

  } else if (splicing_iit != NULL) {
    /* Handle known splicing, intron level */
    if (watsonp == true) {
      if (cdna_direction > 0) {
	for (cL = 0; cL < length2L - 1; cL++) {
	  splicesitepos = leftoffset + cL;
	  if (IIT_low_exists_signed_p(splicing_iit,splicing_divint_crosstable[chrnum],
				      splicesitepos,/*sign*/+1) == true) {
	    debug5(printf("Found known donor at %u\n",splicesitepos));
	    left_known[cL] = KNOWN_SPLICESITE_REWARD;
	  }
	}

	for (cR = 0; cR < length2R - 1; cR++) {
	  splicesitepos = rightoffset - cR + 1;
	  if (IIT_high_exists_signed_p(splicing_iit,splicing_divint_crosstable[chrnum],
				       splicesitepos+1U,/*sign*/+1) == true) {
	    debug5(printf("Found known acceptor at %u\n",splicesitepos+1U));
	    right_known[cR] = KNOWN_SPLICESITE_REWARD;
	  }
	}

      } else {
	for (cL = 0; cL < length2L - 1; cL++) {
	  splicesitepos = leftoffset + cL;
	  if (IIT_low_exists_signed_p(splicing_iit,splicing_divint_crosstable[chrnum],
				      splicesitepos,/*sign*/-1) == true) {
	    debug5(printf("Found known antiacceptor at %u\n",splicesitepos));
	    left_known[cL] = KNOWN_SPLICESITE_REWARD;
	  }
	}

	for (cR = 0; cR < length2R - 1; cR++) {
	  splicesitepos = rightoffset - cR + 1;
	  if (IIT_high_exists_signed_p(splicing_iit,splicing_divint_crosstable[chrnum],
				       splicesitepos+1U,/*sign*/-1) == true) {
	    debug5(printf("Found known antidonor at %u\n",splicesitepos+1U));
	    right_known[cR] = KNOWN_SPLICESITE_REWARD;
	  }
	}
      }

    } else {
      if (cdna_direction > 0) {
	for (cL = 0; cL < length2L - 1; cL++) {
	  splicesitepos = (chrhigh - chroffset) - leftoffset - cL + 1;
	  if (IIT_high_exists_signed_p(splicing_iit,splicing_divint_crosstable[chrnum],
				       splicesitepos+1U,/*sign*/-1) == true) {
	    debug5(printf("Found known antidonor at %u\n",splicesitepos+1U));
	    left_known[cL] = KNOWN_SPLICESITE_REWARD;
	  }
	}

	for (cR = 0; cR < length2R - 1; cR++) {
	  splicesitepos = (chrhigh - chroffset) - rightoffset + cR;
	  if (IIT_low_exists_signed_p(splicing_iit,splicing_divint_crosstable[chrnum],
				      splicesitepos,/*sign*/-1) == true) {
	    debug5(printf("Found known antiacceptor at %u\n",splicesitepos));
	    right_known[cR] = KNOWN_SPLICESITE_REWARD;
	  }
	}

      } else {
	for (cL = 0; cL < length2L - 1; cL++) {
	  splicesitepos = (chrhigh - chroffset) - leftoffset - cL + 1;
	  if (IIT_high_exists_signed_p(splicing_iit,splicing_divint_crosstable[chrnum],
				       splicesitepos+1U,/*sign*/+1) == true) {
	    debug5(printf("Found known acceptor at %u\n",splicesitepos));
	    left_known[cL] = KNOWN_SPLICESITE_REWARD;
	  }
	}

	for (cR = 0; cR < length2R - 1; cR++) {
	  splicesitepos = (chrhigh - chroffset) - rightoffset + cR;
	  if (IIT_low_exists_signed_p(splicing_iit,splicing_divint_crosstable[chrnum],
				      splicesitepos,/*sign*/+1) == true) {
	    debug5(printf("Found known donor at %u\n",splicesitepos));
	    right_known[cR] = KNOWN_SPLICESITE_REWARD;
	  }
	}
      }
    }
  }

  /* Perform computations */
#if 0
  /* Allows unlimited indel lengths */
  rbandL = length2L - length1 + extraband_paired;
  lbandL = extraband_paired;

  rbandR = length2R - length1 + extraband_paired;
  lbandR = extraband_paired;
#else
  /* Limit indels to 3 bp around splice sites */
  rbandL = 3;
  lbandL = 3;

  rbandR = 3;
  lbandR = 3;
#endif


  if (novelsplicingp == false && splicing_iit != NULL && (donor_typeint < 0 || acceptor_typeint < 0)) {
    /* Constrain to given introns */
    for (rL = 1, rR = length1-1; rL < length1; rL++, rR--) {
      debug3(printf("\nGenomic insert: At row %d on left and %d on right\n",rL,rR));
      if ((cloL = rL - lbandL) < 1) {
	cloL = 1;
      }
      if ((chighL = rL + rbandL) > length2L-1) {
	chighL = length2L-1;
      }

      if ((cloR = rR - lbandR) < 1) {
	cloR = 1;
      }
      if ((chighR = rR + rbandR) > length2R-1) {
	chighR = length2R-1;
      }

      /* Test indel on left... */
      for (cL = cloL; cL <= chighL; cL++) {
	/* The following check limits genomic inserts (horizontal) and
	   multiple cDNA inserts (vertical). */
	if (left_known[cL] > 0) {
	  scoreL = matrixL[rL][cL].nogap;

	  if (directionsL[rL][cL].nogap == HORIZ || directionsL[rL][cL].nogap == VERT) {
	    /* Favor gaps away from intron if possible */
	    scoreL -= 1;
	  }

	  /* ...but not on right */
	  cR = rR;

	  /* Disallow leftoffset + cL >= rightoffset - cR, or cR >= rightoffset - leftoffset - cL */
	  if (cR < rightoffset-leftoffset-cL) {
	    if (right_known[cR] > 0) {
	      scoreR = matrixR[rR][cR].nogap;

#if 0
	      /* since we are on diagonal */
	      if (directionsR[rR][cR].nogap == HORIZ || directionsR[rR][cR].nogap == VERT) {
		/* Favor gaps away from intron if possible */
		scoreR -= 1;
	      }
#endif

	      if ((score = scoreL + scoreR) > bestscore || (score == bestscore && jump_late_p)) {
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

      /* Test indel on right... */
      for (cR = cloR; cR <= chighR; cR++) {
	if (right_known[cR] > 0) {
	  scoreR = matrixR[rR][cR].nogap;

	  if (directionsR[rR][cR].nogap == HORIZ || directionsR[rR][cR].nogap == VERT) {
	    /* Favor gaps away from intron if possible */
	    scoreR -= 1;
	  }

	  /* ...but not on left */
	  cL = rL;

	  /* Disallow leftoffset + cL >= rightoffset - cR, or cL >= rightoffset - leftoffset - cR */
	  if (cL < rightoffset-leftoffset-cR) {
	    if (left_known[cL] > 0) {
	      scoreL = matrixL[rL][cL].nogap;
	      
#if 0
	      /* since we are on diagonal */
	      if (directionsL[rL][cL].nogap == HORIZ || directionsL[rL][cL].nogap == VERT) {
		/* Favor gaps away from intron if possible */
		scoreL -= 1;
	      }
#endif

	      if ((score = scoreL + scoreR) > bestscore || (score == bestscore && jump_late_p)) {
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

    *finalscore = bestscore;
    *best_introntype = NONINTRON;

  } else if (use_probabilities_p == false) {
    /* Check for genomic insertion */
    for (rL = 1, rR = length1-1; rL < length1; rL++, rR--) {
      debug3(printf("\nGenomic insert: At row %d on left and %d on right\n",rL,rR));
      if ((cloL = rL - lbandL) < 1) {
	cloL = 1;
      }
      if ((chighL = rL + rbandL) > length2L-1) {
	chighL = length2L-1;
      }

      if ((cloR = rR - lbandR) < 1) {
	cloR = 1;
      }
      if ((chighR = rR + rbandR) > length2R-1) {
	chighR = length2R-1;
      }

      /* Test indel on left... */
      for (cL = cloL; cL <= chighL; cL++) {
	/* The following check limits genomic inserts (horizontal) and
	   multiple cDNA inserts (vertical). */
	if (1) {
	  scoreL = matrixL[rL][cL].nogap;
	  scoreL += left_known[cL];

	  if (directionsL[rL][cL].nogap == HORIZ || directionsL[rL][cL].nogap == VERT) {
	    /* Favor gaps away from intron if possible */
	    scoreL -= 1;
	  }

	  /* ...but not on right */
	  cR = rR;

	  /* Disallow leftoffset + cL >= rightoffset - cR, or cR >= rightoffset - leftoffset - cL */
	  if (cR < rightoffset-leftoffset-cL) {
	    if (1) {
	      scoreR = matrixR[rR][cR].nogap;
	      scoreR += right_known[cR];

#if 0
	      /* since we are on diagonal */
	      if (directionsR[rR][cR].nogap == HORIZ || directionsR[rR][cR].nogap == VERT) {
		/* Favor gaps away from intron if possible */
		scoreR -= 1;
	      }
#endif

	      scoreI = intron_score(&introntype,leftdi[cL],rightdi[cR],
				    cdna_direction,canonical_reward,finalp);

	      if ((score = scoreL + scoreI + scoreR) > bestscore || (score == bestscore && jump_late_p)) {
		debug3(printf("At %d left to %d right, score is (%d)+(%d)+(%d) = %d (bestscore)\n",
			      cL,cR,scoreL,scoreI,scoreR,scoreL+scoreI+scoreR));
		bestscore = score;
		bestscoreI = scoreI;
		*bestrL = rL;
		*bestrR = rR;
		*bestcL = cL;
		*bestcR = cR;
		*best_introntype = introntype;
	      } else {
		debug3a(printf("At %d left to %d right, score is (%d)+(%d)+(%d) = %d\n",
			       cL,cR,scoreL,scoreI,scoreR,score));
	      }
	    }
	  }
	}
      }

      /* Test indel on right...*/
      for (cR = cloR; cR <= chighR; cR++) {
	if (1) {
	  scoreR = matrixR[rR][cR].nogap;
	  scoreR += right_known[cR];

	  if (directionsR[rR][cR].nogap == HORIZ || directionsR[rR][cR].nogap == VERT) {
	    /* Favor gaps away from intron if possible */
	    scoreR -= 1;
	  }

	  /* ...but not on left */
	  cL = rL;

	  /* Disallow leftoffset + cL >= rightoffset - cR, or cL >= rightoffset - leftoffset - cR */
	  if (cL < rightoffset-leftoffset-cR) {
	    if (1) {
	      scoreL = matrixL[rL][cL].nogap;
	      scoreL += left_known[cL];

#if 0
	      /* since we are on diagonal */
	      if (directionsL[rL][cL].nogap == HORIZ || directionsL[rL][cL].nogap == VERT) {
		/* Favor gaps away from intron if possible */
		scoreL -= 1;
	      }
#endif

	      scoreI = intron_score(&introntype,leftdi[cL],rightdi[cR],
				    cdna_direction,canonical_reward,finalp);

	      if ((score = scoreL + scoreI + scoreR) > bestscore || (score == bestscore && jump_late_p)) {
		debug3(printf("At %d left to %d right, score is (%d)+(%d)+(%d) = %d (bestscore)\n",
			      cL,cR,scoreL,scoreI,scoreR,scoreL+scoreI+scoreR));
		bestscore = score;
		bestscoreI = scoreI;
		*bestrL = rL;
		*bestrR = rR;
		*bestcL = cL;
		*bestcR = cR;
		*best_introntype = introntype;
	      } else {
		debug3a(printf("At %d left to %d right, score is (%d)+(%d)+(%d) = %d\n",
			       cL,cR,scoreL,scoreI,scoreR,score));
	      }
	    }
	  }
	}
      }
    }

#if 0
    /* Previously had code here to allow cDNA insert */
#endif

    if (halfp == true) {
      *finalscore = bestscore - bestscoreI/2;
    } else {
      *finalscore = bestscore;
    }

  } else {
    /* Do search based on splice site probabilities */
    left_probabilities = (double *) CALLOC(length2L+1,sizeof(double));
    right_probabilities = (double *) CALLOC(length2R+1,sizeof(double));

    if (watsonp == true) {
      if (cdna_direction > 0) {
	for (cL = 0; cL < length2L - 1; cL++) {
	  splicesitepos = chroffset + leftoffset + cL;
	  left_probabilities[cL] = left_known[cL] ? 1.0 : Maxent_hr_donor_prob(splicesitepos,chroffset);
	}

	for (cR = 0; cR < length2R - 1; cR++) {
	  splicesitepos = chroffset + rightoffset - cR + 1;
	  right_probabilities[cR] = right_known[cR] ? 1.0 : Maxent_hr_acceptor_prob(splicesitepos,chroffset);
	}

      } else {
	for (cL = 0; cL < length2L - 1; cL++) {
	  splicesitepos = chroffset + leftoffset + cL;
	  left_probabilities[cL] = left_known[cL] ? 1.0 : Maxent_hr_antiacceptor_prob(splicesitepos,chroffset);
	}

	for (cR = 0; cR < length2R - 1; cR++) {
	  splicesitepos = chroffset + rightoffset - cR + 1;
	  right_probabilities[cR] = right_known[cR] ? 1.0 : Maxent_hr_antidonor_prob(splicesitepos,chroffset);
	}
      }

    } else {
      if (cdna_direction > 0) {
	for (cL = 0; cL < length2L - 1; cL++) {
	  splicesitepos = chrhigh - leftoffset - cL + 1;
	  left_probabilities[cL] = left_known[cL] ? 1.0 : Maxent_hr_antidonor_prob(splicesitepos,chroffset);
	}

	for (cR = 0; cR < length2R - 1; cR++) {
	  splicesitepos = chrhigh - rightoffset + cR;
	  right_probabilities[cR] = right_known[cR] ? 1.0 : Maxent_hr_antiacceptor_prob(splicesitepos,chroffset);
	}

      } else {
	for (cL = 0; cL < length2L - 1; cL++) {
	  splicesitepos = chrhigh - leftoffset - cL + 1;
	  left_probabilities[cL] = left_known[cL] ? 1.0 : Maxent_hr_acceptor_prob(splicesitepos,chroffset);
	}

	for (cR = 0; cR < length2R - 1; cR++) {
	  splicesitepos = chrhigh - rightoffset + cR;
	  right_probabilities[cR] = right_known[cR] ? 1.0 : Maxent_hr_donor_prob(splicesitepos,chroffset);
	}
      }
    }

    bestprob = 0.0;

    /* Check for genomic insert */
    for (rL = 1, rR = length1-1; rL < length1; rL++, rR--) {
      debug3(printf("\nAt row %d on left and %d on right\n",rL,rR));
      if ((cloL = rL - lbandL) < 1) {
	cloL = 1;
      }
      if ((chighL = rL + rbandL) > length2L-1) {
	chighL = length2L-1;
      }

      if ((cloR = rR - lbandR) < 1) {
	cloR = 1;
      }
      if ((chighR = rR + rbandR) > length2R-1) {
	chighR = length2R-1;
      }

      /* Test indel on left... */
      for (cL = cloL; cL <= chighL; cL++) {
	/* The following check limits genomic inserts (horizontal) and
	   multiple cDNA inserts (vertical). */
	if (1) {
	  probL = left_probabilities[cL];

	  /* ...but not on right */
	  cR = rR;

	  /* Disallow leftoffset + cL >= rightoffset - cR, or cR >= rightoffset - leftoffset - cL */
	  if (cR < rightoffset-leftoffset-cL) {
	    if (1) {
	      probR = right_probabilities[cR];

	      if (probL + probR <= bestprob) {
		debug3a(printf("At %d left to %d right, prob is %f + %f = %f\n",
			       cL,cR,probL,probR,probL+probR));
	      } else {
		debug3(printf("At %d left to %d right, prob is %f + %f = %f\n",
			      cL,cR,probL,probR,probL+probR));

		scoreL = matrixL[rL][cL].nogap;
		scoreL += left_known[cL];

		if (directionsL[rL][cL].nogap == HORIZ || directionsL[rL][cL].nogap == VERT) {
		  /* Favor gaps away from intron if possible */
		  scoreL -= 1;
		}

		scoreR = matrixR[rR][cR].nogap;
		scoreR += right_known[cR];
	    
#if 0
		/* since we are on diagonal */
		if (directionsR[rR][cR].nogap == HORIZ || directionsR[rR][cR].nogap == VERT) {
		  /* Favor gaps away from intron if possible */
		  scoreR -= 1;
		}
#endif

		scoreI = intron_score(&introntype,leftdi[cL],rightdi[cR],
				      cdna_direction,canonical_reward,finalp);

		if (scoreL + scoreI + scoreR >= score_threshold) {
		  debug3(printf(" (bestscore where score %d >= threshold %d)\n",
				scoreL+scoreI+scoreR,score_threshold));
		  bestprob = probL + probR;
		  *bestrL = rL;
		  *bestrR = rR;
		  *bestcL = cL;
		  *bestcR = cR;
		} else {
		  debug3a(printf(" (bestscore, but score %d < threshold %d)\n",
				 scoreL+scoreI+scoreR,score_threshold));
		}
	      }
	    }
	  }
	}
      }

      /* Test indel on right... */
      for (cR = cloR; cR <= chighR; cR++) {
	if (1) {
	  probR = right_probabilities[cR];

	  /* ...but not on left */
	  cL = rL;

	  /* Disallow leftoffset + cL >= rightoffset - cR, or cL >= rightoffset - leftoffset - cR */
	  if (cL < rightoffset-leftoffset-cR) {
	    if (1) {
	      probL = left_probabilities[cL];

	      if (probL + probR <= bestprob) {
		debug3a(printf("At %d left to %d right, prob is %f + %f = %f\n",
			       cL,cR,probL,probR,probL+probR));
	      } else {
		debug3(printf("At %d left to %d right, prob is %f + %f = %f\n",
			      cL,cR,probL,probR,probL+probR));

		scoreL = matrixL[rL][cL].nogap;
		scoreL += left_known[cL];

#if 0
		/* since we are on diagonal */
		if (directionsL[rL][cL] == HORIZ || directionsL[rL][cL] == VERT) {
		  /* Favor gaps away from intron if possible */
		  scoreL -= 1;
		}
#endif

		scoreR = matrixR[rR][cR].nogap;
		scoreR += right_known[cR];
	    
		if (directionsR[rR][cR].nogap == HORIZ || directionsR[rR][cR].nogap == VERT) {
		  /* Favor gaps away from intron if possible */
		  scoreR -= 1;
		}

		scoreI = intron_score(&introntype,leftdi[cL],rightdi[cR],
				      cdna_direction,canonical_reward,finalp);

		if (scoreL + scoreI + scoreR >= score_threshold) {
		  debug3(printf(" (bestscore where score %d >= threshold %d)\n",
				scoreL+scoreI+scoreR,score_threshold));
		  bestprob = probL + probR;
		  *bestrL = rL;
		  *bestrR = rR;
		  *bestcL = cL;
		  *bestcR = cR;
		} else {
		  debug3(printf(" (bestscore, but score %d < threshold %d)\n",
				scoreL+scoreI+scoreR,score_threshold));
		}
	      }
	    }
	  }
	}
      }
    }

#if 0
    /* Previously had code to check for cDNA insert of up to 3 */
#endif

    FREE(left_probabilities);
    FREE(right_probabilities);

    /* Compute finalscore */
    scoreL = matrixL[*bestrL][*bestcL].nogap;
    scoreL += left_known[*bestcL];

    if (directionsL[*bestrL][*bestcL].nogap == HORIZ || directionsL[*bestrL][*bestcL].nogap == VERT) {
      /* Favor gaps away from intron if possible */
      scoreL -= 1;
    }

    scoreR = matrixR[*bestrR][*bestcR].nogap;
    scoreR += right_known[*bestcR];

    if (directionsR[*bestrR][*bestcR].nogap == HORIZ || directionsR[*bestrR][*bestcR].nogap == VERT) {
      /* Favor gaps away from intron if possible */
      scoreR -= 1;
    }

    scoreI = intron_score(&introntype,leftdi[*bestcL],rightdi[*bestcR],
			  cdna_direction,canonical_reward,finalp);

    if (halfp == true) {
      *finalscore = scoreL + scoreI + scoreR - scoreI/2;
    } else {
      debug3(printf("At %d left to %d right, score is (%d)+(%d)+(%d) = %d\n",
		    *bestcL,*bestcR,scoreL,scoreI,scoreR,scoreL+scoreI+scoreR));
      *finalscore = scoreL + scoreI + scoreR;
    }
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


  if (finalp == true && result == true) {
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
			 struct Int3_T **matrixL, struct Int3_T **matrixR,
			 int diaglength, int offset2L, int revoffset2R,
			 Genomicpos_T chroffset, Genomicpos_T chrhigh,
			 int cdna_direction, bool watsonp, bool jump_late_p) {
  int bestscore = NEG_INFINITY, score, scoreL, scoreR, scoreI;
  int rcL, rcR;
  char left1, left2, right2, right1, left1_alt, left2_alt, right2_alt, right1_alt;
  int introntype, leftdi, rightdi;

  *bestrcL = *bestrcR = 0;

  for (rcL = 1; rcL <= diaglength; rcL++) {
    rcR = diaglength - rcL;
    left1 = get_genomic_nt(&left1_alt,offset2L+rcL,chroffset,chrhigh,watsonp);
    left2 = get_genomic_nt(&left2_alt,offset2L+rcL+1,chroffset,chrhigh,watsonp);
#ifdef EXTRACT_GENOMICSEG
    assert(left1 == sequenceuc2L[rcL]);
    assert(left2 == sequenceuc2L[rcL+1]);
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

    right2 = get_genomic_nt(&right2_alt,revoffset2R-rcR-1,chroffset,chrhigh,watsonp);
    right1 = get_genomic_nt(&right1_alt,revoffset2R-rcR,chroffset,chrhigh,watsonp);
#ifdef EXTRACT_GENOMICSEG
    assert(right2 == revsequenceuc2R[-1-rcR]);
    assert(right1 == revsequenceuc2R[-rcR]);
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
    
    scoreL = matrixL[rcL][rcL].nogap;
    scoreI = intron_score(&introntype,leftdi,rightdi,
			  cdna_direction,/*canonical_reward*/FINAL_CANONICAL_INTRON_HIGHQ,/*finalp*/true);
    scoreR = matrixR[rcR][rcR].nogap;

    if ((score = scoreL + scoreI + scoreR) > bestscore || (score == bestscore && jump_late_p)) {
      *bestrcL = rcL;
      *bestrcR = rcR;
      bestscore = score;
    }
  }

  *finalscore = bestscore;
  return;
}
#endif



#if 0
static void
bridge_dual_break_fwd (int *finalscore, int *bestrcL, int *bestrcR,
		       struct Int3_T **matrixL, struct Int3_T **matrixR,
		       int diaglength, int offset2L, int revoffset2R,
		       Genomicpos_T chroffset, Genomicpos_T chrhigh, bool watsonp) {
  int bestscore, scoreL, scoreR, scoreI;
  int rcL, rcR;
  int bestscoreL_GT, bestscoreL_GC, bestscoreL_AT, bestscoreL_XX;
  int bestrcL_GT, bestrcL_GC, bestrcL_AT, bestrcL_XX;
  int bestscoreR_AG, bestscoreR_AC, bestscoreR_XX;
  int bestrcR_AG, bestrcR_AC, bestrcR_XX;
  char left1, left2, right2, right1, left1_alt, left2_alt, right2_alt, right1_alt;
  int introntype;

  *bestrcL = *bestrcR = 0;

  bestrcL_GT = bestrcL_GC = bestrcL_AT = bestrcL_XX = 0;
  bestscoreL_GT = bestscoreL_GC = bestscoreL_AT = bestscoreL_XX = NEG_INFINITY;

  for (rcL = 1; rcL <= diaglength; rcL++) {
    if ((scoreL = matrixL[rcL][rcL].nogap) >= bestscoreL_XX) {
      bestscoreL_XX = scoreL;
      bestrcL_XX = rcL;
    }

    left1 = get_genomic_nt(&left1_alt,offset2L+rcL,chroffset,chrhigh,watsonp);
    left2 = get_genomic_nt(&left2_alt,offset2L+rcL+1,chroffset,chrhigh,watsonp);
#ifdef EXTRACT_GENOMICSEG
    assert(left1 == sequenceuc2L[rcL]);
    assert(left2 == sequenceuc2L[rcL+1]);
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
  bestscoreR_AG = bestscoreR_AC = bestscoreR_XX = NEG_INFINITY;

  for (rcR = 1; rcR <= diaglength; rcR++) {
    if ((scoreR = matrixR[rcR][rcR].nogap) >= bestscoreR_XX) {
      bestscoreR_XX = scoreR;
      bestrcR_XX = rcR;
    }
    right2 = get_genomic_nt(&right2_alt,revoffset2R-rcR-1,chroffset,chrhigh,watsonp);
    right1 = get_genomic_nt(&right1_alt,revoffset2R-rcR,chroffset,chrhigh,watsonp);
#ifdef EXTRACT_GENOMICSEG
    assert(right2 == revsequenceuc2R[-1-rcR]);
    assert(right1 == revsequenceuc2R[-rcR]);
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

  *finalscore = bestscore;
  /*
  fprintf(stderr,"%d + %d >=? %d\n",*bestrcL,*bestrcR,diaglength1);
  fprintf(stderr,"%d + %d >=? %d\n\n",*bestrcL,*bestrcR,diaglength2);
  */
  if (*bestrcL + *bestrcR >= diaglength) {
    bridge_dual_break_nogap(&(*finalscore),&(*bestrcL),&(*bestrcR),matrixL,matrixR,
			    diaglength,offset2L,revoffset2R,chroffset,chrhigh,
			    /*cdna_direction*/+1,watsonp,/*jump_late_p*/false);
  }

  return;
}
#endif


#if 0
static void
bridge_dual_break_rev (int *finalscore, int *bestrcL, int *bestrcR,
		       struct Int3_T **matrixL, struct Int3_T **matrixR,
		       int diaglength, int offset2L, int revoffset2R,
		       Genomicpos_T chroffset, Genomicpos_T chrhigh, bool watsonp) {
  int bestscore, scoreL, scoreR, scoreI;
  int rcL, rcR;
  int bestscoreL_CT, bestscoreL_GT, bestscoreL_XX;
  int bestrcL_CT, bestrcL_GT, bestrcL_XX;
  int bestscoreR_AC, bestscoreR_GC, bestscoreR_AT, bestscoreR_XX;
  int bestrcR_AC, bestrcR_GC, bestrcR_AT, bestrcR_XX;
  char left1, left2, right2, right1, left1_alt, left2_alt, right2_alt, right1_alt;
  int introntype;

  *bestrcL = *bestrcR = 0;

  bestrcL_CT = bestrcL_GT = bestrcL_XX = 0;
  bestscoreL_CT = bestscoreL_GT = bestscoreL_XX = NEG_INFINITY;

  for (rcL = 1; rcL <= diaglength; rcL++) {
    if ((scoreL = matrixL[rcL][rcL].nogap) >= bestscoreL_XX) {
      bestscoreL_XX = scoreL;
      bestrcL_XX = rcL;
    }

    left1 = get_genomic_nt(&left1_alt,offset2L+rcL,chroffset,chrhigh,watsonp);
    left2 = get_genomic_nt(&left2_alt,offset2L+rcL+1,chroffset,chrhigh,watsonp);
#ifdef EXTRACT_GENOMICSEG
    assert(left1 == sequenceuc2L[rcL]);
    assert(left2 == sequenceuc2L[rcL+1]);
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
  bestscoreR_AC = bestscoreR_GC = bestscoreR_AT = bestscoreR_XX = NEG_INFINITY;

  for (rcR = 1; rcR <= diaglength; rcR++) {
    if ((scoreR = matrixR[rcR][rcR].nogap) >= bestscoreR_XX) {
      bestscoreR_XX = scoreR;
      bestrcR_XX = rcR;
    }

    right2 = get_genomic_nt(&right2_alt,revoffset2R-rcR-1,chroffset,chrhigh,watsonp);
    right1 = get_genomic_nt(&right1_alt,revoffset2R-rcR,chroffset,chrhigh,watsonp);
#ifdef EXTRACT_GENOMICSEG
    assert(right2 == revsequenceuc2R[-1-rcR]);
    assert(right1 == revsequenceuc2R[-rcR]);
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

  *finalscore = bestscore;
  /*
  fprintf(stderr,"%d + %d >=? %d\n",*bestrcL,*bestrcR,diaglength1);
  fprintf(stderr,"%d + %d >=? %d\n\n",*bestrcL,*bestrcR,diaglength2);
  */
  if (*bestrcL + *bestrcR >= diaglength) {
    bridge_dual_break_nogap(&(*finalscore),&(*bestrcL),&(*bestrcR),matrixL,matrixR,
			    diaglength,offset2L,revoffset2R,chroffset,chrhigh,
			    /*cdna_direction*/-1,watsonp,/*jump_late_p*/true);
  }

  return;
}
#endif



/************************************************************************/

List_T
Dynprog_single_gap (int *dynprogindex, int *finalscore, int *nmatches, int *nmismatches, int *nopens, int *nindels,
		    T dynprog, char *sequence1, char *sequenceuc1,
		    int length1, int length2, int offset1, int offset2,
		    Genomicpos_T chroffset, Genomicpos_T chrhigh,
#ifdef PMAP
		    char *queryaaseq,
#endif
		    int cdna_direction, bool watsonp, bool jump_late_p, Pairpool_T pairpool,
		    int extraband_single, double defect_rate, int close_indels_mode, bool widebandp) {
  List_T pairs = NULL;
#ifdef PMAP
  char *instsequence1;
#endif
  Mismatchtype_T mismatchtype;
  int open, extend;
  struct Int3_T **matrix;
  struct Direction3_T **directions;
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

  /* Length1: maxlookback+MAXPEELBACK.  Length2 +EXTRAMATERIAL */
  debug(
	printf("%c:  ",*dynprogindex > 0 ? (*dynprogindex-1)%26+'a' : (-(*dynprogindex)-1)%26+'A');
	printf("Aligning single gap middle with wideband = %d and extraband %d\n",widebandp,extraband_single);
	printf("At query offset %d-%d, %.*s\n",offset1,offset1+length1-1,length1,sequence1));
#ifdef EXTRACT_GENOMICSEG
  debug(printf("At genomic offset %d-%d, %.*s\n",offset2,offset2+length2-1,length2,sequence2));
#endif
  debug(printf("\n"));

  if (length1 > dynprog->maxlength1 || length2 > dynprog->maxlength2) {
    debug(printf("length1 %d or length2 %d is too long.  Returning NULL\n",length1,length2));
    *finalscore = -10000;
    *nmatches = *nmismatches = *nopens = *nindels = 0;
#if 0
    /* Don't push a gapholder for single gap, because gapholder already exists */
    pairs = Pairpool_push_gapholder(NULL,pairpool,length1,length2,/*knownp*/false);
#endif
    *dynprogindex += (*dynprogindex > 0 ? +1 : -1);
    return (List_T) NULL;
  }

#ifdef PMAP
  instsequence1 = instantiate_codons(sequenceuc1,queryaaseq,offset1,length1);
#endif
  
    /* If extraband_single is too large, then gaps may be inserted on
       both sides, like this

       CACCC   AGAGCAGGCACTGCCT
       |||||--- ||| ||---||||| 
       CACCCCAGGGAGGAG   CTGCCC

    */

  matrix = compute_scores_lookup_fwd(&directions,dynprog,
#ifdef PMAP
				     instsequence1,
#else
				     sequence1,
#endif
				     offset2,length1,length2,
				     chroffset,chrhigh,watsonp,
				     mismatchtype,open,extend,
				     extraband_single,widebandp,jump_late_p);
  *finalscore = matrix[length1][length2].nogap;

  *nmatches = *nmismatches = *nopens = *nindels = 0;
  pairs = traceback(NULL,&(*nmatches),&(*nmismatches),&(*nopens),&(*nindels),
		    directions,length1,length2,
#ifdef PMAP
		    instsequence1,instsequence1,
#else
		    sequence1,sequenceuc1,
#endif
		    offset1,offset2,pairpool,/*revp*/false,
		    chroffset,chrhigh,cdna_direction,watsonp,*dynprogindex);

#ifdef PMAP
  FREE(instsequence1);
#endif

  /*
  Directions_free(directions);
  Matrix_free(matrix);
  */

  debug(printf("End of dynprog single gap\n"));

  *dynprogindex += (*dynprogindex > 0 ? +1 : -1);
  return List_reverse(pairs);
}




/* Sequence 1L and 1R represent the two ends of the cDNA insertion */
List_T
Dynprog_cdna_gap (int *dynprogindex, int *finalscore, bool *incompletep,
		  T dynprogL, T dynprogR, char *sequence1L, char *sequenceuc1L, 
		  char *revsequence1R, char *revsequenceuc1R,
#if 0
		  char *sequence2, char *sequenceuc2,
#endif
		  int length1L, int length1R, int length2,
		  int offset1L, int revoffset1R, int offset2,
		  Genomicpos_T chroffset, Genomicpos_T chrhigh,
#ifdef PMAP
		  char *queryaaseq,
#endif
		  int cdna_direction, bool watsonp, bool jump_late_p, Pairpool_T pairpool,
		  int extraband_paired, double defect_rate) {
  List_T pairs = NULL;
#ifdef PMAP
  char *instsequence1, *instsequence1L, *instrevsequence1R;
#endif
  Mismatchtype_T mismatchtype;
  struct Int3_T **matrixL, **matrixR;
  int mismatch, open, extend;
  struct Direction3_T **directionsL, **directionsR;
  int revoffset2, bestrL, bestrR, bestcL, bestcR, k;
  int nmatches, nmismatches, nopens, nindels;
  int queryjump, genomejump;
  char c2, c2_alt;


  if (length2 <= 1) {
    return NULL;
  }

  debug(
	printf("\n");
	printf("%c:  ",*dynprogindex > 0 ? (*dynprogindex-1)%26+'a' : (-(*dynprogindex)-1)%26+'A');
	printf("Aligning cdna gap\n");
	printf("At query offset %d-%d, %.*s\n",offset1L,offset1L+length1L-1,length1L,sequence1L);
	printf("At query offset %d-%d, %.*s\n",revoffset1R-length1R+1,revoffset1R,length1R,&(revsequence1R[-length1R+1]));
	printf("Whole piece at query offset %d-%d, %.*s\n",offset1L,revoffset1R,revoffset1R-offset1L+1,sequence1L));
#ifdef EXTRACT_GENOMICSEG
  debug(printf("At genomic offset %d-%d, %.*s\n",offset2,offset2+length2-1,length2,sequence2));
#endif
  debug(printf("\n"));

  /* ?check if offsets are too close.  But this eliminates a segment
     of the cDNA.  Should check in stage 3, and do single gap instead. */
  /*
  if (offset1L+length1L-1 >= revoffset1R-length1R+1) {
    debug(printf("Bounds don't make sense\n"));
    *finalscore = -100000.0;
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

  if (length2 > dynprogR->maxlength1 || length1R > dynprogR->maxlength2) {
    debug(printf("length2 %d or length1R %d is too long.  Returning NULL\n",length2,length1R));
#if 0
    revoffset2 = offset2 + length2 - 1;
    queryjump = revoffset1R - offset1L + 1;
    genomejump = revoffset2 - offset2 + 1;
    pairs = Pairpool_push_gapholder(NULL,pairpool,queryjump,genomejump,/*knownp*/false);
#endif
    *dynprogindex += (*dynprogindex > 0 ? +1 : -1);
    return (List_T) NULL;
  }

  if (length2 > dynprogL->maxlength1 || length1L > dynprogL->maxlength2) {
    debug(printf("length2 %d or length1L %d is too long.  Returning NULL\n",length2,length1L));
#if 0
    revoffset2 = offset2 + length2 - 1;
    queryjump = revoffset1R - offset1L + 1;
    genomejump = revoffset2 - offset2 + 1;
    pairs = Pairpool_push_gapholder(NULL,pairpool,queryjump,genomejump,/*knownp*/false);
#endif
    *dynprogindex += (*dynprogindex > 0 ? +1 : -1);
    return (List_T) NULL;
  }

#if 0
  /* Right side looks like 5' end */
  /* Note: sequence 1 and 2 flipped, because 1 has extramaterial */
  revsequence2 = &(sequence2[length2-1]);
  revsequenceuc2 = &(sequenceuc2[length2-1]);
#endif
  revoffset2 = offset2+length2-1;

#ifdef PMAP
  instsequence1 = instantiate_codons(sequence1L,queryaaseq,offset1L,revoffset1R-offset1L+1);
  instsequence1L = instsequence1;
  instrevsequence1R = &(instsequence1[revoffset1R-offset1L]);
#endif

  matrixR = compute_scores_lookup_rev_12(&directionsR,dynprogR,
#ifdef PMAP
					 instrevsequence1R,
#else
					 revsequence1R,
#endif
					 revoffset2,length2,length1R,
					 chroffset,chrhigh,watsonp,
					 mismatchtype,open,extend,
					 extraband_paired,/*widebandp*/true,
					 /*for revp true*/!jump_late_p);

  /* ? previously used compute_scores, not compute_scores_lookup */
  matrixL = compute_scores_lookup_fwd_12(&directionsL,dynprogL,
#ifdef PMAP
					 instsequence1L,
#else
					 sequence1L,
#endif
					 offset2,length2,length1L,
					 chroffset,chrhigh,watsonp,
					 mismatchtype,open,extend,
					 extraband_paired,/*widebandp*/true,
					 jump_late_p);

  bridge_cdna_gap(&(*finalscore),&bestrL,&bestrR,&bestcL,&bestcR,matrixL,matrixR,
		  length2,length1L,length1R,extraband_paired,
		  open,extend,offset1L,revoffset1R,jump_late_p);

  nmatches = nmismatches = nopens = nindels = 0;
  pairs = traceback_cdna(NULL,&nmatches,&nmismatches,&nopens,&nindels,
			 directionsR,bestrR,bestcR,
#ifdef PMAP
			 instrevsequence1R,instrevsequence1R,
#else
			 revsequence1R,revsequenceuc1R,
#endif
			 revoffset1R,revoffset2,pairpool,/*revp*/true,
			 chroffset,chrhigh,cdna_direction,watsonp,*dynprogindex);
  pairs = List_reverse(pairs);

  queryjump = (revoffset1R-bestcR) - (offset1L+bestcL) + 1;
  genomejump = (revoffset2-bestrR) - (offset2+bestrL) + 1;
  /* No need to revise queryjump or genomejump, because the above
     coordinates are internal to the gap. */

  if (queryjump == INSERT_PAIRS && genomejump == INSERT_PAIRS) {
    /* Add cDNA insertion, if any */
    for (k = revoffset1R-bestcR; k >= offset1L+bestcL; k--) {
#ifdef PMAP
      debug(printf("cDNA insertion, Pushing [%d,%d] (%c,-)\n",k,revoffset2-bestrR+1,instsequence1L[k-offset1L]));
      pairs = Pairpool_push(pairs,pairpool,k,revoffset2-bestrR+1,instsequence1L[k-offset1L],SHORTGAP_COMP,
			    /*genome*/' ',/*genomealt*/' ',*dynprogindex);
#else
      debug(printf("cDNA insertion, Pushing [%d,%d] (%c,-)\n",k,revoffset2-bestrR+1,sequence1L[k-offset1L]));
      pairs = Pairpool_push(pairs,pairpool,k,revoffset2-bestrR+1,sequence1L[k-offset1L],SHORTGAP_COMP,
			    /*genome*/' ',/*genomealt*/' ',*dynprogindex);
#endif
    }
    debug(printf("\n"));

    /* This loop not yet checked for get_genomic_nt giving correct answer */
    for (k = revoffset2-bestrR; k >= offset2+bestrL; k--) {
      c2 = get_genomic_nt(&c2_alt,k,chroffset,chrhigh,watsonp);
      debug(printf("genome insertion, Pushing [%d,%d] (-,%c)\n",offset1L+bestcL,k,c2));
#if 0
      assert(c2 == sequence2[k-offset2]);
      pairs = Pairpool_push(pairs,pairpool,offset1L+bestcL,k,' ',SHORTGAP_COMP,
			    sequence2[k-offset2],/*genomealt*/GENOMEALT_DEFERRED,*dynprogindex);
#else
      pairs = Pairpool_push(pairs,pairpool,offset1L+bestcL,k,' ',SHORTGAP_COMP,c2,c2_alt,*dynprogindex);
#endif
    }
    debug(printf("\n"));

  } else {

    /* Add gapholder to be solved in the future */
#ifndef NOGAPHOLDER
    debug(printf("Pushing a gap with queryjump = %d, genomejump = %d\n",queryjump,genomejump));
    pairs = Pairpool_push_gapholder(pairs,pairpool,/*queryjump*/UNKNOWNJUMP,/*genomejump*/UNKNOWNJUMP,/*knownp*/false);
#endif
    *incompletep = true;
  }

  pairs = traceback_cdna(pairs,&nmatches,&nmismatches,&nopens,&nindels,
			 directionsL,bestrL,bestcL,
#ifdef PMAP
			 instsequence1L,instsequence1L,
#else
			 sequence1L,sequenceuc1L,
#endif
			 offset1L,offset2,pairpool,/*revp*/false,
			 chroffset,chrhigh,cdna_direction,watsonp,*dynprogindex);

  /*
  Directions_free(directionsR);
  Directions_free(directionsL);
  Matrix_free(matrixR);
  Matrix_free(matrixL);
  */

#ifdef PMAP
  FREE(instsequence1);
#endif

  if (List_length(pairs) == 1) {
    /* Only a gap added */
    pairs = (List_T) NULL;
  }

  debug(printf("End of dynprog cDNA gap\n"));

  *dynprogindex += (*dynprogindex > 0 ? +1 : -1);
  return List_reverse(pairs);
}


/* A genome gap is usually an intron.  Sequence 2L and 2R represent
   the two genomic ends of the intron. */
List_T
Dynprog_genome_gap (int *dynprogindex, int *finalscore, int *new_leftgenomepos, int *new_rightgenomepos,
		    double *left_prob, double *right_prob,
		    int *nmatches, int *nmismatches, int *nopens, int *nindels,
		    int *exonhead, int *introntype, T dynprogL, T dynprogR, 
		    char *sequence1, char *sequenceuc1,
		    int length1, int length2L, int length2R, 
		    int offset1, int offset2L, int revoffset2R, 
		    Chrnum_T chrnum, Genomicpos_T chroffset, Genomicpos_T chrhigh,
#ifdef PMAP
		    char *queryaaseq,
#endif
		    int cdna_direction, bool watsonp, bool jump_late_p, Pairpool_T pairpool, int extraband_paired,
		    double defect_rate, int maxpeelback, bool halfp, bool finalp, bool use_probabilities_p,
		    int score_threshold, bool splicingp) {
  List_T pairs = NULL;
#ifdef PMAP
  char *instsequence1, *instrevsequence1;
#else
  char *revsequence1, *revsequenceuc1;
#endif
  Mismatchtype_T mismatchtype;
  int open, extend;
  struct Int3_T **matrixL, **matrixR;
  int canonical_reward;
  struct Direction3_T **directionsL, **directionsR;
  int revoffset1, bestrL, bestrR, bestcL, bestcR;
  int maxhorizjump, maxvertjump;
  /* int queryjump, genomejump; */

  debug(
	printf("\n");
	printf("%c:  ",*dynprogindex > 0 ? (*dynprogindex-1)%26+'a' : (-(*dynprogindex)-1)%26+'A');
	printf("Aligning genome gap with cdna_direction %d\n",cdna_direction);
	printf("At query offset %d-%d, %.*s\n",offset1,offset1+length1-1,length1,sequence1));
#ifdef EXTRACT_GENOMICSEG
  debug(printf("At genomic offset %d-%d, %.*s\n",offset2L,offset2L+length2L-1,length2L,sequence2L));
  debug(printf("At genomic offset %d-%d, %.*s\n",revoffset2R-length2R+1,revoffset2R,length2R,&(revsequence2R[-length2R+1])));
#endif
  debug(printf("\n"));

  /* ?check if offsets are too close.  But this eliminates a segment
     of the cDNA.  Should check in stage 3, and do single gap instead. */
  /*
  if (offset2L+length2L-1 >= revoffset2R-length2R+1) {
    debug(printf("Bounds don't make sense\n"));
    *finalscore = -100000.0;
    return NULL;
  }
  */

  *nmatches = *nmismatches = *nopens = *nindels = 0;
  *left_prob = *right_prob = 0.0;
  if (length1 <= 1) {
    *finalscore = NEG_INFINITY;
    *introntype = NONINTRON;
    return (List_T) NULL;
  }

  if (defect_rate < DEFECT_HIGHQ) {
    mismatchtype = HIGHQ;
    if (length1 > maxpeelback * 4) {
      debug(printf("length1 %d is greater than maxpeelback %d * 4, so using single gap penalties\n",
		   length1,maxpeelback));
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
    if (length1 > maxpeelback * 4) {
      debug(printf("length1 %d is greater than maxpeelback %d * 4, so using single gap penalties\n",
		   length1,maxpeelback));
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
    if (length1 > maxpeelback * 4) {
      debug(printf("length1 %d is greater than maxpeelback %d * 4, so using single gap penalties\n",
		   length1,maxpeelback));
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

  if (length1 > dynprogL->maxlength1 || length2L > dynprogL->maxlength2) {
    debug(printf("length1 %d or length2L %d is too long.  Returning NULL\n",length1,length2L));
    *new_leftgenomepos = offset2L-1;
    *new_rightgenomepos = revoffset2R+1;
    *exonhead = revoffset1 = offset1+length1-1;
#ifndef NOGAPHOLDER
    /*
    queryjump = revoffset1 - offset1 + 1;
    genomejump = revoffset2R - offset2L + 1;
    pairs = Pairpool_push_gapholder(NULL,pairpool,queryjump,genomejump,false);
    */
#endif
    *dynprogindex += (*dynprogindex > 0 ? +1 : -1);
    *finalscore = NEG_INFINITY;
    *introntype = NONINTRON;
    return (List_T) NULL;
  }

  if (length1 > dynprogR->maxlength1 || length2R > dynprogR->maxlength2) {
    debug(printf("length1 %d or length2R %d is too long.  Returning NULL\n",length1,length2R));
    *new_leftgenomepos = offset2L-1;
    *new_rightgenomepos = revoffset2R+1;
    *exonhead = revoffset1 = offset1+length1-1;
#ifndef NOGAPHOLDER
    /*
    queryjump = revoffset1 - offset1 + 1;
    genomejump = revoffset2R - offset2L + 1;
    pairs = Pairpool_push_gapholder(NULL,pairpool,queryjump,genomejump,false);
    */
#endif
    *dynprogindex += (*dynprogindex > 0 ? +1 : -1);
    *finalscore = NEG_INFINITY;
    *introntype = NONINTRON;
    return (List_T) NULL;
  }

#ifdef PMAP
  instsequence1 = instantiate_codons(sequence1,queryaaseq,offset1,length1);
  instrevsequence1 = &(instsequence1[length1-1]);
#else
  revsequence1 = &(sequence1[length1-1]);
  revsequenceuc1 = &(sequenceuc1[length1-1]);
#endif
  revoffset1 = offset1+length1-1;

  matrixL = compute_scores_lookup_fwd(&directionsL,dynprogL,
#ifdef PMAP
				      instsequence1,
#else
				      sequence1,
#endif
				      offset2L,length1,length2L,
				      chroffset,chrhigh,watsonp,
				      mismatchtype,open,extend,
				      extraband_paired,/*widebandp*/true,
				      jump_late_p);

  matrixR = compute_scores_lookup_rev(&directionsR,dynprogR,
#ifdef PMAP
				      instrevsequence1,
#else
				      revsequence1,
#endif
				      revoffset2R,length1,length2R,
				      chroffset,chrhigh,watsonp,
				      mismatchtype,open,extend,
				      extraband_paired,/*widebandp*/true,
				      /*for revp true*/!jump_late_p);

  if (bridge_intron_gap(&(*finalscore),&bestrL,&bestrR,&bestcL,&bestcR,
			&(*introntype),&(*left_prob),&(*right_prob),
			matrixL,matrixR,directionsL,directionsR,
			offset2L,revoffset2R,length1,length2L,length2R,
			cdna_direction,watsonp,extraband_paired,
			canonical_reward,maxhorizjump,maxvertjump,offset2L,revoffset2R,
			chrnum,chroffset,chrhigh,halfp,finalp,use_probabilities_p,
			score_threshold,jump_late_p) == false) {
    return (List_T) NULL;
  } else {
    *new_leftgenomepos = offset2L+(bestcL-1);
    *new_rightgenomepos = revoffset2R-(bestcR-1);
    debug(printf("New leftgenomepos = %d, New rightgenomepos = %d\n",*new_leftgenomepos,*new_rightgenomepos));

    *exonhead = revoffset1-(bestrR-1);

    pairs = traceback(NULL,&(*nmatches),&(*nmismatches),&(*nopens),&(*nindels),
		      directionsR,bestrR,bestcR,
#ifdef PMAP
		      instrevsequence1,instrevsequence1,
#else
		      revsequence1,revsequenceuc1,
#endif
		      revoffset1,revoffset2R,pairpool,/*revp*/true,
		      chroffset,chrhigh,cdna_direction,watsonp,*dynprogindex);
    pairs = List_reverse(pairs);

    /* queryjump = (revoffset1-bestrR) - (offset1+bestrL) + 1; */
    /* genomejump = (revoffset2R-bestcR) - (offset2L+bestcL) + 1; */
    /* No need to revise queryjump or genomejump, because the above
       coordinates are internal to the gap. */

    debug(printf("Pushing a gap\n"));
#ifndef NOGAPHOLDER
    /* pairs = Pairpool_push_gapholder(pairs,pairpool,queryjump,genomejump,false); */
    pairs = Pairpool_push_gapholder(pairs,pairpool,/*queryjump*/UNKNOWNJUMP,/*genomejump*/UNKNOWNJUMP,/*knownp*/false);
#endif

    pairs = traceback(pairs,&(*nmatches),&(*nmismatches),&(*nopens),&(*nindels),
		      directionsL,bestrL,bestcL,
#ifdef PMAP
		      instsequence1,instsequence1,
#else
		      sequence1,sequenceuc1,
#endif
		      offset1,offset2L,pairpool,/*revp*/false,
		      chroffset,chrhigh,cdna_direction,watsonp,*dynprogindex);

    /*
      Directions_free(directionsR);
      Directions_free(directionsL);
      Matrix_free(matrixR);
      Matrix_free(matrixL);
    */

#ifdef PMAP
    FREE(instsequence1);
#endif

    if (List_length(pairs) == 1) {
      /* Only a gap inserted */
      pairs = (List_T) NULL;
    }

    debug(printf("End of dynprog genome gap\n"));

    *dynprogindex += (*dynprogindex > 0 ? +1 : -1);
    return List_reverse(pairs);
  }
}





static int
binary_search (int lowi, int highi, Genomicpos_T *positions, Genomicpos_T goal) {
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
		  char *revsequence1, char *revsequenceuc1,
		  int length1, int length2, int revoffset1, int revoffset2, 
		  Genomicpos_T chroffset, Genomicpos_T chrhigh,
#ifdef PMAP
		  char *queryaaseq,
#endif
		  int cdna_direction, bool watsonp, bool jump_late_p, Pairpool_T pairpool,
		  int extraband_end, double defect_rate, Endalign_T endalign) {
  List_T pairs = NULL;
#ifdef PMAP
  char *instsequence1, *instrevsequence1;
#endif
  Pair_T pair;
  Mismatchtype_T mismatchtype;
  struct Int3_T **matrix;
  int open, extend;
  int bestr, bestc;
  struct Direction3_T **directions;
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
  if (length1 <= 0) {
    /* Needed to avoid abort by Matrix_alloc */
    debug6(printf("length1 %d <= 0, so returning NULL\n",length1));
    *nmatches = *nmismatches = *nopens = *nindels = 0;
    *finalscore = 0;
    return (List_T) NULL;
  } else if (endalign == QUERYEND_NOGAPS) {
    /* Don't shorten length1 */
  } else if (length1 > dynprog->maxlength1) {
    debug6(printf("length1 %d is too long.  Chopping to %d\n",length1,dynprog->maxlength1));
    length1 = dynprog->maxlength1;
  }
  if (length2 <= 0) {
    /* Needed to avoid abort by Matrix_alloc */
    debug6(printf("length2 %d <= 0, so returning NULL\n",length2));
    *nmatches = *nmismatches = *nopens = *nindels = 0;
    *finalscore = 0;
    return (List_T) NULL;
  } else if (endalign == QUERYEND_NOGAPS) {
    /* Don't shorten length2 */
  } else if (length2 > dynprog->maxlength2) {
    debug6(printf("length2 %d is too long.  Chopping to %d\n",length2,dynprog->maxlength2));
    length2 = dynprog->maxlength2;
  }

#ifdef PMAP
  instsequence1 = instantiate_codons(&(revsequence1[-length1+1]),queryaaseq,revoffset1-length1+1,length1);
  instrevsequence1 = &(instsequence1[length1-1]);
  debug6(printf("At query offset %d-%d, %.*s\n",revoffset1-length1+1,revoffset1,length1,&(instrevsequence1[-length1+1])));
#else
  debug6(printf("At query offset %d-%d, %.*s\n",revoffset1-length1+1,revoffset1,length1,&(revsequence1[-length1+1])));
#endif

#ifdef EXTRACT_GENOMICSEG
  debug6(printf("At genomic offset %d-%d, %.*s\n",
		revoffset2-length2+1,revoffset2,length2,&(revsequence2[-length2+1])));
#endif


  if (endalign == QUERYEND_GAP || endalign == BEST_LOCAL) {
    matrix = compute_scores_lookup_rev(&directions,dynprog,
#ifdef PMAP
				       instrevsequence1,
#else
				       revsequence1,
#endif
				       revoffset2,length1,length2,
				       chroffset,chrhigh,watsonp,
				       mismatchtype,open,extend,
				       extraband_end,/*widebandp*/true,
				       /*for revp true*/!jump_late_p);
    find_best_endpoint(&(*finalscore),&bestr,&bestc,matrix,length1,length2,extraband_end,
		       !jump_late_p);

  } else if (endalign == QUERYEND_INDELS) {
    matrix = compute_scores_lookup_rev(&directions,dynprog,
#ifdef PMAP
				       instrevsequence1,
#else
				       revsequence1,
#endif
				       revoffset2,length1,length2,
				       chroffset,chrhigh,watsonp,
				       mismatchtype,open,extend,
				       extraband_end,/*widebandp*/true,
				       /*for revp true*/!jump_late_p);
    find_best_endpoint_to_queryend_indels(&(*finalscore),&bestr,&bestc,matrix,length1,length2,extraband_end,
					  !jump_late_p);
    /* *finalscore = 0 -- Splicetrie procedures need to know finalscore */

  } else if (endalign == QUERYEND_NOGAPS) {
    find_best_endpoint_to_queryend_nogaps(&bestr,&bestc,length1,length2);
    /* *finalscore = 0;	-- Splicetrie procedures need to know finalscore */

  } else {
    fprintf(stderr,"Unexpected endalign value %d\n",endalign);
    abort();
  }


#ifdef PMAP
  initpos = revoffset1-(bestc-1);
  debug6(printf("Initial query pos is %d\n",initpos));
  if ((initmod = initpos % 3) > 0) {
    if (bestr + initmod < length1 && bestc + initmod < length2) {
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
			     instrevsequence1,instrevsequence1,
#else
			     revsequence1,revsequenceuc1,
#endif
			     revoffset1,revoffset2,pairpool,chroffset,chrhigh,
			     /*revp*/true,watsonp,*dynprogindex);
    *finalscore = (*nmatches)*FULLMATCH + (*nmismatches)*MISMATCH_ENDQ;

  } else {
    pairs = traceback(NULL,&(*nmatches),&(*nmismatches),&(*nopens),&(*nindels),
		      directions,bestr,bestc,
#ifdef PMAP
		      instrevsequence1,instrevsequence1,
#else
		      revsequence1,revsequenceuc1,
#endif
		      revoffset1,revoffset2,pairpool,/*revp*/true,
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
  FREE(instsequence1);
#endif

  debug6(printf("End of dynprog end5 gap\n\n"));

  *dynprogindex += (*dynprogindex > 0 ? +1 : -1);
  return List_reverse(pairs);
}



/* revsequence2 is the splicejunction */
List_T
Dynprog_end5_splicejunction (int *dynprogindex, int *finalscore, int *nmatches, int *nmismatches, 
			     int *nopens, int *nindels, T dynprog, 
			     char *revsequence1, char *revsequenceuc1,
			     char *revsequence2, char *revsequenceuc2, char *revsequencealt2,
			     int length1, int length2, int revoffset1, int revoffset2_anchor, int revoffset2_far,
			     Genomicpos_T chroffset, Genomicpos_T chrhigh,
#ifdef PMAP
			     char *queryaaseq,
#endif
			     int cdna_direction, bool watsonp, bool jump_late_p, Pairpool_T pairpool,
			     int extraband_end, double defect_rate, int contlength) {
  List_T pairs = NULL;
#ifdef PMAP
  char *instsequence1, *instrevsequence1;
#endif
  Pair_T pair;
  Mismatchtype_T mismatchtype;
  struct Int3_T **matrix;
  int open, extend;
  int bestr, bestc;
  struct Direction3_T **directions;
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
  if (length1 <= 0 || length1 > dynprog->maxlength1) {
    /* Needed to avoid abort by Matrix_alloc */
    debug6(printf("length1 %d <= 0, so returning NULL\n",length1));
    *nmatches = *nmismatches = *nopens = *nindels = 0;
    *finalscore = 0;
    return (List_T) NULL;
  }
  if (length2 <= 0 || length2 > dynprog->maxlength2) {
    /* Needed to avoid abort by Matrix_alloc */
    debug6(printf("length2 %d <= 0, so returning NULL\n",length2));
    *nmatches = *nmismatches = *nopens = *nindels = 0;
    *finalscore = 0;
    return (List_T) NULL;
  }

#ifdef PMAP
  instsequence1 = instantiate_codons(&(revsequence1[-length1+1]),queryaaseq,revoffset1-length1+1,length1);
  instrevsequence1 = &(instsequence1[length1-1]);
  debug6(printf("At query offset %d-%d, %.*s\n",revoffset1-length1+1,revoffset1,length1,&(instrevsequence1[-length1+1])));
#else
  debug6(printf("At query offset %d-%d, %.*s\n",revoffset1-length1+1,revoffset1,length1,&(revsequence1[-length1+1])));
#endif

  debug6(printf("At genomic offset %d-%d, %.*s\n",
		revoffset2_anchor-length2+1,revoffset2_anchor,length2,&(revsequence2[-length2+1])));
  
#ifdef PMAP
  initpos = revoffset1-(bestc-1);
  debug6(printf("Initial query pos is %d\n",initpos));
  if ((initmod = initpos % 3) > 0) {
    if (bestr + initmod < length1 && bestc + initmod < length2) {
      debug6(printf("Rounding down by %d\n",initmod));
      bestr += initmod;
      bestc += initmod;
    }
  }
#endif

  matrix = compute_scores_genomicseg_rev(&directions,dynprog,
#ifdef PMAP
					 instrevsequence1,revsequence2,revsequencealt2,
#else
					 revsequence1,revsequence2,revsequencealt2,
#endif
					 length1,length2,
					 mismatchtype,open,extend,
					 extraband_end,/*widebandp*/true,
					 /*for revp true*/!jump_late_p);
  find_best_endpoint_to_queryend_indels(&(*finalscore),&bestr,&bestc,matrix,length1,length2,extraband_end,
					!jump_late_p);

  *nmatches = *nmismatches = *nopens = *nindels = 0;
  pairs = traceback_local(NULL,&(*nmatches),&(*nmismatches),&(*nopens),&(*nindels),
			  directions,&bestr,&bestc,/*endc*/contlength,
#ifdef PMAP
			  instrevsequence1,instrevsequence1,
#else
			  revsequence1,revsequenceuc1,
#endif
			  revsequence2,revsequenceuc2,revsequencealt2,
			  revoffset1,revoffset2_far,pairpool,/*revp*/true,
			  chroffset,chrhigh,cdna_direction,watsonp,*dynprogindex);

  pairs = Pairpool_push_gapholder(pairs,pairpool,/*queryjump*/0,/*genomejump*/revoffset2_anchor - revoffset2_far,
				  /*knownp*/true);

  pairs = traceback_local(pairs,&(*nmatches),&(*nmismatches),&(*nopens),&(*nindels),
			  directions,&bestr,&bestc,/*endc*/0,
#ifdef PMAP
			  instrevsequence1,instrevsequence1,
#else
			  revsequence1,revsequenceuc1,
#endif
			  revsequence2,revsequenceuc2,revsequencealt2,
			  revoffset1,revoffset2_anchor,pairpool,/*revp*/true,
			  chroffset,chrhigh,cdna_direction,watsonp,*dynprogindex);

  /* Score compared with perfect score, so heavy weight on mismatches may not be necessary */
  *finalscore = (*nmatches)*FULLMATCH + (*nmismatches)*MISMATCH_ENDQ + (*nopens)*open + (*nindels)*extend;

  /* Add 1 to count the match already in the alignment */
  pairs = List_reverse(pairs); /* Look at 5' end to remove excess gaps */
  while (pairs != NULL && (pair = List_head(pairs)) && pair->comp == INDEL_COMP) {
    pairs = List_next(pairs);
  }

#ifdef PMAP
  FREE(instsequence1);
#endif

  debug6(Pair_dump_list(pairs,true));
  debug6(printf("End of dynprog end5 gap splicejunction\n\n"));

  *dynprogindex += (*dynprogindex > 0 ? +1 : -1);
  return List_reverse(pairs);
}



List_T
Dynprog_end3_gap (int *dynprogindex, int *finalscore, int *nmatches, int *nmismatches, 
		  int *nopens, int *nindels, T dynprog, 
		  char *sequence1, char *sequenceuc1,
		  int length1, int length2, int offset1, int offset2, 
		  Genomicpos_T chroffset, Genomicpos_T chrhigh,
#ifdef PMAP
		  char *queryaaseq,
#endif
		  int cdna_direction, bool watsonp, bool jump_late_p, Pairpool_T pairpool,
		  int extraband_end, double defect_rate, Endalign_T endalign) {
  List_T pairs = NULL;
#ifdef PMAP
  char *instsequence1;
#endif
  Pair_T pair;
  Mismatchtype_T mismatchtype;
  struct Int3_T **matrix;
  int open, extend;
  int bestr, bestc;
  struct Direction3_T **directions;
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
  if (length1 <= 0) {
    /* Needed to avoid abort by Matrix_alloc */
    *nmatches = *nmismatches = *nopens = *nindels = 0;
    *finalscore = 0;
    return (List_T) NULL;
  } else if (endalign == QUERYEND_NOGAPS) {
    /* Don't shorten length1 */
  } else if (length1 > dynprog->maxlength1) {
    debug6(printf("length1 %d is too long.  Chopping to %d\n",length1,dynprog->maxlength1));
    length1 = dynprog->maxlength1;
  }
  if (length2 <= 0) {
    /* Needed to avoid abort by Matrix_alloc */
    *nmatches = *nmismatches = *nopens = *nindels = 0;
    *finalscore = 0;
    return (List_T) NULL;
  } else if (endalign == QUERYEND_NOGAPS) {
    /* Don't shorten length2 */
  } else if (length2 > dynprog->maxlength2) {
    debug6(printf("length2 %d is too long.  Chopping to %d\n",length2,dynprog->maxlength2));
    length2 = dynprog->maxlength2;
  }

#ifdef PMAP
  instsequence1 = instantiate_codons(sequence1,queryaaseq,offset1,length1);
  debug6(printf("At query offset %d-%d, %.*s\n",offset1,offset1+length1-1,length1,instsequence1));
#else
  debug6(printf("At query offset %d-%d, %.*s\n",offset1,offset1+length1-1,length1,sequence1));
#endif
	
#ifdef EXTRACT_GENOMICSEG
  debug6(printf("At genomic offset %d-%d, %.*s\n",offset2,offset2+length2-1,length2,sequence2));
#endif


  if (endalign == QUERYEND_GAP || endalign == BEST_LOCAL) {
    matrix = compute_scores_lookup_fwd(&directions,dynprog,
#ifdef PMAP
				       instsequence1,
#else
				       sequence1,
#endif
				       offset2,length1,length2,
				       chroffset,chrhigh,watsonp,
				       mismatchtype,open,extend,
				       extraband_end,/*widebandp*/true,
				       jump_late_p);
    find_best_endpoint(&(*finalscore),&bestr,&bestc,matrix,length1,length2,extraband_end,
		       jump_late_p);

  } else if (endalign == QUERYEND_INDELS) {
    matrix = compute_scores_lookup_fwd(&directions,dynprog,
#ifdef PMAP
				       instsequence1,
#else
				       sequence1,
#endif
				       offset2,length1,length2,
				       chroffset,chrhigh,watsonp,
				       mismatchtype,open,extend,
				       extraband_end,/*widebandp*/true,
				       jump_late_p);
    find_best_endpoint_to_queryend_indels(&(*finalscore),&bestr,&bestc,matrix,length1,length2,extraband_end,
					  jump_late_p);
    /* *finalscore = 0; -- Splicetrie procedures need to know finalscore */

  } else if (endalign == QUERYEND_NOGAPS) {
    find_best_endpoint_to_queryend_nogaps(&bestr,&bestc,length1,length2);
    /* *finalscore = 0; -- Splicetrie procedures need to know finalscore */

  } else {
    fprintf(stderr,"Unexpected endalign value %d\n",endalign);
    abort();
  }

#ifdef PMAP
  termpos = offset1+(bestc-1);
  debug6(printf("Final query pos is %d\n",termpos));
  if ((termmod = termpos % 3) < 2) {
    if (bestr + (2 - termmod) < length1 && bestc + (2 - termmod) < length2) {
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
			     instsequence1,instsequence1,
#else
			     sequence1,sequenceuc1,
#endif
			     offset1,offset2,pairpool,chroffset,chrhigh,
			     /*revp*/false,watsonp,*dynprogindex);
    *finalscore = (*nmatches)*FULLMATCH + (*nmismatches)*MISMATCH_ENDQ;

  } else {
    pairs = traceback(NULL,&(*nmatches),&(*nmismatches),&(*nopens),&(*nindels),
		      directions,bestr,bestc,
#ifdef PMAP
		      instsequence1,instsequence1,
#else
		      sequence1,sequenceuc1,
#endif
		      offset1,offset2,pairpool,/*revp*/false,
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
  FREE(instsequence1);
#endif

  debug6(printf("End of dynprog end3 gap\n\n"));

  *dynprogindex += (*dynprogindex > 0 ? +1 : -1);
  return pairs;			/* not List_reverse(pairs) */
}


/* sequence2 is the splicejunction */
List_T
Dynprog_end3_splicejunction (int *dynprogindex, int *finalscore, int *nmatches, int *nmismatches, 
			     int *nopens, int *nindels, T dynprog, 
			     char *sequence1, char *sequenceuc1,
			     char *sequence2, char *sequenceuc2, char *sequencealt2,
			     int length1, int length2, int offset1, int offset2_anchor, int offset2_far,
			     Genomicpos_T chroffset, Genomicpos_T chrhigh,
#ifdef PMAP
			     char *queryaaseq,
#endif
			     int cdna_direction, bool watsonp, bool jump_late_p, Pairpool_T pairpool,
			     int extraband_end, double defect_rate, int contlength) {
  List_T pairs = NULL;
#ifdef PMAP
  char *instsequence1;
#endif
  Pair_T pair;
  Mismatchtype_T mismatchtype;
  struct Int3_T **matrix;
  int open, extend;
  int bestr, bestc;
  struct Direction3_T **directions;
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
  if (length1 <= 0 || length1 > dynprog->maxlength1) {
    /* Needed to avoid abort by Matrix_alloc */
    *nmatches = *nmismatches = *nopens = *nindels = 0;
    *finalscore = 0;
    return (List_T) NULL;
  }
  if (length2 <= 0 || length2 > dynprog->maxlength2) {
    /* Needed to avoid abort by Matrix_alloc */
    *nmatches = *nmismatches = *nopens = *nindels = 0;
    *finalscore = 0;
    return (List_T) NULL;
  }

#ifdef PMAP
  instsequence1 = instantiate_codons(sequence1,queryaaseq,offset1,length1);
  debug6(printf("At query offset %d-%d, %.*s\n",offset1,offset1+length1-1,length1,instsequence1));
#else
  debug6(printf("At query offset %d-%d, %.*s\n",offset1,offset1+length1-1,length1,sequence1));
#endif
	
  debug6(printf("At genomic offset %d-%d, %.*s\n",
		offset2_anchor,offset2_anchor+length2-1,length2,sequence2));


  /* find_best_endpoint_to_queryend_nogaps(&bestr,&bestc,length1,length2); */
  /* bestr = bestc = length1; */
  /* *finalscore = 0; -- Splicetrie procedures need to know finalscore */

#ifdef PMAP
  termpos = offset1+(bestc-1);
  debug6(printf("Final query pos is %d\n",termpos));
  if ((termmod = termpos % 3) < 2) {
    if (bestr + (2 - termmod) < length1 && bestc + (2 - termmod) < length2) {
      debug6(printf("Rounding up by %d\n",2 - termmod));
      bestr += 2 - termmod;
      bestc += 2 - termmod;
    }
  }
#endif

  matrix = compute_scores_genomicseg_fwd(&directions,dynprog,
#ifdef PMAP
					 instsequence1,sequence2,sequencealt2,
#else
					 sequence1,sequence2,sequencealt2,
#endif
					 length1,length2,
					 mismatchtype,open,extend,
					 extraband_end,/*widebandp*/true,
					 jump_late_p);
  find_best_endpoint_to_queryend_indels(&(*finalscore),&bestr,&bestc,matrix,length1,length2,extraband_end,
					jump_late_p);

  *nmatches = *nmismatches = *nopens = *nindels = 0;
  pairs = traceback_local(NULL,&(*nmatches),&(*nmismatches),&(*nopens),&(*nindels),
			  directions,&bestr,&bestc,/*endc*/contlength,
#ifdef PMAP
			  instsequence1,instsequence1,
#else
			  sequence1,sequenceuc1,
#endif
			  sequence2,sequenceuc2,sequencealt2,
			  offset1,offset2_far,pairpool,/*revp*/false,
			  chroffset,chrhigh,cdna_direction,watsonp,*dynprogindex);

  pairs = Pairpool_push_gapholder(pairs,pairpool,/*queryjump*/0,/*genomejump*/offset2_far - offset2_anchor,
				  /*knownp*/true);

  pairs = traceback_local(pairs,&(*nmatches),&(*nmismatches),&(*nopens),&(*nindels),
			  directions,&bestr,&bestc,/*endc*/0,
#ifdef PMAP
			  instsequence1,instsequence1,
#else
			  sequence1,sequenceuc1,
#endif
			  sequence2,sequenceuc2,sequencealt2,
			  offset1,offset2_anchor,pairpool,/*revp*/false,
			  chroffset,chrhigh,cdna_direction,watsonp,*dynprogindex);

  /* Score compared with perfect score, so heavy weight on mismatches may not be necessary */
  *finalscore = (*nmatches)*FULLMATCH + (*nmismatches)*MISMATCH_ENDQ + (*nopens)*open + (*nindels)*extend;

  /* Add 1 to count the match already in the alignment */
  pairs = List_reverse(pairs); /* Look at 3' end to remove excess gaps */
  while (pairs != NULL && (pair = List_head(pairs)) && pair->comp == INDEL_COMP) {
    pairs = List_next(pairs);
  }

#ifdef PMAP
  FREE(instsequence1);
#endif

  debug6(Pair_dump_list(pairs,true));
  debug6(printf("End of dynprog end3 gap splicejunction\n\n"));

  *dynprogindex += (*dynprogindex > 0 ? +1 : -1);
  return pairs;			/* not List_reverse(pairs) */
}


static void
make_contjunction_5 (char *splicejunction, char *splicejunction_alt, Genomicpos_T splicecoord,
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
    printf(" (fwd)  contjunction: %.*s alt:%.*s\n",contlength,proximal,contlength,proximal_alt);
  } else {
    printf(" (rev)  contjunction: %.*s alt:%.*s\n",contlength,proximal,contlength,proximal_alt);
  }
#endif

  return;
}



/* Fills in just the distal part, keeping the proximal part same for contlength */
void
Dynprog_make_splicejunction_5 (char *splicejunction, char *splicejunction_alt, Genomicpos_T splicecoord,
			       int splicelength, int contlength, Splicetype_T far_splicetype,
			       bool watsonp) {
  char *distal, *distal_alt;

  debug7(printf("make_splicejunction_5 at %u, splice (%s), contlength %d, splicelength %d:",
		splicecoord,Splicetype_string(far_splicetype),contlength, splicelength));

  distal = &(splicejunction[0]);
  distal_alt = &(splicejunction_alt[0]);

  if (far_splicetype == ACCEPTOR) {
    Genome_fill_buffer_blocks_noterm(splicecoord,splicelength,distal,distal_alt);

  } else if (far_splicetype == ANTIDONOR) {
    Genome_fill_buffer_blocks_noterm(splicecoord,splicelength,distal,distal_alt);

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
    printf(" (fwd)  splicejunction: %s\n",splicejunction);
  } else {
    printf(" (rev)  splicejunction: %s\n",splicejunction);
  }
#endif

  return;
}


static void
make_contjunction_3 (char *splicejunction, char *splicejunction_alt, Genomicpos_T splicecoord,
		     int splicelength, int contlength, Splicetype_T anchor_splicetype,
		     bool watsonp) {
  char *proximal, *proximal_alt;

  debug7(printf("make_contjunction_3 at %u, splice (%s), contlength %d, splicelength %d:",
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
    printf(" (fwd)  contjunction: %.*s\n",contlength,proximal);
  } else {
    printf(" (rev)  contjunction: %.*s\n",contlength,proximal);
  }
#endif

  return;
}



/* Fills in just the distal part, keeping the proximal part same for contlength */
void
Dynprog_make_splicejunction_3 (char *splicejunction, char *splicejunction_alt, Genomicpos_T splicecoord,
			       int splicelength, int contlength, Splicetype_T far_splicetype,
			       bool watsonp) {
  char *distal, *distal_alt;

  debug7(printf("make_splicejunction_3 at %u, splice (%s), contlength %d, splicelength %d:",
		splicecoord,Splicetype_string(far_splicetype),contlength,splicelength));

  distal = &(splicejunction[contlength]);
  distal_alt = &(splicejunction_alt[contlength]);

  if (far_splicetype == DONOR) {
    Genome_fill_buffer_blocks_noterm(splicecoord-splicelength,splicelength,distal,distal_alt);

  } else if (far_splicetype == ANTIACCEPTOR) {
    Genome_fill_buffer_blocks_noterm(splicecoord-splicelength,splicelength,distal,distal_alt);

  } else if (far_splicetype == ANTIDONOR) {
    Genome_fill_buffer_blocks_noterm(splicecoord,splicelength,distal,distal_alt);

  } else if (far_splicetype == ACCEPTOR) {
    Genome_fill_buffer_blocks_noterm(splicecoord,splicelength,distal,distal_alt);
    
  } else {
    fprintf(stderr,"Unexpected far_splicetype value %d\n",far_splicetype);
    abort();
  }

  if (watsonp == false) {
    make_complement_inplace(distal,splicelength);
  }

#ifdef DEBUG7
  if (watsonp == true) {
    printf(" (fwd)  splicejunction: %s\n",splicejunction);
  } else {
    printf(" (rev)  splicejunction: %s\n",splicejunction);
  }
#endif

  return;
}


List_T
Dynprog_end5_known (bool *knownsplicep, int *dynprogindex, int *finalscore,
		    int *ambig_end_length, Splicetype_T *ambig_splicetype,
		    int *nmatches, int *nmismatches, int *nopens, int *nindels, T dynprog, 
		    char *revsequence1, char *revsequenceuc1,
		    int length1, int length2, int revoffset1, int revoffset2, 
		    Genomicpos_T chroffset, Genomicpos_T chrhigh,
		    Genomicpos_T knownsplice_limit_low, Genomicpos_T knownsplice_limit_high,
#ifdef PMAP
		    char *queryaaseq,
#endif
		    int cdna_direction, bool watsonp, bool jump_late_p, Pairpool_T pairpool,
		    int extraband_end, double defect_rate) {
  List_T best_pairs = NULL, orig_pairs;
  Pair_T pair;
  Genomicpos_T low, high, far_limit_low, far_limit_high;
  Splicetype_T anchor_splicetype, far_splicetype;
  int contlength, splicelength, endlength;
  char *splicejunction, *splicejunction_alt;
#ifdef EXTRACT_GENOMICSEG
  char *splicejunction_test;
#endif
  int jstart, j;

  int orig_score, threshold_miss_score, perfect_score;
  int obsmax_penalty;


  assert(length2 >= length1);

  debug7(
	printf("%c:  ",*dynprogindex > 0 ? (*dynprogindex-1)%26+'a' : (-(*dynprogindex)-1)%26+'A');
	printf("Aligning 5' end gap, known\n")
	);

  *ambig_end_length = 0;

  /* We can just chop lengths to work, since we're not constrained on 5' end */
  if (length1 <= 0) {
    /* Needed to avoid abort by Matrix_alloc */
    *finalscore = 0;
    *knownsplicep = false;
    return (List_T) NULL;
  }
  if (length2 <= 0) {
    /* Needed to avoid abort by Matrix_alloc */
    *finalscore = 0;
    *knownsplicep = false;
    return (List_T) NULL;
  }

  debug7(printf("At query offset %d-%d, %.*s\n",revoffset1-length1+1,revoffset1,length1,&(revsequence1[-length1+1])));
#ifdef EXTRACT_GENOMICSEG
  debug7(printf("At genomic offset %d-%d, %.*s\n",revoffset2-length2+1,revoffset2,length2,&(revsequence2[-length2+1])));
#endif

  perfect_score = length1*FULLMATCH;

  /* Try without splicing, all the way to query end */
  best_pairs = Dynprog_end5_gap(&(*dynprogindex),&(*finalscore),&(*nmatches),&(*nmismatches),
				&(*nopens),&(*nindels),dynprog,revsequence1,revsequenceuc1,
				length1,length2,revoffset1,revoffset2,chroffset,chrhigh,
#ifdef PMAP
				queryaaseq,
#endif
				cdna_direction,watsonp,jump_late_p,pairpool,
				extraband_end,defect_rate,/*endalign*/QUERYEND_NOGAPS);
  orig_score = *finalscore;
  orig_pairs = best_pairs;
  threshold_miss_score = perfect_score - orig_score;
  debug7(printf("perfect score %d - score %d = threshold %d\n",
		perfect_score,orig_score,threshold_miss_score));
  *knownsplicep = false;


  if (threshold_miss_score > 0 && length2 > 0) {
    /* Try known splicing */
    splicejunction = (char *) CALLOC(length2+1,sizeof(char));
    splicejunction_alt = (char *) CALLOC(length2+1,sizeof(char));
#ifdef EXTRACT_GENOMICSEG
    splicejunction_test = (char *) CALLOC(length2+1,sizeof(char));
#endif

    endlength = length1;
    if (watsonp == true) {
      low = chroffset + revoffset2-endlength + 2;
      high = chroffset + revoffset2 + 1;
      debug7(printf("5' watson\n"));
      debug7(printf("Calculating low %u (%u) = %u + %d-%d + 2\n",
		    low,low-chroffset,chroffset,revoffset2,endlength));
      debug7(printf("Calculating high %u (%u) = %u + %d + 1\n",
		    high,high-chroffset,chroffset,revoffset2));
      if (cdna_direction > 0) {
	anchor_splicetype = ACCEPTOR;
	far_splicetype = DONOR;
      } else {
	anchor_splicetype = ANTIDONOR;
	far_splicetype = ANTIACCEPTOR;
      }
    } else {
      low = chrhigh - revoffset2;
      high = chrhigh - (revoffset2-endlength) - 1;
      debug7(printf("5' crick\n"));
      debug7(printf("Calculating low %u (%u) = %u - %d\n",
		    low,low-chroffset,chrhigh,revoffset2));
      debug7(printf("Calculating high %u (%u) = %u - (%d-%d) - 1\n",
		    high,high-chroffset,chrhigh,revoffset2,endlength));
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
	debug7(printf("contlength %d, splicelength %d, length1 %d, length2 %d\n",
		      contlength,length2-contlength,length1,length2));
	assert(contlength >= 0 && contlength < length1);

#ifdef EXTRACT_GENOMICSEG
	debug7(printf("cont: %.*s\n",contlength,&(revsequence2[-contlength+1])));
#endif
	splicelength = length2 - contlength;
	assert(splicelength > 0);
	debug7(printf("  Saw %u (%u) of type %s (cont length %d, splice length %d)\n",
		      splicesites[j],splicesites[j]-chroffset,Splicetype_string(splicetypes[j]),contlength,splicelength));

	make_contjunction_5(splicejunction,splicejunction_alt,splicesites[j],splicelength,contlength,anchor_splicetype,watsonp);
#ifdef EXTRACT_GENOMICSEG
	strncpy(&(splicejunction_test[splicelength]),&(revsequence2[-contlength+1]),contlength);
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
	  debug7(printf("  Running Splicetrie_solve_end5 on observed splice sites with revoffset2 %d\n",revoffset2));
	  best_pairs = Splicetrie_solve_end5(best_pairs,triecontents_obs,trieoffsets_obs,j,
					     far_limit_low,far_limit_high,
					     &(*finalscore),&(*nmatches),&(*nmismatches),
					     &(*nopens),&(*nindels),&(*knownsplicep),&(*ambig_end_length),
					     &threshold_miss_score,/*obsmax_penalty*/0,perfect_score,
					     /*anchor_splicesite*/splicesites[j],splicejunction,splicejunction_alt,
					     splicelength,contlength,far_splicetype,
					     chroffset,chrhigh,&(*dynprogindex),dynprog,
					     revsequence1,revsequenceuc1,length1,length2,revoffset1,revoffset2,
#ifdef PMAP
					     queryaaseq,
#endif
					     cdna_direction,watsonp,jump_late_p,pairpool,extraband_end,defect_rate);
	  debug7(printf("  Result\n"));
	  debug7(Pair_dump_list(best_pairs,/*zerobasedp*/true));
	  obsmax_penalty += FULLMATCH;
	}

	if (threshold_miss_score - obsmax_penalty > 0 && trieoffsets_max != NULL) {
	  debug7(printf("  Running Splicetrie_solve_end5 on maxdistance splice sites with revoffset2 %d\n",revoffset2));
	  best_pairs = Splicetrie_solve_end5(best_pairs,triecontents_max,trieoffsets_max,j,
					     far_limit_low,far_limit_high,
					     &(*finalscore),&(*nmatches),&(*nmismatches),
					     &(*nopens),&(*nindels),&(*knownsplicep),&(*ambig_end_length),
					     &threshold_miss_score,obsmax_penalty,perfect_score,
					     /*anchor_splicesite*/splicesites[j],splicejunction,splicejunction_alt,
					     splicelength,contlength,far_splicetype,
					     chroffset,chrhigh,&(*dynprogindex),dynprog,
					     revsequence1,revsequenceuc1,length1,length2,revoffset1,revoffset2,
#ifdef PMAP
					     queryaaseq,
#endif
					     cdna_direction,watsonp,jump_late_p,pairpool,extraband_end,defect_rate);
	  debug7(printf("  Result\n"));
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
      if (length1 > dynprog->maxlength1) {
	debug7(printf("length1 %d is too long.  Chopping to %d\n",length1,dynprog->maxlength1));
	length1 = dynprog->maxlength1;
      }
      if (length2 > dynprog->maxlength2) {
	debug7(printf("length2 %d is too long.  Chopping to %d\n",length2,dynprog->maxlength2));
	length2 = dynprog->maxlength2;
      }
      orig_pairs = Dynprog_end5_gap(&(*dynprogindex),&(*finalscore),&(*nmatches),&(*nmismatches),
				    &(*nopens),&(*nindels),dynprog,revsequence1,revsequenceuc1,
				    length1,length2,revoffset1,revoffset2,chroffset,chrhigh,
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
		    char *sequence1, char *sequenceuc1,
		    int length1, int length2, int offset1, int offset2, int querylength,
		    Genomicpos_T chroffset, Genomicpos_T chrhigh,
		    Genomicpos_T knownsplice_limit_low, Genomicpos_T knownsplice_limit_high,
#ifdef PMAP
		    char *queryaaseq,
#endif
		    int cdna_direction, bool watsonp, bool jump_late_p, Pairpool_T pairpool,
		    int extraband_end, double defect_rate) {
  List_T best_pairs = NULL, orig_pairs;
  Pair_T pair;
  Genomicpos_T low, high, far_limit_low, far_limit_high;
  Splicetype_T anchor_splicetype, far_splicetype;
  int contlength, splicelength, endlength;
  char *splicejunction, *splicejunction_alt;
#ifdef EXTRACT_GENOMICSEG
  char *splicejunction_test;
#endif
  int jstart, j;

  int orig_score, threshold_miss_score, perfect_score;
  int obsmax_penalty;


  assert(length2 >= length1);

  debug7(
	printf("%c:  ",*dynprogindex > 0 ? (*dynprogindex-1)%26+'a' : (-(*dynprogindex)-1)%26+'A');
	printf("Aligning 3' end gap, known\n")
	);

  *ambig_end_length = 0;

  /* We can just chop lengths to work, since we're not constrained on 3' end */
  if (length1 <= 0) {
    /* Needed to avoid abort by Matrix_alloc */
    *finalscore = 0;
    *knownsplicep = false;
    return (List_T) NULL;
  }
  if (length2 <= 0) {
    /* Needed to avoid abort by Matrix_alloc */
    *finalscore = 0;
    *knownsplicep = false;
    return (List_T) NULL;
  }


  debug7(printf("At query offset %d-%d, %.*s\n",offset1,offset1+length1-1,length1,sequence1));
#ifdef EXTRACT_GENOMICSEG
  debug7(printf("At genomic offset %d-%d, %.*s\n",offset2,offset2+length2-1,length2,sequence2));
#endif

  perfect_score = length1*FULLMATCH;

  /* Try without splicing, all the way to query end */
  best_pairs = Dynprog_end3_gap(&(*dynprogindex),&(*finalscore),&(*nmatches),&(*nmismatches),
				&(*nopens),&(*nindels),dynprog,sequence1,sequenceuc1,
				length1,length2,offset1,offset2,chroffset,chrhigh,
#ifdef PMAP
				queryaaseq,
#endif
				cdna_direction,watsonp,jump_late_p,pairpool,
				extraband_end,defect_rate,/*endalign*/QUERYEND_NOGAPS);
  orig_score = *finalscore;
  orig_pairs = best_pairs;
  threshold_miss_score = perfect_score - orig_score;
  debug7(printf("perfect score %d - score %d = threshold %d\n",
		perfect_score,orig_score,threshold_miss_score));
  *knownsplicep = false;


  if (threshold_miss_score > 0 && length2 > 0) {
    /* Try known splicing */
    splicejunction = (char *) CALLOC(length2+1,sizeof(char));
    splicejunction_alt = (char *) CALLOC(length2+1,sizeof(char));
#ifdef EXTRACT_GENOMICSEG
    splicejunction_test = (char *) CALLOC(length2+1,sizeof(char));
#endif

    endlength = length1;
    if (watsonp == true) {
      low = chroffset + offset2;
      high = chroffset + offset2+endlength - 1;
      debug7(printf("3' watson\n"));
      debug7(printf("Calculating low %u (%u) = %u + %d\n",
		    low,low-chroffset,chroffset,offset2));
      debug7(printf("Calculating high %u (%u) = %u + %d+%d - 1\n",
		    high,high-chroffset,chroffset,offset2,endlength));
      if (cdna_direction > 0) {
	anchor_splicetype = DONOR;
	far_splicetype = ACCEPTOR;
      } else {
	anchor_splicetype = ANTIACCEPTOR;
	far_splicetype = ANTIDONOR;
      }
    } else {
      low = chrhigh - (offset2+endlength) + 2;
      high = chrhigh - offset2 + 1;
      debug7(printf("3' crick\n"));
      debug7(printf("Calculating low %u (%u) = %u - (%d+%d)\n",
		    low,low-chroffset,chrhigh,offset2,endlength));
      debug7(printf("Calculating high %u (%u) = %u - %d + 1\n",
		    high,high-chroffset,chrhigh,offset2));
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
	debug7(printf("contlength %d, splicelength %d, length1 %d, length2 %d\n",
		      contlength,length2-contlength,length1,length2));
	assert(contlength >= 0 && contlength < length1);

#ifdef EXTRACT_GENOMICSEG
	debug7(printf("cont: %.*s\n",contlength,sequence2));
#endif
	splicelength = length2 - contlength;
	assert(splicelength > 0);
	debug7(printf("  Saw %u (%u) of type %s (cont length %d, splice length %d)\n",
		      splicesites[j],splicesites[j]-chroffset,Splicetype_string(splicetypes[j]),contlength,splicelength));

	make_contjunction_3(splicejunction,splicejunction_alt,splicesites[j],splicelength,contlength,anchor_splicetype,watsonp);
#ifdef EXTRACT_GENOMICSEG
	strncpy(splicejunction_test,sequence2,contlength);
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
	  debug7(printf("  Running Splicetrie_solve_end3 on observed splice sites with offset2 %d\n",offset2));
	  best_pairs = Splicetrie_solve_end3(best_pairs,triecontents_obs,trieoffsets_obs,j,
					     far_limit_low,far_limit_high,
					     &(*finalscore),&(*nmatches),&(*nmismatches),
					     &(*nopens),&(*nindels),&(*knownsplicep),&(*ambig_end_length),
					     &threshold_miss_score,/*obsmax_penalty*/0,perfect_score,
					     /*anchor_splicesite*/splicesites[j],splicejunction,splicejunction_alt,
					     splicelength,contlength,far_splicetype,
					     chroffset,chrhigh,&(*dynprogindex),dynprog,
					     sequence1,sequenceuc1,length1,length2,offset1,offset2,
#ifdef PMAP
					     queryaaseq,
#endif
					     cdna_direction,watsonp,jump_late_p,pairpool,extraband_end,defect_rate);
	  debug7(printf("  Result\n"));
	  debug7(Pair_dump_list(best_pairs,/*zerobasedp*/true));
	  obsmax_penalty += FULLMATCH;
	}

	if (threshold_miss_score - obsmax_penalty > 0 && trieoffsets_max != NULL) {
	  debug7(printf("  Running Splicetrie_solve_end3 on maxdistance splice sites with offset2 %d\n",offset2));
	  best_pairs = Splicetrie_solve_end3(best_pairs,triecontents_max,trieoffsets_max,j,
					     far_limit_low,far_limit_high,
					     &(*finalscore),&(*nmatches),&(*nmismatches),
					     &(*nopens),&(*nindels),&(*knownsplicep),&(*ambig_end_length),
					     &threshold_miss_score,obsmax_penalty,perfect_score,
					     /*anchor_splicesite*/splicesites[j],splicejunction,splicejunction_alt,
					     splicelength,contlength,far_splicetype,
					     chroffset,chrhigh,&(*dynprogindex),dynprog,
					     sequence1,sequenceuc1,length1,length2,offset1,offset2,
#ifdef PMAP
					     queryaaseq,
#endif
					     cdna_direction,watsonp,jump_late_p,pairpool,extraband_end,defect_rate);
	  debug7(printf("  Result\n"));
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
      if (length1 > dynprog->maxlength1) {
	debug7(printf("length1 %d is too long.  Chopping to %d\n",length1,dynprog->maxlength1));
	length1 = dynprog->maxlength1;
      }
      if (length2 > dynprog->maxlength2) {
	debug7(printf("length2 %d is too long.  Chopping to %d\n",length2,dynprog->maxlength2));
	length2 = dynprog->maxlength2;
      }
      orig_pairs = Dynprog_end3_gap(&(*dynprogindex),&(*finalscore),&(*nmatches),&(*nmismatches),
				    &(*nopens),&(*nindels),dynprog,sequence1,sequenceuc1,
				    length1,length2,offset1,offset2,chroffset,chrhigh,
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
make_microexon_pairs_double (int offset1L, int offset1M, int offset1R,
			     int offset2L, int offset2M, int offset2R,
			     int lengthL, int lengthM, int lengthR,
			     char *queryseq, char *queryuc,
			     Genomicpos_T chroffset, Genomicpos_T chrhigh, bool watsonp,
			     Pairpool_T pairpool, char gapchar, int dynprogindex) {
  List_T pairs = NULL;
  Pair_T gappair;
  char c1, c2, c2_alt;
  int i;

  /* Left segment */
  for (i = 0; i < lengthL; i++) {
    c1 = queryseq[offset1L+i];

    c2 = get_genomic_nt(&c2_alt,offset2L+i,chroffset,chrhigh,watsonp);
#ifdef EXTRACT_GENOMICSEG
    assert(c2 == genomicseg[offset2L+i]);
#endif

    if (queryuc[offset1L+i] == c2) {
      pairs = Pairpool_push(pairs,pairpool,offset1L+i,offset2L+i,c1,DYNPROG_MATCH_COMP,c2,c2_alt,
			    dynprogindex);
    } else if (consistent_array[(int) c1][(int) c2] == true) {
      pairs = Pairpool_push(pairs,pairpool,offset1L+i,offset2L+i,c1,AMBIGUOUS_COMP,c2,c2_alt,
			    dynprogindex);
    } else {
      pairs = Pairpool_push(pairs,pairpool,offset1L+i,offset2L+i,c1,MISMATCH_COMP,c2,c2_alt,
			    dynprogindex);
    }
  }

  /* First gap */
  /* Don't have to adjust querypos/genomepos, since no cdna/genome skips allowed */
#if 0
  pairs = Pairpool_push_gapholder(pairs,pairpool,/*queryjump*/0,offset2M-(offset2L+lengthL),/*knownp*/false);
#endif
  pairs = Pairpool_push_gapholder(pairs,pairpool,/*queryjump*/UNKNOWNJUMP,/*genomejump*/UNKNOWNJUMP,/*knownp*/false);
  
  /* Assign pair->comp, because might occur after assign_gap_types */
  gappair = (Pair_T) List_head(pairs);
  gappair->comp = gapchar;


  /* Microexon */
  for (i = 0; i < lengthM; i++) {
    c1 = queryseq[offset1M+i];

    c2 = get_genomic_nt(&c2_alt,offset2M+i,chroffset,chrhigh,watsonp);
#ifdef EXTRACT_GENOMICSEG
    assert(c2 == genomicseg[offset2M+i]);
#endif

    if (queryuc[offset1M+i] == c2) {
      pairs = Pairpool_push(pairs,pairpool,offset1M+i,offset2M+i,c1,DYNPROG_MATCH_COMP,c2,c2_alt,
			    dynprogindex);
    } else if (consistent_array[(int) c1][(int) c2] == true) {
      pairs = Pairpool_push(pairs,pairpool,offset1M+i,offset2M+i,c1,AMBIGUOUS_COMP,c2,c2_alt,
			    dynprogindex);
    } else {
      pairs = Pairpool_push(pairs,pairpool,offset1M+i,offset2M+i,c1,MISMATCH_COMP,c2,c2_alt,
			    dynprogindex);
    }
  }

  /* Second gap */
  /* Don't have to adjust querypos/genomepos, since no cdna/genome skips allowed */
  if (lengthR == 0) {
    /* If lengthR is zero, then we will have a gap after a gap */
    Except_raise(&microexon_error,__FILE__,__LINE__);
  } else {
#if 0
    pairs = Pairpool_push_gapholder(pairs,pairpool,/*queryjump*/0,offset2R-(offset2M+lengthM),/*knownp*/false);
#endif
    pairs = Pairpool_push_gapholder(pairs,pairpool,/*queryjump*/UNKNOWNJUMP,/*genomejump*/UNKNOWNJUMP,/*knownp*/false);
  }

  /* Assign pair->comp, because might occur after assign_gap_types */
  gappair = (Pair_T) List_head(pairs);
  gappair->comp = gapchar;

  
  /* Right segment */
  for (i = 0; i < lengthR; i++) {
    c1 = queryseq[offset1R+i];

    c2 = get_genomic_nt(&c2_alt,offset2R+i,chroffset,chrhigh,watsonp);
#ifdef EXTRACT_GENOMICSEG
    assert(c2 == genomicseg[offset2R+i]);
#endif

    if (queryuc[offset1R+i] == c2) {
      pairs = Pairpool_push(pairs,pairpool,offset1R+i,offset2R+i,c1,DYNPROG_MATCH_COMP,c2,c2_alt,
			    dynprogindex);
    } else if (consistent_array[(int) c1][(int) c2] == true) {
      pairs = Pairpool_push(pairs,pairpool,offset1R+i,offset2R+i,c1,AMBIGUOUS_COMP,c2,c2_alt,
			    dynprogindex);
    } else {
      pairs = Pairpool_push(pairs,pairpool,offset1R+i,offset2R+i,c1,MISMATCH_COMP,c2,c2_alt,
			    dynprogindex);
    }
  }

  return pairs;
}


#if 0
static List_T
make_microexon_pairs_single (int offset1L, int offset1R,
			     int offset2L, int offset2R,
			     int lengthL, int lengthR, char *queryseq, char *queryuc,
			     Genomicpos_T chroffset, Genomicpos_T chrhigh, bool watsonp,
			     Pairpool_T pairpool, char gapchar, int dynprogindex) {
  List_T pairs = NULL;
  Pair_T gappair;
  char c1, c2, c2_alt;
  int i;

  /* Microexon */
  for (i = 0; i < lengthL; i++) {
    c1 = queryseq[offset1L+i];

    c2 = get_genomic_nt(&c2_alt,offset2L+i,chroffset,chrhigh,watsonp);
#ifdef EXTRACT_GENOMICSEG
    assert(c2 == genomicseg[offset2L+i]);
#endif

    if (queryuc[offset1L+i] == c2) {
      pairs = Pairpool_push(pairs,pairpool,offset1L+i,offset2L+i,c1,DYNPROG_MATCH_COMP,c2,c2_alt,
			    dynprogindex);
    } else if (consistent_array[(int) c1][(int) c2] == true) {
      pairs = Pairpool_push(pairs,pairpool,offset1L+i,offset2L+i,c1,AMBIGUOUS_COMP,c2,c2_alt,
			    dynprogindex);
    } else {
      pairs = Pairpool_push(pairs,pairpool,offset1L+i,offset2L+i,c1,MISMATCH_COMP,c2,c2_alt,
			    dynprogindex);
    }
  }

  /* Gap */
  /* Don't have to adjust querypos/genomepos, since no cdna/genome skips allowed */
#if 0
  pairs = Pairpool_push_gapholder(pairs,pairpool,/*queryjump*/0,offset2R-(offset2L+lengthL),/*knownp*/false);
#endif
  pairs = Pairpool_push_gapholder(pairs,pairpool,/*queryjump*/UNKNOWNJUMP,/*genomejump*/UNKNOWNJUMP,/*knownp*/false);

  /* Assign pair->comp, because might occur after assign_gap_types */
  gappair = (Pair_T) List_head(pairs);
  gappair->comp = gapchar;
  

  /* Right segment */
  for (i = 0; i < lengthR; i++) {
    c1 = queryseq[offset1R+i];

    c2 = get_genomic_nt(&c2_alt,offset2R+i,chroffset,chrhigh,watsonp);
#ifdef EXTRACT_GENOMICSEG
    assert(c2 == genomicseg[offset2R+i]);
#endif

    if (queryuc[offset1R+i] == c2) {
      pairs = Pairpool_push(pairs,pairpool,offset1R+i,offset2R+i,c1,DYNPROG_MATCH_COMP,c2,c2_alt,
			    dynprogindex);
    } else if (consistent_array[(int) c1][(int) c2] == true) {
      pairs = Pairpool_push(pairs,pairpool,offset1R+i,offset2R+i,c1,AMBIGUOUS_COMP,c2,c2_alt,
			    dynprogindex);
    } else {
      pairs = Pairpool_push(pairs,pairpool,offset1R+i,offset2R+i,c1,MISMATCH_COMP,c2,c2_alt,
			    dynprogindex);
    }
  }

  return pairs;
}
#endif


List_T
Dynprog_microexon_int (double *bestprob2, double *bestprob3, int *dynprogindex, int *microintrontype,
		       char *sequence1, char *sequenceuc1,
		       int length1, int length2L, int length2R,
		       int offset1, int offset2L, int revoffset2R, int cdna_direction,
#ifdef PMAP
		       char *queryaaseq, char *genomicuc,
#endif
		       char *queryseq, char *queryuc,
		       Genomicpos_T chroffset, Genomicpos_T chrhigh, bool watsonp,
		       Pairpool_T pairpool, double defect_rate) {
  List_T pairs = NULL;
#ifdef PMAP
  char *instsequence1;
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
  double pvalue, bestprob = 0.0, prob2, prob3;
  Genomicpos_T splicesitepos;


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
	       length2L,sequence2L,length2R,&(revsequence2R[-length2R+1])));
#else
  debug(printf("Begin microexon search\n"));
#endif

#ifdef PMAP
  instsequence1 = instantiate_codons(sequence1,queryaaseq,offset1,length1);
  debug(printf("  Query sequence is %.*s\n",length1,instsequence1));
#else
  debug(printf("  Query sequence is %.*s\n",length1,sequence1));
#endif
  span = revoffset2R-offset2L;
  debug(printf("  Genomic span is of length %d\n",span));

  if (span <= 0) {
    fprintf(stderr,"Bug in Dynprog_microexon_int.  span %d <= 0.  Please report to twu@gene.com\n",span);
    abort();
  } else {
    min_microexon_length = ceil(-log(1.0-pow(1.0-pvalue,1.0/(double) span))/log(4));
  }
  min_microexon_length -= 8;	/* Two donor-acceptor pairs */
  debug(printf("  Min microexon length is %d\n",min_microexon_length));
  if (min_microexon_length > MAX_MICROEXON_LENGTH) {
#ifdef PMAP
    FREE(instsequence1);
#endif
    *microintrontype = NONINTRON;
    return NULL;
  } else if (min_microexon_length < MIN_MICROEXON_LENGTH) {
    min_microexon_length = MIN_MICROEXON_LENGTH;
  }

  debug(printf("\nFinding starting boundary on left\n"));
  leftbound = 0;
  nmismatches = 0;
  while (leftbound < length1 - 1 && nmismatches <= 1) {
    debug(printf("  leftbound = %d, nmismatches = %d.",leftbound,nmismatches));
#ifdef PMAP
    c = get_genomic_nt(&c_alt,offset2L+leftbound,chroffset,chrhigh,watsonp);
    debug(printf("  Comparing %c with %c\n",instsequence1[leftbound],c));
    if (matchtable[instsequence1[leftbound]-'A'][c-'A'] == false) {
      nmismatches++;
    }
#else
    c = get_genomic_nt(&c_alt,offset2L+leftbound,chroffset,chrhigh,watsonp);
#ifdef EXTRACT_GENOMICSEG
    assert(c == sequenceuc2L[leftbound]);
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
  i = length1-1;
  nmismatches = 0;
  while (i >= 0 && nmismatches <= 1) {
    debug(printf("  rightbound = %d, nmismatches = %d.",rightbound,nmismatches));
#ifdef PMAP
    c = get_genomic_nt(&c_alt,revoffset2R-rightbound,chroffset,chrhigh,watsonp);
    debug(printf("  Comparing %c with %c\n",instsequence1[i],c));
    if (matchtable[instsequence1[i]-'A'][c-'A'] == false) {
      nmismatches++;
    }
#else
    c = get_genomic_nt(&c_alt,revoffset2R-rightbound,chroffset,chrhigh,watsonp);
#ifdef EXTRACT_GENOMICSEG
    assert(c == revsequenceuc2R[-rightbound]);
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
    left1 = get_genomic_nt(&left1_alt,offset2L+cL,chroffset,chrhigh,watsonp);
    left2 = get_genomic_nt(&left2_alt,offset2L+cL+1,chroffset,chrhigh,watsonp);
#ifdef EXTRACT_GENOMICSEG
    assert(left1 == sequenceuc2L[cL]);
    assert(left2 == sequenceuc2L[cL+1]);
#endif

    debug(printf("  %d: %c%c\n",cL,left1,left2));
    if (left1 == intron1 && left2 == intron2) {
      mincR = length1 - MAX_MICROEXON_LENGTH - cL;
      if (mincR < 1) {
	mincR = 1;
      }
      maxcR = length1 - min_microexon_length - cL;
      if (maxcR > rightbound) {
	maxcR = rightbound;
	}
      debug(printf("  Found left GT at %d.  Scanning from %d - cL - (1-7), or %d to %d\n",
		   cL,length1,mincR,maxcR));
      for (cR = mincR; cR <= maxcR; cR++) {
	right2 = get_genomic_nt(&right2_alt,revoffset2R-cR-1,chroffset,chrhigh,watsonp);
	right1 = get_genomic_nt(&right1_alt,revoffset2R-cR,chroffset,chrhigh,watsonp);
#ifdef EXTRACT_GENOMICSEG
	assert(right2 == revsequenceuc2R[-cR-1]);
	assert(right1 == revsequenceuc2R[-cR]);
#endif
	debug(printf("   Checking %d: %c%c\n",cR,right2,right1));
	if (right2 == intron3 && right1 == intron4) {
	  middlelength = length1 - cL - cR;
#ifdef PMAP
	  debug(printf("  Found pair at %d to %d, length %d.  Middle sequence is %.*s\n",
		       cL,cR,middlelength,middlelength,&(instsequence1[cL])));
#else
	  debug(printf("  Found pair at %d to %d, length %d.  Middle sequence is %.*s\n",
		       cL,cR,middlelength,middlelength,&(sequence1[cL])));
#endif
	  
	  textleft = offset2L + cL + MICROINTRON_LENGTH;
	  textright = revoffset2R - cR - MICROINTRON_LENGTH;
#ifdef PMAP
	  hits = BoyerMoore(&(instsequence1[cL]),middlelength,&(genomicuc[textleft]),textright-textleft);
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
	      debug(printf("  Successful microexon at %d >>> %d..%d >>> %d\n",offset2L+cL,candidate,candidate+middlelength,revoffset2R-cR));

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

  if (bestcL < 0 || bestcR < 0) {
    debug(printf("End of dynprog microexon int\n"));

#ifdef PMAP
    FREE(instsequence1);
#endif

    *microintrontype = NONINTRON;
    return NULL;

  } else {
    pairs = make_microexon_pairs_double(offset1,offset1+bestcL,offset1+bestcL+best_middlelength,
					offset2L,candidate,revoffset2R-bestcR+1,
					/*lengthL*/bestcL,/*lengthM*/best_middlelength,/*lengthR*/bestcR,
#ifdef PMAP
					&(instsequence1[-offset1]),&(instsequence1[-offset1]),
#else
					queryseq,queryuc,
#endif
					chroffset,chrhigh,watsonp,pairpool,gapchar,*dynprogindex);
#ifdef PMAP
    FREE(instsequence1);
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
    result = (int) (log(1.0-ENDSEQUENCE_PVALUE)/log(1-p1));
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
		     char *revsequence1, char *revsequenceuc1, char *revsequence2, char *revsequenceuc2,
		     int length1, int length2, int revoffset1, int revoffset2, int cdna_direction,
#ifdef PMAP
		     char *queryaaseq,
#endif
		     char *queryseq, char *queryuc, char *genomicseg, char *genomicuc,
		     Pairpool_T pairpool, bool end_microexons_p) {
  List_T pairs = NULL;
#ifdef PMAP
  char *instsequence1, *instrevsequence1;
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
	       length2,&(revsequence2[-length2+1])));
#else
  debug(printf("Begin microexon search at 5'\n"));
#endif

#ifdef PMAP
  instsequence1 = instantiate_codons(&(revsequence1[-length1+1]),queryaaseq,revoffset1-length1+1,length1);
  instrevsequence1 = &(instsequence1[length1-1]);
  debug(printf("  Query sequence is %.*s\n",length1,&(instrevsequence1[-length1+1])));
#else
  debug(printf("  Query sequence is %.*s\n",length1,&(revsequence1[-length1+1])));
#endif

  *microexonlength = 0;
  if (length2 < length1) {
    maxc = length2 - MIN_MICROEXON_LENGTH;
  } else {
    maxc = length1 - MIN_MICROEXON_LENGTH;
  }
  for (c = 0; c < maxc; c++) {
    right2 = revsequenceuc2[-c-1];
    right1 = revsequenceuc2[-c];
    debug(printf("   Checking %c%c\n",right2,right1));
#ifdef PMAP
    if (c > 0 && matchtable[instrevsequence1[-c+1]-'A'][revsequenceuc2[-c+1]-'A'] == false) {
      nmismatches++;
    }
#else
    if (c > 0 && revsequenceuc1[-c+1] != revsequenceuc2[-c+1]) {
      nmismatches++;
    }
#endif
    if (nmismatches > 1) {
#ifdef PMAP
      debug(printf("   Aborting at %c !~ %c\n",instrevsequence1[-c+1],revsequence2[-c+1]));
      FREE(instsequence1);
#else
      debug(printf("   Aborting at %c != %c\n",revsequence1[-c+1],revsequence2[-c+1]));
#endif
      *microintrontype = NONINTRON;
      return NULL;
    }
    if (right2 == intron3 && right1 == intron4) {
      endlength = length1 - c;
#ifdef PMAP
      debug(printf("  Found acceptor at %d, length %d.  End sequence is %.*s\n",
		       c,endlength,endlength,&(instrevsequence1[-endlength+1])));
#else
      debug(printf("  Found acceptor at %d, length %d.  End sequence is %.*s\n",
		       c,endlength,endlength,&(revsequence1[-endlength+1])));
#endif

      textright = revoffset2 - c - MICROINTRON_LENGTH;
      textleft = textright - search_length(endlength,textright,end_microexons_p) + MICROINTRON_LENGTH;
#ifdef PMAP
      hits = BoyerMoore(&(instrevsequence1[-c-endlength+1]),endlength,&(genomicuc[textleft]),textright-textleft);
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
	  pairs = make_microexon_pairs_single(revoffset1-c-endlength+1,revoffset1-c+1,
					      candidate,revoffset2-c+1,endlength,c,
#ifdef PMAP
					      &(instrevsequence1[-revoffset1]),&(instrevsequence1[-revoffset1]),
#else
					      queryseq,queryuc,
#endif
					      chroffset,watsonp,pairpool,gapchar,*dynprogindex);
#ifdef PMAP
	  FREE(instsequence1);
#endif
	  *dynprogindex += (*dynprogindex > 0 ? +1 : -1);
	  return pairs;
	}
      }
      
      Intlist_free(&hits);
    }
  }

#ifdef PMAP
  FREE(instsequence1);
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
		     char *sequence1, char *sequenceuc1, char *sequence2, char *sequenceuc2,
		     int length1, int length2, int offset1, int offset2, int cdna_direction,
#ifdef PMAP
		     char *queryaaseq,
#endif
		     char *queryseq, char *queryuc, char *genomicseg, char *genomicuc,
		     int genomiclength, Pairpool_T pairpool, bool end_microexons_p) {
  List_T pairs = NULL;
#ifdef PMAP
  char *instsequence1;
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
  debug(printf("Begin microexon search at 3' for %.*s\n",length2,sequence2));
#else
  debug(printf("Begin microexon search at 3'\n"));
#endif

#ifdef PMAP
  instsequence1 = instantiate_codons(sequence1,queryaaseq,offset1,length1);
  debug(printf("  Query sequence is %.*s\n",length1,instsequence1));
#else
  debug(printf("  Query sequence is %.*s\n",length1,sequence1));
#endif

  *microexonlength = 0;
  if (length2 < length1) {
    maxc = length2 - MIN_MICROEXON_LENGTH;
  } else {
    maxc = length1 - MIN_MICROEXON_LENGTH;
  }
  for (c = 0; c < maxc; c++) {
    left1 = sequenceuc2[c];
    left2 = sequenceuc2[c+1];
    debug(printf("   Checking %c%c\n",left1,left2));
#ifdef PMAP
    if (c > 0 && matchtable[instsequence1[c-1]-'A'][sequenceuc2[c-1]-'A'] == false) {
      nmismatches++;
    }
#else
    if (c > 0 && sequenceuc1[c-1] != sequenceuc2[c-1]) {
      nmismatches++;
    }
#endif
    if (nmismatches > 1) {
#ifdef PMAP
      debug(printf("   Aborting at %c !~ %c\n",instsequence1[c-1],sequence2[c-1]));
      FREE(instsequence1);
#else
      debug(printf("   Aborting at %c != %c\n",sequence1[c-1],sequence2[c-1]));
#endif
      *microintrontype = NONINTRON;
      return NULL;
    }
    if (left1 == intron1 && left2 == intron2) {
      endlength = length1 - c;
#ifdef PMAP
      debug(printf("  Found donor at %d, length %d.  End sequence is %.*s\n",
		   c,endlength,endlength,&(instsequence1[c])));
#else
      debug(printf("  Found donor at %d, length %d.  End sequence is %.*s\n",
		   c,endlength,endlength,&(sequence1[c])));
#endif

      textleft = offset2 + c;
      textright = textleft + search_length(endlength,genomiclength-textleft,end_microexons_p);
#ifdef PMAP
      hits = BoyerMoore(&(instsequence1[c]),endlength,&(genomicuc[textleft]),textright-textleft);
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
	  pairs = make_microexon_pairs_single(offset1,offset1+c,
					      offset2,candidate,c,endlength,
#ifdef PMAP
					      &(instsequence1[-offset1]),&(instsequence1[-offset1]),
#else
					      queryseq,queryuc,
#endif
					      genomicseg,genomicuc,
					      pairpool,gapchar,*dynprogindex);
#ifdef PMAP
	  FREE(instsequence1);
#endif
	  *dynprogindex += (*dynprogindex > 0 ? +1 : -1);
	  return pairs;
	}
      }
      
      Intlist_free(&hits);
    }
  }

#ifdef PMAP
  FREE(instsequence1);
#endif

  debug(printf("End of dynprog microexon 3\n"));

  *microintrontype = NONINTRON;
  return NULL;
}
#endif

