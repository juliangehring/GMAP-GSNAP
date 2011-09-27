static char rcsid[] = "$Id: dynprog.c 27450 2010-08-05 19:02:48Z twu $";

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
#include "mem.h"
#include "comp.h"
#include "pair.h"
#include "pairdef.h"
#include "boyer-moore.h"
#include "intron.h"

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

/* Prints all bridge scores */
#ifdef DEBUG3
#define debug3(x) x
#else
#define debug3(x)
#endif

/* Codon instantiation */
#ifdef DEBUG4
#define debug4(x) x
#else
#define debug4(x)
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

typedef enum {HIGHQ, MEDQ, LOWQ, END} Mismatchtype_T;
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
#define MISMATCH_HIGHQ -5
#define MISMATCH_MEDQ -4
#define MISMATCH_LOWQ -3

/* Previously allowed lower mismatch scores on end to allow more
   complete alignments to the end, and because ends are typically of
   lower quality.  Previously made equal to FULLMATCH, because
   criterion is nmatches >= nmismatches.  However, extensions at ends
   appear to defeat purpose of trimming, so increase mismatch at end
   from -3 to -4. */
#define MISMATCH_END -4


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
#define END_OPEN_HIGHQ -10
#define END_OPEN_MEDQ -10
#define END_OPEN_LOWQ -10

#define END_EXTEND_HIGHQ -3
#define END_EXTEND_MEDQ -3
#define END_EXTEND_LOWQ -3


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
#define CANONICAL_INTRON_HIGHQ 20 /* GT-AG */
#define CANONICAL_INTRON_MEDQ  26
#define CANONICAL_INTRON_LOWQ  32

#define FINAL_CANONICAL_INTRON_HIGHQ 50 /* GT-AG */
#define FINAL_CANONICAL_INTRON_MEDQ  56
#define FINAL_CANONICAL_INTRON_LOWQ  62

/* Prefer alternate intron to other non-canonicals, but don't
   introduce mismatches or gaps to identify */
#define ALTERNATE_INTRON 15	/* GC-AG or AT-AC */
#define FINAL_ALTERNATE_INTRON 40  /* GC-AG or AT-AC.  Amount above
				      regular should approximately
				      match FINAL_CANONICAL_INTRON -
				      CANONICAL_INTRON */

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


/************************************************************************
 * Matrix
 ************************************************************************/

/* Makes a matrix of dimensions 0..length1 x 0..length2 inclusive */
static int **
Matrix_alloc (int length1, int length2, int **ptrs, int *space) {
  int **matrix, i;

  if (length1 < 0 || length2 < 0) {
    fprintf(stderr,"lengths are negative: %d %d\n",length1,length2);
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

static void
Matrix_print (int **matrix, int length1, int length2, char *sequence1, char *sequence2,
	      bool revp) {
  int i, j;

  printf("  ");
  for (j = 0; j <= length2; ++j) {
    if (j == 0) {
      printf("    ");
    } else {
      printf("  %c ",revp ? sequence2[-j+1] : sequence2[j-1]);
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
      printf("%3d ",matrix[i][j]);
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
		  char *sequence1, char *sequence2, bool revp) {
  int i, j;
  char buffer[4];

  printf("  ");
  for (j = 0; j <= length2; ++j) {
    if (j == 0) {
      printf("    ");
    } else {
      printf("  %c ",revp ? sequence2[-j+1] : sequence2[j-1]);
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


#define QUERY_MAXLENGTH 500
#define GENOMIC_MAXLENGTH 2000


#define T Dynprog_T
struct T {
  int maxlength1;
  int maxlength2;

  int **matrix_ptrs, *matrix_space;
  Direction_T **directions_ptrs, *directions_space;
  int **jump_ptrs, *jump_space;
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

  new->matrix_ptrs = (int **) CALLOC(maxlength1+1,sizeof(int *));
  new->matrix_space = (int *) CALLOC((maxlength1+1)*(maxlength2+1),sizeof(int));
  new->directions_ptrs = (Direction_T **) CALLOC(maxlength1+1,sizeof(Direction_T *));
  new->directions_space = (Direction_T *) CALLOC((maxlength1+1)*(maxlength2+1),sizeof(Direction_T));
  new->jump_ptrs = (int **) CALLOC(maxlength1+1,sizeof(int *));
  new->jump_space = (int *) CALLOC((maxlength1+1)*(maxlength2+1),sizeof(int));

  return new;
}


void
Dynprog_free (T *old) {
  if (*old) {
    FREE((*old)->matrix_ptrs);
    FREE((*old)->matrix_space);
    FREE((*old)->directions_ptrs);
    FREE((*old)->directions_space);
    FREE((*old)->jump_ptrs);
    FREE((*old)->jump_space);

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

  if ((index = aa_index_table[aa]) < 0) {
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
      index = aa_index_table[aa];
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
pairdistance_init () {
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
      pairdistance_array[END][c1][c2] = MISMATCH_END;
    }
  }
#else
  for (c1 = 'A'; c1 <= 'Z'; c1++) {
    for (c2 = 'A'; c2 < 'Z'; c2++) {
      pairdistance_array[HIGHQ][c1][c2] = MISMATCH_HIGHQ;
      pairdistance_array[MEDQ][c1][c2] = MISMATCH_MEDQ;
      pairdistance_array[LOWQ][c1][c2] = MISMATCH_LOWQ;
      pairdistance_array[END][c1][c2] = MISMATCH_END;
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

  for (c = 'A'; c < 'Z'; c++) {
    permute_cases(c,c,FULLMATCH);
  }

  return;
}


/************************************************************************/

/* static int max_jump_penalty_lookup; */
static int *jump_penalty_array[NJUMPTYPES];


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
  return open + extend*(length-1);
#endif
}


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
	      int extramaterial_end, int extramaterial_paired) {
  pairdistance_init();
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

static int **
compute_scores_lookup (Direction_T ***directions, int ***jump, T this, 
		       char *sequence1, char *sequence2, bool revp, 
		       int length1, int length2, 
		       Mismatchtype_T mismatchtype, Jumptype_T jumptype, bool init_jump_penalty_p,
		       int extraband, bool widebandp, bool onesidegapp) {
  int **matrix;
  int r, c, r1, c1, na1, na2;
  int bestscore, score, bestjump, j;
  Direction_T bestdir;
  int lband, rband, clo, chigh, rlo;
  int **pairdistance_array_type;
  int *jump_penalty_array_type;

  pairdistance_array_type = pairdistance_array[mismatchtype];
  jump_penalty_array_type = jump_penalty_array[jumptype];

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
  debug(printf("Lengths are %d and %d, so bands are %d on left and %d on right\n",length1,length2,lband,rband));

  matrix = Matrix_alloc(length1,length2,this->matrix_ptrs,this->matrix_space);
  *directions = Directions_alloc(length1,length2,this->directions_ptrs,this->directions_space);
  *jump = Matrix_alloc(length1,length2,this->jump_ptrs,this->jump_space);

  matrix[0][0] = 0;
  (*directions)[0][0] = STOP;
  (*jump)[0][0] = 0;

  /* Row 0 initialization */
  if (init_jump_penalty_p == true) {
    for (c = 1; c <= rband && c <= length2; c++) {
      matrix[0][c] = jump_penalty_array_type[c];
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
      matrix[r][0] = jump_penalty_array_type[r];
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
	    score = matrix[r][c1] + jump_penalty_array_type[j];
	    if (score > bestscore) {
	      bestscore = score;
	      bestdir = HORIZ;
	      bestjump = j;
	    }
	  }
	}
      } else {
	for (c1 = c-1, j = 1; c1 >= clo; c1--, j++) {
	  if ((*directions)[r][c1] == DIAG) {
	    score = matrix[r][c1] + jump_penalty_array_type[j];
	    if (score > bestscore) {
	      bestscore = score;
	      bestdir = HORIZ;
	      bestjump = j;
	    }
#ifdef RIGHTANGLE
	  } else if ((*directions)[r][c1] == VERT && (*jump)[r][c1] >= 3) {
	    score = matrix[r][c1] + jump_penalty_array_type[j] - open;
	    if (score > bestscore) {
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
	    score = matrix[r1][c] + jump_penalty_array_type[j];
	    if (score > bestscore) {
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
	    score = matrix[r1][c] + jump_penalty_array_type[j];
	    if (score > bestscore) {
	      bestscore = score;
	      bestdir = VERT;
	      bestjump = j;
	    }
#ifdef RIGHTANGLE
	  } else if ((*directions)[r1][c] == HORIZ && (*jump)[r1][c] >= 3) {
	    score = matrix[r1][c] + jump_penalty_array_type[j] - open;
	    if (score > bestscore) {
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

static int **
compute_scores (Direction_T ***directions, int ***jump, T this, 
		char *sequence1, char *sequence2, bool revp, 
		int length1, int length2, 
		Mismatchtype_T mismatchtype,
		int mismatch, int open, int extend, bool init_jump_penalty_p,
		int extraband, bool widebandp, bool onesidegapp) {
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
	    if (score > bestscore) {
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
	    if (score > bestscore) {
	      bestscore = score;
	      bestdir = HORIZ;
	      bestjump = j;
	    }
#ifdef RIGHTANGLE
	  } else if ((*directions)[r][c1] == VERT && (*jump)[r][c1] >= 3) {
	    score = matrix[r][c1] + jump_penalty(j,open,extend) - open;
	    if (score > bestscore) {
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
	    if (score > bestscore) {
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
	    if (score > bestscore) {
	      bestscore = score;
	      bestdir = VERT;
	      bestjump = j;
	    }
#ifdef RIGHTANGLE
	  } else if ((*directions)[r1][c] == HORIZ && (*jump)[r1][c] >= 3) {
	    score = matrix[r1][c] + jump_penalty(j,open,extend) - open;
	    if (score > bestscore) {
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

static void
find_best_endpoint (int *finalscore, int *bestr, int *bestc, int **matrix, 
		    int length1, int length2, int extraband_end_or_paired) {
  int bestscore = 0;
  int r, c;
  int rband, lband, clo, chigh;
  /* No need for loffset or cmid because they apply only for cdnaend
     == FIVE, which doesn't require searching */

  *bestr = *bestc = 0;

  rband = length2 - length1 + extraband_end_or_paired;
  lband = extraband_end_or_paired;

  for (r = 1; r <= length1; r++) {
    if ((clo = r - lband) < 1) {
      clo = 1;
    }
    if ((chigh = r + rband) > length2) {
      chigh = length2;
    }
    for (c = clo; c <= chigh; c++) {
      if (matrix[r][c] > bestscore) {
	*bestr = r;
	*bestc = c;
	bestscore = matrix[r][c];
      }
    }
  }
  *finalscore = bestscore;
  return;
}


static void
find_best_endpoint_onegap (int *finalscore, int *bestr, int *bestc, int **matrix, 
			   int length1, int length2, int extraband_end,
			   bool extend_mismatch_p) {
  int bestscore = -1000000;
  int r, c;
  int rband, lband, clo, chigh;
  /* No need for loffset or cmid because they apply only for cdnaend
     == FIVE, which doesn't require searching */

  *bestr = *bestc = 0;

  rband = lband = 1;

  for (r = 1; r <= length1; r++) {
    if ((clo = r - lband) < 1) {
      clo = 1;
    }
    if ((chigh = r + rband) > length2) {
      chigh = length2;
    }
    for (c = clo; c <= chigh; c++) {
      if (matrix[r][c] >= bestscore) {
	*bestr = r;
	*bestc = c;
	bestscore = matrix[r][c];
      }
    }
  }

  *finalscore = bestscore;

  if (extend_mismatch_p == true) {
    if (*bestr == length1 - 1) {
      if ((chigh = *bestr + rband) > length2) {
	chigh = length2;
      }
      if (*bestc < chigh) {
	debug(printf("Extending by 1 just to be complete\n"));
	*bestr += 1;
	*bestc += 1;
      }
    }
  }
  return;
}


/* Finds best score along diagonal */
static void
find_best_endpoint_nogap_fwd (int *finalscore, int *bestr, int *bestc, int **matrix, int diaglength,
			      char *sequenceuc2, char *dinucleotide) {
  int bestscore = -1000000, score;
  int rc;
  /* No need for loffset or cmid because they apply only for cdnaend
     == FIVE, which doesn't require searching */

  *bestr = 0;

  for (rc = 1; rc <= diaglength; rc++) {
    score = matrix[rc][rc];
    /* printf("At rc = %d, got score %d plus %c%c\n",rc,score,sequenceuc2[rc],sequenceuc2[rc+1]); */
    if (sequenceuc2[rc] == dinucleotide[0] && sequenceuc2[rc+1] == dinucleotide[1]) {
      score += 1000;
    }
    if (score >= bestscore) {
      /* printf("At rc = %d, got best score %d\n",rc,score); */
      *bestr = rc;
      bestscore = score;
    }
  }
  *bestc = *bestr;
  *finalscore = bestscore;
  return;
}

static void
find_best_endpoint_nogap_rev (int *finalscore, int *bestr, int *bestc, int **matrix, int diaglength,
			      char *revsequenceuc2, char *dinucleotide) {
  int bestscore = -1000000, score;
  int rc;
  /* No need for loffset or cmid because they apply only for cdnaend
     == FIVE, which doesn't require searching */

  *bestr = 0;

  for (rc = 1; rc <= diaglength; rc++) {
    score = matrix[rc][rc];
    /* printf("At rc = %d, got score %d plus %c%c\n",rc,score,revsequenceuc2[-1-rc],revsequenceuc2[-rc]); */
    if (revsequenceuc2[-1-rc] == dinucleotide[0] && revsequenceuc2[-rc] == dinucleotide[1]) {
      score += FINAL_CANONICAL_INTRON_HIGHQ/2;
    }
    if (score >= bestscore) {
      /* printf("At rc = %d, got best score %d\n",rc,score); */
      *bestr = rc;
      bestscore = score;
    }
  }
  *bestc = *bestr;
  *finalscore = bestscore;
  return;
}


static void
find_best_endpoint_length2 (int *finalscore, int *bestr, int **matrix, Direction_T **directions,
			    int length1, int length2, int extraband_end) {
  int bestscore = -1000000;
  int r;
  /* No need for loffset or cmid because they apply only for cdnaend
     == FIVE, which doesn't require searching */

  *bestr = 0;

  for (r = 1; r <= length1; r++) {
    if (directions[r][length2] != STOP && matrix[r][length2] >= bestscore) {
      *bestr = r;
      bestscore = matrix[r][length2];
    }
  }

  *finalscore = bestscore;

  return;
}


static List_T
add_queryskip (List_T pairs, int r, int c, int dist, char *querysequence, char *querysequenceuc, 
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
			  c1,INDEL_COMP,' ',dynprogindex);
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
		int queryoffset, int genomeoffset, Pairpool_T pairpool, bool revp, bool cdna_gap_p,
		int cdna_direction, int dynprogindex) {
  int j;
  char left1, left2, right2, right1, c2;
  int introntype;
  int querycoord, leftgenomecoord, rightgenomecoord, genomecoord, temp, step;

  if (cdna_gap_p == false) {
    querycoord = r-1;
    leftgenomecoord = c-dist;
    rightgenomecoord = c-1;
  } else {
    querycoord = c-1;
    leftgenomecoord = r-dist;
    rightgenomecoord = r-1;
  }
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
    left1 = genomesequenceuc[leftgenomecoord];
    left2 = genomesequenceuc[leftgenomecoord+1];
    right2 = genomesequenceuc[rightgenomecoord-1];
    right1 = genomesequenceuc[rightgenomecoord];
	
#ifdef PMAP
    introntype = Intron_type(left1,left2,right2,right1,/*cdna_direction*/+1);
#else
    introntype = Intron_type(left1,left2,right2,right1,cdna_direction);
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
      c2 = genomesequence[genomecoord];
      debug(printf("Pushing %d,%d [%d,%d] (-,%c),",r,c,queryoffset+querycoord,genomeoffset+genomecoord,c2));
      pairs = Pairpool_push(pairs,pairpool,queryoffset+querycoord,genomeoffset+genomecoord,
			    ' ',INDEL_COMP,c2,dynprogindex);
      debug(
	    if (cdna_gap_p == false) {
	      c--;
	    } else {
	      r--;
	    }
	    );
      genomecoord += step;
    }
  } else {
    debug(printf("Large gap %c%c..%c%c.  Adding gap of type %d.\n",
		 left1,left2,right2,right1,introntype));
#ifndef NOGAPHOLDER
#if 0
    pairs = Pairpool_push_gapholder(pairs,pairpool,/*queryjump*/0,/*genomejump*/dist);
#endif
    pairs = Pairpool_push_gapholder(pairs,pairpool,/*queryjump*/UNKNOWNJUMP,/*genomejump*/UNKNOWNJUMP);
#endif
  }
  
  return pairs;
}


/* Preferences: Continue in same direction if possible.  This has the
   effect of extending gaps to maximum size. Then, take diagonal if
   possible. Finally, take vertical if possible, because this will use
   up sequence 1, which is the query (cDNA) sequence. */

/* revp means both rev1p and rev2p, which must have equal values */
/* Iterative version */
static List_T
traceback (List_T pairs, int *nmatches, int *nmismatches, int *nopens, int *nindels,
	   Direction_T **directions, int **jump, int r, int c, 
	   char *querysequence, char *querysequenceuc, char *genomesequence, char *genomesequenceuc,
	   int queryoffset, int genomeoffset, Pairpool_T pairpool, 
	   bool revp, bool cdna_gap_p, int cdna_direction, int dynprogindex) {
  char c1, c2;
  int dist;
  bool add_dashes_p;
  int querycoord, genomecoord;

  debug(printf("Starting traceback at r=%d,c=%d (offset1=%d, offset2=%d)\n",r,c,queryoffset,genomeoffset));

  while (directions[r][c] != STOP) {
    if (directions[r][c] == DIAG) {
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
      }

      c1 = querysequence[querycoord];
      c2 = genomesequence[genomecoord];

      if (querysequenceuc[querycoord] == genomesequenceuc[genomecoord]) {
	debug(printf("D%d: Pushing %d,%d [%d,%d] (%c,%c) - match\n",
		     jump[r][c],r,c,queryoffset+querycoord,genomeoffset+genomecoord,c1,c2));
	*nmatches += 1;
	pairs = Pairpool_push(pairs,pairpool,queryoffset+querycoord,genomeoffset+genomecoord,c1,DYNPROG_MATCH_COMP,c2,
			      dynprogindex);

      } else if (consistent_array[(int) c1][(int) c2] == true) {
	debug(printf("D%d: Pushing %d,%d [%d,%d] (%c,%c) - ambiguous\n",
		     jump[r][c],r,c,queryoffset+querycoord,genomeoffset+genomecoord,c1,c2));
	*nmatches += 1;
	pairs = Pairpool_push(pairs,pairpool,queryoffset+querycoord,genomeoffset+genomecoord,c1,AMBIGUOUS_COMP,c2,
			      dynprogindex);

      } else {
	debug(printf("D%d: Pushing %d,%d [%d,%d] (%c,%c) - mismatch\n",
		     jump[r][c],r,c,queryoffset+querycoord,genomeoffset+genomecoord,c1,c2));
	*nmismatches += 1;
	pairs = Pairpool_push(pairs,pairpool,queryoffset+querycoord,genomeoffset+genomecoord,c1,MISMATCH_COMP,c2,
			      dynprogindex);
      }
      r--; c--;

    } else if (directions[r][c] == HORIZ) {
      dist = jump[r][c];
      debug(printf("H%d: ",dist));
      if (cdna_gap_p == false) {
	pairs = add_genomeskip(&add_dashes_p,pairs,r,c,dist,genomesequence,genomesequenceuc,
			       queryoffset,genomeoffset,pairpool,revp,cdna_gap_p,
			       cdna_direction,dynprogindex);
	if (add_dashes_p == true) {
	  *nopens += 1;
	  *nindels += dist;
	}
      } else {
	pairs = add_queryskip(pairs,r,c,dist,querysequence,querysequenceuc,
			      queryoffset,genomeoffset,pairpool,revp,cdna_gap_p,dynprogindex);
	*nopens += 1;
	*nindels += dist;
      }
      c -= dist;
      debug(printf("\n"));

    } else if (directions[r][c] == VERT) {
      dist = jump[r][c];
      debug(printf("V%d: ",dist));
      if (cdna_gap_p == false) {
	pairs = add_queryskip(pairs,r,c,dist,querysequence,querysequenceuc,
			      queryoffset,genomeoffset,pairpool,revp,cdna_gap_p,dynprogindex);
	*nopens += 1;
	*nindels += dist;
      } else {
	pairs = add_genomeskip(&add_dashes_p,pairs,r,c,dist,genomesequence,genomesequenceuc,
			       queryoffset,genomeoffset,pairpool,revp,cdna_gap_p,
			       cdna_direction,dynprogindex);
	if (add_dashes_p == true) {
	  *nopens += 1;
	  *nindels += dist;
	}
      }
      r -= dist;
      debug(printf("\n"));

    } else {
      abort();
    }
  }
  return pairs;
}


static void
countback (int *nmatches, int *nmismatches, int *nopens, int *nindels,
	   Direction_T **directions, int **jump, int r, int c, 
	   char *sequenceuc1, char *sequenceuc2, int offset1, int offset2, 
	   Pairpool_T pairpool, bool revp) {
  char cuc1, cuc2;
  int dist;

  while (directions[r][c] != STOP) {
    if (directions[r][c] == DIAG) {
      cuc1 = revp ? sequenceuc1[-r+1] : sequenceuc1[r-1];
      cuc2 = revp ? sequenceuc2[-c+1] : sequenceuc2[c-1];
      if (cuc1 == cuc2) {
	*nmatches += 1;
      } else {
	*nmismatches += 1;
      }
      r--; c--;
    } else if (directions[r][c] == HORIZ) {
      debug(printf("H%d: ",jump[r][c]));
      *nopens += 1;
      *nindels += (dist = jump[r][c]);
      c -= dist;
    } else if (directions[r][c] == VERT) {
      debug(printf("V%d: ",jump[r][c]));
      *nopens += 1;
      *nindels += (dist = jump[r][c]);
      r -= dist;
    } else {
      abort();
    }
  }

  return;
}


/* We have switched length2 for columns, and length 1L and 1R for rows */
static void
bridge_cdna_gap (int *finalscore, int *bestrL, int *bestrR, int *bestcL, int *bestcR,
		 int **matrixL, int **matrixR, Direction_T **directionsL, Direction_T **directionsR, 
		 int length2, int length1L, int length1R, int extraband_paired,
		 int open, int extend, int leftoffset, int rightoffset) {
  int bestscore = -100000, scoreL, scoreR, pen, end_reward = 0;
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
	scoreL = matrixL[rL][cL];
	
	/* Disallow leftoffset + cL >= rightoffset - cR, or cR >= rightoffset - leftoffset - cL */
	/* debug3(printf("  Disallowing cR to be >= %d\n",rightoffset-leftoffset-cL)); */
	for (cR = cloR; cR <= chighR && cR < rightoffset-leftoffset-cL; cR++) {
	  scoreR = matrixR[rR][cR];

	  if (scoreL + scoreR + pen + end_reward > bestscore) {
	    /*
	    debug3(printf("At %d left to %d right, score is (%d)+(%d)+(%d)+(%d) = %d (bestscore)\n",
			  cL,cR,scoreL,scoreR,pen,end_reward,scoreL+scoreR+pen+end_reward));
	    */

	    bestscore = scoreL + scoreR + pen + end_reward;
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
    case GCAG_FWD: case ATAC_FWD: 
      scoreI = finalp == true ? FINAL_ALTERNATE_INTRON : ALTERNATE_INTRON;
      break;
    default: *introntype = NONINTRON; scoreI = 0.0;
    }
  }
#else
  if ((*introntype = leftdi & rightdi) == NONINTRON) {
    scoreI = 0.0;
  } else if (cdna_direction > 0) {
    switch (*introntype) {
    case GTAG_FWD: scoreI = canonical_reward; break;
    case GCAG_FWD: case ATAC_FWD: 
      scoreI = finalp == true ? FINAL_ALTERNATE_INTRON : ALTERNATE_INTRON;
      break;
    default: *introntype = NONINTRON; scoreI = 0.0;
    }
  } else if (cdna_direction < 0) {
    switch (*introntype) {
    case GTAG_REV: scoreI = canonical_reward; break;
    case GCAG_REV: case ATAC_REV:
      scoreI = finalp == true ? FINAL_ALTERNATE_INTRON : ALTERNATE_INTRON;
      break;
    default: *introntype = NONINTRON; scoreI = 0.0;
    }
  } else {
    switch (*introntype) {
    case GTAG_FWD: case GTAG_REV: scoreI = canonical_reward; break;
    case GCAG_FWD: case GCAG_REV: case ATAC_FWD: case ATAC_REV: 
      scoreI = finalp == true ? FINAL_ALTERNATE_INTRON : ALTERNATE_INTRON;
      break;
    default: *introntype = NONINTRON; scoreI = 0.0;
    }
  }
#endif

  return scoreI;
}


static void
bridge_intron_gap (int *finalscore, int *bestrL, int *bestrR, int *bestcL, int *bestcR, int *best_introntype, 
		   int **matrixL, int **matrixR, Direction_T **directionsL, Direction_T **directionsR, 
		   int **jumpL, int **jumpR,  char *sequenceuc2L, char *revsequenceuc2R, 
		   int length1, int length2L, int length2R,
		   int cdna_direction, int extraband_paired, int canonical_reward,
		   int maxhorizjump, int maxvertjump, int leftoffset, int rightoffset, bool halfp,
		   bool finalp) {
  int bestscore = -100000, scoreL, scoreR, scoreI, bestscoreI = -100000;
  int rL, rR, cL, cR;
  int lbandL, rbandL, cloL, chighL;
  int lbandR, rbandR, cloR, chighR;
  char left1, left2, right2, right1;
  int *leftdi, *rightdi, introntype;

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
    left1 = sequenceuc2L[cL];
    left2 = sequenceuc2L[cL+1];
    if (left1 == 'G' && left2 == 'T') {
      leftdi[cL] = LEFT_GT;
    } else if (left1 == 'G' && left2 == 'C') {
      leftdi[cL] = LEFT_GC;
    } else if (left1 == 'A' && left2 == 'T') {
      leftdi[cL] = LEFT_AT;
    } else if (left1 == 'C' && left2 == 'T') {
      leftdi[cL] = LEFT_CT;
    } else {
      leftdi[cL] = 0x00;
    }
  }
  leftdi[length2L-1] = leftdi[length2L] = 0x00;

  for (cR = 0; cR < length2R - 1; cR++) {
    right2 = revsequenceuc2R[-cR-1];
    right1 = revsequenceuc2R[-cR];
    if (right2 == 'A' && right1 == 'G') {
      rightdi[cR] = RIGHT_AG;
    } else if (right2 == 'A' && right1 == 'C') {
      rightdi[cR] = RIGHT_AC;
    } else if (right2 == 'G' && right1 == 'C') {
      rightdi[cR] = RIGHT_GC;
    } else if (right2 == 'A' && right1 == 'T') {
      rightdi[cR] = RIGHT_AT;
    } else {
      rightdi[cR] = 0x00;
    }
  }
  rightdi[length2R-1] = rightdi[length2R] = 0x00;

  /* Perform computations */
  rbandL = length2L - length1 + extraband_paired;
  lbandL = extraband_paired;

  rbandR = length2R - length1 + extraband_paired;
  lbandR = extraband_paired;

  for (rL = 1, rR = length1-1; rL < length1; rL++, rR--) {
    /* debug3(printf("\nAt row %d on left and %d on right\n",rL,rR)); */
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

    for (cL = cloL; cL <= chighL; cL++) {
      /* The following check limits genomic inserts (horizontal) and
         multiple cDNA inserts (vertical). */
      if (directionsL[rL][cL] == DIAG ||
	  (directionsL[rL][cL] == HORIZ && jumpL[rL][cL] <= maxhorizjump) ||
	  (directionsL[rL][cL] == VERT && jumpL[rL][cL] <= maxvertjump) ||
	  directionsL[rL][cL] == STOP) {
	scoreL = matrixL[rL][cL];
	
	if (directionsL[rL][cL] == HORIZ || directionsL[rL][cL] == VERT) {
	  /* Favor gaps away from intron if possible */
	  scoreL -= 1;
	}

	/* Disallow leftoffset + cL >= rightoffset - cR, or cR >= rightoffset - leftoffset - cL */
	/* debug3(printf("  Disallowing cR to be >= %d\n",rightoffset-leftoffset-cL)); */
	for (cR = cloR; cR <= chighR && cR < rightoffset-leftoffset-cL; cR++) {
	  if (directionsR[rR][cR] == DIAG ||
	      (directionsR[rR][cR] == HORIZ && jumpR[rR][cR] <= maxhorizjump) ||
	      (directionsR[rR][cR] == VERT && jumpR[rR][cR] <= maxvertjump) ||
	      directionsR[rR][cR] == STOP) {
	    scoreR = matrixR[rR][cR];

	    if (directionsR[rR][cR] == HORIZ || directionsR[rR][cR] == VERT) {
	      /* Favor gaps away from intron if possible */
	      scoreR -= 1;
	    }
	    
	    scoreI = intron_score(&introntype,leftdi[cL],rightdi[cR],
				  cdna_direction,canonical_reward,finalp);
	    
	    if (scoreL + scoreI + scoreR > bestscore) {
	      /*
		debug3(printf("At %d left to %d right, score is (%d)+(%d)+(%d) = %d (bestscore)\n",
		cL,cR,scoreL,scoreI,scoreR,scoreL+scoreI+scoreR));
	      */
	      bestscore = scoreL + scoreI + scoreR;
	      bestscoreI = scoreI;
	      *bestrL = rL;
	      *bestrR = rR;
	      *bestcL = cL;
	      *bestcR = cR;
	      *best_introntype = introntype;

	    } else {
	      /*
		debug3(printf("At %d left to %d right, score is (%d)+(%d)+(%d) = %d\n",
		cL,cR,scoreL,scoreI,scoreR,scoreL+scoreI+scoreR));
	      */
	    }
	  }
	}
      }
    }
  }
      
  if (halfp == true) {
    *finalscore = bestscore - bestscoreI/2;
  } else {
    *finalscore = bestscore;
  }
  debug3(printf("Returning final score of %d at (%d,%d) left to (%d,%d) right\n",
		*finalscore,*bestrL,*bestcL,*bestrR,*bestcR));
  FREE(rightdi);
  FREE(leftdi);

  return;
}


static void
bridge_dual_break_nogap (int *finalscore, int *bestrcL, int *bestrcR, int **matrixL, int **matrixR,
			 int diaglength, char *sequenceuc2L, char *revsequenceuc2R, int cdna_direction) {
  int bestscore = -1000000, scoreL, scoreR, scoreI;
  int rcL, rcR;
  char left1, left2, right2, right1;
  int introntype, leftdi, rightdi;

  *bestrcL = *bestrcR = 0;

  for (rcL = 1; rcL <= diaglength; rcL++) {
    rcR = diaglength - rcL;
    left1 = sequenceuc2L[rcL];
    left2 = sequenceuc2L[rcL+1];
    if (left1 == 'G' && left2 == 'T') {
      leftdi = LEFT_GT;
    } else if (left1 == 'G' && left2 == 'C') {
      leftdi = LEFT_GC;
    } else if (left1 == 'A' && left2 == 'T') {
      leftdi = LEFT_AT;
    } else if (left1 == 'C' && left2 == 'T') {
      leftdi = LEFT_CT;
    } else {
      leftdi = 0x00;
    }

    right2 = revsequenceuc2R[-1-rcR];
    right1 = revsequenceuc2R[-rcR];
    if (right2 == 'A' && right1 == 'G') {
      rightdi = RIGHT_AG;
    } else if (right2 == 'A' && right1 == 'C') {
      rightdi = RIGHT_AC;
    } else if (right2 == 'G' && right1 == 'C') {
      rightdi = RIGHT_GC;
    } else if (right2 == 'A' && right1 == 'T') {
      rightdi = RIGHT_AT;
    } else {
      rightdi = 0x00;
    }
    
    scoreL = matrixL[rcL][rcL];
    scoreI = intron_score(&introntype,leftdi,rightdi,
			  cdna_direction,/*canonical_reward*/FINAL_CANONICAL_INTRON_HIGHQ,/*finalp*/true);
    scoreR = matrixR[rcR][rcR];

    if (scoreL + scoreI + scoreR > bestscore) {
      *bestrcL = rcL;
      *bestrcR = rcR;
      bestscore = scoreL + scoreI + scoreR;
    }
  }

  *finalscore = bestscore;
  return;
}



static void
bridge_dual_break_fwd (int *finalscore, int *bestrcL, int *bestrcR, int **matrixL, int **matrixR,
		       int diaglength, char *sequenceuc2L, char *revsequenceuc2R) {
  int bestscore, scoreL, scoreR, scoreI;
  int rcL, rcR;
  int bestscoreL_GT, bestscoreL_GC, bestscoreL_AT, bestscoreL_XX;
  int bestrcL_GT, bestrcL_GC, bestrcL_AT, bestrcL_XX;
  int bestscoreR_AG, bestscoreR_AC, bestscoreR_XX;
  int bestrcR_AG, bestrcR_AC, bestrcR_XX;
  char left1, left2, right2, right1;
  int introntype;

  *bestrcL = *bestrcR = 0;

  bestrcL_GT = bestrcL_GC = bestrcL_AT = bestrcL_XX = 0;
  bestscoreL_GT = bestscoreL_GC = bestscoreL_AT = bestscoreL_XX = -1000000;

  for (rcL = 1; rcL <= diaglength; rcL++) {
    if ((scoreL = matrixL[rcL][rcL]) >= bestscoreL_XX) {
      bestscoreL_XX = scoreL;
      bestrcL_XX = rcL;
    }
    left1 = sequenceuc2L[rcL];
    left2 = sequenceuc2L[rcL+1];
    if (left1 == 'G' && left2 == 'T') {
      if (scoreL >= bestscoreL_GT) {
	bestscoreL_GT = scoreL;
	bestrcL_GT = rcL;
      }
    } else if (left1 == 'G' && left2 == 'C') {
      if (scoreL >= bestscoreL_GC) {
	bestscoreL_GC = scoreL;
	bestrcL_GC = rcL;
      }
    } else if (left1 == 'A' && left2 == 'T') {
      if (scoreL >= bestscoreL_AT) {
	bestscoreL_AT = scoreL;
	bestrcL_AT = rcL;
      }
    }
  }


  bestrcR_AG = bestrcR_AC = bestrcR_XX = 0;
  bestscoreR_AG = bestscoreR_AC = bestscoreR_XX = -1000000;

  for (rcR = 1; rcR <= diaglength; rcR++) {
    if ((scoreR = matrixR[rcR][rcR]) >= bestscoreR_XX) {
      bestscoreR_XX = scoreR;
      bestrcR_XX = rcR;
    }
    right2 = revsequenceuc2R[-1-rcR];
    right1 = revsequenceuc2R[-rcR];
    if (right2 == 'A' && right1 == 'G') {
      if (scoreR >= bestscoreR_AG) {
	bestscoreR_AG = scoreR;
	bestrcR_AG = rcR;
      }
    } else if (right2 == 'A' && right1 == 'C') {
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
			    diaglength,sequenceuc2L,revsequenceuc2R,/*cdna_direction*/+1);
  }

  return;
}

static void
bridge_dual_break_rev (int *finalscore, int *bestrcL, int *bestrcR, int **matrixL, int **matrixR,
		       int diaglength, char *sequenceuc2L, char *revsequenceuc2R) {
  int bestscore, scoreL, scoreR, scoreI;
  int rcL, rcR;
  int bestscoreL_CT, bestscoreL_GT, bestscoreL_XX;
  int bestrcL_CT, bestrcL_GT, bestrcL_XX;
  int bestscoreR_AC, bestscoreR_GC, bestscoreR_AT, bestscoreR_XX;
  int bestrcR_AC, bestrcR_GC, bestrcR_AT, bestrcR_XX;
  char left1, left2, right2, right1;
  int introntype;

  *bestrcL = *bestrcR = 0;

  bestrcL_CT = bestrcL_GT = bestrcL_XX = 0;
  bestscoreL_CT = bestscoreL_GT = bestscoreL_XX = -1000000;

  for (rcL = 1; rcL <= diaglength; rcL++) {
    if ((scoreL = matrixL[rcL][rcL]) >= bestscoreL_XX) {
      bestscoreL_XX = scoreL;
      bestrcL_XX = rcL;
    }
    left1 = sequenceuc2L[rcL];
    left2 = sequenceuc2L[rcL+1];
    if (left1 == 'C' && left2 == 'T') {
      if (scoreL >= bestscoreL_CT) {
	bestscoreL_CT = scoreL;
	bestrcL_CT = rcL;
      }
    } else if (left1 == 'G' && left2 == 'T') {
      if (scoreL >= bestscoreL_GT) {
	bestscoreL_GT = scoreL;
	bestrcL_GT = rcL;
      }
    }
  }


  bestrcR_AC = bestrcR_GC = bestrcR_AT = bestrcR_XX = 0;
  bestscoreR_AC = bestscoreR_GC = bestscoreR_AT = bestscoreR_XX = -1000000;

  for (rcR = 1; rcR <= diaglength; rcR++) {
    if ((scoreR = matrixR[rcR][rcR]) >= bestscoreR_XX) {
      bestscoreR_XX = scoreR;
      bestrcR_XX = rcR;
    }
    right2 = revsequenceuc2R[-1-rcR];
    right1 = revsequenceuc2R[-rcR];
    if (right2 == 'A' && right1 == 'C') {
      if (scoreR >= bestscoreR_AC) {
	bestscoreR_AC = scoreR;
	bestrcR_AC = rcR;
      }
    } else if (right2 == 'G' && right1 == 'C') {
      if (scoreR >= bestscoreR_GC) {
	bestscoreR_GC = scoreR;
	bestrcR_GC = rcR;
      }
    } else if (right2 == 'A' && right1 == 'T') {
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
			    diaglength,sequenceuc2L,revsequenceuc2R,/*cdna_direction*/-1);
  }

  return;
}



/************************************************************************/

List_T
Dynprog_single_gap (int *dynprogindex, int *finalscore, int *nmatches, int *nmismatches, int *nopens, int *nindels,
		    T dynprog, char *sequence1, char *sequenceuc1, char *sequence2, char *sequenceuc2,
		    int length1, int length2, int offset1, int offset2,
#ifdef PMAP
		    char *queryaaseq,
#endif
		    int cdna_direction, Pairpool_T pairpool, int extraband_single, double defect_rate,
		    bool widebandp) {
  List_T pairs = NULL;
#ifdef PMAP
  char *instsequence1;
#endif
  Mismatchtype_T mismatchtype;
  Jumptype_T jumptype;
  int **matrix, **jump;
  Direction_T **directions;

  if (defect_rate < DEFECT_HIGHQ) {
    mismatchtype = HIGHQ;
    jumptype = SINGLE_HIGHQ;
  } else if (defect_rate < DEFECT_MEDQ) {
    mismatchtype = MEDQ;
    jumptype = SINGLE_MEDQ;
  } else {
    mismatchtype = LOWQ;
    jumptype = SINGLE_LOWQ;
  }


  /* Length1: maxlookback+MAXPEELBACK.  Length2 +EXTRAMATERIAL */

  debug(
	printf("%c:  ",*dynprogindex > 0 ? (*dynprogindex-1)%26+'a' : (-(*dynprogindex)-1)%26+'A');
	printf("Aligning single gap middle with wideband = %d and extraband %d\n",widebandp,extraband_single);
	printf("At query offset %d-%d, %.*s\n",offset1,offset1+length1-1,length1,sequence1);
	printf("At genomic offset %d-%d, %.*s\n",offset2,offset2+length2-1,length2,sequence2);
	printf("\n");
	);

  if (length1 > dynprog->maxlength1 || length2 > dynprog->maxlength2) {
    debug(printf("length1 %d or length2 %d is too long.  Returning NULL\n",length1,length2));
    *finalscore = -10000;
    *nmatches = *nmismatches = *nopens = *nindels = 0;
    /* Don't push a gapholder for single gap, because gapholder already exists */
    /* pairs = Pairpool_push_gapholder(NULL,pairpool,length1,length2); */
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

  matrix = compute_scores_lookup(&directions,&jump,dynprog,
#ifdef PMAP
				 instsequence1,sequence2,
#else
				 sequence1,sequence2,
#endif
				 /*revp*/false,length1,length2,
				 mismatchtype,jumptype,/*init_jump_penalty_p*/true,
				 extraband_single,widebandp,/*onesidegapp*/true);
  *finalscore = matrix[length1][length2];

  *nmatches = *nmismatches = *nopens = *nindels = 0;
  pairs = traceback(NULL,&(*nmatches),&(*nmismatches),&(*nopens),&(*nindels),
		    directions,jump,length1,length2,
#ifdef PMAP
		    instsequence1,instsequence1,sequence2,sequenceuc2,
#else
		    sequence1,sequenceuc1,sequence2,sequenceuc2,
#endif
		    offset1,offset2,pairpool,
		    /*revp*/false,/*cdna_gap_p*/false,cdna_direction,*dynprogindex);

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
		  char *sequence2, char *sequenceuc2,
		  int length1L, int length1R, int length2,
		  int offset1L, int revoffset1R, int offset2,
#ifdef PMAP
		  char *queryaaseq,
#endif
		  int cdna_direction, Pairpool_T pairpool, int extraband_paired, double defect_rate,
		  bool forcep) {
  List_T pairs = NULL;
  char *revsequence2, *revsequenceuc2;
#ifdef PMAP
  char *instsequence1, *instsequence1L, *instrevsequence1R;
#endif
  Mismatchtype_T mismatchtype;
  Jumptype_T jumptype;
  int **matrixL, **matrixR, **jumpL, **jumpR, mismatch, open, extend;
  Direction_T **directionsL, **directionsR;
  int revoffset2, bestrL, bestrR, bestcL, bestcR, k;
  int nmatches, nmismatches, nopens, nindels;
  int queryjump, genomejump;

  if (length2 <= 1) {
    return NULL;
  }

  debug(
	printf("\n");
	printf("%c:  ",*dynprogindex > 0 ? (*dynprogindex-1)%26+'a' : (-(*dynprogindex)-1)%26+'A');
	printf("Aligning cdna gap\n");
	printf("At query offset %d-%d, %.*s\n",offset1L,offset1L+length1L-1,length1L,sequence1L);
	printf("At query offset %d-%d, %.*s\n",revoffset1R-length1R+1,revoffset1R,length1R,&(revsequence1R[-length1R+1]));
	printf("Whole piece at query offset %d-%d, %.*s\n",offset1L,revoffset1R,revoffset1R-offset1L+1,sequence1L);
	printf("At genomic offset %d-%d, %.*s\n",offset2,offset2+length2-1,length2,sequence2);
	printf("\n");
	);

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
    jumptype = CDNA_HIGHQ;
    open = CDNA_OPEN_HIGHQ;
    extend = CDNA_EXTEND_HIGHQ;
  } else if (defect_rate < DEFECT_MEDQ) {
    mismatchtype = MEDQ;
    mismatch = MISMATCH_MEDQ;
    jumptype = CDNA_MEDQ;
    open = CDNA_OPEN_MEDQ;
    extend = CDNA_EXTEND_MEDQ;
  } else {
    mismatchtype = LOWQ;
    mismatch = MISMATCH_LOWQ;
    jumptype = CDNA_LOWQ;
    open = CDNA_OPEN_LOWQ;
    extend = CDNA_EXTEND_LOWQ;
  }

  if (length2 > dynprogR->maxlength1 || length1R > dynprogR->maxlength2) {
    debug(printf("length2 %d or length1R %d is too long.  Returning NULL\n",length2,length1R));
    /*
    revoffset2 = offset2 + length2 - 1;
    queryjump = revoffset1R - offset1L + 1;
    genomejump = revoffset2 - offset2 + 1;
    pairs = Pairpool_push_gapholder(NULL,pairpool,queryjump,genomejump);
    */
    *dynprogindex += (*dynprogindex > 0 ? +1 : -1);
    return (List_T) NULL;
  }

  if (length2 > dynprogL->maxlength1 || length1L > dynprogL->maxlength2) {
    debug(printf("length2 %d or length1L %d is too long.  Returning NULL\n",length2,length1L));
    /*
    revoffset2 = offset2 + length2 - 1;
    queryjump = revoffset1R - offset1L + 1;
    genomejump = revoffset2 - offset2 + 1;
    pairs = Pairpool_push_gapholder(NULL,pairpool,queryjump,genomejump);
    */
    *dynprogindex += (*dynprogindex > 0 ? +1 : -1);
    return (List_T) NULL;
  }

  /* Right side looks like 5' end */
  /* Note: sequence 1 and 2 flipped, because 1 has extramaterial */
  revsequence2 = &(sequence2[length2-1]);
  revsequenceuc2 = &(sequenceuc2[length2-1]);
  revoffset2 = offset2+length2-1;

#ifdef PMAP
  instsequence1 = instantiate_codons(sequence1L,queryaaseq,offset1L,revoffset1R-offset1L+1);
  instsequence1L = instsequence1;
  instrevsequence1R = &(instsequence1[revoffset1R-offset1L]);
#endif

  matrixR = compute_scores_lookup(&directionsR,&jumpR,dynprogR,
#ifdef PMAP
				  revsequence2,instrevsequence1R,
#else
				  revsequence2,revsequence1R,
#endif
				  /*revp*/true,length2,length1R,mismatchtype,jumptype,
				  /*init_jump_penalty_p*/true,
				  extraband_paired,/*widebandp*/true,/*onesidegapp*/false);

  matrixL = compute_scores(&directionsL,&jumpL,dynprogL,
#ifdef PMAP
			   sequence2,instsequence1L,
#else
			   sequence2,sequence1L,
#endif
			   /*revp*/false,length2,length1L,mismatchtype,
			   mismatch,open,extend,/*init_jump_penalty_p*/true,
			   extraband_paired,/*widebandp*/true,/*onesidegapp*/false);

  bridge_cdna_gap(&(*finalscore),&bestrL,&bestrR,&bestcL,&bestcR,matrixL,matrixR,
		  directionsL,directionsR,length2,length1L,length1R,extraband_paired,
		  open,extend,offset1L,revoffset1R);

  nmatches = nmismatches = nopens = nindels = 0;
  pairs = traceback(NULL,&nmatches,&nmismatches,&nopens,&nindels,
		    directionsR,jumpR,bestrR,bestcR,
#ifdef PMAP
		    instrevsequence1R,instrevsequence1R,revsequence2,revsequenceuc2,
#else
		    revsequence1R,revsequenceuc1R,revsequence2,revsequenceuc2,
#endif
		    revoffset1R,revoffset2,pairpool,/*revp*/true,/*cdna_gap_p*/true,
		    cdna_direction,*dynprogindex);

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
      pairs = Pairpool_push(pairs,pairpool,k,revoffset2-bestrR+1,instsequence1L[k-offset1L],SHORTGAP_COMP,' ',
			    *dynprogindex);
#else
      debug(printf("cDNA insertion, Pushing [%d,%d] (%c,-)\n",k,revoffset2-bestrR+1,sequence1L[k-offset1L]));
      pairs = Pairpool_push(pairs,pairpool,k,revoffset2-bestrR+1,sequence1L[k-offset1L],SHORTGAP_COMP,' ',
			    *dynprogindex);
#endif
    }
    debug(printf("\n"));

    /* Add genome insertion, if any */
    for (k = revoffset2-bestrR; k >= offset2+bestrL; k--) {
      debug(printf("genome insertion, Pushing [%d,%d] (-,%c)\n",offset1L+bestcL,k,sequence2[k-offset2]));
      pairs = Pairpool_push(pairs,pairpool,offset1L+bestcL,k,' ',SHORTGAP_COMP,sequence2[k-offset2],
			    *dynprogindex);
    }
    debug(printf("\n"));

  } else {

    /* Add gapholder to be solved in the future */
#ifndef NOGAPHOLDER
    debug(printf("Pushing a gap with queryjump = %d, genomejump = %d\n",queryjump,genomejump));
    pairs = Pairpool_push_gapholder(pairs,pairpool,/*queryjump*/UNKNOWNJUMP,/*genomejump*/UNKNOWNJUMP);
#endif
    *incompletep = true;
  }

  pairs = traceback(pairs,&nmatches,&nmismatches,&nopens,&nindels,
		    directionsL,jumpL,bestrL,bestcL,
#ifdef PMAP
		    instsequence1L,instsequence1L,sequence2,sequenceuc2,
#else
		    sequence1L,sequenceuc1L,sequence2,sequenceuc2,
#endif
		    offset1L,offset2,pairpool,/*revp*/false,/*cdna_gap_p*/true,
		    cdna_direction,*dynprogindex);

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
		    double defect_rate, int maxpeelback, bool halfp, bool finalp) {
  List_T pairs = NULL;
#ifdef PMAP
  char *instsequence1, *instrevsequence1;
#else
  char *revsequence1, *revsequenceuc1;
#endif
  Mismatchtype_T mismatchtype;
  Jumptype_T jumptype;
  int **matrixL, **matrixR, **jumpL, **jumpR;
  int canonical_reward;
  Direction_T **directionsL, **directionsR;
  int revoffset1, bestrL, bestrR, bestcL, bestcR;
  int maxhorizjump, maxvertjump;
  /* int queryjump, genomejump; */

  debug(
	printf("\n");
	printf("%c:  ",*dynprogindex > 0 ? (*dynprogindex-1)%26+'a' : (-(*dynprogindex)-1)%26+'A');
	printf("Aligning genome gap with cdna_direction %d\n",cdna_direction);
	printf("At query offset %d-%d, %.*s\n",offset1,offset1+length1-1,length1,sequence1);
	printf("At genomic offset %d-%d, %.*s\n",offset2L,offset2L+length2L-1,length2L,sequence2L);
	printf("At genomic offset %d-%d, %.*s\n",revoffset2R-length2R+1,revoffset2R,length2R,&(revsequence2R[-length2R+1]));
	printf("\n");
	);

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
  if (length1 <= 1) {
    *finalscore = -100000;
    return NULL;
  }

  if (defect_rate < DEFECT_HIGHQ) {
    mismatchtype = HIGHQ;
    if (length1 > maxpeelback * 4) {
      debug(printf("length1 %d is greater than maxpeelback %d * 4, so using single gap penalties\n",
		   length1,maxpeelback));
      jumptype = SINGLE_HIGHQ;
    } else {
      jumptype = PAIRED_HIGHQ;
    }
    if (finalp == true) {
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
      jumptype = SINGLE_MEDQ;
    } else {
      jumptype = PAIRED_MEDQ;
    }
    if (finalp == true) {
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
      jumptype = SINGLE_LOWQ;
    } else {
      jumptype = PAIRED_LOWQ;
    }
    if (finalp == true) {
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
    pairs = Pairpool_push_gapholder(NULL,pairpool,queryjump,genomejump);
    */
#endif
    *dynprogindex += (*dynprogindex > 0 ? +1 : -1);
    *finalscore = -100000;
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
    pairs = Pairpool_push_gapholder(NULL,pairpool,queryjump,genomejump);
    */
#endif
    *dynprogindex += (*dynprogindex > 0 ? +1 : -1);
    *finalscore = -1000000;
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

  matrixL = compute_scores_lookup(&directionsL,&jumpL,dynprogL,
#ifdef PMAP
				  instsequence1,sequence2L,
#else
				  sequence1,sequence2L,
#endif
				  /*revp*/false,length1,length2L,mismatchtype,jumptype,
				  /*init_jump_penalty_p*/true,
				  extraband_paired,/*widebandp*/true,/*onesidegapp*/false);

  matrixR = compute_scores_lookup(&directionsR,&jumpR,dynprogR,
#ifdef PMAP
				  instrevsequence1,revsequence2R,
#else
				  revsequence1,revsequence2R,
#endif
				  /*revp*/true,length1,length2R,mismatchtype,jumptype,
				  /*init_jump_penalty_p*/true,
				  extraband_paired,/*widebandp*/true,/*onesidegapp*/false);

  bridge_intron_gap(&(*finalscore),&bestrL,&bestrR,&bestcL,&bestcR,&(*introntype),matrixL,matrixR,
		    directionsL,directionsR,jumpL,jumpR,sequenceuc2L,revsequenceuc2R,
		    length1,length2L,length2R,cdna_direction,extraband_paired,
		    canonical_reward,maxhorizjump,maxvertjump,offset2L,revoffset2R,halfp,finalp);

  *new_leftgenomepos = offset2L+(bestcL-1);
  *new_rightgenomepos = revoffset2R-(bestcR-1);
  debug(printf("New leftgenomepos = %d, New rightgenomepos = %d\n",*new_leftgenomepos,*new_rightgenomepos));

  *exonhead = revoffset1-(bestrR-1);

  pairs = traceback(NULL,&(*nmatches),&(*nmismatches),&(*nopens),&(*nindels),
		    directionsR,jumpR,bestrR,bestcR,
#ifdef PMAP
		    instrevsequence1,instrevsequence1,revsequence2R,revsequenceuc2R,
#else
		    revsequence1,revsequenceuc1,revsequence2R,revsequenceuc2R,
#endif
		    revoffset1,revoffset2R,pairpool,/*revp*/true,/*cdna_gap_p*/false,
		    cdna_direction,*dynprogindex);
  pairs = List_reverse(pairs);

  /* queryjump = (revoffset1-bestrR) - (offset1+bestrL) + 1; */
  /* genomejump = (revoffset2R-bestcR) - (offset2L+bestcL) + 1; */
  /* No need to revise queryjump or genomejump, because the above
     coordinates are internal to the gap. */

  debug(printf("Pushing a gap\n"));
#ifndef NOGAPHOLDER
  /* pairs = Pairpool_push_gapholder(pairs,pairpool,queryjump,genomejump); */
  pairs = Pairpool_push_gapholder(pairs,pairpool,/*queryjump*/UNKNOWNJUMP,/*genomejump*/UNKNOWNJUMP);
#endif

  pairs = traceback(pairs,&(*nmatches),&(*nmismatches),&(*nopens),&(*nindels),
		    directionsL,jumpL,bestrL,bestcL,
#ifdef PMAP
		    instsequence1,instsequence1,sequence2L,sequenceuc2L,
#else
		    sequence1,sequenceuc1,sequence2L,sequenceuc2L,
#endif
		    offset1,offset2L,pairpool,/*revp*/false,/*cdna_gap_p*/false,
		    cdna_direction,*dynprogindex);

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


List_T
Dynprog_end5_gap (int *dynprogindex, int *finalscore, int *nmatches, int *nmismatches, 
		  int *nopens, int *nindels, T dynprog, 
		  char *revsequence1, char *revsequenceuc1,
		  char *revsequence2, char *revsequenceuc2,
		  int length1, int length2, int revoffset1, int revoffset2, 
#ifdef PMAP
		  char *queryaaseq,
#endif
		  int cdna_direction, Pairpool_T pairpool, int extraband_end, double defect_rate) {
  List_T pairs = NULL;
#ifdef PMAP
  char *instsequence1, *instrevsequence1;
#endif
  Pair_T pair;
  Mismatchtype_T mismatchtype;
  Jumptype_T jumptype;
  int **matrix, **jump;
  int bestr, bestc;
  Direction_T **directions;
#ifdef PMAP
  int initpos, initmod;
#endif

  debug(
	printf("%c:  ",*dynprogindex > 0 ? (*dynprogindex-1)%26+'a' : (-(*dynprogindex)-1)%26+'A');
	printf("Aligning 5' end gap\n")
	);

  mismatchtype = END;
  if (defect_rate < DEFECT_HIGHQ) {
    jumptype = END_HIGHQ; 
  } else if (defect_rate < DEFECT_MEDQ) {
    jumptype = END_MEDQ;
  } else {
    jumptype = END_LOWQ;
  }

  /* We can just chop lengths to work, since we're not constrained on 5' end */
  if (length1 > dynprog->maxlength1) {
    debug(printf("length1 %d is too long.  Chopping to %d\n",length1,dynprog->maxlength1));
    length1 = dynprog->maxlength1;
  }
  if (length2 > dynprog->maxlength2) {
    debug(printf("length2 %d is too long.  Chopping to %d\n",length2,dynprog->maxlength2));
    length2 = dynprog->maxlength2;
  }

#ifdef PMAP
  instsequence1 = instantiate_codons(&(revsequence1[-length1+1]),queryaaseq,revoffset1-length1+1,length1);
  instrevsequence1 = &(instsequence1[length1-1]);
  debug(printf("At query offset %d-%d, %.*s\n",revoffset1-length1+1,revoffset1,length1,&(instrevsequence1[-length1+1])));
#else
  debug(printf("At query offset %d-%d, %.*s\n",revoffset1-length1+1,revoffset1,length1,&(revsequence1[-length1+1])));
#endif

  debug(printf("At genomic offset %d-%d, %.*s\n\n",revoffset2-length2+1,revoffset2,length2,&(revsequence2[-length2+1])));


  /* Can set extraband_end to zero if we are not allowing gaps */
  matrix = compute_scores_lookup(&directions,&jump,dynprog,
#ifdef PMAP
				 instrevsequence1,revsequence2,
#else
				 revsequence1,revsequence2,
#endif
				 /*revp*/true,length1,length2,mismatchtype,jumptype,
				 /*init_jump_penalty_p*/false,
				 extraband_end,/*widebandp*/true,/*onesidegapp*/true);

  find_best_endpoint(&(*finalscore),&bestr,&bestc,matrix,length1,length2,extraband_end);

#ifdef PMAP
  initpos = revoffset1-(bestc-1);
  debug(printf("Initial query pos is %d\n",initpos));
  if ((initmod = initpos % 3) > 0) {
    if (bestr + initmod < length1 && bestc + initmod < length2) {
      debug(printf("Rounding down by %d\n",initmod));
      bestr += initmod;
      bestc += initmod;
    }
  }
#endif

  *nmatches = *nmismatches = *nopens = *nindels = 0;
  pairs = traceback(NULL,&(*nmatches),&(*nmismatches),&(*nopens),&(*nindels),
		    directions,jump,bestr,bestc,
#ifdef PMAP
		    instrevsequence1,instrevsequence1,revsequence2,revsequenceuc2,
#else
		    revsequence1,revsequenceuc1,revsequence2,revsequenceuc2,
#endif
		    revoffset1,revoffset2,pairpool,/*revp*/true,/*cdna_gap_p*/false,
		    cdna_direction,*dynprogindex);

  if ((*nmatches + 1) >= *nmismatches) {
    /* Add 1 to count the match already in the alignment */
    pairs = List_reverse(pairs); /* Look at 5' end to remove excess gaps */
    while (pairs != NULL && (pair = List_head(pairs)) && pair->comp == INDEL_COMP) {
      pairs = List_next(pairs);
    }
  } else {
    *finalscore = 0;
    /* No need to free pairs */
    pairs = NULL;
  }

  /*
    Directions_free(directions);
    Matrix_free(matrix);
  */
  
#ifdef PMAP
  FREE(instsequence1);
#endif

  debug(printf("End of dynprog end5 gap\n"));

  *dynprogindex += (*dynprogindex > 0 ? +1 : -1);
  return List_reverse(pairs);
}


List_T
Dynprog_end3_gap (int *dynprogindex, int *finalscore, int *nmatches, int *nmismatches, 
		  int *nopens, int *nindels, T dynprog, 
		  char *sequence1, char *sequenceuc1,
		  char *sequence2, char *sequenceuc2,
		  int length1, int length2, int offset1, int offset2, 
#ifdef PMAP
		  char *queryaaseq,
#endif
		  int cdna_direction, Pairpool_T pairpool, int extraband_end, double defect_rate) {
  List_T pairs = NULL;
#ifdef PMAP
  char *instsequence1;
#endif
  Pair_T pair;
  Mismatchtype_T mismatchtype;
  Jumptype_T jumptype;
  int **matrix, **jump;
  int bestr, bestc;
  Direction_T **directions;
#ifdef PMAP
  int termpos, termmod;
#endif

  debug(
	printf("%c:  ",*dynprogindex > 0 ? (*dynprogindex-1)%26+'a' : (-(*dynprogindex)-1)%26+'A');
	printf("Aligning 3' end gap\n")
	);

  mismatchtype = END;
  if (defect_rate < DEFECT_HIGHQ) {
    jumptype = END_HIGHQ;
  } else if (defect_rate < DEFECT_MEDQ) {
    jumptype = END_MEDQ;
  } else {
    jumptype = END_LOWQ;
  }

  /* We can just chop lengths to work, since we're not constrained on 3' end */
  if (length1 > dynprog->maxlength1) {
    debug(printf("length1 %d is too long.  Chopping to %d\n",length1,dynprog->maxlength1));
    length1 = dynprog->maxlength1;
  }
  if (length2 > dynprog->maxlength2) {
    debug(printf("length2 %d is too long.  Chopping to %d\n",length2,dynprog->maxlength2));
    length2 = dynprog->maxlength2;
  }

#ifdef PMAP
  instsequence1 = instantiate_codons(sequence1,queryaaseq,offset1,length1);
  debug(printf("At query offset %d-%d, %.*s\n",offset1,offset1+length1-1,length1,instsequence1));
#else
  debug(printf("At query offset %d-%d, %.*s\n",offset1,offset1+length1-1,length1,sequence1));
#endif
	
  debug(printf("At genomic offset %d-%d, %.*s\n\n",offset2,offset2+length2-1,length2,sequence2));


  /* Can set extraband_end to zero if we are not allowing gaps */
  matrix = compute_scores_lookup(&directions,&jump,dynprog,
#ifdef PMAP
				 instsequence1,sequence2,
#else
				 sequence1,sequence2,
#endif
				 /*revp*/false,length1,length2,mismatchtype,jumptype,
				 /*init_jump_penalty_p*/false,
				 extraband_end,/*widebandp*/true,/*onesidegapp*/true);

  find_best_endpoint(&(*finalscore),&bestr,&bestc,matrix,length1,length2,extraband_end);

#ifdef PMAP
  termpos = offset1+(bestc-1);
  debug(printf("Final query pos is %d\n",termpos));
  if ((termmod = termpos % 3) < 2) {
    if (bestr + (2 - termmod) < length1 && bestc + (2 - termmod) < length2) {
      debug(printf("Rounding up by %d\n",2 - termmod));
      bestr += 2 - termmod;
      bestc += 2 - termmod;
    }
  }
#endif

  *nmatches = *nmismatches = *nopens = *nindels = 0;
  pairs = traceback(NULL,&(*nmatches),&(*nmismatches),&(*nopens),&(*nindels),
		    directions,jump,bestr,bestc,
#ifdef PMAP
		    instsequence1,instsequence1,sequence2,sequenceuc2,
#else
		    sequence1,sequenceuc1,sequence2,sequenceuc2,
#endif
		    offset1,offset2,pairpool,/*revp*/false,/*cdna_gap_p*/false,
		    cdna_direction,*dynprogindex);

  if ((*nmatches + 1) >= *nmismatches) {
    /* Add 1 to count the match already in the alignment */
    pairs = List_reverse(pairs); /* Look at 3' end to remove excess gaps */
    while (pairs != NULL && (pair = List_head(pairs)) && pair->comp == INDEL_COMP) {
      pairs = List_next(pairs);
    }
  } else {
    *finalscore = 0;
    /* No need to free pairs */
    pairs = NULL;
  }

  /*
    Directions_free(directions);
    Matrix_free(matrix);
  */

#ifdef PMAP
  FREE(instsequence1);
#endif

  debug(printf("End of dynprog end3 gap\n"));

  *dynprogindex += (*dynprogindex > 0 ? +1 : -1);
  return pairs;			/* not List_reverse(pairs) */
}


int
Dynprog_internal_gap_stats (T dynprog, char *sequenceuc1, char *sequenceuc2,
			    int length1, int length2, int offset1, int offset2, 
#ifdef PMAP
			    char *queryaaseq,
#endif
			    int cdna_direction, int extraband_end, double defect_rate) {
  int finalscore;
#ifdef PMAP
  char *instsequence1;
#endif
  Mismatchtype_T mismatchtype;
  Jumptype_T jumptype;
  int **matrix, **jump;
  int bestr;
  Direction_T **directions;
#ifdef PMAP
  int termpos, termmod;
#endif

  debug(
	printf("Aligning 3' end gap\n")
	);

  mismatchtype = END;
  if (defect_rate < DEFECT_HIGHQ) {
    jumptype = END_HIGHQ;
  } else if (defect_rate < DEFECT_MEDQ) {
    jumptype = END_MEDQ;
  } else {
    jumptype = END_LOWQ;
  }

  /* We can just chop lengths to work, since we're not constrained on 3' end */
  if (length1 > dynprog->maxlength1) {
    debug(printf("length1 %d is too long.  Chopping to %d\n",length1,dynprog->maxlength1));
    length1 = dynprog->maxlength1;
  }
  if (length2 > dynprog->maxlength2) {
    debug(printf("length2 %d is too long.  Chopping to %d\n",length2,dynprog->maxlength2));
    length2 = dynprog->maxlength2;
  }

#ifdef PMAP
  instsequence1 = instantiate_codons(sequenceuc1,queryaaseq,offset1,length1);
  debug(printf("At query offset %d-%d, %.*s\n",offset1,offset1+length1-1,length1,instsequence1));
#else
  debug(printf("At query offset %d-%d, %.*s\n",offset1,offset1+length1-1,length1,sequenceuc1));
#endif
	
  debug(printf("At genomic offset %d-%d, %.*s\n\n",offset2,offset2+length2-1,length2,sequenceuc2));


  /* Can set extraband_end to zero if we are not allowing gaps */
  matrix = compute_scores_lookup(&directions,&jump,dynprog,
#ifdef PMAP
				 instsequence1,sequenceuc2,
#else
				 sequenceuc1,sequenceuc2,
#endif
				 /*revp*/false,length1,length2,mismatchtype,
				 jumptype,/*init_jump_penalty_p*/false,
				 extraband_end,/*widebandp*/true,/*onesidegapp*/true);

  find_best_endpoint_length2(&finalscore,&bestr,matrix,directions,length1,length2,extraband_end);

  if (offset2 == 3945 && offset2+length2-1 == 4110) {
    Matrix_print(matrix,length1,length2,sequenceuc1,sequenceuc2,/*revp*/false);
  }


#if 0
  /* May also need to take out of companion program that produces pairs */
#ifdef PMAP
  termpos = offset1+(bestc-1);
  debug(printf("Final query pos is %d\n",termpos));
  if ((termmod = termpos % 3) < 2) {
    if (bestr + (2 - termmod) < length1 && bestc + (2 - termmod) < length2) {
      debug(printf("Rounding up by %d\n",2 - termmod));
      bestr += 2 - termmod;
      bestc += 2 - termmod;
    }
  }
#endif
#endif

  /*
    Directions_free(directions);
    Matrix_free(matrix);
  */

#ifdef PMAP
  FREE(instsequence1);
#endif

  debug(printf("End of dynprog internal gap stats, score = %d\n",finalscore));

  return finalscore;
}


List_T
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
		    int cdna_direction, Pairpool_T pairpool, int extraband_end, double defect_rate) {
  List_T pairs = NULL;
#ifdef PMAP
  char *instsequence1L, *instsequence1R, *instrevsequence1R;
#endif
  Mismatchtype_T mismatchtype;
  Jumptype_T jumptype;
  int **matrixL, **matrixR, **jumpL, **jumpR;
  int bestrcL, bestrcR;
  int diaglength;
  Direction_T **directionsL, **directionsR;
  int nmatches, nmismatches, nopens, nindels;
  /* int queryjump, genomejump; */

  debug(
	printf("%c:  ",*dynprogindex > 0 ? (*dynprogindex-1)%26+'a' : (-(*dynprogindex)-1)%26+'A');
	printf("Aligning dual break\n")
	);

  mismatchtype = END;
  if (defect_rate < DEFECT_HIGHQ) {
    jumptype = END_HIGHQ;
  } else if (defect_rate < DEFECT_MEDQ) {
    jumptype = END_MEDQ;
  } else {
    jumptype = END_LOWQ;
  }

  if (length1 < length2) {
    diaglength = length1;
  } else {
    diaglength = length2;
  }

  if (diaglength > dynprogL->maxlength1 || diaglength > dynprogL->maxlength2) {
    debug(printf("diaglength %d (vs max %d or max %d) is too long.  Returning NULL\n",
		 diaglength,dynprogL->maxlength1,dynprogL->maxlength2));
#ifndef NOGAPHOLDER
    /*
    queryjump = revoffset1R - offset1L + 1;
    genomejump = revoffset2R - offset2L + 1;
    pairs = Pairpool_push_gapholder(NULL,pairpool,queryjump,genomejump);
    */
#endif
    *dynprogindex += (*dynprogindex > 0 ? +1 : -1);
    return (List_T) NULL;
  }
  if (diaglength > dynprogR->maxlength1 || diaglength > dynprogR->maxlength2) {
    debug(printf("diaglength %d (vs max %d or max %d) is too long.  Returning NULL\n",
		 diaglength,dynprogR->maxlength1,dynprogR->maxlength2));
    /*
    queryjump = revoffset1R - offset1L + 1;
    genomejump = revoffset2R - offset2L + 1;
    pairs = Pairpool_push_gapholder(NULL,pairpool,queryjump,genomejump);
    */
    *dynprogindex += (*dynprogindex > 0 ? +1 : -1);
    return (List_T) NULL;
  }

#ifdef PMAP
  instsequence1L = instantiate_codons(sequence1L,queryaaseq,offset1L,length1);
  instsequence1R = instantiate_codons(&(revsequence1R[-length1+1]),queryaaseq,revoffset1R-length1+1,length1);
  instrevsequence1R = &(instsequence1R[length1-1]);

  debug(printf("At query offset %d-%d, %.*s\n",offset1L,offset1L+length1-1,length1,instsequence1L));
  debug(printf("At query offset %d-%d, %.*s\n",revoffset1R-length1+1,revoffset1R,length1,&(instrevsequence1R[-length1+1])));
#else
  debug(printf("At query offset %d-%d, %.*s\n",offset1L,offset1L+length1-1,length1,sequence1L));
  debug(printf("At query offset %d-%d, %.*s\n",revoffset1R-length1+1,revoffset1R,length1,&(revsequence1R[-length1+1])));
#endif
	
  debug(printf("At genomic offset %d-%d, %.*s\n\n",offset2L,offset2L+length2-1,length2,sequence2L));
  debug(printf("At genomic offset %d-%d, %.*s\n\n",revoffset2R-length2+1,revoffset2R,length2,&(revsequence2R[-length2+1])));

  /* Can set extraband_end to zero if we are not allowing gaps */
  matrixL = compute_scores_lookup(&directionsL,&jumpL,dynprogL,
#ifdef PMAP
				  instsequence1L,sequence2L,
#else
				  sequence1L,sequence2L,
#endif
				  /*revp*/false,diaglength,diaglength,mismatchtype,jumptype,
				  /*init_jump_penalty_p*/true,
				  extraband_end,/*widebandp*/false,/*onesidegapp*/true);

  matrixR = compute_scores_lookup(&directionsR,&jumpR,dynprogR,
#ifdef PMAP
				  instrevsequence1R,revsequence2R,
#else
				  revsequence1R,revsequence2R,
#endif
				  /*revp*/true,diaglength,diaglength,mismatchtype,jumptype,
				  /*init_jump_penalty_p*/true,
				  extraband_end,/*widebandp*/false,/*onesidegapp*/true);

  if (cdna_direction >= 0) {
    bridge_dual_break_fwd(&(*finalscore),&bestrcL,&bestrcR,matrixL,matrixR,diaglength,
			  sequenceuc2L,revsequenceuc2R);
  } else {
    bridge_dual_break_rev(&(*finalscore),&bestrcL,&bestrcR,matrixL,matrixR,diaglength,
			  sequenceuc2L,revsequenceuc2R);
  }

  nmatches = nmismatches = nopens = nindels = 0;
  pairs = traceback(NULL,&nmatches,&nmismatches,&nopens,&nindels,
		    directionsR,jumpR,bestrcR,bestrcR,
#ifdef PMAP
		    instrevsequence1R,instrevsequence1R,revsequence2R,revsequenceuc2R,
#else
		    revsequence1R,revsequenceuc1R,revsequence2R,revsequenceuc2R,
#endif
		    revoffset1R,revoffset2R,pairpool,/*revp*/true,/*cdna_gap_p*/false,
		    cdna_direction,*dynprogindex);
  pairs = List_reverse(pairs);

  /* queryjump = (revoffset1R-bestrcR) - (offset1L+bestrcL) + 1; */
  /* genomejump = (revoffset2R-bestrcR) - (offset2L+bestrcL) + 1; */
  /* No need to revise queryjump or genomejump, because the above
     coordinates are internal to the gap. */

  debug(printf("Pushing a gap\n"));
  pairs = Pairpool_push_gapholder(pairs,pairpool,/*queryjump*/UNKNOWNJUMP,/*genomejump*/UNKNOWNJUMP);
  /*
  pair = List_head(pairs);
  pair->comp = DUALBREAK_COMP;
  */

  pairs = traceback(pairs,&nmatches,&nmismatches,&nopens,&nindels,
		    directionsL,jumpL,bestrcL,bestrcL,
#ifdef PMAP
		    instsequence1L,instsequence1L,sequence2L,sequenceuc2L,
#else
		    sequence1L,sequenceuc1L,sequence2L,sequenceuc2L,
#endif
		    offset1L,offset2L,pairpool,/*revp*/false,/*cdna_gap_p*/false,
		    cdna_direction,*dynprogindex);

  if (List_length(pairs) == 1) {
    /* Only a gapholder was added */
    pairs = (List_T) NULL;
  }

  /*
  Directions_free(directionsR);
  Directions_free(directionsL);
  Matrix_free(matrixR);
  Matrix_free(matrixL);
  */

#ifdef PMAP
  FREE(instsequence1R);
  FREE(instsequence1L);
#endif

  debug(printf("End of dynprog dual break\n"));

  *dynprogindex += (*dynprogindex > 0 ? +1 : -1);
  return List_reverse(pairs);
}




static const Except_T microexon_error = {"Microexon error"};

static List_T
make_microexon_pairs_double (int offset1L, int offset1M, int offset1R,
			     int offset2L, int offset2M, int offset2R,
			     int lengthL, int lengthM, int lengthR,
			     char *queryseq, char *queryuc, char *genomicseg, char *genomicuc,
			     Pairpool_T pairpool, char gapchar, int dynprogindex) {
  List_T pairs = NULL;
  char c1, c2;
  int i;

  /* Left segment */
  for (i = 0; i < lengthL; i++) {
    c1 = queryseq[offset1L+i];
    c2 = genomicseg[offset2L+i];
    if (queryuc[offset1L+i] == genomicuc[offset2L+i]) {
      pairs = Pairpool_push(pairs,pairpool,offset1L+i,offset2L+i,c1,DYNPROG_MATCH_COMP,c2,
			    dynprogindex);
#ifdef PMAP
    } else if (consistent_array[(int) c1][(int) c2] == true) {
      pairs = Pairpool_push(pairs,pairpool,offset1L+i,offset2L+i,c1,AMBIGUOUS_COMP,c2,
			    dynprogindex);
#endif
    } else {
      pairs = Pairpool_push(pairs,pairpool,offset1L+i,offset2L+i,c1,MISMATCH_COMP,c2,
			    dynprogindex);
    }
  }

  /* First gap */
  /* Don't have to adjust querypos/genomepos, since no cdna/genome skips allowed */
#if 0
  pairs = Pairpool_push_gapholder(pairs,pairpool,/*queryjump*/0,offset2M-(offset2L+lengthL));
#endif
  pairs = Pairpool_push_gapholder(pairs,pairpool,/*queryjump*/UNKNOWNJUMP,/*genomejump*/UNKNOWNJUMP);
  
  /* Microexon */
  for (i = 0; i < lengthM; i++) {
    c1 = queryseq[offset1M+i];
    c2 = genomicseg[offset2M+i];
    if (queryuc[offset1M+i] == genomicuc[offset2M+i]) {
      pairs = Pairpool_push(pairs,pairpool,offset1M+i,offset2M+i,c1,DYNPROG_MATCH_COMP,c2,
			    dynprogindex);
#ifdef PMAP
    } else if (consistent_array[(int) c1][(int) c2] == true) {
      pairs = Pairpool_push(pairs,pairpool,offset1M+i,offset2M+i,c1,AMBIGUOUS_COMP,c2,
			    dynprogindex);
#endif
    } else {
      pairs = Pairpool_push(pairs,pairpool,offset1M+i,offset2M+i,c1,MISMATCH_COMP,c2,
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
    pairs = Pairpool_push_gapholder(pairs,pairpool,/*queryjump*/0,offset2R-(offset2M+lengthM));
#endif
    pairs = Pairpool_push_gapholder(pairs,pairpool,/*queryjump*/UNKNOWNJUMP,/*genomejump*/UNKNOWNJUMP);
  }
  
  /* Right segment */
  for (i = 0; i < lengthR; i++) {
    c1 = queryseq[offset1R+i];
    c2 = genomicseg[offset2R+i];
    if (queryuc[offset1R+i] == genomicuc[offset2R+i]) {
      pairs = Pairpool_push(pairs,pairpool,offset1R+i,offset2R+i,c1,DYNPROG_MATCH_COMP,c2,
			    dynprogindex);
#ifdef PMAP
    } else if (consistent_array[(int) c1][(int) c2] == true) {
      pairs = Pairpool_push(pairs,pairpool,offset1R+i,offset2R+i,c1,AMBIGUOUS_COMP,c2,
			    dynprogindex);
#endif
    } else {
      pairs = Pairpool_push(pairs,pairpool,offset1R+i,offset2R+i,c1,MISMATCH_COMP,c2,
			    dynprogindex);
    }
  }

  return pairs;
}


static List_T
make_microexon_pairs_single (int offset1L, int offset1R,
			     int offset2L, int offset2R,
			     int lengthL, int lengthR,
			     char *queryseq, char *queryuc, char *genomicseg, char *genomicuc,
			     Pairpool_T pairpool, char gapchar, int dynprogindex) {
  List_T pairs = NULL;
  Pair_T gappair;
  char c1, c2;
  int i;

  /* Microexon */
  for (i = 0; i < lengthL; i++) {
    c1 = queryseq[offset1L+i];
    c2 = genomicseg[offset2L+i];
    if (queryuc[offset1L+i] == genomicuc[offset2L+i]) {
      pairs = Pairpool_push(pairs,pairpool,offset1L+i,offset2L+i,c1,DYNPROG_MATCH_COMP,c2,
			    dynprogindex);
#ifdef PMAP
    } else if (consistent_array[(int) c1][(int) c2] == true) {
      pairs = Pairpool_push(pairs,pairpool,offset1L+i,offset2L+i,c1,AMBIGUOUS_COMP,c2,
			    dynprogindex);
#endif
    } else {
      pairs = Pairpool_push(pairs,pairpool,offset1L+i,offset2L+i,c1,MISMATCH_COMP,c2,
			    dynprogindex);
    }
  }

  /* Gap */
  /* Don't have to adjust querypos/genomepos, since no cdna/genome skips allowed */
#if 0
  pairs = Pairpool_push_gapholder(pairs,pairpool,/*queryjump*/0,offset2R-(offset2L+lengthL));
#endif
  pairs = Pairpool_push_gapholder(pairs,pairpool,/*queryjump*/UNKNOWNJUMP,/*genomejump*/UNKNOWNJUMP);


  /* Assign pair->comp, because might occur after assign_gap_types */
  gappair = (Pair_T) List_head(pairs);
  gappair->comp = gapchar;
  
  /* Right segment */
  for (i = 0; i < lengthR; i++) {
    c1 = queryseq[offset1R+i];
    c2 = genomicseg[offset2R+i];
    if (queryuc[offset1R+i] == genomicuc[offset2R+i]) {
      pairs = Pairpool_push(pairs,pairpool,offset1R+i,offset2R+i,c1,DYNPROG_MATCH_COMP,c2,
			    dynprogindex);
#ifdef PMAP
    } else if (consistent_array[(int) c1][(int) c2] == true) {
      pairs = Pairpool_push(pairs,pairpool,offset1R+i,offset2R+i,c1,DYNPROG_MATCH_COMP,c2,
			    dynprogindex);
#endif
    } else {
      pairs = Pairpool_push(pairs,pairpool,offset1R+i,offset2R+i,c1,MISMATCH_COMP,c2,
			    dynprogindex);
    }
  }

  return pairs;
}

List_T
Dynprog_microexon_int (int *dynprogindex, int *microintrontype, char *sequence1, char *sequenceuc1,
		       char *sequence2L, char *sequenceuc2L,
		       char *revsequence2R, char *revsequenceuc2R,
		       int length1, int length2L, int length2R,
		       int offset1, int offset2L, int revoffset2R, int cdna_direction,
#ifdef PMAP
		       char *queryaaseq,
#endif
		       char *queryseq, char *queryuc, char *genomicseg, char *genomicuc,
		       Pairpool_T pairpool, double defect_rate) {
  List_T pairs = NULL;
#ifdef PMAP
  char *instsequence1;
#endif
  Intlist_T hits = NULL, p;
  int middlelength, cL, cR, mincR, maxcR, leftbound, rightbound, textleft, textright, candidate, i;
  int min_microexon_length, span, nmismatches;
  char left1, left2, right2, right1;
  char intron1, intron2, intron3, intron4, gapchar;
  double pvalue;

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
    abort();
  }
#endif

  debug(printf("Begin microexon search for %.*s and %.*s\n",
	       length2L,sequence2L,length2R,&(revsequence2R[-length2R+1])));
#ifdef PMAP
  instsequence1 = instantiate_codons(sequence1,queryaaseq,offset1,length1);
  debug(printf("  Query sequence is %.*s\n",length1,instsequence1));
#else
  debug(printf("  Query sequence is %.*s\n",length1,sequence1));
#endif
  span = revoffset2R-offset2L;
  debug(printf("  Genomic span is of length %d\n",span));

  if (span <= 0) {
    fprintf(stderr,"Bug in Dynprog_microexon_int.  span = %d.  Please report to twu@gene.com\n",span);
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
    debug(printf("  Comparing %c with %c\n",instsequence1[leftbound],sequenceuc2L[leftbound]));
    if (matchtable[instsequence1[leftbound]-'A'][sequenceuc2L[leftbound]-'A'] == false) {
      nmismatches;
    }
#else
    debug(printf("  Comparing %c with %c\n",sequenceuc1[leftbound],sequenceuc2L[leftbound]));
    if (sequenceuc1[leftbound] != sequenceuc2L[leftbound]) {
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
    debug(printf("  Comparing %c with %c\n",instsequence1[i],revsequenceuc2R[-rightbound]));
    if (matchtable[instsequence1[i]-'A'][revsequenceuc2R[-rightbound]-'A'] == false) {
      nmismatches++;
    }
#else
    debug(printf("  Comparing %c with %c\n",sequenceuc1[i],revsequenceuc2R[-rightbound]));
    if (sequenceuc1[i] != revsequenceuc2R[-rightbound]) {
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
    left1 = sequenceuc2L[cL];
    left2 = sequenceuc2L[cL+1];
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
	right2 = revsequenceuc2R[-cR-1];
	right1 = revsequenceuc2R[-cR];
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
	  hits = BoyerMoore(&(sequenceuc1[cL]),middlelength,&(genomicuc[textleft]),textright-textleft);
#endif
	  for (p = hits; p != NULL; p = Intlist_next(p)) {
	    candidate = textleft + Intlist_head(p);
	    if (genomicuc[candidate - 2] == intron3 &&
		genomicuc[candidate - 1] == intron4 &&
		genomicuc[candidate + middlelength] == intron1 &&
		genomicuc[candidate + middlelength + 1] == intron2) {
	      debug(printf("  Successful microexon at %d >>> %d..%d >>> %d\n",offset2L+cL,candidate,candidate+middlelength,revoffset2R-cR));
	      
	      Intlist_free(&hits);
	      pairs = make_microexon_pairs_double(offset1,offset1+cL,offset1+cL+middlelength,
						  offset2L,candidate,revoffset2R-cR+1,
						  /*lengthL*/cL,/*lengthM*/middlelength,/*lengthR*/cR,
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
    }
  }

#ifdef PMAP
  FREE(instsequence1);
#endif

  debug(printf("End of dynprog microexon int\n"));

  *microintrontype = NONINTRON;
  return NULL;
}


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
    abort();
  }
#endif

  debug(printf("Begin microexon search at 5' for %.*s\n",
	       length2,&(revsequence2[-length2+1])));
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

  debug(printf("End of dynprog microexon 5\n"));

  *microintrontype = NONINTRON;
  return NULL;
}


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
    abort();
  }
#endif

  debug(printf("Begin microexon search at 3' for %.*s\n",length2,sequence2));

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
