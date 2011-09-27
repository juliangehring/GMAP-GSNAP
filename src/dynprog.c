static char rcsid[] = "$Id: dynprog.c,v 1.124 2005/10/14 19:34:29 twu Exp $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "dynprog.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>		/* For ceil, log, pow */
#include "bool.h"
#include "mem.h"
#include "comp.h"
#include "pair.h"
#include "pairdef.h"
#include "boyer-moore.h"

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

/*
#define RIGHTANGLE 1
*/

#define MICROEXON_PVALUE 0.01
#define ENDSEQUENCE_PVALUE 0.001 /* Have stricter threshold for making end exons */

#define MIN_MICROEXON_LENGTH 3
#define MAX_MICROEXON_LENGTH 12	/* Should be oligomer length - 1 plus peelback */
#define MICROINTRON_LENGTH 9

#define MIN_NONINTRON 50

#define FULLMATCH 3
#define HALFMATCH 1
#define AMBIGUOUS 0

typedef enum {HIGHQ, MEDQ, LOWQ, END} Mismatchtype_T;
#define NMISMATCHTYPES 4

#define MISMATCH_HIGHQ -10
#define MISMATCH_MEDQ -9
#define MISMATCH_LOWQ -8

/* Allow lower mismatch scores on end to allow more complete
   alignments to the end, and because ends are typically of lower
   quality.  Make equal to FULLMATCH, because criterion is nmatches >=
   nmismatches. */
#define MISMATCH_END -3


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

#define PAIRED_OPEN_HIGHQ -17
#define PAIRED_OPEN_MEDQ -15
#define PAIRED_OPEN_LOWQ -13

#define PAIRED_EXTEND_HIGHQ -7
#define PAIRED_EXTEND_MEDQ -6
#define PAIRED_EXTEND_LOWQ -5

#define SINGLE_OPEN_HIGHQ -17
#define SINGLE_OPEN_MEDQ -15
#define SINGLE_OPEN_LOWQ -13

#define SINGLE_EXTEND_HIGHQ -7
#define SINGLE_EXTEND_MEDQ -6
#define SINGLE_EXTEND_LOWQ -5


/* cDNA insertions are biologically not meaningful, so look for a good
   gap opening somewhere */
#define CDNA_OPEN_HIGHQ -10
#define CDNA_OPEN_MEDQ -9
#define CDNA_OPEN_LOWQ -8

#define CDNA_EXTEND_HIGHQ -7
#define CDNA_EXTEND_MEDQ -6
#define CDNA_EXTEND_LOWQ -5

/* Ends tend to be of lower quality, so we drop the scores a notch */
#define END_OPEN_HIGHQ -9
#define END_OPEN_MEDQ -7
#define END_OPEN_LOWQ -5

#define END_EXTEND_HIGHQ -6
#define END_EXTEND_MEDQ -5
#define END_EXTEND_LOWQ -4


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

#define CANONICAL_INTRON_HIGHQ 28 /* GT-AG */
#define CANONICAL_INTRON_MEDQ  31
#define CANONICAL_INTRON_LOWQ  34

/* Prefer alternate intron to other non-canonicals, but don't
   introduce mismatches or gaps to identify */
#define ALTERNATE_INTRON 10	/* GC-AG or AT-AC */

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

  matrix = ptrs;
  matrix[0] = space;
  for (i = 1; i <= length1; i++) {
    matrix[i] = matrix[i-1] + (length2 + 1);
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
  memset(space,0,(length1+1)*(length2+1)*sizeof(int));

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
    directions[i] = directions[i-1] + (length2 + 1);
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
  memset(space,0,(length1+1)*(length2+1)*sizeof(Direction_T));

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
  *maxlength2 = *maxlength1 + extraquerygap;
  if (extramaterial_end > extramaterial_paired) {
    *maxlength2 += extramaterial_end;
  } else {
    *maxlength2 += extramaterial_paired;
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

static T
Dynprog_temp (int maxlength1, int maxlength2) {
  T new = (T) MALLOC(sizeof(*new));

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


static int **pairdistance_array[NMISMATCHTYPES];
static bool **consistent_array;

int
Dynprog_pairdistance (int c1, int c2) {
  return pairdistance_array[HIGHQ][c1][c2];
}

static void
permute_cases (int NA1, int NA2, int score) {
  int i;
  int na1, na2;

  na1 = tolower(NA1);
  na2 = tolower(NA2);

  consistent_array[na1][na2] = true;
  consistent_array[na1][NA2] = true;
  consistent_array[NA1][na2] = true;
  consistent_array[NA1][NA2] = true;

  consistent_array[na2][na1] = true;
  consistent_array[na2][NA1] = true;
  consistent_array[NA2][na1] = true;
  consistent_array[NA2][NA1] = true;
  for (i = 0; i < NMISMATCHTYPES; i++) {
    pairdistance_array[i][na1][na2] = score;
    pairdistance_array[i][na1][NA2] = score;
    pairdistance_array[i][NA1][na2] = score;
    pairdistance_array[i][NA1][NA2] = score;

    pairdistance_array[i][na2][na1] = score;
    pairdistance_array[i][na2][NA1] = score;
    pairdistance_array[i][NA2][na1] = score;
    pairdistance_array[i][NA2][NA1] = score;
  }

  return;
}


static void
pairdistance_init () {
  int i, j;
  int c, c1, c2;

  consistent_array = (bool **) CALLOC(128,sizeof(bool *));
  for (j = 0; j < 128; j++) {
    consistent_array[j] = (bool *) CALLOC(128,sizeof(bool));
  }
  for (i = 0; i < NMISMATCHTYPES; i++) {
    pairdistance_array[i] = (int **) CALLOC(128,sizeof(int *));
    for (j = 0; j < 128; j++) {
      pairdistance_array[i][j] = (int *) CALLOC(128,sizeof(int));
    }
  }

  for (c1 = 'A'; c1 <= 'z'; c1++) {
    for (c2 = 'A'; c2 < 'z'; c2++) {
      pairdistance_array[HIGHQ][c1][c2] = MISMATCH_HIGHQ;
      pairdistance_array[MEDQ][c1][c2] = MISMATCH_MEDQ;
      pairdistance_array[LOWQ][c1][c2] = MISMATCH_LOWQ;
      pairdistance_array[END][c1][c2] = MISMATCH_END;
    }
  }

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
  int ncodons;

  ncodons = (length - 1)/3;
  return open + extend + ncodons*3*extend;
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
  int i, j;

  for (i = 0; i < NJUMPTYPES; i++) {
    FREE(jump_penalty_array[i]);
  }

  for (i = 0; i < NMISMATCHTYPES; i++) {
    for (j = 0; j < 128; j++) {
      FREE(pairdistance_array[i][j]);
    }
    FREE(pairdistance_array[i]);
  }
  for (j = 0; j < 128; j++) {
    FREE(consistent_array[j]);
  }
  FREE(consistent_array);

  return;
}

  /************************************************************************/

static int **
compute_scores_lookup (Direction_T ***directions, int ***jump, T this, 
		       char *sequence1, char *sequence2, bool revp, 
		       int length1, int length2, 
		       Mismatchtype_T mismatchtype, Jumptype_T jumptype,
		       int extraband) {
  int **matrix;
  int r, c, r1, c1, na1, na2;
  int bestscore, score, bestjump, j;
  Direction_T bestdir;
  int lband, rband, clo, chigh, rlo;
  int **pairdistance_array_type;
  int *jump_penalty_array_type;

  pairdistance_array_type = pairdistance_array[mismatchtype];
  jump_penalty_array_type = jump_penalty_array[jumptype];

  if (length2 >= length1) {
    rband = length2 - length1 + extraband;
    lband = extraband;
  } else {
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
  for (c = 1; c <= rband && c <= length2; c++) {
    matrix[0][c] = jump_penalty_array_type[c];
    (*directions)[0][c] = HORIZ;
    (*jump)[0][c] = c;
  }

  /* Column 0 initialization */
  for (r = 1; r <= lband && r <= length1; r++) {
    matrix[r][0] = jump_penalty_array_type[r];
    (*directions)[r][0] = VERT;
    (*jump)[r][0] = r;
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

      /* Vertical case */
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
		Mismatchtype_T mismatchtype, Jumptype_T jumptype,
		int mismatch, int open, int extend, int extraband) {
  int **matrix;
  int r, c, r1, c1, na1, na2;
  int bestscore, score, bestjump, j;
  Direction_T bestdir;
  int lband, rband, clo, chigh, rlo;
  int **pairdistance_array_type;

  pairdistance_array_type = pairdistance_array[mismatchtype];

  if (length2 >= length1) {
    rband = length2 - length1 + extraband;
    lband = extraband;
  } else {
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
  for (c = 1; c <= rband && c <= length2; c++) {
    matrix[0][c] = jump_penalty(c,open,extend);
    (*directions)[0][c] = HORIZ;
    (*jump)[0][c] = c;
  }

  /* Column 0 initialization */
  for (r = 1; r <= lband && r <= length1; r++) {
    matrix[r][0] = jump_penalty(r,open,extend);
    (*directions)[r][0] = VERT;
    (*jump)[r][0] = r;
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

      /* Vertical case */
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

#if 0
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
#endif

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


#if 0
/* Finds best score along diagonal.  Not used anymore. */
static void
find_best_endpoint_nogap (int *finalscore, int *bestr, int *bestc, int **matrix, 
			  int length1, int length2, int extraband_end_or_paired) {
  int bestscore = -1000000;
  int r, c;
  /* No need for loffset or cmid because they apply only for cdnaend
     == FIVE, which doesn't require searching */

  *bestr = *bestc = 0;

  for (r = 1; r <= length1; r++) {
    c = r;
    if (matrix[r][c] >= bestscore) {
      *bestr = r;
      *bestc = c;
      bestscore = matrix[r][c];
    }
  }
  *finalscore = bestscore;
  return;
}
#endif



/************************************************************************
 *  Scheme for matching dinucleotide pairs
 ************************************************************************/

/* Pieces for logical AND */
#define LEFT_GT  0x21		/* 100001 */
#define LEFT_GC	 0x10		/* 010000 */
#define LEFT_AT  0x08		/* 001000 */
#define LEFT_CT  0x06		/* 000110 */

#define RIGHT_AG 0x30		/* 110000 */
#define RIGHT_AC 0x0C		/* 001100 */
#define RIGHT_GC 0x02		/* 000010 */
#define RIGHT_AT 0x01		/* 000001 */


static List_T
add_gap (List_T pairs, int bestr, int bestcL, int bestcR,
	 char *sequence1, char *sequence2L, char *revsequence2R,
	 int length1, int length2L, int length2R, 
	 int offset1, int offset2L, int revoffset2R, 
	 int introntype, int ngap, Pairpool_T pairpool) {
  char comp, c2;
  int intronlength, genomicpos;
  int i;

  intronlength = (revoffset2R-bestcR) - (offset2L+bestcL) + 1;
  switch (introntype) {
  case GTAG_FWD: comp = FWD_CANONICAL_INTRON_COMP; break;
  case GCAG_FWD: comp = FWD_GCAG_INTRON_COMP; break;
  case ATAC_FWD: comp = FWD_ATAC_INTRON_COMP; break;
  case NONINTRON:
    if (intronlength < MIN_NONINTRON) {
      comp = SHORTGAP_COMP;	/* Will be printed as INDEL_COMP, but need to score as NONINTRON_COMP */
    } else {
      comp = NONINTRON_COMP;
    }
    break;
  case ATAC_REV: comp = REV_ATAC_INTRON_COMP; break;
  case GCAG_REV: comp = REV_GCAG_INTRON_COMP; break;
  case GTAG_REV: comp = REV_CANONICAL_INTRON_COMP; break;
  default: fprintf(stderr,"Unexpected intron type %d\n",introntype);
    exit(9);
  }
  debug(printf("Adding gap of type %c\n",comp));

  if (comp == SHORTGAP_COMP) {
    /* Treat as an insertion rather than an intron */
    for (i = 0; i < intronlength; i++) {
      /* c2 = (bestcR+i < length2R) ? revsequence2R[-bestcR-i] : '?'; */
      c2 = revsequence2R[-bestcR-i];
      pairs = Pairpool_push(pairs,pairpool,offset1+bestr,revoffset2R-(bestcR+i),' ',comp,c2);
    }
  } else {
    genomicpos = revoffset2R-bestcR;
    if (ngap + ngap + 3 > intronlength) {
      for (i = 0; i < intronlength; i++) {
	/* c2 = (bestcR+i < length2R) ? revsequence2R[-bestcR-i] : '?'; */
	c2 = revsequence2R[-bestcR-i];
	pairs = Pairpool_push(pairs,pairpool,offset1+bestr,genomicpos,' ',comp,c2);
      }
    } else {
      for (i = 0; i < ngap; i++) {
	/* c2 = (bestcR+i < length2R) ? revsequence2R[-bestcR-i] : '?'; */
	c2 = revsequence2R[-bestcR-i];
	pairs = Pairpool_push(pairs,pairpool,offset1+bestr,genomicpos,' ',comp,c2);
      }
    
      pairs = Pairpool_push(pairs,pairpool,offset1+bestr,genomicpos,' ',INTRONGAP_COMP,INTRONGAP_CHAR);
      pairs = Pairpool_push(pairs,pairpool,offset1+bestr,genomicpos,' ',INTRONGAP_COMP,INTRONGAP_CHAR);
      pairs = Pairpool_push(pairs,pairpool,offset1+bestr,genomicpos,' ',INTRONGAP_COMP,INTRONGAP_CHAR);
    
      for (i = ngap-1; i >= 0; i--) {
	/* c2 = (bestcL+i < length2L) ? sequence2L[bestcL+i] : '?'; */
	c2 = sequence2L[bestcL+i];
	pairs = Pairpool_push(pairs,pairpool,offset1+bestr,genomicpos,' ',comp,c2);
      }
    }
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
	   char *sequence1, char *sequenceuc1, char *sequence2, char *sequenceuc2,
	   int offset1, int offset2, Pairpool_T pairpool, 
	   bool revp, int ngap, int cdna_direction) {
  char c1, c2;
  int dist, j, genomicpos;
  char left1, left2, right2, right1, comp;
  int introntype, leftdi, rightdi;
  bool add_dashes_p;

  while (directions[r][c] != STOP) {
    if (directions[r][c] == DIAG) {
      if (revp == true) {
	c1 = sequence1[-r+1];
	c2 = sequence2[-c+1];

	if (sequenceuc1[-r+1] == sequenceuc2[-c+1]) {
	  debug(printf("D%d: Pushing %d,%d (%c,%c) - match\n",jump[r][c],r,c,c1,c2));
	  *nmatches += 1;
	  pairs = Pairpool_push(pairs,pairpool,offset1-(r-1),offset2-(c-1),c1,DYNPROG_MATCH_COMP,c2);

	} else if (consistent_array[(int) c1][(int) c2] == true) {
	  debug(printf("D%d: Pushing %d,%d (%c,%c) - ambiguous\n",jump[r][c],r,c,c1,c2));
	  *nmatches += 1;
	  pairs = Pairpool_push(pairs,pairpool,offset1-(r-1),offset2-(c-1),c1,AMBIGUOUS_COMP,c2);

	} else {
	  debug(printf("D%d: Pushing %d,%d (%c,%c) - mismatch\n",jump[r][c],r,c,c1,c2));
	  *nmismatches += 1;
	  pairs = Pairpool_push(pairs,pairpool,offset1-(r-1),offset2-(c-1),c1,MISMATCH_COMP,c2);
	}

      } else {
	c1 = sequence1[r-1];
	c2 = sequence2[c-1];

	if (sequenceuc1[r-1] == sequenceuc2[c-1]) {
	  debug(printf("D%d: Pushing %d,%d (%c,%c) - match\n",jump[r][c],r,c,c1,c2));
	  *nmatches += 1;
	  pairs = Pairpool_push(pairs,pairpool,offset1+r-1,offset2+c-1,c1,DYNPROG_MATCH_COMP,c2);

	} else if (consistent_array[(int) c1][(int) c2] == true) {
	  debug(printf("D%d: Pushing %d,%d (%c,%c) - ambiguous\n",jump[r][c],r,c,c1,c2));
	  *nmatches += 1;
	  pairs = Pairpool_push(pairs,pairpool,offset1+r-1,offset2+c-1,c1,AMBIGUOUS_COMP,c2);

	} else {
	  debug(printf("D%d: Pushing %d,%d (%c,%c) - mismatch\n",jump[r][c],r,c,c1,c2));
	  *nmismatches += 1;
	  pairs = Pairpool_push(pairs,pairpool,offset1+r-1,offset2+c-1,c1,MISMATCH_COMP,c2);
	}
      }
      r--; c--;

    } else if (directions[r][c] == HORIZ) {
      debug(printf("H%d: ",jump[r][c]));
      if (ngap == 0 || (dist = jump[r][c]) < MICROINTRON_LENGTH) {
	add_dashes_p = true;
      } else {
	/* Check for intron */
	if (revp == 1) {
	  left1 = sequenceuc2[-c+1];
	  left2 = sequenceuc2[-(c-1)+1];
	  right2 = sequenceuc2[-(c-(dist-2))+1];
	  right1 = sequenceuc2[-(c-(dist-1))+1];

	} else {
	  left1 = sequenceuc2[(c-(dist-1))-1];
	  left2 = sequenceuc2[(c-(dist-2))-1];
	  right2 = sequenceuc2[(c-1)-1];
	  right1 = sequenceuc2[c-1];
	}
	
	if (left1 == 'G' && left2 == 'T') {
	  leftdi = LEFT_GT;
	} else if (left1 == 'G' && left2 == 'C') {
	  leftdi = LEFT_GC;
	} else if (left1 == 'A' && left2 == 'T') {
	  leftdi = LEFT_AT;
#ifndef PMAP
	} else if (left1 == 'C' && left2 == 'T') {
	  leftdi = LEFT_CT;
#endif
	} else {
	  leftdi = 0x00;
	}

	if (right2 == 'A' && right1 == 'G') {
	  rightdi = RIGHT_AG;
	} else if (right2 == 'A' && right1 == 'C') {
	  rightdi = RIGHT_AC;
#ifndef PMAP
	} else if (right2 == 'G' && right1 == 'C') {
	  rightdi = RIGHT_GC;
	} else if (right2 == 'A' && right1 == 'T') {
	  rightdi = RIGHT_AT;
#endif
	} else {
	  rightdi = 0x00;
	}

	introntype = leftdi & rightdi;
#ifdef PMAP
	switch (introntype) {
	case GTAG_FWD: comp = FWD_CANONICAL_INTRON_COMP; break;
	case GCAG_FWD: comp = FWD_GCAG_INTRON_COMP; break;
	case ATAC_FWD: comp = FWD_ATAC_INTRON_COMP; break;
	default: comp = NONINTRON_COMP;
	}
#else
	if (cdna_direction > 0) {
	  switch (introntype) {
	  case GTAG_FWD: comp = FWD_CANONICAL_INTRON_COMP; break;
	  case GCAG_FWD: comp = FWD_GCAG_INTRON_COMP; break;
	  case ATAC_FWD: comp = FWD_ATAC_INTRON_COMP; break;
	  default: comp = NONINTRON_COMP;
	  }
	} else if (cdna_direction < 0) {
	  switch(introntype) {
	  case ATAC_REV: comp = REV_ATAC_INTRON_COMP; break;
	  case GCAG_REV: comp = REV_GCAG_INTRON_COMP; break;
	  case GTAG_REV: comp = REV_CANONICAL_INTRON_COMP; break;
	  default: comp = NONINTRON_COMP;
	  }
	} else {
	  comp = NONINTRON_COMP;
	}
#endif

	if (comp == NONINTRON_COMP) {
	  add_dashes_p = true;
	} else {
	  add_dashes_p = false;
	}
      }

      if (add_dashes_p == true) {
	/* Add dashes */
	*nopens += 1;
	*nindels += (dist = jump[r][c]);
	for (j = 0; j < dist; j++) {
	  c2 = revp ? sequence2[-c+1] : sequence2[c-1];
	  debug(printf("Pushing %d,%d (-,%c),",r,c,c2));
	  if (revp) {
	    pairs = Pairpool_push(pairs,pairpool,offset1-(r-1),offset2-(c-1),' ',INDEL_COMP,c2);
	  } else {
	    pairs = Pairpool_push(pairs,pairpool,offset1+r-1,offset2+c-1,' ',INDEL_COMP,c2);
	  }
	  c--;
	}
      } else {
	debug(printf("Large gap %c%c..%c%c.  Adding gap of type %c.\n",
		     left1,left2,right2,right1,comp));

	if (revp) {
	  genomicpos = offset2-(c-1);
	  for (j = 0; j < ngap; j++) {
	    c2 = sequence2[-(c-j)+1];
	    pairs = Pairpool_push(pairs,pairpool,offset1-(r-1),genomicpos,' ',comp,c2);
	  }
	  pairs = Pairpool_push(pairs,pairpool,offset1-(r-1),genomicpos,' ',INTRONGAP_COMP,INTRONGAP_CHAR);
	  pairs = Pairpool_push(pairs,pairpool,offset1-(r-1),genomicpos,' ',INTRONGAP_COMP,INTRONGAP_CHAR);
	  pairs = Pairpool_push(pairs,pairpool,offset1-(r-1),genomicpos,' ',INTRONGAP_COMP,INTRONGAP_CHAR);

	  c -= dist - 1;
	  for (j = ngap-1; j >= 0; j--) {
	    c2 = sequence2[-(c+j)+1];
	    pairs = Pairpool_push(pairs,pairpool,offset1-(r-1),genomicpos,' ',comp,c2);
	  }
	  c--;

	} else {
	  genomicpos = offset2+(c-1);
	  for (j = 0; j < ngap; j++) {
	    c2 = sequence2[c-j-1];
	    pairs = Pairpool_push(pairs,pairpool,offset1+(r-1),genomicpos,' ',comp,c2);
	  }
	  pairs = Pairpool_push(pairs,pairpool,offset1+(r-1),genomicpos,' ',INTRONGAP_COMP,INTRONGAP_CHAR);
	  pairs = Pairpool_push(pairs,pairpool,offset1+(r-1),genomicpos,' ',INTRONGAP_COMP,INTRONGAP_CHAR);
	  pairs = Pairpool_push(pairs,pairpool,offset1+(r-1),genomicpos,' ',INTRONGAP_COMP,INTRONGAP_CHAR);

	  c -= dist - 1;
	  for (j = ngap-1; j >= 0; j--) {
	    c2 = sequence2[c+j-1];
	    pairs = Pairpool_push(pairs,pairpool,offset1+(r-1),genomicpos,' ',comp,c2);
	  }
	  c--;

	}
      }

      debug(printf("\n"));

    } else if (directions[r][c] == VERT) {
      debug(printf("V%d: ",jump[r][c]));
      *nopens += 1;
      *nindels += (dist = jump[r][c]);
      for (j = 0; j < dist; j++) {
	c1 = revp ? sequence1[-r+1] : sequence1[r-1];
	debug(printf("Pushing %d,%d (%c,-), ",r,c,c1));
	if (revp) {
	  pairs = Pairpool_push(pairs,pairpool,offset1-(r-1),offset2-(c-1),c1,INDEL_COMP,' ');
	} else {
	  pairs = Pairpool_push(pairs,pairpool,offset1+r-1,offset2+c-1,c1,INDEL_COMP,' ');
	}
	r--;
      }
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
  int bestscore = -100000, scoreL, scoreR;
  int rL, rR, cL, cR;
  int lbandL, rbandL, cloL, chighL;
  int lbandR, rbandR, cloR, chighR;
  int pen, end_reward;

  /* Perform computations */
  rbandL = length1L - length2 + extraband_paired;
  lbandL = extraband_paired;

  rbandR = length1R - length2 + extraband_paired;
  lbandR = extraband_paired;

  /* Need a double loop on rows here, in contrast with a single loop
     for introns, because we allow a genomic insertion that doesn't
     match the cDNA.  So we need to add a penalty for a genomic
     insertion */
  for (rL = 0; rL <= length2; rL++) {
    if (rL == 0 || rL == length2) {
      end_reward = 1;
    } else {
      end_reward = 0;
    }

    /* Note: opening penalty is added at the bottom of the loop */
    for (rR = length2-rL, pen = 0; rR >= 0; rR--, pen += extend) {
      debug3(printf("\nAt row %d on left and %d on right\n",rL,rR));
      if ((cloL = rL - lbandL) < 0) {
	cloL = 0;
      }
      if ((chighL = rL + rbandL) > length1L) {
	chighL = length1L;
      }

      if ((cloR = rR - lbandR) < 0) {
	cloR = 0;
      }
      if ((chighR = rR + rbandR) > length1R) {
	chighR = length1R;
      }

      for (cL = cloL; cL <= chighL; cL++) {
	scoreL = matrixL[rL][cL];
	
	for (cR = cloR; cR <= chighR; cR++) {
	  scoreR = matrixR[rR][cR];

	  if (scoreL + scoreR + pen + end_reward > bestscore) {
	    debug3(printf("At %d left to %d right, score is (%d)+(%d)+(%d)+(%d) = %d (bestscore)\n",
			  cL,cR,scoreL,scoreR,pen,end_reward,scoreL+scoreR+pen+end_reward));
#if 0
	    if (leftoffset + cL >= rightoffset - cR) {
	      debug3(printf("  Disallowed because %d >= %d\n",leftoffset+cL,rightoffset-cR));
	    } else {
#endif
	      bestscore = scoreL + scoreR + pen + end_reward;
	      *bestrL = rL;
	      *bestrR = rR;
	      *bestcL = cL;
	      *bestcR = cR;
#if 0
	    }
#endif
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
intron_score (int *introntype, int leftdi, int rightdi, int cdna_direction, int canonical_reward, bool endp) {
  int scoreI;

#ifdef PMAP
  if ((*introntype = leftdi & rightdi) == NONINTRON) {
    scoreI = 0.0;
  } else {
    switch (*introntype) {
    case GTAG_FWD: scoreI = canonical_reward; break;
    case GCAG_FWD: case ATAC_FWD: 
      if (endp == true) {
	*introntype = NONINTRON; scoreI = 0.0;
      } else {
	scoreI = ALTERNATE_INTRON;
      }
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
      if (endp == true) {
	*introntype = NONINTRON; scoreI = 0.0;
      } else {
	scoreI = ALTERNATE_INTRON;
      }
      break;
    default: *introntype = NONINTRON; scoreI = 0.0;
    }
  } else if (cdna_direction < 0) {
    switch (*introntype) {
    case GTAG_REV: scoreI = canonical_reward; break;
    case GCAG_REV: case ATAC_REV:
      if (endp == true) {
	*introntype = NONINTRON; 
	scoreI = 0.0;
      } else {
	scoreI = ALTERNATE_INTRON;
      }
      break;
    default: *introntype = NONINTRON; scoreI = 0.0;
    }
  } else {
    switch (*introntype) {
    case GTAG_FWD: case GTAG_REV: scoreI = canonical_reward; break;
    case GCAG_FWD: case GCAG_REV: case ATAC_FWD: case ATAC_REV: 
      if (endp == true) {
	*introntype = NONINTRON;	      
	scoreI = 0.0;
      } else {
	scoreI = ALTERNATE_INTRON;
      }
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
		   int cdna_direction, int extraband_paired, bool endp, int canonical_reward,
		   int maxhorizjump, int maxvertjump, int leftoffset, int rightoffset) {
  int bestscore = -100000, scoreL, scoreR, scoreI;
  int rL, rR, cL, cR;
  int lbandL, rbandL, cloL, chighL;
  int lbandR, rbandR, cloR, chighR;
  char left1, left2, right2, right1;
  int *leftdi, *rightdi, introntype;

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

  for (rL = 0, rR = length1; rL <= length1; rL++, rR--) {
    debug3(printf("\nAt row %d on left and %d on right\n",rL,rR));
    if ((cloL = rL - lbandL) < 0) {
      cloL = 0;
    }
    if ((chighL = rL + rbandL) > length2L) {
      chighL = length2L;
    }

    if ((cloR = rR - lbandR) < 0) {
      cloR = 0;
    }
    if ((chighR = rR + rbandR) > length2R) {
      chighR = length2R;
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

	for (cR = cloR; cR <= chighR; cR++) {
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
				  cdna_direction,canonical_reward,endp);
	    
	    if (scoreL + scoreI + scoreR > bestscore) {
	      debug3(printf("At %d left to %d right, score is (%d)+(%d)+(%d) = %d (bestscore)\n",
			    cL,cR,scoreL,scoreI,scoreR,scoreL+scoreI+scoreR));
#if 0
	      if (leftoffset + cL >= rightoffset - cR) {
		debug3(printf("  Disallowed because %d >= %d\n",leftoffset+cL,rightoffset-cR));
	      } else {
#endif
		debug3(printf("  Allowed because %d < %d\n",leftoffset+cL,rightoffset-cR));
		bestscore = scoreL + scoreI + scoreR;
		*bestrL = rL;
		*bestrR = rR;
		*bestcL = cL;
		*bestcR = cR;
		*best_introntype = introntype;
#if 0
	      }
#endif
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
      
  *finalscore = bestscore;
  debug3(printf("Returning final score of %d at (%d,%d) left to (%d,%d) right\n",
		*finalscore,*bestrL,*bestcL,*bestrR,*bestcR));
  FREE(rightdi);
  FREE(leftdi);

  return;
}


/************************************************************************/

List_T
Dynprog_single_gap (int *finalscore, int *nmatches, int *nmismatches, int *nopens, int *nindels,
		    T dynprog, char *sequence1, char *sequenceuc1, char *sequence2, char *sequenceuc2,
		    int length1, int length2, int offset1, int offset2,
#ifdef PMAP
		    char *queryaaseq,
#endif
		    Pairpool_T pairpool, int extraband_single, double defect_rate) {
  List_T pairs = NULL;
#ifdef PMAP
  char *instsequence1;
#endif
  Mismatchtype_T mismatchtype;
  Jumptype_T jumptype;
  int **matrix, **jump, mismatch, open, extend;
  Direction_T **directions;
  bool allocp = false;

  /* Length1: maxlookback+MAXPEELBACK.  Length2 +EXTRAMATERIAL */

  debug(
	printf("Aligning single gap middle\n");
	printf("At query offset %d-%d, %.*s\n",offset1,offset1+length1-1,length1,sequence1);
	printf("At genomic offset %d-%d, %.*s\n",offset2,offset2+length2-1,length2,sequence2);
	printf("\n");
	);

  if (defect_rate < DEFECT_HIGHQ) {
    mismatchtype = HIGHQ;
    mismatch = MISMATCH_HIGHQ;
    jumptype = SINGLE_HIGHQ;
    open = SINGLE_OPEN_HIGHQ;
    extend = SINGLE_EXTEND_HIGHQ;
  } else if (defect_rate < DEFECT_MEDQ) {
    mismatchtype = MEDQ;
    mismatch = MISMATCH_MEDQ;
    jumptype = SINGLE_MEDQ;
    open = SINGLE_OPEN_MEDQ;
    extend = SINGLE_EXTEND_MEDQ;
  } else {
    mismatchtype = LOWQ;
    mismatch = MISMATCH_LOWQ;
    jumptype = SINGLE_LOWQ;
    open = SINGLE_OPEN_LOWQ;
    extend = SINGLE_EXTEND_LOWQ;
  }

  /* endp is false */
  if (length1 > dynprog->maxlength1 || length2 > dynprog->maxlength2) {
    dynprog = Dynprog_temp(length1,length2);
    allocp = true;
  }

#ifdef PMAP
  instsequence1 = instantiate_codons(sequence1,queryaaseq,offset1,length1);
#endif
  
  if (allocp == true) {
    matrix = compute_scores(&directions,&jump,dynprog,
#ifdef PMAP
			    instsequence1,sequence2,
#else
			    sequence1,sequence2,
#endif
			    /*revp*/false,length1,length2,
			    mismatchtype,jumptype,
			    mismatch,open,extend,extraband_single);
  } else {
    matrix = compute_scores_lookup(&directions,&jump,dynprog,
#ifdef PMAP
				   instsequence1,sequence2,
#else
				   sequence1,sequence2,
#endif
				   /*revp*/false,length1,length2,
				   mismatchtype,jumptype,extraband_single);
  }
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
		    /*revp*/false,/*ngap*/0,/*cdna_direction*/0);


#ifdef PMAP
  FREE(instsequence1);
#endif

  /*
  Directions_free(directions);
  Matrix_free(matrix);
  */

  if (allocp == true) {
    Dynprog_free(&dynprog);
  }

  return List_reverse(pairs);
}




/* Sequence 1L and 1R represent the two ends of the cDNA insertion */
List_T
Dynprog_cdna_gap (int *finalscore, T dynprogL, T dynprogR, 
		  char *sequence1L, char *sequenceuc1L, 
		  char *revsequence1R, char *revsequenceuc1R,
		  char *sequence2, char *sequenceuc2,
		  int length1L, int length1R, int length2,
		  int offset1L, int revoffset1R, int offset2,
#ifdef PMAP
		  char *queryaaseq,
#endif
		  Pairpool_T pairpool, int extraband_paired, double defect_rate) {
  List_T pairs = NULL, p;
  Pair_T pair;
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
  bool allocRp = false, allocLp = false;

  if (length2 <= 0) {
    return NULL;
  }

  debug(
	printf("\n");
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
    dynprogR = Dynprog_temp(length2,length1R);
    allocRp = true;
  }
  if (length2 > dynprogL->maxlength1 || length1L > dynprogL->maxlength2) {
    dynprogL = Dynprog_temp(length2,length1L);
    allocLp = true;
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

  if (allocRp == true) {
    matrixR = compute_scores(&directionsR,&jumpR,dynprogR,
#ifdef PMAP
			     revsequence2,instrevsequence1R,
#else
			     revsequence2,revsequence1R,
#endif
			     /*revp*/true,length2,length1R,mismatchtype,jumptype,
			     mismatch,open,extend,extraband_paired);
  } else {
    matrixR = compute_scores_lookup(&directionsR,&jumpR,dynprogR,
#ifdef PMAP
				    revsequence2,instrevsequence1R,
#else
				    revsequence2,revsequence1R,
#endif
				    /*revp*/true,length2,length1R,mismatchtype,jumptype,
				    extraband_paired);
  }

  if (allocLp == true) {
    matrixL = compute_scores(&directionsL,&jumpL,dynprogL,
#ifdef PMAP
			     sequence2,instsequence1L,
#else
			     sequence2,sequence1L,
#endif
			     /*revp*/false,length2,length1L,mismatchtype,jumptype,
			     mismatch,open,extend,extraband_paired);
  } else {
    matrixL = compute_scores_lookup(&directionsL,&jumpL,dynprogL,
#ifdef PMAP
				    sequence2,instsequence1L,
#else
				    sequence2,sequence1L,
#endif
				    /*revp*/false,length2,length1L,mismatchtype,jumptype,
				    extraband_paired);
  }

  bridge_cdna_gap(&(*finalscore),&bestrL,&bestrR,&bestcL,&bestcR,matrixL,matrixR,
		  directionsL,directionsR,length2,length1L,length1R,extraband_paired,
		  open,extend,offset1L,revoffset1R);

  nmatches = nmismatches = nopens = nindels = 0;
  pairs = traceback(NULL,&nmatches,&nmismatches,&nopens,&nindels,
		    directionsR,jumpR,bestrR,bestcR,
#ifdef PMAP
		    revsequence2,revsequenceuc2,instrevsequence1R,instrevsequence1R,
#else
		    revsequence2,revsequenceuc2,revsequence1R,revsequenceuc1R,
#endif
		    revoffset2,revoffset1R,pairpool,/*revp*/true,/*ngap*/0,/*cdna_direction*/0);
  pairs = List_reverse(pairs);

  /* cDNA and genome are flipped here.  Use c instead of r.  Unflip later. */
  /* Add genome insertion, if any */
  debug(printf("\n"));
  for (k = revoffset2-bestrR; k >= offset2+bestrL; k--) {
    debug(printf("genome insertion, Pushing %d (%c)\n",k+1,sequence2[k-offset2]));
    pairs = Pairpool_push(pairs,pairpool,k,revoffset1R-bestcR+1,sequence2[k-offset2],INDEL_COMP,' ');
  }
  debug(printf("\n"));

  /* Add cDNA insertion */
  for (k = revoffset1R-bestcR; k >= offset1L+bestcL; k--) {
#ifdef PMAP
    debug(printf("cDNA insertion, Pushing %d (%c)\n",k+1,instsequence1L[k-offset1L]));
    pairs = Pairpool_push(pairs,pairpool,offset2+bestrL,k,' ',INDEL_COMP,instsequence1L[k-offset1L]);
#else
    debug(printf("cDNA insertion, Pushing %d (%c)\n",k+1,sequence1L[k-offset1L]));
    pairs = Pairpool_push(pairs,pairpool,offset2+bestrL,k,' ',INDEL_COMP,sequence1L[k-offset1L]);
#endif
  }
  debug(printf("\n"));

  pairs = traceback(pairs,&nmatches,&nmismatches,&nopens,&nindels,
		    directionsL,jumpL,bestrL,bestcL,
#ifdef PMAP
		    sequence2,sequenceuc2,instsequence1L,instsequence1L,
#else
		    sequence2,sequenceuc2,sequence1L,sequenceuc1L,
#endif
		    offset2,offset1L,pairpool,/*revp*/false,/*ngap*/0,/*cdna_direction*/0);

  /* Unflip cDNA and genome */
  for (p = pairs; p != NULL; p = List_next(p)) {
    pair = (Pair_T) List_head(p);
    Pair_flip(pair);
  }

  /*
  Directions_free(directionsR);
  Directions_free(directionsL);
  Matrix_free(matrixR);
  Matrix_free(matrixL);
  */

#ifdef PMAP
  FREE(instsequence1);
#endif

  if (allocLp == true) {
    Dynprog_free(&dynprogL);
  }
  if (allocRp == true) {
    Dynprog_free(&dynprogR);
  }

  return List_reverse(pairs);
}


/* A genome gap is usually an intron.  Sequence 2L and 2R represent
   the two genomic ends of the intron. */
List_T
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
		    bool endp, double defect_rate, bool returnpairsp, bool addgapp) {
  List_T pairs = NULL;
#ifdef PMAP
  char *instsequence1, *instrevsequence1;
#else
  char *revsequence1, *revsequenceuc1;
#endif
  Mismatchtype_T mismatchtype;
  Jumptype_T jumptype;
  int **matrixL, **matrixR, **jumpL, **jumpR, mismatch, open, extend;
  int canonical_reward;
  Direction_T **directionsL, **directionsR;
  int revoffset1, bestrL, bestrR, bestcL, bestcR;
  int maxhorizjump, maxvertjump;
  bool allocRp = false, allocLp = false;

  debug(
	printf("\n");
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

  if (defect_rate < DEFECT_HIGHQ) {
    mismatchtype = HIGHQ;
    mismatch = MISMATCH_HIGHQ;
    jumptype = PAIRED_HIGHQ;
    open = PAIRED_OPEN_HIGHQ;
    extend = PAIRED_EXTEND_HIGHQ;
    canonical_reward = CANONICAL_INTRON_HIGHQ;
    maxhorizjump = MAXHORIZJUMP_HIGHQ;
    maxvertjump = MAXVERTJUMP_HIGHQ;
  } else if (defect_rate < DEFECT_MEDQ) {
    mismatchtype = MEDQ;
    mismatch = MISMATCH_MEDQ;
    jumptype = PAIRED_MEDQ;
    open = PAIRED_OPEN_MEDQ;
    extend = PAIRED_EXTEND_MEDQ;
    canonical_reward = CANONICAL_INTRON_MEDQ;
    maxhorizjump = MAXHORIZJUMP_MEDQ;
    maxvertjump = MAXVERTJUMP_MEDQ;
  } else {
    mismatchtype = LOWQ;
    mismatch = MISMATCH_LOWQ;
    jumptype = PAIRED_LOWQ;
    open = PAIRED_OPEN_LOWQ;
    extend = PAIRED_EXTEND_LOWQ;
    canonical_reward = CANONICAL_INTRON_LOWQ;
    maxhorizjump = MAXHORIZJUMP_LOWQ;
    maxvertjump = MAXVERTJUMP_LOWQ;
  }

  if (length1 > dynprogL->maxlength1 || length2L > dynprogL->maxlength2) {
    dynprogL = Dynprog_temp(length1,length2L);
    allocLp = true;
  }
  if (length1 > dynprogR->maxlength1 || length2R > dynprogR->maxlength2) {
    dynprogR = Dynprog_temp(length1,length2R);
    allocRp = true;
  }

#ifdef PMAP
  instsequence1 = instantiate_codons(sequence1,queryaaseq,offset1,length1);
  instrevsequence1 = &(instsequence1[length1-1]);
#else
  revsequence1 = &(sequence1[length1-1]);
  revsequenceuc1 = &(sequenceuc1[length1-1]);
#endif
  revoffset1 = offset1+length1-1;

  if (allocLp == true) {
    matrixL = compute_scores(&directionsL,&jumpL,dynprogL,
#ifdef PMAP
			     instsequence1,sequence2L,
#else
			     sequence1,sequence2L,
#endif
			     /*revp*/false,length1,length2L,mismatchtype,jumptype,
			     mismatch,open,extend,extraband_paired);
  } else {
    matrixL = compute_scores_lookup(&directionsL,&jumpL,dynprogL,
#ifdef PMAP
				    instsequence1,sequence2L,
#else
				    sequence1,sequence2L,
#endif
				    /*revp*/false,length1,length2L,mismatchtype,jumptype,
				    extraband_paired);
  }

  if (allocRp == true) {
    matrixR = compute_scores(&directionsR,&jumpR,dynprogR,
#ifdef PMAP
			     instrevsequence1,revsequence2R,
#else
			     revsequence1,revsequence2R,
#endif
			     /*revp*/true,length1,length2R,mismatchtype,jumptype,
			     mismatch,open,extend,extraband_paired);
  } else {
    matrixR = compute_scores_lookup(&directionsR,&jumpR,dynprogR,
#ifdef PMAP
				    instrevsequence1,revsequence2R,
#else
				    revsequence1,revsequence2R,
#endif
				    /*revp*/true,length1,length2R,mismatchtype,jumptype,
				    extraband_paired);
  }

  bridge_intron_gap(&(*finalscore),&bestrL,&bestrR,&bestcL,&bestcR,&(*introntype),matrixL,matrixR,
		    directionsL,directionsR,jumpL,jumpR,sequenceuc2L,revsequenceuc2R,
		    length1,length2L,length2R,cdna_direction,extraband_paired,endp,
		    canonical_reward,maxhorizjump,maxvertjump,offset2L,revoffset2R);


  *exonhead = revoffset1-(bestrR-1);

  *nmatches = *nmismatches = *nopens = *nindels = 0;
  if (returnpairsp == false) {
    countback(&(*nmatches),&(*nmismatches),&(*nopens),&(*nindels),
	      directionsR,jumpR,bestrR,bestcR,
#ifdef PMAP
	      instrevsequence1,revsequenceuc2R,
#else
	      revsequenceuc1,revsequenceuc2R,
#endif
	      revoffset1,revoffset2R,pairpool,/*revp*/true);
    countback(&(*nmatches),&(*nmismatches),&(*nopens),&(*nindels),
	      directionsL,jumpL,bestrL,bestcL,
#ifdef PMAP
	      instsequence1,sequenceuc2L,
#else
	      sequenceuc1,sequenceuc2L,
#endif
	      offset1,offset2L,pairpool,/*revp*/false);

  } else {
    pairs = traceback(NULL,&(*nmatches),&(*nmismatches),&(*nopens),&(*nindels),
		      directionsR,jumpR,bestrR,bestcR,
#ifdef PMAP
		      instrevsequence1,instrevsequence1,revsequence2R,revsequenceuc2R,
#else
		      revsequence1,revsequenceuc1,revsequence2R,revsequenceuc2R,
#endif
		      revoffset1,revoffset2R,pairpool,/*revp*/true,ngap,cdna_direction);
    pairs = List_reverse(pairs);
    if (addgapp) {
      pairs = add_gap(pairs,bestrL,bestcL,bestcR,
#ifdef PMAP
		      instsequence1,sequence2L,revsequence2R,
#else
		      sequence1,sequence2L,revsequence2R,
#endif
		      length1,length2L,length2R,offset1,offset2L,revoffset2R,
		      *introntype,ngap,pairpool);
    }
    pairs = traceback(pairs,&(*nmatches),&(*nmismatches),&(*nopens),&(*nindels),
		      directionsL,jumpL,bestrL,bestcL,
#ifdef PMAP
		      instsequence1,instsequence1,sequence2L,sequenceuc2L,
#else
		      sequence1,sequenceuc1,sequence2L,sequenceuc2L,
#endif
		      offset1,offset2L,pairpool,/*revp*/false,ngap,cdna_direction);
  }


  /*
  Directions_free(directionsR);
  Directions_free(directionsL);
  Matrix_free(matrixR);
  Matrix_free(matrixL);
  */

#ifdef PMAP
  FREE(instsequence1);
#endif

  if (allocRp == true) {
    Dynprog_free(&dynprogR);
  }
  if (allocLp == true) {
    Dynprog_free(&dynprogL);
  }

  return List_reverse(pairs);
}


List_T
Dynprog_end5_gap (int *finalscore, int *nmatches, int *nmismatches, 
		  int *nopens, int *nindels, T dynprog, 
		  char *revsequence1, char *revsequenceuc1,
		  char *revsequence2, char *revsequenceuc2,
		  int length1, int length2, int revoffset1, int revoffset2, 
#ifdef PMAP
		  char *queryaaseq,
#endif
		  Pairpool_T pairpool, int extraband_end, double defect_rate,
		  int cdna_direction, int ngap, bool extend_mismatch_p) {
  List_T pairs = NULL;
#ifdef PMAP
  char *instsequence1, *instrevsequence1;
#endif
  Pair_T pair;
  Mismatchtype_T mismatchtype;
  Jumptype_T jumptype;
  int **matrix, **jump, mismatch, open, extend;
  int bestr, bestc;
  Direction_T **directions;
  bool allocp = false;

  debug(
	printf("\n");
	printf("Aligning 5' end gap\n");
	printf("At query offset %d-%d, %.*s\n",revoffset1-length1+1,revoffset1,length1,&(revsequence1[-length1+1]));
	printf("At genomic offset %d-%d, %.*s\n",revoffset2-length2+1,revoffset2,length2,&(revsequence2[-length2+1]));
	printf("\n")
	);

  mismatchtype = END;
  mismatch = MISMATCH_END;
  if (defect_rate < DEFECT_HIGHQ) {
    jumptype = END_HIGHQ; 
    open = END_OPEN_HIGHQ;
    extend = END_EXTEND_HIGHQ;
  } else if (defect_rate < DEFECT_MEDQ) {
    jumptype = END_MEDQ;
    open = END_OPEN_MEDQ;
    extend = END_EXTEND_MEDQ;
  } else {
    jumptype = END_LOWQ;
    open = END_OPEN_LOWQ;
    extend = END_EXTEND_LOWQ;
  }

  if (length1 > dynprog->maxlength1 || length2 > dynprog->maxlength2) {
    dynprog = Dynprog_temp(length1,length2);
    allocp = true;
  }

#ifdef PMAP
  instsequence1 = instantiate_codons(&(revsequence1[-length1+1]),queryaaseq,revoffset1-length1+1,length1);
  instrevsequence1 = &(instsequence1[length1-1]);
#endif

  if (allocp == true) {
    matrix = compute_scores(&directions,&jump,dynprog,
#ifdef PMAP
			    instrevsequence1,revsequence2,
#else
			    revsequence1,revsequence2,
#endif
			    /*revp*/true,length1,length2,mismatchtype,jumptype,
			    mismatch,open,extend,extraband_end);
  } else {
    matrix = compute_scores_lookup(&directions,&jump,dynprog,
#ifdef PMAP
				   instrevsequence1,revsequence2,
#else
				   revsequence1,revsequence2,
#endif
				   /*revp*/true,length1,length2,mismatchtype,jumptype,
				   extraband_end);
  }
  find_best_endpoint_onegap(&(*finalscore),&bestr,&bestc,matrix,length1,length2,
			    extraband_end,extend_mismatch_p);


  *nmatches = *nmismatches = *nopens = *nindels = 0;
  pairs = traceback(NULL,&(*nmatches),&(*nmismatches),&(*nopens),&(*nindels),
		    directions,jump,bestr,bestc,
#ifdef PMAP
		    instrevsequence1,instrevsequence1,revsequence2,revsequenceuc2,
#else
		    revsequence1,revsequenceuc1,revsequence2,revsequenceuc2,
#endif
		    revoffset1,revoffset2,pairpool,/*revp*/true,ngap,cdna_direction);

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

  if (allocp == true) {
    Dynprog_free(&dynprog);
  }

  return List_reverse(pairs);
}


List_T
Dynprog_end3_gap (int *finalscore, int *nmatches, int *nmismatches, 
		  int *nopens, int *nindels, T dynprog, 
		  char *sequence1, char *sequenceuc1,
		  char *sequence2, char *sequenceuc2,
		  int length1, int length2, int offset1, int offset2, 
#ifdef PMAP
		  char *queryaaseq,
#endif
		  Pairpool_T pairpool, int extraband_end, double defect_rate,
		  int cdna_direction, int ngap, bool extend_mismatch_p) {
  List_T pairs = NULL;
#ifdef PMAP
  char *instsequence1;
#endif
  Pair_T pair;
  Mismatchtype_T mismatchtype;
  Jumptype_T jumptype;
  int **matrix, **jump, mismatch, open, extend;
  int bestr, bestc;
  Direction_T **directions;
  bool allocp = false;

  debug(
	printf("\n");
	printf("Aligning 3' end gap\n");
	printf("At query offset %d-%d, %.*s\n",offset1,offset1+length1-1,length1,sequence1);
	printf("At genomic offset %d-%d, %.*s\n",offset2,offset2+length2-1,length2,sequence2);
	printf("\n")
	);

  mismatchtype = END;
  mismatch = MISMATCH_END;
  if (defect_rate < DEFECT_HIGHQ) {
    jumptype = END_HIGHQ;
    open = END_OPEN_HIGHQ;
    extend = END_EXTEND_HIGHQ;
  } else if (defect_rate < DEFECT_MEDQ) {
    jumptype = END_MEDQ;
    open = END_OPEN_MEDQ;
    extend = END_EXTEND_MEDQ;
  } else {
    jumptype = END_LOWQ;
    open = END_OPEN_LOWQ;
    extend = END_EXTEND_LOWQ;
  }

  /* endp is true */
  if (length1 > dynprog->maxlength1 || length2 > dynprog->maxlength2) {
    dynprog = Dynprog_temp(length1,length2);
    allocp = true;
  }
  
#ifdef PMAP
  instsequence1 = instantiate_codons(sequence1,queryaaseq,offset1,length1);
#endif

  if (allocp == true) {
    matrix = compute_scores(&directions,&jump,dynprog,
#ifdef PMAP
			    instsequence1,sequence2,
#else
			    sequence1,sequence2,
#endif
			    /*revp*/false,length1,length2,mismatchtype,jumptype,
			    mismatch,open,extend,extraband_end);
  } else {
    matrix = compute_scores_lookup(&directions,&jump,dynprog,
#ifdef PMAP
				   instsequence1,sequence2,
#else
				   sequence1,sequence2,
#endif
				   /*revp*/false,length1,length2,mismatchtype,jumptype,
				   extraband_end);
  }
  find_best_endpoint_onegap(&(*finalscore),&bestr,&bestc,matrix,length1,length2,
			    extraband_end,extend_mismatch_p);


  *nmatches = *nmismatches = *nopens = *nindels = 0;
  pairs = traceback(NULL,&(*nmatches),&(*nmismatches),&(*nopens),&(*nindels),
		    directions,jump,bestr,bestc,
#ifdef PMAP
		    instsequence1,instsequence1,sequence2,sequenceuc2,
#else
		    sequence1,sequenceuc1,sequence2,sequenceuc2,
#endif
		    offset1,offset2,pairpool,/*revp*/false,ngap,cdna_direction);

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

  if (allocp == true) {
    Dynprog_free(&dynprog);
  }

  return pairs;			/* not List_reverse(pairs) */
}



static List_T
make_microexon_pairs_double (int offset1L, int offset1M, int offset1R,
			     int offset2L, int offset2M, int offset2R,
			     int lengthL, int lengthM, int lengthR,
#ifdef PMAP
			     char *queryaaseq,
#endif
			     char *queryseq, char *queryuc, char *genomicseg, char *genomicuc,
			     Pairpool_T pairpool, int ngap, char gapchar) {
  List_T pairs = NULL;
  char c1, c2;
  int i, j, k, genomicpos;

  /* Left segment */
  for (i = 0; i < lengthL; i++) {
    c1 = queryseq[offset1L+i];
    c2 = genomicseg[offset2L+i];
    if (queryuc[offset1L+i] == genomicuc[offset2L+i]) {
      pairs = Pairpool_push(pairs,pairpool,offset1L+i,offset2L+i,c1,DYNPROG_MATCH_COMP,c2);
    } else {
      pairs = Pairpool_push(pairs,pairpool,offset1L+i,offset2L+i,c1,MISMATCH_COMP,c2);
    }
  }

  /* First gap */
  j = i-1;
  genomicpos = offset2L+i;
  for (k = 0; k < ngap; k++) {
    c2 = genomicseg[offset2L+(++j)];
    pairs = Pairpool_push(pairs,pairpool,offset1L+i,genomicpos--,' ',gapchar,c2);
  }

  pairs = Pairpool_push(pairs,pairpool,offset1L+i,genomicpos,' ',INTRONGAP_COMP,INTRONGAP_CHAR);
  pairs = Pairpool_push(pairs,pairpool,offset1L+i,genomicpos,' ',INTRONGAP_COMP,INTRONGAP_CHAR);
  pairs = Pairpool_push(pairs,pairpool,offset1L+i,genomicpos,' ',INTRONGAP_COMP,INTRONGAP_CHAR);

  for (j = ngap; j >= 1; j--) {
    c2 = genomicseg[offset2M-j];
    pairs = Pairpool_push(pairs,pairpool,offset1L+i,genomicpos,' ',gapchar,c2);
  }
  
  /* Microexon */
  for (i = 0; i < lengthM; i++) {
    c1 = queryseq[offset1M+i];
    c2 = genomicseg[offset2M+i];
    if (queryuc[offset1M+i] == genomicuc[offset2M+i]) {
      pairs = Pairpool_push(pairs,pairpool,offset1M+i,offset2M+i,c1,DYNPROG_MATCH_COMP,c2);
    } else {
      pairs = Pairpool_push(pairs,pairpool,offset1M+i,offset2M+i,c1,MISMATCH_COMP,c2);
    }
  }

  /* Second gap */
  j = i-1;
  genomicpos = offset2M+i;
  for (k = 0; k < ngap; k++) {
    c2 = genomicseg[offset2M+(++j)];
    pairs = Pairpool_push(pairs,pairpool,offset1M+i,genomicpos--,' ',gapchar,c2);
  }

  pairs = Pairpool_push(pairs,pairpool,offset1M+i,genomicpos,' ',INTRONGAP_COMP,INTRONGAP_CHAR);
  pairs = Pairpool_push(pairs,pairpool,offset1M+i,genomicpos,' ',INTRONGAP_COMP,INTRONGAP_CHAR);
  pairs = Pairpool_push(pairs,pairpool,offset1M+i,genomicpos,' ',INTRONGAP_COMP,INTRONGAP_CHAR);

  for (j = ngap; j >= 1; j--) {
    c2 = genomicseg[offset2R-j];
    pairs = Pairpool_push(pairs,pairpool,offset1M+i,genomicpos,' ',gapchar,c2);
  }
  
  /* Right segment */
  for (i = 0; i < lengthR; i++) {
    c1 = queryseq[offset1R+i];
    c2 = genomicseg[offset2R+i];
    if (queryuc[offset1R+i] == genomicuc[offset2R+i]) {
      pairs = Pairpool_push(pairs,pairpool,offset1R+i,offset2R+i,c1,DYNPROG_MATCH_COMP,c2);
    } else {
      pairs = Pairpool_push(pairs,pairpool,offset1R+i,offset2R+i,c1,MISMATCH_COMP,c2);
    }
  }

  return pairs;
}


static List_T
make_microexon_pairs_single (int offset1L, int offset1R,
			     int offset2L, int offset2R,
			     int lengthL, int lengthR,
#ifdef PMAP
			     char *queryaaseq,
#endif
			     char *queryseq, char *queryuc, char *genomicseg, char *genomicuc,
			     Pairpool_T pairpool, int ngap, char gapchar) {
  List_T pairs = NULL;
  char c1, c2;
  int i, j, k, genomicpos;

  /* Microexon */
  for (i = 0; i < lengthL; i++) {
    c1 = queryseq[offset1L+i];
    c2 = genomicseg[offset2L+i];
    if (queryuc[offset1L+i] == genomicuc[offset2L+i]) {
      pairs = Pairpool_push(pairs,pairpool,offset1L+i,offset2L+i,c1,DYNPROG_MATCH_COMP,c2);
    } else {
      pairs = Pairpool_push(pairs,pairpool,offset1L+i,offset2L+i,c1,MISMATCH_COMP,c2);
    }
  }

  /* Gap */
  j = i-1;
  genomicpos = offset2L+i;
  for (k = 0; k < ngap; k++) {
    c2 = genomicseg[offset2L+(++j)];
    pairs = Pairpool_push(pairs,pairpool,offset1L+i,genomicpos--,' ',gapchar,c2);
  }

  pairs = Pairpool_push(pairs,pairpool,offset1L+i,genomicpos,' ',INTRONGAP_COMP,INTRONGAP_CHAR);
  pairs = Pairpool_push(pairs,pairpool,offset1L+i,genomicpos,' ',INTRONGAP_COMP,INTRONGAP_CHAR);
  pairs = Pairpool_push(pairs,pairpool,offset1L+i,genomicpos,' ',INTRONGAP_COMP,INTRONGAP_CHAR);

  for (j = ngap; j >= 1; j--) {
    c2 = genomicseg[offset2R-j];
    pairs = Pairpool_push(pairs,pairpool,offset1L+i,genomicpos,' ',gapchar,c2);
  }
  
  /* Right segment */
  for (i = 0; i < lengthR; i++) {
    c1 = queryseq[offset1R+i];
    c2 = genomicseg[offset2R+i];
    if (queryuc[offset1R+i] == genomicuc[offset2R+i]) {
      pairs = Pairpool_push(pairs,pairpool,offset1R+i,offset2R+i,c1,DYNPROG_MATCH_COMP,c2);
    } else {
      pairs = Pairpool_push(pairs,pairpool,offset1R+i,offset2R+i,c1,MISMATCH_COMP,c2);
    }
  }

  return pairs;
}

List_T
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
		       Pairpool_T pairpool, int ngap) {
  Intlist_T hits = NULL, p;
  int middlelength, cL, cR, mincR, maxcR, leftbound, rightbound, textleft, textright, candidate, i;
  int min_microexon_length, span, nmismatches;
  char left1, left2, right2, right1;
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

  debug(printf("Begin microexon search for %.*s and %.*s\n",
	       length2L,sequence2L,length2R,&(revsequence2R[-length2R+1])));
  debug(printf("  Query sequence is %.*s\n",length1,sequence1));
  span = revoffset2R-offset2L;
  debug(printf("  Genomic span is of length %d\n",span));

  if (span <= 0) {
    fprintf(stderr,"Bug in Dynprog_microexon_int.  span = %d.  Please report to twu@gene.com\n",span);
    abort();
  } else {
    min_microexon_length = ceil(-log(1.0-pow(1.0-MICROEXON_PVALUE,1.0/(double) span))/log(4));
  }
  min_microexon_length -= 8;	/* Two donor-acceptor pairs */
  debug(printf("  Min microexon length is %d\n",min_microexon_length));
  if (min_microexon_length > MAX_MICROEXON_LENGTH) {
    return NULL;
  } else if (min_microexon_length < MIN_MICROEXON_LENGTH) {
    min_microexon_length = MIN_MICROEXON_LENGTH;
  }

  leftbound = 0;
  nmismatches = 0;
  while (leftbound < length1 - 1 && nmismatches <= 1) {
    if (sequenceuc1[leftbound] != sequenceuc2L[leftbound]) {
      nmismatches++;
    }
    leftbound++;
  }
  leftbound--;			/* This is where the last mismatch occurred */

  rightbound = 0;
  i = length1-1;
  nmismatches = 0;
  while (i >= 0 && nmismatches <= 1) {
    if (sequenceuc1[i] != revsequenceuc2R[-rightbound]) {
      nmismatches++;
    }
    rightbound++;
    i--;
  }
  rightbound--;			/* This is where the last mismatch occurred */

  debug(printf("  Left must start before %d from left end.  Right must start after %d from right end\n",
	       leftbound,rightbound));

  for (cL = 0; cL <= leftbound; cL++) {
    left1 = sequenceuc2L[cL];
    left2 = sequenceuc2L[cL+1];
    debug(printf("  %d: %c%c\n",cL,left1,left2));
    if (left1 == intron1 && left2 == intron2) {
      mincR = length1 - MAX_MICROEXON_LENGTH - cL;
      if (mincR < 0) {
	mincR = 0;
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
	  debug(printf("  Found pair at %d to %d, length %d.  Middle sequence is %.*s\n",
		       cL,cR,middlelength,middlelength,&(sequence1[cL])));
	  
	  textleft = offset2L + cL + MICROINTRON_LENGTH;
	  textright = revoffset2R - cR - MICROINTRON_LENGTH;
	  hits = BoyerMoore(&(sequenceuc1[cL]),middlelength,&(genomicuc[textleft]),textright-textleft);
	  for (p = hits; p != NULL; p = Intlist_next(p)) {
	    candidate = textleft + Intlist_head(p);
	    if (genomicuc[candidate - 2] == intron3 &&
		genomicuc[candidate - 1] == intron4 &&
		genomicuc[candidate + middlelength] == intron1 &&
		genomicuc[candidate + middlelength + 1] == intron2) {
	      debug(printf("  Successful microexon at %d,%d,%d\n",offset2L,candidate,revoffset2R-cR));
	      
	      Intlist_free(&hits);
	      return make_microexon_pairs_double(offset1,offset1+cL,offset1+cL+middlelength,
						 offset2L,candidate,revoffset2R-cR+1,
						 cL,middlelength,cR,
#ifdef PMAP
						 queryaaseq,
#endif
						 queryseq,queryuc,genomicseg,genomicuc,
						 pairpool,ngap,gapchar);
	    }
	  }
	  Intlist_free(&hits);
	}
      }
    }
  }
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
Dynprog_microexon_5 (int *microintrontype, int *microexonlength,
		     char *revsequence1, char *revsequenceuc1, char *revsequence2, char *revsequenceuc2,
		     int length1, int length2, int revoffset1, int revoffset2, int cdna_direction,
#ifdef PMAP
		     char *queryaaseq,
#endif
		     char *queryseq, char *queryuc, char *genomicseg, char *genomicuc,
		     Pairpool_T pairpool, int ngap, bool end_microexons_p) {
  Intlist_T hits = NULL, p;
  int endlength, c, textleft, textright, candidate, nmismatches = 0;
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
  debug(printf("  Query sequence is %.*s\n",length1,&(revsequence1[-length1+1])));

  *microexonlength = 0;
  for (c = 0; c < length1 - MIN_MICROEXON_LENGTH; c++) {
    right2 = revsequenceuc2[-c-1];
    right1 = revsequenceuc2[-c];
    debug(printf("   Checking %c%c\n",right2,right1));
    if (c > 0 && revsequenceuc1[-c+1] != revsequenceuc2[-c+1]) {
      nmismatches++;
    }
    if (nmismatches > 1) {
      debug(printf("   Aborting at %c != %c\n",revsequence1[-c+1],revsequence2[-c+1]));
      *microintrontype = NONINTRON;
      return NULL;
    }
    if (right2 == intron3 && right1 == intron4) {
      endlength = length1 - c;
      debug(printf("  Found acceptor at %d, length %d.  End sequence is %.*s\n",
		       c,endlength,endlength,&(revsequence1[-endlength+1])));

      textright = revoffset2 - c - MICROINTRON_LENGTH;
      textleft = textright - search_length(endlength,textright,end_microexons_p) + MICROINTRON_LENGTH;
      hits = BoyerMoore(&(revsequenceuc1[-c-endlength+1]),endlength,&(genomicuc[textleft]),textright-textleft);
      for (p = hits; p != NULL; p = Intlist_next(p)) {
	candidate = textleft + Intlist_head(p);
	if (genomicseg[candidate + endlength] == intron1 &&
	    genomicseg[candidate + endlength + 1] == intron2) {
	  debug(printf("  Successful microexon at %d\n",candidate));

	  Intlist_free(&hits);
	  *microexonlength = endlength;
	  return make_microexon_pairs_single(revoffset1-c-endlength+1,revoffset1-c+1,
					     candidate,revoffset2-c+1,endlength,c,
#ifdef PMAP
					     queryaaseq,
#endif
					     queryseq,queryuc,genomicseg,genomicuc,
					     pairpool,ngap,gapchar);
	}
      }
      
      Intlist_free(&hits);
    }
  }
  *microintrontype = NONINTRON;
  return NULL;
}


List_T
Dynprog_microexon_3 (int *microintrontype, int *microexonlength, 
		     char *sequence1, char *sequenceuc1, char *sequence2, char *sequenceuc2,
		     int length1, int length2, int offset1, int offset2, int cdna_direction,
#ifdef PMAP
		     char *queryaaseq,
#endif
		     char *queryseq, char *queryuc, char *genomicseg, char *genomicuc,
		     int genomiclength, Pairpool_T pairpool, int ngap, bool end_microexons_p) {
  Intlist_T hits = NULL, p;
  int endlength, c, textleft, textright, candidate, nmismatches = 0;
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
  debug(printf("  Query sequence is %.*s\n",length1,sequence1));

  *microexonlength = 0;
  for (c = 0; c < length1 - MIN_MICROEXON_LENGTH; c++) {
    left1 = sequenceuc2[c];
    left2 = sequenceuc2[c+1];
    debug(printf("   Checking %c%c\n",left1,left2));
    if (c > 0 && sequenceuc1[c-1] != sequenceuc2[c-1]) {
      nmismatches++;
    }
    if (nmismatches > 1) {
      debug(printf("   Aborting at %c != %c\n",sequence1[c-1],sequence2[c-1]));
      *microintrontype = NONINTRON;
      return NULL;
    }
    if (left1 == intron1 && left2 == intron2) {
      endlength = length1 - c;
      debug(printf("  Found donor at %d, length %d.  End sequence is %.*s\n",
		   c,endlength,endlength,&(sequence1[c])));

      textleft = offset2 + c;
      textright = textleft + search_length(endlength,genomiclength-textleft,end_microexons_p);
      hits = BoyerMoore(&(sequenceuc1[c]),endlength,&(genomicuc[textleft]),textright-textleft);
      for (p = hits; p != NULL; p = Intlist_next(p)) {
	candidate = textleft + Intlist_head(p);
	if (genomicseg[candidate - 2] == intron3 &&
	    genomicseg[candidate - 1] == intron4) {
	  debug(printf("  Successful microexon at %d\n",candidate));

	  Intlist_free(&hits);
	  *microexonlength = endlength;
	  return make_microexon_pairs_single(offset1,offset1+c,
					     offset2,candidate,c,endlength,
#ifdef PMAP
					     queryaaseq,
#endif
					     queryseq,queryuc,genomicseg,genomicuc,
					     pairpool,ngap,gapchar);
	}
      }
      
      Intlist_free(&hits);
    }
  }
  *microintrontype = NONINTRON;
  return NULL;
}
