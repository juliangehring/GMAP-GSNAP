static char rcsid[] = "$Id: translation.c,v 1.36 2005/10/11 22:10:38 twu Exp $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "translation.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>		/* For toupper */
#include "mem.h"
#include "comp.h"
#include "pairdef.h"
#include "complement.h"
#include "list.h"
#include "mutation.h"
#ifdef PMAP
#include "dynprog.h"
#endif


#define IGNORE_MARGIN 20

#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/* Translation via reference */
#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif

/* translate_est_forward and translate_est_backward */
#ifdef DEBUG2
#define debug2(x) x
#else
#define debug2(x)
#endif

/* Mutations */
#ifdef DEBUG3
#define debug3(x) x
#else
#define debug3(x)
#endif

/* Nucleotide consistency */
#ifdef DEBUG4
#define debug4(x) x
#else
#define debug4(x)
#endif


typedef enum {FRAME0, FRAME1, FRAME2, NOFRAME} Frame_T;

#define T Translation_T
struct T {
  char aa;
  Frame_T frame;
};


static struct T *
Translation_array_new (int translationlen) {
  struct T *new;
  int i;

  new = (struct T *) CALLOC(translationlen,sizeof(struct T));
  for (i = 0; i < translationlen; i++) {
    new[i].aa = ' ';
    new[i].frame = NOFRAME;
  }

  return new;
}

static void
Translation_dump (struct Pair_T *pairs, struct T *translation, int translationlen) {
  int i;

  for (i = 0; i < translationlen; i++) {
    if (pairs[i].aamarker == true) {
      printf("=> %c %c ",pairs[i].aa_g,pairs[i].aa_e);
    } else {
      printf("       ");
    }
    printf("%d: %d %d ",i,pairs[i].querypos,pairs[i].aapos);
    switch (translation[i].frame) {
    case NOFRAME: printf("%c %c %c",' ',' ',' '); break;
    case FRAME0: printf("%c %c %c",translation[i].aa,' ',' '); break;
    case FRAME1: printf("%c %c %c",' ',translation[i].aa,' '); break;
    case FRAME2: printf("%c %c %c",' ',' ',translation[i].aa); break;
    }
    printf("\n");
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
    -1, 17, 18, -1, 19, -1,

    -1, -1, -1, -1, -1, -1,

  /* a,  b,  c,  d,  e,  f,  g,  h,  i,  j, */
     0, -1,  1,  2,  3,  4,  5,  6,  7, -1,

  /* k,  l,  m,  n,  o,  p,  q,  r,  s,  t  */
     8,  9, 10, 11, -1, 12, 13, 14, 15, 16,

  /* u,  v,  w,  x,  y,  z  */
    -1, 17, 18, -1, 19, -1,

    -1, -1, -1, -1, -1};


static char *full_iupac_table[21][2] =
  {{"GCN",""},			/* A */
   {"TGY",""},			/* C */
   {"GAY",""},			/* D */
   {"GAR",""},			/* E */
   {"TTY",""},			/* F */
   {"GGN",""},			/* G */
   {"CAY",""},			/* H */
   {"ATH",""},			/* I */
   {"AAR",""},			/* K */
   {"TTR","CTN"},		/* L */
   {"ATG",""},			/* M */
   {"AAY",""},			/* N */
   {"CCN",""},			/* P */
   {"CAR",""},			/* Q */
   {"CGN","AGR"},		/* R */
   {"TCN","AGY"},		/* S */
   {"ACN",""},			/* T */
   {"GTN",""},			/* V */
   {"TGG",""},			/* W */
   {"TAY",""},			/* Y */
   {"TAR","TGA"}};		/* STOP */

static int n_full_iupac_table[21] = {1,1,1,1,1,1,1,1,1,2,1,1,1,1,2,2,1,1,1,1,2};
#endif

/************************************************************************/


static char
get_codon (char a, char b, char c) {
  switch (b) {
  case 'T':
    switch (a) {
    case 'T':
      switch (c) {
      case 'T': case 'C': case 'Y': return 'F';
      case 'A': case 'G': case 'R': return 'L';
      default: 	return 'X';
      }
    case 'C': return 'L';
    case 'A': 
      switch (c) {
      case 'G': return 'M';
      case 'T': case 'A': case 'C': case 'H': return 'I';
      default: return 'X';
      }
    case 'G': return 'V';
    default: return 'X';
  }
  case 'C':
    switch (a) {
    case 'T': return 'S';
    case 'C': return 'P';
    case 'A': return 'T';
    case 'G': return 'A';
    default: return 'X';
    }
  case 'A':
    switch (a) {
    case 'T':
      switch (c) {
      case 'T': case 'C': case 'Y': return 'Y';
      case 'A': case 'G': case 'R': return '*';
      default: return 'X';
      }
    case 'C':
      switch (c) {
      case 'T': case 'C': case 'Y': return 'H';
      case 'A': case 'G': case 'R': return 'Q';
      default: return 'X';
      }
    case 'A':
      switch (c) {
      case 'T': case 'C': case 'Y': return 'N';
      case 'A': case 'G': case 'R': return 'K';
      default: return 'X';
      }
    case 'G':
      switch (c) {
      case 'T': case 'C': case 'Y': return 'D';
      case 'A': case 'G': case 'R': return 'E';
      default: return 'X';
      }
    default: return 'X';
    }
  case 'G':
    switch (a) {
    case 'T':
      switch (c) {
      case 'T': case 'C': case 'Y': return 'C';
      case 'A': return '*';
      case 'G': return 'W';
      default: return 'X';
      }
    case 'C': return 'R';
    case 'A':
      switch (c) {
      case 'T': case 'C': case 'Y': return 'S';
      case 'A': case 'G': case 'R': return 'R';
      default: return 'X';
      }
    case 'G': return 'G';
    default: return 'X';
    }
  default: return 'X';
  }
  return 'X';
}


static void
find_bounds_forward (int *translation_frame, int *translation_start, 
		     int *translation_end, int *translation_length,
		     bool *endstopp, struct T *translation, 
		     int translationlen, bool fulllengthp) {
  int beststart0, beststart1, beststart2, bestend0, bestend1, bestend2;
  int bestorf0 = 0, bestorf1 = 0, bestorf2 = 0, orf0 = 0, orf1 = 0, orf2 = 0;
  int start0 = 0, start1 = 0, start2 = 0;
  bool needmet0p, needmet1p, needmet2p;
  bool endstop0p = false, endstop1p = false, endstop2p = false;
  char codon;
  int i, frame;

  if (fulllengthp == true) {
    needmet0p = needmet1p = needmet2p = true;
  } else {
    needmet0p = needmet1p = needmet2p = false;
  }

  for (i = 0; i < translationlen; i++) {
    debug(printf("%d %c: %d %d %d\n",i,translation[i].aa,orf0,orf1,orf2));
    frame = translation[i].frame;
    if ((codon = translation[i].aa) != ' ') {
      if (frame == FRAME0) {
	if (needmet0p) {
	  if (codon == 'M') {
	    orf0 = 1;
	    start0 = i;
	    needmet0p = false;
	  }
	} else if (codon == '*') {
	  orf0++;
	  if (orf0 > bestorf0) {
	    debug(printf("Frame 0: Best orf is %d\n",orf0));
	    bestorf0 = orf0;
	    beststart0 = start0;
	    bestend0 = i;
	    endstop0p = true;
	  }
	  needmet0p = true;
	} else {
	  debug(printf("Incrementing orf0\n"));
	  orf0++;
	}
      } else if (frame == FRAME1) {
	if (needmet1p) {
	  if (codon == 'M') {
	    orf1 = 1;
	    start1 = i;
	    needmet1p = false;
	  }
	} else if (codon == '*') {
	  orf1++;
	  if (orf1 > bestorf1) {
	    debug(printf("Frame 1: Best orf is %d\n",orf1));
	    bestorf1 = orf1;
	    beststart1 = start1;
	    bestend1 = i;
	    endstop1p = true;
	  }
	  needmet1p = true;
	} else {
	  debug(printf("Incrementing orf1\n"));
	  orf1++;
	}
      } else if (frame == FRAME2) {
	if (needmet2p) {
	  if (codon == 'M') {
	    orf2 = 1;
	    start2 = i;
	    needmet2p = false;
	  }
	} else if (codon == '*') {
	  orf2++;
	  if (orf2 > bestorf2) {
	    debug(printf("Frame 2: Best orf is %d\n",orf2));
	    bestorf2 = orf2;
	    beststart2 = start2;
	    bestend2 = i;
	    endstop2p = true;
	  }
	  needmet2p = true;
	} else {
	  debug(printf("Incrementing orf2\n"));
	  orf2++;
	}
      } else {
	fprintf(stderr,"No frame at %d\n",i);
      }
    }
  }
  
  /* Handle last segments */
  if (orf0 > bestorf0) {
    debug(printf("Frame 0: Best orf is %d\n",orf0));
    bestorf0 = orf0;
    beststart0 = start0;
    bestend0 = translationlen-1;
    endstop0p = false;
  }
  if (orf1 > bestorf1) {
    debug(printf("Frame 1: Best orf is %d\n",orf1));
    bestorf1 = orf1;
    beststart1 = start1;
    bestend1 = translationlen-1;
    endstop1p = false;
  }
  if (orf2 > bestorf2) {
    debug(printf("Frame 2: Best orf is %d\n",orf2));
    bestorf2 = orf2;
    beststart2 = start2;
    bestend2 = translationlen-1;
    endstop2p = false;
  }

  /* Find overall best */
  *translation_length = bestorf0;
  *endstopp = endstop0p;
  if (bestorf1 > *translation_length) {
    *translation_length = bestorf1;
    *endstopp = endstop1p;
  }
  if (bestorf2 > *translation_length) {
    *translation_length = bestorf2;
    *endstopp = endstop2p;
  }

  if (bestorf2 == *translation_length) {
    debug(printf("Assigning frame 2\n"));
    *translation_frame = FRAME2;
    *translation_start = beststart2;
    *translation_end = bestend2;
  } else if (bestorf1 == *translation_length) {
    debug(printf("Assigning frame 1\n"));
    *translation_frame = FRAME1;
    *translation_start = beststart1;
    *translation_end = bestend1;
  } else if (bestorf0 == *translation_length) {
    debug(printf("Assigning frame 0\n"));
    *translation_frame = FRAME0;
    *translation_start = beststart0;
    *translation_end = bestend0;
  } else {
    abort();
  }

  debug(printf("Frame0: %d, Frame1: %d, Frame2: %d\n",bestorf0,bestorf1,bestorf2));
  debug(printf("Best orf is %d codons\n",*translation_length));
  debug(printf("Frame0: %d %d\n",beststart0,bestend0));
  debug(printf("Frame1: %d %d\n",beststart1,bestend1));
  debug(printf("Frame2: %d %d\n",beststart2,bestend2));
  debug(printf("Value of endstopp is %d\n",*endstopp));

  return;
}

static void
find_bounds_backward (int *translation_frame, int *translation_start,
		      int *translation_end, int *translation_length,
		      bool *endstopp, struct T *translation, int translationlen, bool fulllengthp) {
  int beststart0, beststart1, beststart2, bestend0, bestend1, bestend2;
  int bestorf0 = 0, bestorf1 = 0, bestorf2 = 0, orf0 = 0, orf1 = 0, orf2 = 0;
  int start0 = translationlen-1, start1 = translationlen-1, start2 = translationlen-1;
  bool needmet0p, needmet1p, needmet2p;
  bool endstop0p = false, endstop1p = false, endstop2p = false;
  char codon;
  int i, frame;

  if (fulllengthp == true) {
    needmet0p = needmet1p = needmet2p = true;
  } else {
    needmet0p = needmet1p = needmet2p = false;
  }

  for (i = translationlen-1; i >= 0; --i) {
    frame = translation[i].frame;
    if ((codon = translation[i].aa) != ' ') {
      if (frame == FRAME0) {
	if (needmet0p) {
	  if (codon == 'M') {
	    orf0 = 1;
	    start0 = i;
	    needmet0p = false;
	  }
	} else if (codon == '*') {
	  orf0++;
	  if (orf0 > bestorf0) {
	    debug(printf("Frame 0: Best orf is %d\n",orf0));
	    bestorf0 = orf0;
	    bestend0 = i;
	    beststart0 = start0;
	    endstop0p = true;
	  }
	  needmet0p = true;
	} else {
	  orf0++;
	}
      } else if (frame == FRAME1) {
	if (needmet1p) {
	  if (codon == 'M') {
	    orf1 = 1;
	    start1 = i;
	    needmet1p = false;
	  }
	} else if (codon == '*') {
	  orf1++;
	  if (orf1 > bestorf1) {
	    debug(printf("Frame 1: Best orf is %d\n",orf1));
	    bestorf1 = orf1;
	    bestend1 = i;
	    beststart1 = start1;
	    endstop1p = true;
	  }
	  needmet1p = true;
	} else {
	  orf1++;
	}
      } else if (frame == FRAME2) {
	if (needmet2p) {
	  if (codon == 'M') {
	    orf2 = 1;
	    start2 = i;
	    needmet2p = false;
	  }
	} else if (codon == '*') {
	  orf2++;
	  if (orf2 > bestorf2) {
	    debug(printf("Frame 2: Best orf is %d\n",orf2));
	    bestorf2 = orf2;
	    bestend2 = i;
	    beststart2 = start2;
	    endstop2p = true;
	  }
	  needmet2p = true;
	} else {
	  orf2++;
	}
      }
    }
  }
  
  /* Handle last segments */
  if (orf0 > bestorf0) {
    debug(printf("Frame 0: Best orf is %d\n",orf0));
    bestorf0 = orf0;
    bestend0 = 0;
    beststart0 = start0;
    endstop0p = false;
  }
  if (orf1 > bestorf1) {
    debug(printf("Frame 1: Best orf is %d\n",orf1));
    bestorf1 = orf1;
    bestend1 = 0;
    beststart1 = start1;
    endstop1p = false;
  }
  if (orf2 > bestorf2) {
    debug(printf("Frame 2: Best orf is %d\n",orf2));
    bestorf2 = orf2;
    bestend2 = 0;
    beststart2 = start2;
    endstop2p = false;
  }

  /* Find overall best */
  *translation_length = bestorf0;
  *endstopp = endstop0p;
  if (bestorf1 > *translation_length) {
    *translation_length = bestorf1;
    *endstopp = endstop1p;
  }
  if (bestorf2 > *translation_length) {
    *translation_length = bestorf2;
    *endstopp = endstop2p;
  }

  if (bestorf0 == *translation_length) {
    debug(printf("Assigning frame 0\n"));
    *translation_frame = FRAME0;
    *translation_start = beststart0;
    *translation_end = bestend0;
  } else if (bestorf1 == *translation_length) {
    debug(printf("Assigning frame 1\n"));
    *translation_frame = FRAME1;
    *translation_start = beststart1;
    *translation_end = bestend1;
  } else if (bestorf2 == *translation_length) {
    debug(printf("Assigning frame 2\n"));
    *translation_frame = FRAME2;
    *translation_start = beststart2;
    *translation_end = bestend2;
  } else {
    abort();
  }

  debug(printf("Frame0: %d, Frame1: %d, Frame2: %d\n",bestorf0,bestorf1,bestorf2));
  debug(printf("Best orf is %d codons\n",*translation_length));
  debug(printf("Frame0: %d %d\n",beststart0,bestend0));
  debug(printf("Frame1: %d %d\n",beststart1,bestend1));
  debug(printf("Frame2: %d %d\n",beststart2,bestend2));
  debug(printf("Value of endstopp is %d\n",*endstopp));

  return;
}


static void
translate_pairs_cdna (int *translation_start, int *translation_length, struct Pair_T *pairs, int npairs,
		      char *queryaaseq_ptr) {
  struct Pair_T *ptr, *pair;
  int i;

  *translation_length = 0;

  /* Go backward so we put aa at beginning of codon */
  ptr = &(pairs[npairs]);
  for (i = npairs-1; i >= 0; --i) {
    pair = --ptr;
    pair->aapos = pair->querypos/3 + 1;
    if (pair->gapp == true) {
      /* pair->aa_e = ' '; */
    } else if (pair->cdna == ' ') {
      /* pair->aa_e = ' '; */
    } else {
      if (pair->querypos % 3 == 0) {
	pair->aa_e = queryaaseq_ptr[pair->querypos/3];
	pair->aamarker = true;
	*translation_start = i;
	(*translation_length) += 1;
      }
    }
  }

  return;
}  


static struct Translation_T *
translate_pairs_forward (struct Pair_T *pairs, int npairs, bool revcompp) {
  struct T *translation;
  struct Pair_T *ptr, *pair;
  int i, gpos = 0;
  char codon, nt2 = 'X', nt1 = 'X', nt0 = 'X';
  char complCode[128] = COMPLEMENT;

  translation = Translation_array_new(npairs);

  /* Go backward so we put aa at beginning of codon */
  ptr = &(pairs[npairs]);
  for (i = npairs-1; i >= 0; --i) {
    pair = --ptr;
    if (pair->gapp == true) {
      /* translation[i].aa = ' '; */
    } else if (pair->genome == ' ') {
      /* translation[i].aa = ' '; */
    } else {
      nt2 = nt1;
      nt1 = nt0;
      nt0 = revcompp ? complCode[toupper((int) pair->genome)] : toupper((int) pair->genome);

      codon = get_codon(nt0,nt1,nt2);
      if (gpos < 2 && codon == 'X') {
	/* translation[i].aa = ' '; */
      } else {
	translation[i].aa = codon;
	switch (gpos % 3) {
	case 0: translation[i].frame = FRAME0; break;
	case 1: translation[i].frame = FRAME1; break;
	case 2: translation[i].frame = FRAME2; break;
	}
      }
      gpos++;
    }
  }

  return translation;
}  

static struct Translation_T *
translate_pairs_backward (struct Pair_T *pairs, int npairs, bool revcompp) {
  struct T *translation;
  struct Pair_T *ptr, *pair;
  int i, gpos = 0;
  char codon, nt2 = 'X', nt1 = 'X', nt0 = 'X';
  char complCode[128] = COMPLEMENT;

  translation = Translation_array_new(npairs);

  /* Go forwards so we put aa at beginning of codon */
  ptr = pairs;
  for (i = 0; i < npairs; i++) {
    pair = ptr++;
    if (pair->gapp == true) {
      /* translation[i].aa = ' '; */
    } else if (pair->genome == ' ') {
      /* translation[i].aa = ' '; */
    } else {
      nt2 = nt1;
      nt1 = nt0;
      nt0 = revcompp ? complCode[toupper((int) pair->genome)] : toupper((int) pair->genome);

      codon = get_codon(nt0,nt1,nt2);
      if (gpos < 2 && codon == 'X') {
	/* translation[i].aa = ' '; */
      } else {
	translation[i].aa = codon;
	switch (gpos % 3) {
	case 0: translation[i].frame = FRAME0; break;
	case 1: translation[i].frame = FRAME1; break;
	case 2: translation[i].frame = FRAME2; break;
	}
      }
      gpos++;
    }
  }

  return translation;
}  


/* The only difference between standard and alt is the precedence of
   aamarker and cdna == ' ' */
static int
count_aa_genomic_standard (int aapos, struct Pair_T *pairs, 
			   int npairs, int starti, int end) {
  int number_aa = 0, j;

  j = starti;
  if (pairs[j].gapp == true) {
  } else if (pairs[j].genome == ' ') {
  } else {
    number_aa++;
  }
  j++;

  while (j <= end) {
    if (pairs[j].gapp == true) {
    } else if (pairs[j].aamarker == true) {
      return number_aa;
    } else if (pairs[j].genome == ' ') {
    } else {
      number_aa++;
    }
    j++;
  }

  while (j < npairs && number_aa < 3) {
    if (pairs[j].gapp == true) {
    } else if (pairs[j].aamarker == true) {
      return number_aa;
    } else if (pairs[j].genome == ' ') {
    } else {
      number_aa++;
    }
    j++;
  }

  return number_aa;
}

static int
count_aa_genomic_alt (int *j, int aapos, struct Pair_T *pairs, 
		      int npairs, int starti, int end) {
  int number_aa = 0;

  *j = starti;
  if (pairs[*j].gapp == true) {
  } else if (pairs[*j].genome == ' ') {
  } else {
    number_aa++;
  }
  (*j)++;

  while (*j <= end) {
    if (pairs[*j].gapp == true) {
    } else if (pairs[*j].genome == ' ') {
    } else if (pairs[*j].aamarker == true) {
      return number_aa;
    } else {
      number_aa++;
    }
    (*j)++;
  }

  while (*j < npairs && number_aa < 3) {
    if (pairs[*j].gapp == true) {
    } else if (pairs[*j].genome == ' ') {
    } else if (pairs[*j].aamarker == true) {
      return number_aa;
    } else {
      number_aa++;
    }
    (*j)++;
  }

  return number_aa;
}


#ifdef PMAP
static int
get_codon_genomic (int *j, int aapos, char aa_e, 
		   struct Pair_T *pairs, int npairs, int starti, int end) {
  char nt_g2 = 'X', nt_g1 = 'X', nt_g0 = 'X';
  char *codon, *bestcodon;
  int ncorrespondences = 0, number_aa = 0, index, ncodons, choice, bestscore, score, i;

  *j = starti;
  while (number_aa < 3) {
    if (pairs[*j].gapp == true) {
    } else if (pairs[*j].genome == ' ') {
    } else {
      nt_g0 = nt_g1;
      nt_g1 = nt_g2;
      nt_g2 = toupper(pairs[*j].genome);

      if (pairs[*j].cdna != ' ') {
	ncorrespondences++;
      }

      number_aa++;
    }
    (*j)++;
  }

  /* This avoids putting aa in an intron */
  while (*j < npairs-1 && pairs[*j].gapp == true) {
    (*j)++;
  }

  /* Make cDNA nucleotides consistent */
  if (ncorrespondences == 3) {
    index = aa_index_table[aa_e];
    if ((ncodons = n_full_iupac_table[index]) > 1) {
      debug4(printf("Making amino acid %c at %d consistent with %c%c%c\n",
		    aa_e,starti,nt_g0,nt_g1,nt_g2));
      bestscore = -1000000;
      for (choice = 0; choice < ncodons; choice++) {
	codon = full_iupac_table[index][choice];
	score = Dynprog_pairdistance(codon[0],nt_g0) + Dynprog_pairdistance(codon[1],nt_g1) +
	  Dynprog_pairdistance(codon[2],nt_g2);
	debug4(printf("  Score of codon %s is %d+%d+%d=%d\n",
		      codon,Dynprog_pairdistance(codon[0],nt_g0),Dynprog_pairdistance(codon[1],nt_g1),
		      Dynprog_pairdistance(codon[2],nt_g2),score));
	if (score > bestscore) {
	  bestcodon = codon;
	  bestscore = score;
	}
      }

      i = starti;
      number_aa = 0;
      while (number_aa < 3) {
	if (pairs[i].gapp == true) {
	} else if (pairs[i].cdna == ' ') {
	} else {
	  if (pairs[i].comp == AMBIGUOUS_COMP || pairs[i].comp == MISMATCH_COMP) {
	    debug4(printf("  Position %d: %c => %c\n",number_aa,pairs[i].cdna,bestcodon[number_aa]));
	    pairs[i].cdna = bestcodon[number_aa];
	    if (number_aa == 0) {
	      debug4(printf("    Comparing %c with %c\n",pairs[i].cdna,nt_g0));
	      if (Dynprog_pairdistance(pairs[i].cdna,nt_g0) < 0) {
		debug4(printf("    Changed comp to space\n"));
		pairs[i].comp = MISMATCH_COMP;
	      }
	    } else if (number_aa == 1) {
	      debug4(printf("    Comparing %c with %c\n",pairs[i].cdna,nt_g1));
	      if (Dynprog_pairdistance(pairs[i].cdna,nt_g1) < 0) {
		debug4(printf("    Changed comp to space\n"));
		pairs[i].comp = MISMATCH_COMP;
	      }
	    } else {
	      debug4(printf("    Comparing %c with %c\n",pairs[i].cdna,nt_g2));
	      if (Dynprog_pairdistance(pairs[i].cdna,nt_g2) < 0) {
		debug4(printf("    Changed comp to space\n"));
		pairs[i].comp = MISMATCH_COMP;
	      }
	    }
	  }
	  number_aa++;
	}
	i++;
      }
    }
  }

  return get_codon(nt_g0,nt_g1,nt_g2);
}
#endif


/* The only difference between standard and alt is the precedence of
   aamarker and cdna == ' ' */
static int
count_aa_forward_standard (int aapos, struct Pair_T *pairs, 
			   int npairs, int starti, int end) {
  int number_aa = 0, j;

  j = starti;
  if (pairs[j].gapp == true) {
  } else if (pairs[j].cdna == ' ') {
  } else {
    number_aa++;
  }
  j++;

  while (j <= end) {
    if (pairs[j].gapp == true) {
    } else if (pairs[j].aamarker == true) {
      return number_aa;
    } else if (pairs[j].cdna == ' ') {
    } else {
      number_aa++;
    }
    j++;
  }

  while (j < npairs && number_aa < 3) {
    if (pairs[j].gapp == true) {
    } else if (pairs[j].aamarker == true) {
      return number_aa;
    } else if (pairs[j].cdna == ' ') {
    } else {
      number_aa++;
    }
    j++;
  }

  return number_aa;
}

static int
count_aa_forward_alt (int *j, int aapos, struct Pair_T *pairs, 
		      int npairs, int starti, int end) {
  int number_aa = 0;

  *j = starti;
  if (pairs[*j].gapp == true) {
  } else if (pairs[*j].cdna == ' ') {
  } else {
    number_aa++;
  }
  (*j)++;

  while (*j <= end) {
    if (pairs[*j].gapp == true) {
    } else if (pairs[*j].cdna == ' ') {
    } else if (pairs[*j].aamarker == true) {
      return number_aa;
    } else {
      number_aa++;
    }
    (*j)++;
  }

  while (*j < npairs && number_aa < 3) {
    if (pairs[*j].gapp == true) {
    } else if (pairs[*j].cdna == ' ') {
    } else if (pairs[*j].aamarker == true) {
      return number_aa;
    } else {
      number_aa++;
    }
    (*j)++;
  }

  return number_aa;
}

static int
get_codon_forward (int *j, int aapos, struct Pair_T *pairs, int npairs, int starti, int end, 
		   char *complCode, bool revcompp) {
  char nt2 = 'X', nt1 = 'X', nt0 = 'X';
  int number_aa = 0;

  *j = starti;
  while (number_aa < 3) {
    if (pairs[*j].gapp == true) {
    } else if (pairs[*j].cdna == ' ') {
    } else {
      nt0 = nt1;
      nt1 = nt2;
      nt2 = revcompp ? complCode[toupper((int) pairs[*j].cdna)] : toupper((int) pairs[*j].cdna);
      number_aa++;
    }
    (*j)++;
  }

  /* This avoids putting aa in an intron */
  while (*j < npairs-1 && pairs[*j].gapp == true) {
    (*j)++;
  }

  return get_codon(nt0,nt1,nt2);
}

/* The only difference between standard and alt is the precedence of
   aamarker and cdna == ' ' */
static int
count_aa_backward_standard (int aapos, struct Pair_T *pairs, 
			    int npairs, int start, int endi) {
  int number_aa = 0, j;

  j = endi;
  if (pairs[j].gapp == true) {
  } else if (pairs[j].cdna == ' ') {
  } else {
    number_aa++;
  }
  j--;

  while (j >= start) {
    if (pairs[j].gapp == true) {
    } else if (pairs[j].aamarker == true) {
      return number_aa;
    } else if (pairs[j].cdna == ' ') {
    } else {
      number_aa++;
    }
    j--;
  }

  while (j >= 0 && number_aa < 3) {
    if (pairs[j].gapp == true) {
    } else if (pairs[j].aamarker == true) {
      return number_aa;
    } else if (pairs[j].cdna == ' ') {
    } else {
      number_aa++;
    }
    j--;
  }

  return number_aa;
}

static int
count_aa_backward_alt (int *j, int aapos, struct Pair_T *pairs, 
		       int npairs, int start, int endi) {
  int number_aa = 0;

  *j = endi;
  if (pairs[*j].gapp == true) {
  } else if (pairs[*j].cdna == ' ') {
  } else {
    number_aa++;
  }
  (*j)--;

  while (*j >= start) {
    if (pairs[*j].gapp == true) {
    } else if (pairs[*j].cdna == ' ') {
    } else if (pairs[*j].aamarker == true) {
      return number_aa;
    } else {
      number_aa++;
    }
    (*j)--;
  }

  while (*j >= 0 && number_aa < 3) {
    if (pairs[*j].gapp == true) {
    } else if (pairs[*j].cdna == ' ') {
    } else if (pairs[*j].aamarker == true) {
      return number_aa;
    } else {
      number_aa++;
    }
    (*j)--;
  }

  return number_aa;
}

static int
get_codon_backward (int *j, int aapos, struct Pair_T *pairs, int start, int endi, 
		    char *complCode, bool revcompp) {
  char nt2 = 'X', nt1 = 'X', nt0 = 'X';
  int number_aa = 0;

  *j = endi;
  while (number_aa < 3) {
    if (pairs[*j].gapp == true) {
    } else if (pairs[*j].cdna == ' ') {
    } else {
      nt0 = nt1;
      nt1 = nt2;
      nt2 = revcompp ? complCode[toupper((int) pairs[*j].cdna)] : toupper((int) pairs[*j].cdna);
      number_aa++;
    }
    (*j)--;
  }

  /* This avoids putting aa in an intron */
  while (*j > 0 && pairs[*j].gapp == true) {
    (*j)--;
  }

  return get_codon(nt0,nt1,nt2);
}



#ifdef PMAP
static void
translate_genomic_via_marker (struct Pair_T *pairs, int npairs, int start, int end) {
  struct Pair_T *pair;
  int i, j, number_aa;

  debug2(printf("translate_genomic called with %d %d\n",start,end));

  i = start;
  while (i <= end) {
    pair = &(pairs[i]);
    if (pair->aamarker == true) {
      number_aa = count_aa_genomic_standard(pair->aapos,pairs,npairs,i,end);
      if (number_aa == 3) {
	pair->aa_g = get_codon_genomic(&j,pair->aapos,pair->aa_e,pairs,npairs,i,end);
	i += 3;

      } else if (number_aa < 3) {
	debug2(printf("At %d, saw %d < 3 aa standard\n",pair->aapos,number_aa));
	
	number_aa = count_aa_genomic_alt(&j,pair->aapos,pairs,npairs,i,end);
	debug2(printf("At %d, saw %d aa alt\n",pair->aapos,number_aa));

	if (number_aa == 3) {
	  /* Interpret */
	  pair->aa_g = get_codon_genomic(&j,pair->aapos,pair->aa_e,pairs,npairs,i,end);
	  i = j;
	} else if (number_aa < 2) {
	  /* Skip */
	  debug2(printf("Skipping\n"));
	  i = j;
	} else {
	  /* Substitute */
	  pair->aa_g = pair->aa_e;
	  i++;
	}

      } else if (number_aa < 6) {
	/* Insertion of less than a codon */
	debug2(printf("At %d, saw %d < 6 aa\n",pair->aapos,number_aa));
	pair->aa_g = pair->aa_e;
	count_aa_genomic_alt(&j,pair->aapos,pairs,npairs,i,end);
	i = j;

      } else {
	/* Insertion of one or more codons */
	debug2(printf("At %d, saw %d >= 6 aa\n",pair->aapos,number_aa));
	pair->aa_g = get_codon_genomic(&j,pair->aapos,pair->aa_e,pairs,npairs,i,end);

	number_aa -= 3;
	while (number_aa >= 3) {
	  i = j;
	  pairs[i].aa_g = get_codon_genomic(&j,pair->aapos,pair->aa_e,pairs,npairs,i,end);

	  number_aa -= 3;
	}
      }
    } else {
      i++;
    }
  }

  return;
}  
#endif


static void
translate_est_forward_via_marker (struct Pair_T *pairs, int npairs, bool revcompp,
				  int start, int end) {
  struct Pair_T *pair;
  int i, j, number_aa;
  char complCode[128] = COMPLEMENT;

  debug2(printf("translate_est_forward called with %d %d\n",start,end));

  i = start;
  while (i <= end) {
    pair = &(pairs[i]);
    if (pair->aamarker == true) {
      number_aa = count_aa_forward_standard(pair->aapos,pairs,npairs,i,end);
      if (number_aa == 3) {
	pair->aa_e = get_codon_forward(&j,pair->aapos,pairs,npairs,i,end,complCode,
				       revcompp);
	i += 3;

      } else if (number_aa < 3) {
	debug2(printf("At %d, saw %d < 3 aa standard\n",pair->aapos,number_aa));
	
	number_aa = count_aa_forward_alt(&j,pair->aapos,pairs,npairs,i,end);
	debug2(printf("At %d, saw %d aa alt\n",pair->aapos,number_aa));

	if (number_aa == 3) {
	  /* Interpret */
	  pair->aa_e = get_codon_forward(&j,pair->aapos,pairs,npairs,i,end,complCode,
					 revcompp);
	  i = j;
	} else if (number_aa < 2) {
	  /* Skip */
	  debug2(printf("Skipping\n"));
	  i = j;
	} else {
	  /* Substitute */
	  pair->aa_e = pair->aa_g;
	  i++;
	}

      } else if (number_aa < 6) {
	/* Insertion of less than a codon */
	debug2(printf("At %d, saw %d < 6 aa\n",pair->aapos,number_aa));
	pair->aa_e = pair->aa_g;
	count_aa_forward_alt(&j,pair->aapos,pairs,npairs,i,end);
	i = j;

      } else {
	/* Insertion of one or more codons */
	debug2(printf("At %d, saw %d >= 6 aa\n",pair->aapos,number_aa));
	pair->aa_e = get_codon_forward(&j,pair->aapos,pairs,npairs,i,end,complCode,
				       revcompp);
	number_aa -= 3;
	while (number_aa >= 3) {
	  i = j;
	  pairs[i].aa_e = get_codon_forward(&j,pair->aapos,pairs,npairs,i,end,complCode,
					    revcompp);
	  number_aa -= 3;
	}
      }
    } else {
      i++;
    }
  }

  return;
}  


static void
translate_est_backward_via_marker (struct Pair_T *pairs, int npairs, bool revcompp,
				   int start, int end) {
  struct Pair_T *pair;
  int i, j, number_aa;
  char complCode[128] = COMPLEMENT;

  debug2(printf("translate_est_backward called with %d %d\n",start,end));

  i = end;
  while (i >= start) {
    pair = &(pairs[i]);
    if (pair->aamarker == true) {
      number_aa = count_aa_backward_standard(pair->aapos,pairs,npairs,start,i);
      if (number_aa == 3) {
	pair->aa_e = get_codon_backward(&j,pair->aapos,pairs,start,i,complCode,
					revcompp);
	i -= 3;

      } else if (number_aa < 3) {
	debug2(printf("At %d, saw %d < 3 aa standard\n",pair->aapos,number_aa));
	
	number_aa = count_aa_backward_alt(&j,pair->aapos,pairs,npairs,start,i);
	debug2(printf("At %d, saw %d aa alt\n",pair->aapos,number_aa));

	if (number_aa == 3) {
	  /* Interpret */
	  pair->aa_e = get_codon_backward(&j,pair->aapos,pairs,start,i,complCode,
					 revcompp);
	  i = j;
	} else if (number_aa < 2) {
	  /* Skip */
	  debug2(printf("Skipping\n"));
	  i = j;
	} else {
	  /* Substitute */
	  pair->aa_e = pair->aa_g;
	  i--;
	}

      } else if (number_aa < 6) {
	/* Insertion of less than a codon */
	debug2(printf("At %d, saw %d < 6 aa\n",pair->aapos,number_aa));
	pair->aa_e = pair->aa_g;
	count_aa_backward_alt(&j,pair->aapos,pairs,npairs,start,i);
	i = j;

      } else {
	/* Insertion of one or more codons */
	debug2(printf("At %d, saw %d >= 6 aa\n",pair->aapos,number_aa));
	pair->aa_e = get_codon_backward(&j,pair->aapos,pairs,start,i,complCode,
				       revcompp);
	number_aa -= 3;
	while (number_aa >= 3) {
	  i = j;
	  pairs[i].aa_e = get_codon_backward(&j,pair->aapos,pairs,start,i,complCode,
					    revcompp);
	  number_aa -= 3;
	}
      }
    } else {
      i--;
    }
  }

  return;
}  


#ifdef PMAP
void
Translation_via_cdna (int *translation_leftpos, int *translation_rightpos, int *translation_length,
		      int *relaastart, int *relaaend,
		      struct Pair_T *pairs, int npairs, char *queryaaseq_ptr) {
  int translation_start, translation_end, i;

  for (i = 0; i < npairs; i++) {
    pairs[i].refquerypos = pairs[i].querypos;
    pairs[i].aa_g = pairs[i].aa_e = ' ';
  }

  translate_pairs_cdna(&translation_start,&(*translation_length),pairs,npairs,queryaaseq_ptr);
  translation_end = npairs-1;

  /* Take care of genomic side */
    
  *relaastart = pairs[translation_start].aapos;
  *relaaend = pairs[translation_end].aapos;

  translate_genomic_via_marker(pairs,npairs,translation_start,translation_end);

  return;
}
#else
void
Translation_via_genomic (int *translation_leftpos, int *translation_rightpos, int *translation_length,
			 int *relaastart, int *relaaend,
			 struct Pair_T *pairs, int npairs, bool backwardp, bool revcompp, bool fulllengthp) {
  char lastaa;
  struct T *translation;
  bool endstopp;
  int i, aapos = 0;
  int translation_frame, translation_start = 0, translation_end = 0;
  int minpos, maxpos;

  if (backwardp == false) {
    translation = translate_pairs_forward(pairs,npairs,revcompp);
    find_bounds_forward(&translation_frame,&translation_start,
			&translation_end,&(*translation_length),&endstopp,
			translation,npairs,fulllengthp);
    if (fulllengthp == true && *translation_length == 0) {
      /* fprintf(stderr,"No full length gene found.  Assuming partial length.\n"); */
      find_bounds_forward(&translation_frame,&translation_start,
			  &translation_end,&(*translation_length),&endstopp,
			  translation,npairs,/*fulllengthp*/false);
    }

  } else {
    translation = translate_pairs_backward(pairs,npairs,revcompp);
    find_bounds_backward(&translation_frame,&translation_start,
			 &translation_end,&(*translation_length),&endstopp,
			 translation,npairs,fulllengthp);
    if (fulllengthp == true && *translation_length == 0) {
      /* fprintf(stderr,"No full length gene found.  Assuming partial length.\n"); */
      find_bounds_backward(&translation_frame,&translation_start,
			   &translation_end,&(*translation_length),&endstopp,
			   translation,npairs,/*fulllengthp*/false);
    }
  }
  /* debug(printf("ref:\n")); */
  debug(Translation_dump(pairs,translation,npairs));

  for (i = 0; i < npairs; i++) {
    pairs[i].refquerypos = pairs[i].querypos;
    pairs[i].aa_g = pairs[i].aa_e = ' ';
  }

  if (translation_start < 0 || translation_end < 0) {
    *translation_leftpos = *translation_rightpos = -1;
    *relaastart = *relaaend = -1;
  } else {
    minpos = npairs;
    maxpos = 0;
    if (backwardp == false) {
      debug(printf("Translation is forward pairs %d..%d\n",translation_start,translation_end));
      for (i = translation_start; i <= translation_end; i++) {
	/* Avoid problems with genome position advancing prematurely */
	if (pairs[i].genome != ' ') {
	  if (translation[i].frame == translation_frame) {
	    if ((pairs[i].aa_g = translation[i].aa) != ' ') {
	      if (pairs[i].querypos < minpos) {
		minpos = pairs[i].querypos;
	      }
	      if (pairs[i].querypos > maxpos) {
		maxpos = pairs[i].querypos;
	      }
	      lastaa = pairs[i].aa_g;
	      aapos++;
	      pairs[i].aamarker = true;
	    }
	  }
	}
	pairs[i].aapos = aapos;
      }
      *translation_leftpos = minpos;
      *translation_rightpos = maxpos;
      if (lastaa == '*') {
	*translation_length -= 1;
      }
      
      if (i < npairs) {
	*translation_rightpos += 1;
	pairs[i++].aapos = aapos;
      }
      if (i < npairs) {
	*translation_rightpos += 1;
	pairs[i].aapos = aapos;
      }
    

    } else {
      debug(printf("Translation is backward pairs %d..%d\n",translation_start,translation_end));
      
      for (i = translation_start; i >= translation_end; --i) {
	/* Avoid problems with genome position advancing prematurely */
	if (pairs[i].genome != ' ') {
	  if (translation[i].frame == translation_frame) {
	    if ((pairs[i].aa_g = translation[i].aa) != ' ') {
	      if (pairs[i].querypos < minpos) {
		minpos = pairs[i].querypos;
	      }
	      if (pairs[i].querypos > maxpos) {
		maxpos = pairs[i].querypos;
	      }
	      lastaa = pairs[i].aa_g;
	      aapos++;
	      pairs[i].aamarker = true;
	    }
	  }
	}
	pairs[i].aapos = aapos;
      }
      *translation_leftpos = minpos;
      *translation_rightpos = maxpos;
      if (lastaa == '*') {
	*translation_length -= 1;
      }
      
      if (i >= 0) {
	*translation_leftpos -= 1;
	pairs[i--].aapos = aapos;
      }
      if (i >= 0) {
	*translation_leftpos -= 1;
	pairs[i].aapos = aapos;
      }
    }
    
    /* Take care of cDNA side */
    
    debug(printf("Translation start = %d\n",translation_start));
    debug(printf("Translation end = %d\n",translation_end));
    
    *relaastart = pairs[translation_start].aapos;
    *relaaend = pairs[translation_end].aapos;

    if (revcompp == false) {
      translate_est_forward_via_marker(pairs,npairs,revcompp,translation_start,translation_end);
    } else {
      translate_est_backward_via_marker(pairs,npairs,revcompp,translation_end,translation_start);
    }
  }

  FREE(translation);

  return;
}
#endif

/* Pairs are always ordered by ascending cDNA position and genomic
   position.  For Crick strand matches, the genomic position needs to
   be standardized. */
static void
bound_via_reference (int *start, int *end, struct Pair_T *pairs, int npairs, bool watsonp, 
		     struct Pair_T *refpairs, int nrefpairs, bool refwatsonp,
		     int genomiclength) {
  int i, j, aapos = 0;
  int refquerypos, genomepos, refgenomepos;

  debug(Pair_dump_array(refpairs,nrefpairs,false));

  *start = *end = -1;
  if (refwatsonp == true && watsonp == true) {
    debug1(printf("refwatsonp == true && watsonp == true\n"));
    i = 0;
    j = 0;
    while (i < nrefpairs && j < npairs) {
      refquerypos = refpairs[i].querypos;
      refgenomepos = refpairs[i].genomepos;
      genomepos = pairs[j].genomepos;
      debug(printf("Comparing ref %d (%c) with %d\n",refgenomepos,refpairs[i].aa_e,genomepos));
      if (pairs[j].genome == ' ') {
	debug(printf("Not incrementing aapos %d\n",aapos));
	pairs[j].refquerypos = refquerypos;
	pairs[j++].aapos = aapos;
      } else if (refgenomepos < genomepos) {
	refquerypos = refpairs[i].querypos;
	aapos = refpairs[i++].aapos;
      } else if (genomepos < refgenomepos) {
	pairs[j].refquerypos = refquerypos;
	pairs[j++].aapos = aapos;
      } else {
	debug(printf("looking at refpairs %d, genomepos %d (%c)\n",
		     i,refgenomepos,refpairs[i].aa_e));
	if (refpairs[i].aa_e != ' ') {
	  if (*start < 0) {
	    *start = j;
	  }
	  *end = j;
	}
	pairs[j].aamarker = refpairs[i].aamarker;
	refquerypos = refpairs[i].querypos;
	aapos = refpairs[i++].aapos;
	pairs[j].refquerypos = refquerypos;
	pairs[j++].aapos = aapos;
      }
    }

    /*
    if (*end < npairs-1) {
      pairs[++*end].aapos = aapos;
    }
    if (*end < npairs-1) {
      pairs[++*end].aapos = aapos;
    }
    */

  } else if (refwatsonp == true && watsonp == false) {
    debug1(printf("refwatsonp == true && watsonp == false\n"));
    i = 0;
    j = npairs-1;
    while (i < nrefpairs && j >= 0) {
      refquerypos = refpairs[i].querypos;
      refgenomepos = refpairs[i].genomepos;
      genomepos = (genomiclength - 1) - pairs[j].genomepos;
      debug(printf("Comparing ref %d (%c) with %d\n",refgenomepos,refpairs[i].aa_e,genomepos));
      if (pairs[j].genome == ' ') {
	pairs[j].refquerypos = refquerypos;
	pairs[j--].aapos = aapos;
      } else if (refgenomepos < genomepos) {
	refquerypos = refpairs[i].querypos;
	aapos = refpairs[i++].aapos;
      } else if (genomepos < refgenomepos) {
	pairs[j].refquerypos = refquerypos;
	pairs[j--].aapos = aapos;
      } else {
	if (refpairs[i].aa_e != ' ') {
	  if (*end < 0) {
	    *end = j;
	    /*
	    aapos = refpairs[i].aapos;
	    if (*end < npairs-1) {
	      pairs[++*end].aapos = aapos;
	    }
	    if (*end < npairs-1) {
	      pairs[++*end].aapos = aapos;
	    }
	    */
	  }
	  *start = j;
	}
	pairs[j].aamarker = refpairs[i].aamarker;
	refquerypos = refpairs[i].querypos;
	aapos = refpairs[i++].aapos;
	pairs[j].refquerypos = refquerypos;
	pairs[j--].aapos = aapos;
      }
    }

  } else if (refwatsonp == false && watsonp == true) {
    debug1(printf("refwatsonp == false && watsonp == true\n"));
    i = nrefpairs-1;
    j = 0;
    while (i >= 0 && j < npairs) {
      refquerypos = refpairs[i].querypos;
      refgenomepos = (genomiclength - 1) - refpairs[i].genomepos;
      genomepos = pairs[j].genomepos;
      debug(printf("Comparing ref %d (%c) with %d\n",refgenomepos,refpairs[i].aa_e,genomepos));
      if (pairs[j].genome == ' ') {
	pairs[j].refquerypos = refquerypos;
	pairs[j++].aapos = aapos;
      } else if (refgenomepos < genomepos) {
	refquerypos = refpairs[i].querypos;
	aapos = refpairs[i--].aapos;
      } else if (genomepos < refgenomepos) {
	pairs[j].refquerypos = refquerypos;
	pairs[j++].aapos = aapos;
      } else {
	if (refpairs[i].aa_e != ' ') {
	  if (*start < 0) {
	    *start = j;
	  }
	  *end = j;
	}
	pairs[j].aamarker = refpairs[i].aamarker;
	refquerypos = refpairs[i].querypos;
	aapos = refpairs[i--].aapos;
	pairs[j].refquerypos = refquerypos;
	pairs[j++].aapos = aapos;
      }
    }

    /*
    if (*end < npairs-1) {
      pairs[++*end].aapos = aapos;
    }
    if (*end < npairs-1) {
      pairs[++*end].aapos = aapos;
    }
    */

  } else {
    debug1(printf("refwatsonp == false && watsonp == false\n"));
    i = nrefpairs-1;
    j = npairs-1;
    while (i >= 0 && j >= 0) {
      refquerypos = refpairs[i].querypos;
      refgenomepos = (genomiclength - 1) - refpairs[i].genomepos;
      genomepos = (genomiclength - 1) - pairs[j].genomepos;
      debug(printf("Comparing ref %d (%c) with %d\n",refgenomepos,refpairs[i].aa_e,genomepos));
      if (pairs[j].genome == ' ') {
	pairs[j].refquerypos = refquerypos;
	pairs[j--].aapos = aapos;
      } else if (refgenomepos < genomepos) {
	refquerypos = refpairs[i].querypos;
	aapos = refpairs[i--].aapos;
      } else if (genomepos < refgenomepos) {
	pairs[j].refquerypos = refquerypos;
	pairs[j--].aapos = aapos;
      } else {
	if (refpairs[i].aa_e != ' ') {
	  if (*end < 0) {
	    *end = j;
	    /*
	    aapos = refpairs[i].aapos;
	    if (*end < npairs-1) {
	      pairs[++*end].aapos = aapos;
	    }
	    if (*end < npairs-1) {
	      pairs[++*end].aapos = aapos;
	    }
	    */
	  }
	  *start = j;
	}
	pairs[j].aamarker = refpairs[i].aamarker;
	refquerypos = refpairs[i].querypos;
	aapos = refpairs[i--].aapos;
	pairs[j].refquerypos = refquerypos;
	pairs[j--].aapos = aapos;
      }
    }
  }    

  debug(printf("Bound by reference: start = %d, end = %d\n",*start,*end));
  debug(Pair_dump_array(pairs,npairs,false));

  return;
}


void
Translation_via_reference (int *relaastart, int *relaaend,
			   struct Pair_T *pairs, int npairs, bool watsonp, bool backwardp, bool revcompp,
			   struct Pair_T *refpairs, int nrefpairs, bool refwatsonp, int genomiclength,
			   bool fixshiftp) {
  struct T *translation;
  int start, end, i;

  for (i = 0; i < npairs; i++) {
    pairs[i].aa_g = pairs[i].aa_e = ' ';
  }

  bound_via_reference(&start,&end,pairs,npairs,watsonp,refpairs,nrefpairs,refwatsonp,
		      genomiclength);

  if (start < 0 || end < 0) {
    *relaastart = *relaaend = -1;
  } else {
    *relaastart = pairs[start].aapos;
    *relaaend = pairs[end].aapos;

    if (backwardp == false) {
      translation = translate_pairs_forward(pairs,npairs,revcompp);
    } else {
      translation = translate_pairs_backward(pairs,npairs,revcompp);
    }

    for (i = 0; i < npairs; i++) {
      if (pairs[i].aamarker == true) {
	pairs[i].aa_g = translation[i].aa;
      } else {
	pairs[i].aa_g = ' ';
      }
      pairs[i].aa_e = ' ';
    }

    debug2(printf("Bounding %d %d to yield %d %d at %c %c\n",
		  start,end,*relaastart,*relaaend,
		  pairs[start].aa_g,pairs[end].aa_g));

    if (backwardp == false) {
      translate_est_forward_via_marker(pairs,npairs,revcompp,start,end);
    } else {
      translate_est_backward_via_marker(pairs,npairs,revcompp,start,end);
    }

    debug1(Translation_dump(pairs,translation,npairs));
    FREE(translation);

  }

  return;
}

static char
lookup_aa (int aapos, struct Pair_T *refpairs, int nrefpairs) {
  int i;

  for (i = 0; i < nrefpairs; i++) {
    if (refpairs[i].aapos == aapos && refpairs[i].aamarker == true) {
      return refpairs[i].aa_g;
    }
  }
  return 'X';
}

static List_T
push_mutations (List_T mutations, struct Pair_T *pair) {
  if (pair->aa_g == ' ') {
    debug3(printf("  Creating insertion mutation at %d\n",pair->aapos));
    return List_push(mutations,Mutation_insertion_new(pair->refquerypos,pair->aapos,pair->aa_e,1));
  } else if (pair->aa_e == ' ') {
    debug3(printf("  Creating deletion mutation at %d\n",pair->aapos));
    return List_push(mutations,Mutation_deletion_new(pair->refquerypos,pair->aapos,pair->aa_g,
						     pair->aapos,pair->aa_g));
  } else if (pair->aa_g == 'X') {
    return mutations;
  } else if (pair->aa_e == 'X') {
    return mutations;
  } else {
    debug3(printf("  Creating substitution mutation at %d\n",pair->aapos));
    return List_push(mutations,Mutation_substitution_new(pair->refquerypos,pair->aapos,pair->aa_g,
							 pair->aa_e));
  }
}


void
Translation_compare (struct Pair_T *pairs, int npairs, struct Pair_T *refpairs, int nrefpairs,
		     int cdna_direction, int relaastart, int relaaend) {
  List_T mutations = NULL, newmutations, p, q;
  Mutation_T x, y, z, *array;
  int i;
  int lastaa, nmutations;
  bool changep = true;
  char aa1, aa2;
  
  if (relaastart < relaaend) {
    for (i = IGNORE_MARGIN; i < npairs - IGNORE_MARGIN; i++) {
      if (toupper((int) pairs[i].aa_g) != toupper((int) pairs[i].aa_e)) {
	mutations = push_mutations(mutations,&(pairs[i]));
      }
    }

    if (refpairs != NULL) {
      lastaa = relaastart - 1;
      for (i = 0; i < npairs; i++) {
	if (pairs[i].gapp == false && pairs[i].aamarker == true) {
	  if (pairs[i].aapos != lastaa + 1) {
	    debug3(printf("  Observed deletion mutation at %d != %d + 1\n",
			  pairs[i].aapos,lastaa));
	    aa1 = lookup_aa(lastaa+1,refpairs,nrefpairs);
	    aa2 = lookup_aa(pairs[i].aapos-1,refpairs,nrefpairs);
	    mutations = List_push(mutations,Mutation_deletion_new(pairs[i].refquerypos,lastaa+1,aa1,
								  pairs[i].aapos-1,aa2));
	  }
	  lastaa = pairs[i].aapos;
	}
      }
    }
  } else {
    for (i = npairs - 1 - IGNORE_MARGIN; i >= IGNORE_MARGIN; i--) {
      if (pairs[i].aa_g == pairs[i].aa_e) {
      } else {
	mutations = push_mutations(mutations,&(pairs[i]));
      }
    }

    if (refpairs != NULL) {
      lastaa = relaaend - 1;
      for (i = npairs - 1; i >= 0; i--) {
	if (pairs[i].gapp == false && pairs[i].aamarker == true) {
	  if (pairs[i].aapos != lastaa + 1) {
	    debug3(printf("  Observed deletion mutation at %d != %d + 1\n",
			  pairs[i].aapos,lastaa));
	    aa1 = lookup_aa(lastaa+1,refpairs,nrefpairs);
	    aa2 = lookup_aa(pairs[i].aapos-1,refpairs,nrefpairs);
	    mutations = List_push(mutations,Mutation_deletion_new(pairs[i].refquerypos,lastaa+1,aa1,
								  pairs[i].aapos-1,aa2));
	  }
	  lastaa = pairs[i].aapos;
	}
      }
    }
  }

  /* Fix multiple insertions and deletions */
  while (changep == true) {
    changep = false;
    newmutations = NULL;

    /* Handle insertions */
    for (p = mutations; p != NULL; p = List_next(p)) {
      x = (Mutation_T) List_head(p);
      if (Mutation_type(x) == INSERTION) {
	for (q = List_next(p); q != NULL; q = List_next(q)) {
	  y = (Mutation_T) List_head(q);
	  if (Mutation_type(y) == INSERTION) {
	    if ((z = Mutation_merge_insertions(x,y)) != NULL) {
	      newmutations = List_push(newmutations,z);
	      Mutation_invalidate(x);
	      Mutation_invalidate(y);
	      changep = true;
	    }
	  }
	}
      }
    }

    /* Handle deletions */
    for (p = mutations; p != NULL; p = List_next(p)) {
      x = (Mutation_T) List_head(p);
      if (Mutation_type(x) == DELETION) {
	for (q = List_next(p); q != NULL; q = List_next(q)) {
	  y = (Mutation_T) List_head(q);
	  if (Mutation_type(y) == DELETION) {
	    if ((z = Mutation_merge_deletions(x,y)) != NULL) {
	      newmutations = List_push(newmutations,z);
	      Mutation_invalidate(x);
	      Mutation_invalidate(y);
	      changep = true;
	    }
	  }
	}
      }
    }
    
    /* Remove invalids */
    while (mutations != NULL) {
      mutations = List_pop(mutations,(void **) &x);
      if (Mutation_validp(x) == false) {
	Mutation_free(&x);
      } else {
	newmutations = List_push(newmutations,x);
      }
    }

    mutations = newmutations;
  }

  /* Handle single events at ends of segmental changes */
  changep = true;
  while (changep == true) {
    changep = false;
    for (p = mutations; p != NULL; p = List_next(p)) {
      x = (Mutation_T) List_head(p);
      if (Mutation_naa(x) > 1) {
	for (q = mutations; q != NULL; q = List_next(q)) {
	  y = (Mutation_T) List_head(q);
	  if (Mutation_naa(y) == 1 && Mutation_at_ends_p(x,y) == true) {
	    debug3(printf("Invalidating "));
	    debug3(Mutation_print(y));
	    debug3(printf("\n"));
	    Mutation_invalidate(y);
	    changep = true;
	  }
	}
      }
    }

    if (changep == true) {
      /* Remove invalids */
      newmutations = NULL;
      while (mutations != NULL) {
	mutations = List_pop(mutations,(void **) &x);
	if (Mutation_validp(x) == false) {
	  Mutation_free(&x);
	} else {
	  newmutations = List_push(newmutations,x);
	}
      }
      mutations = newmutations;
    }
  }

  printf("    Mutations: ");

  nmutations = List_length(mutations);
  if (nmutations > 0) {
    array = (Mutation_T *) List_to_array(mutations,NULL);
    List_free(&mutations);
    qsort(array,nmutations,sizeof(Mutation_T),Mutation_cmp);

    x = array[0];
    Mutation_print(x);

    if (Mutation_type(x) == SUBSTITUTION) {
      printf(" (");
      Pair_dump_aapos(pairs,npairs,Mutation_aapos(x),cdna_direction);
      printf(")");
      printf(" [%d]",Mutation_refquerypos(x)+1+Pair_codon_changepos(pairs,npairs,Mutation_aapos(x),cdna_direction));
    } else {
      printf(" [%d]",Mutation_refquerypos(x)+1);
    }

    Mutation_free(&x);
    for (i = 1; i < nmutations; i++) {
      printf(", ");
      x = array[i];
      Mutation_print(x);

      if (Mutation_type(x) == SUBSTITUTION) {
	printf(" (");
	Pair_dump_aapos(pairs,npairs,Mutation_aapos(x),cdna_direction);
	printf(")");
	printf(" [%d]",Mutation_refquerypos(x)+1+Pair_codon_changepos(pairs,npairs,Mutation_aapos(x),cdna_direction));
      } else {
	printf(" [%d]",Mutation_refquerypos(x)+1);
      }

      Mutation_free(&x);
    }


    FREE(array);
  }

  printf("\n");
  return;
}


