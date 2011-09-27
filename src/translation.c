static char rcsid[] = "$Id: translation.c,v 1.58 2006/03/02 02:20:10 twu Exp $";
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


#define IGNORE_MARGIN 6
#define HORIZON 99
#define MIN_NPAIRS 30


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

/* Mark and assign cdna */
#ifdef DEBUG2
#define debug2(x) x
#else
#define debug2(x)
#endif



typedef enum {FRAME0, FRAME1, FRAME2, NOFRAME} Frame_T;

#define T Translation_T
struct T {
  char aa;
  Frame_T frame;
};


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
    if (pairs[i].aamarker_g == true || pairs[i].aamarker_e == true) {
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


#ifdef PMAP
static void
translate_pairs_cdna (int *translation_start, int *translation_end, int *translation_length,
		      struct Pair_T *pairs, int npairs, char *queryaaseq_ptr) {
  struct Pair_T *ptr, *pair;
  int i;

  /* Go backward so we put aa at beginning of codon */
  /* printf("lastquerypos is %d\n",pairs[npairs-1].querypos); */
  switch (pairs[npairs-1].querypos % 3) {
  case 0: *translation_end = npairs-2; break; /* codon is on last character */
  case 1: *translation_end = npairs-3; break;
  case 2: *translation_end = npairs-1; break;
  }
  if (*translation_end < 0) {
    *translation_end = 0;
  }

  *translation_start = *translation_end;
  *translation_length = 0;

  ptr = &(pairs[(*translation_end)+1]);
  for (i = *translation_end; i >= 0; --i) {
    pair = --ptr;
    pair->aapos = pair->querypos/3 + 1;
    if (pair->gapp == true) {
      /* pair->aa_e = ' '; */
    } else if (pair->cdna == ' ') {
      /* pair->aa_e = ' '; */
    } else {
      if (pair->querypos % 3 == 0) {
	pair->aa_e = queryaaseq_ptr[pair->querypos/3];
	pair->aamarker_e = true;
	*translation_start = i;
	(*translation_length) += 1;
      }
    }
  }

  return;
}  
#endif


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


/************************************************************************/


static int
count_cdna_forward (int *nexti, struct Pair_T *pairs, int npairs, int starti, int end) {
  int ncdna = 0, j;

  j = starti;
  while (j <= end) {
    if (j > starti && pairs[j].aamarker_g == true && pairs[j].cdna != ' ') {
      *nexti = j;
      return ncdna;
    } else if (pairs[j].gapp == false && pairs[j].cdna != ' ') {
      ncdna++;
    }
    j++;
  }

  *nexti = j;
  return ncdna;
}


static int
count_cdna_forward_mod3 (int *nexti, struct Pair_T *pairs, int npairs, int starti, int end) {
  int ncdna = 0, j;
  
  j = starti;
  while (j <= end && ncdna <= HORIZON) {
    if (j > starti && pairs[j].aamarker_g == true && pairs[j].cdna != ' ' &&
	(ncdna % 3) == 0) {
      *nexti = j;
      return ncdna;
    } else if (pairs[j].gapp == false && pairs[j].cdna != ' ') {
      ncdna++;
    }
    j++;
  }

  *nexti = j;
  return 1;			/* any answer that is not mod 0 */
}

static int
count_cdna_backward (int *nexti, struct Pair_T *pairs, int npairs, int starti, int end) {
  int ncdna = 0, j;

  j = starti;
  while (j >= end) {
    if (j < starti && pairs[j].aamarker_g == true && pairs[j].cdna != ' ') {
      *nexti = j;
      return ncdna;
    } else if (pairs[j].gapp == false && pairs[j].cdna != ' ') {
      ncdna++;
    }
    j--;
  }

  *nexti = j;
  return ncdna;
}


static int
count_cdna_backward_mod3 (int *nexti, struct Pair_T *pairs, int npairs, int starti, int end) {
  int ncdna = 0, j;
  
  j = starti;
  while (j >= end && ncdna <= HORIZON) {
    if (j < starti && pairs[j].aamarker_g == true && pairs[j].cdna != ' ' &&
	(ncdna % 3) == 0) {
      *nexti = j;
      return ncdna;
    } else if (pairs[j].gapp == false && pairs[j].cdna != ' ') {
      ncdna++;
    }
    j--;
  }

  *nexti = j;
  return 1;			/* any answer that is not mod 0 */
}


static int
count_genomic (int *nexti, struct Pair_T *pairs, int npairs, int starti, int end) {
  int ngenomic = 0, j;

  j = starti;
  while (j <= end) {
    if (j > starti && pairs[j].aamarker_e == true && pairs[j].genome != ' ') {
      *nexti = j;
      return ngenomic;
    } else if (pairs[j].gapp == false && pairs[j].genome != ' ') {
      ngenomic++;
    }
    j++;
  }

  *nexti = j;
  return ngenomic;
}


static int
count_genomic_mod3 (int *nexti, struct Pair_T *pairs, int npairs, int starti, int end) {
  int ngenomic = 0, j;
  
  j = starti;
  while (j <= end && ngenomic <= HORIZON) {
    if (j > starti && pairs[j].aamarker_e == true && pairs[j].genome != ' ' &&
	(ngenomic % 3) == 0) {
      *nexti = j;
      return ngenomic;
    } else if (pairs[j].gapp == false && pairs[j].genome != ' ') {
      ngenomic++;
    }
    j++;
  }

  *nexti = j;
  return 1;			/* any answer that is not mod 0 */
}

/************************************************************************/

static int
get_codon_forward (int *nexti, struct Pair_T *pairs, int npairs, int starti, 
		   char *complCode, bool revcompp) {
  char nt2 = 'X', nt1 = 'X', nt0 = 'X';
  int ncdna = 0, j;

  j = starti;
  while (j < npairs && ncdna < 3) {
    if (pairs[j].gapp == false && pairs[j].cdna != ' ') {
      nt0 = nt1;
      nt1 = nt2;
      nt2 = revcompp ? complCode[toupper((int) pairs[j].cdna)] : toupper((int) pairs[j].cdna);
      ncdna++;
    }
    j++;
  }

  /* assign function depends on nexti pointing to a valid position */
  while (j <= npairs - 1 && (pairs[j].gapp == true || pairs[j].cdna == ' ')) {
    j++;
  }

  *nexti = j;
  return get_codon(nt0,nt1,nt2);
}

static int
get_codon_backward (int *nexti, struct Pair_T *pairs, int npairs, int starti,
		   char *complCode, bool revcompp) {
  char nt2 = 'X', nt1 = 'X', nt0 = 'X';
  int ncdna = 0, j;

  j = starti;
  while (j >= 0 && ncdna < 3) {
    if (pairs[j].gapp == false && pairs[j].cdna != ' ') {
      nt0 = nt1;
      nt1 = nt2;
      nt2 = revcompp ? complCode[toupper((int) pairs[j].cdna)] : toupper((int) pairs[j].cdna);
      ncdna++;
    }
    --j;
  }

  /* assign function depends on nexti pointing to a valid position */
  while (j >= 0 && (pairs[j].gapp == true || pairs[j].cdna == ' ')) {
    --j;
  }

  *nexti = j;
  return get_codon(nt0,nt1,nt2);
}


#ifdef PMAP
static int
get_codon_genomic (int *nexti, struct Pair_T *pairs, int npairs, int starti) {
  char nt2 = 'X', nt1 = 'X', nt0 = 'X';
  int ngenomic = 0, j;

  j = starti;
  while (ngenomic < 3) {
    if (pairs[j].gapp == false && pairs[j].genome != ' ') {
      nt0 = nt1;
      nt1 = nt2;
      nt2 = toupper(pairs[j].genome);
      ngenomic++;
    }
    j++;
  }

  /* assign function depends on nexti pointing to a valid position */
  while (j <= npairs-1 && (pairs[j].gapp == true || pairs[j].genome == ' ')) {
    j++;
  }

  *nexti = j;
  return get_codon(nt0,nt1,nt2);
}
#endif



static void
assign_cdna_forward (int ncdna, struct Pair_T *pairs, int npairs, char *complCode, bool revcompp, 
		     int starti) {
  struct Pair_T *pair;
  int i, nexti, j = 0;

  i = starti;
  while (j < ncdna) {
    pair = &(pairs[i]);
    pair->aamarker_e = true;
    pair->aa_e = get_codon_forward(&nexti,pairs,npairs,i,complCode,revcompp);
    i = nexti;
    j += 3;
  }
  return;
}

static void
terminate_cdna_forward (struct Pair_T *pairs, int npairs, char *complCode, bool revcompp, 
			int starti) {
  struct Pair_T *pair;
  int i, nexti;
  char lastcodon = ' ';

  i = starti;
  while (i <= npairs-3 && lastcodon != '*') {
    pair = &(pairs[i]);
    pair->aamarker_e = true;
    lastcodon = pair->aa_e = get_codon_forward(&nexti,pairs,npairs,i,complCode,revcompp);
    i = nexti;
  }
  return;
}

static void
assign_cdna_backward (int ncdna, struct Pair_T *pairs, int npairs, char *complCode, bool revcompp,
		      int starti) {
  struct Pair_T *pair;
  int i, nexti, j = 0;

  i = starti;
  while (j < ncdna) {
    pair = &(pairs[i]);
    pair->aamarker_e = true;
    pair->aa_e = get_codon_backward(&nexti,pairs,npairs,i,complCode,revcompp);
    i = nexti;
    j += 3;
  }
  return;
}

static void
terminate_cdna_backward (struct Pair_T *pairs, int npairs, char *complCode, bool revcompp, 
			int starti) {
  struct Pair_T *pair;
  int i, nexti;
  char lastcodon = ' ';

  i = starti;
  /* i > 1 is equivalent to i >= 2 */
  while (i > 1 && lastcodon != '*') {
    pair = &(pairs[i]);
    pair->aamarker_e = true;
    lastcodon = pair->aa_e = get_codon_backward(&nexti,pairs,npairs,i,complCode,revcompp);
    i = nexti;
  }
  return;
}

#ifdef PMAP
static void
assign_genomic (int ngenomic, struct Pair_T *pairs, int npairs, int starti) {
  struct Pair_T *pair;
  int i, nexti, j = 0;

  i = starti;
  while (j < ngenomic) {
    pair = &(pairs[i]);
    pair->aamarker_g = true;
    pair->aa_g = get_codon_genomic(&nexti,pairs,npairs,i);
    i = nexti;
    j += 3;
  }
  return;
}

static void
terminate_genomic (struct Pair_T *pairs, int npairs, int starti) {
  struct Pair_T *pair;
  int i, nexti;
  char lastcodon = ' ';

  i = starti;
  while (i <= npairs - 3 && lastcodon != '*') {
    pair = &(pairs[i]);
    pair->aamarker_g = true;
    lastcodon = pair->aa_g = get_codon_genomic(&nexti,pairs,npairs,i);
    i = nexti;
  }
  return;
}
#endif



static void
mark_cdna_forward (struct Pair_T *pairs, int npairs, bool revcompp, int start, int end) {
  struct Pair_T *pair;
  int i, nexti, nexti_alt, ncdna, ncdna_alt;
  char complCode[128] = COMPLEMENT;

  debug2(printf("mark_cdna_forward called with pairs #%d..%d\n",start,end));

  i = start;
  while (i < end) {
    pair = &(pairs[i]);
    if (pair->aamarker_g == false) {
      i++;
    } else {
      ncdna = count_cdna_forward(&nexti,pairs,npairs,i,end);
      if (ncdna == 3) {
	assign_cdna_forward(ncdna,pairs,npairs,complCode,revcompp,i);
      } else if ((ncdna % 3) == 0) {
	debug2(printf("At %d, saw %d with mod == 0\n",pair->aapos,ncdna));
	assign_cdna_forward(ncdna,pairs,npairs,complCode,revcompp,i);
      } else if (i + 2 > end) {
	debug2(printf("At %d, saw last codon: %d+2 > %d\n",pair->aapos,i,end));
	assign_cdna_forward(ncdna,pairs,npairs,complCode,revcompp,i);
      } else {
	debug2(printf("At %d, saw %d with mod != 0\n",pair->aapos,ncdna));	
	ncdna_alt = count_cdna_forward_mod3(&nexti_alt,pairs,npairs,i,end);

	debug2(printf("  Alternate search yields %d\n",ncdna_alt));
	if ((ncdna_alt % 3) == 0) {
	  debug2(printf("  Using alternate search\n"));
	  assign_cdna_forward(ncdna_alt,pairs,npairs,complCode,revcompp,i);
	  nexti = nexti_alt;

	} else if (ncdna < 3) {
	  debug2(printf("  Skipping\n"));

	} else {
	  debug2(printf("  Using original count - 3 = %d\n",ncdna-3));
	  assign_cdna_forward(ncdna-3,pairs,npairs,complCode,revcompp,i);
	}
      }
      i = nexti;
    }
  }

  terminate_cdna_forward(pairs,npairs,complCode,revcompp,i);

  return;
}  

static void
mark_cdna_backward (struct Pair_T *pairs, int npairs, bool revcompp, int start, int end) {
  struct Pair_T *pair;
  int i, nexti, nexti_alt, ncdna, ncdna_alt;
  char complCode[128] = COMPLEMENT;

  debug2(printf("mark_cdna_backward called with pairs #%d..%d\n",start,end));

  i = start;
  while (i > end) {
    pair = &(pairs[i]);
    if (pair->aamarker_g == false) {
      i--;
    } else {
      ncdna = count_cdna_backward(&nexti,pairs,npairs,i,end);
      if (ncdna == 3) {
	assign_cdna_backward(ncdna,pairs,npairs,complCode,revcompp,i);
      } else if ((ncdna % 3) == 0) {
	debug2(printf("At %d, saw %d with mod == 0\n",pair->aapos,ncdna));
	assign_cdna_backward(ncdna,pairs,npairs,complCode,revcompp,i);
      } else if (i - 2 < end) {
	debug2(printf("At %d, saw last codon: %d-2 < %d\n",pair->aapos,i,end));
	assign_cdna_backward(ncdna,pairs,npairs,complCode,revcompp,i);
      } else {
	debug2(printf("At %d, saw %d with mod != 0\n",pair->aapos,ncdna));	
	ncdna_alt = count_cdna_backward_mod3(&nexti_alt,pairs,npairs,i,end);

	debug2(printf("  Alternate search yields %d\n",ncdna_alt));
	if ((ncdna_alt % 3) == 0) {
	  debug2(printf("  Using alternate search\n"));
	  assign_cdna_backward(ncdna_alt,pairs,npairs,complCode,revcompp,i);
	  nexti = nexti_alt;

	} else if (ncdna < 3) {
	  debug2(printf("Skipping\n"));

	} else {
	  debug2(printf("  Using original count - 3 = %d\n",ncdna-3));
	  assign_cdna_backward(ncdna-3,pairs,npairs,complCode,revcompp,i);
	}
      }
      i = nexti;
    }
  }

  terminate_cdna_backward(pairs,npairs,complCode,revcompp,i);

  return;
}  


#ifdef PMAP
static void
mark_genomic (struct Pair_T *pairs, int npairs, int start, int end) {
  struct Pair_T *pair;
  int i, nexti, nexti_alt, ngenomic, ngenomic_alt;

  debug2(printf("mark_genomic called with %d %d\n",start,end));

  i = start;
  while (i < end) {
    pair = &(pairs[i]);
    if (pair->aamarker_e == false) {
      i++;
    } else {
      ngenomic = count_genomic(&nexti,pairs,npairs,i,end);
      if (ngenomic == 3) {
	assign_genomic(ngenomic,pairs,npairs,i);
      } else if ((ngenomic % 3) == 0) {
	debug2(printf("At %d, saw %d with mod == 0\n",pair->aapos,ngenomic));
	assign_genomic(ngenomic,pairs,npairs,i);
      } else {
	debug2(printf("At %d, saw %d with mod != 0\n",pair->aapos,ngenomic));
	ngenomic_alt = count_genomic_mod3(&nexti_alt,pairs,npairs,i,end);

	debug2(printf("  Alternate search yields %d\n",ngenomic_alt));
	if ((ngenomic_alt % 3) == 0) {
	  assign_genomic(ngenomic_alt,pairs,npairs,i);
	  nexti = nexti_alt;

	} else if (ngenomic < 3) {
	  debug2(printf("Skipping\n"));

	} else {
	  /* Use original count */
	  assign_genomic(ngenomic,pairs,npairs,i);
	}
      }
      i = nexti;
    }
  }

  terminate_genomic(pairs,npairs,i);

  return;
}  
#endif



#ifdef PMAP
void
Translation_via_cdna (int *translation_leftpos, int *translation_rightpos, int *translation_length,
		      int *relaastart, int *relaaend,
		      struct Pair_T *pairs, int npairs, char *queryaaseq_ptr) {
  int translation_start, translation_end, i;

  if (npairs < MIN_NPAIRS) {
    *translation_leftpos = 0;
    *translation_rightpos = 0;
    *translation_length = 0;
    *relaastart = 0;
    *relaaend = 0;
    return;
  }

  for (i = 0; i < npairs; i++) {
    pairs[i].refquerypos = pairs[i].querypos;
    pairs[i].aa_g = pairs[i].aa_e = ' ';
  }

  translate_pairs_cdna(&translation_start,&translation_end,&(*translation_length),
		       pairs,npairs,queryaaseq_ptr);

  *translation_leftpos = translation_start;
  *translation_rightpos = translation_end;

  /* Take care of genomic side */
    
  *relaastart = pairs[translation_start].aapos;
  *relaaend = pairs[translation_end].aapos;

  mark_genomic(pairs,npairs,translation_start,translation_end);

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

  if (npairs < MIN_NPAIRS) {
    return;
  }

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
      /* forward */
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
	      pairs[i].aamarker_g = true;
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
    
      /* Fill in aapos to the end */
      for ( ; i < npairs; i++) {
	pairs[i].aapos = aapos;
      }

    } else {
      /* backward */
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
	      pairs[i].aamarker_g = true;
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

      /* Fill in aapos to the end */
      for ( ; i >= 0; --i) {
	pairs[i].aapos = aapos;
      }
    }
    
    /* Take care of cDNA side */
    
    debug(printf("Translation start = %d\n",translation_start));
    debug(printf("Translation end = %d\n",translation_end));
    
    *relaastart = pairs[translation_start].aapos;
    *relaaend = pairs[translation_end].aapos;

    if (revcompp == false) {
      mark_cdna_forward(pairs,npairs,revcompp,translation_start,translation_end);
    } else {
      mark_cdna_backward(pairs,npairs,revcompp,translation_start,translation_end);
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
	pairs[j].aamarker_g = refpairs[i].aamarker_g;
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
	pairs[j].aamarker_g = refpairs[i].aamarker_g;
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
	pairs[j].aamarker_g = refpairs[i].aamarker_g;
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
	pairs[j].aamarker_g = refpairs[i].aamarker_g;
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

  if (npairs < MIN_NPAIRS) {
    return;
  }

  for (i = 0; i < npairs; i++) {
    pairs[i].aa_g = pairs[i].aa_e = ' ';
  }

  debug2(printf("Translation_via_reference called with backwardp = %d\n",backwardp));
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
      if (pairs[i].aamarker_g == true) {
	pairs[i].aa_g = translation[i].aa;
      } else {
	pairs[i].aa_g = ' ';
      }
      pairs[i].aa_e = ' ';
    }

    debug2(printf("Bounding pairs #%d..%d to yield amino acid coordinates %d (%c) to %d (%c)\n",
		  start,end,*relaastart,pairs[start].aa_g,*relaaend,pairs[end].aa_g));

    if (backwardp == false) {
      mark_cdna_forward(pairs,npairs,revcompp,start,end);
    } else {
      /* Note that we need to flip start and end here */
      mark_cdna_backward(pairs,npairs,revcompp,end,start);
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
    if (refpairs[i].aapos == aapos && refpairs[i].aamarker_g == true) {
      return refpairs[i].aa_g;
    }
  }
  return 'X';
}


static int
next_aapos_fwd (struct Pair_T *pairs, int i, int npairs, int aapos) {
  struct Pair_T *this;

  this = &(pairs[i]);
  while (i < npairs && this->aapos == aapos) {
    this++;
    i++;
  }

  /* Advance to next amino acid on query side, since aapos assures
     us we have reached it on genome side */
  while (i < npairs && this->aa_e == ' ') {
    this++;
    i++;
  }

  return i;
}

static int
next_aapos_rev (struct Pair_T *pairs, int i, int npairs, int aapos) {
  struct Pair_T *this;

  this = &(pairs[i]);
  while (i >= 0 && this->aapos == aapos) {
    --this;
    --i;
  }

  /* Advance to next amino acid on query side, since aapos assures
     us we have reached it on genome side */
  while (i >= 0 && this->aa_e == ' ') {
    --this;
    --i;
  }

  return i;
}


#define MAXMUT 100

static void
fill_aa_fwd (int *strlen_g, int *strlen_c, int *netchars, char *aa_genomicseg, char *aa_queryseq,
	     char *nt_genomicseg, char *nt_queryseq, struct Pair_T *start, struct Pair_T *end) {
  int k1 = 0, k2 = 0, i1 = 0, i2 = 0;
  struct Pair_T *this;

  *netchars = 0;
  for (this = start; this <= end && i1 < MAXMUT && k1 < MAXMUT; this++) {
    if (this->gapp == false) {
      if (this->genome != ' ') {
	nt_genomicseg[i1++] = toupper(this->genome);
      } else {
	*netchars += 1;
      }
      if (this->aa_g != ' ') {
	aa_genomicseg[k1++] = this->aa_g;
      }
    }
  }

  for (this = start; this <= end && i2 < MAXMUT && k2 < MAXMUT; this++) {
    if (this->gapp == false) {
      if (this->cdna != ' ') {
	nt_queryseq[i2++] = toupper(this->cdna);
      } else {
	*netchars -= 1;
      }
      if (this->aa_e != ' ') {
	aa_queryseq[k2++] = this->aa_e;
      }
    }
  }

  if (i1 >= MAXMUT || k1 >= MAXMUT || i2 >= MAXMUT || k2 >= MAXMUT) {
    i1 = i2 = k1 = k2 = 0;
  }
  nt_genomicseg[i1] = '\0';
  aa_genomicseg[k1] = '\0';
  nt_queryseq[i2] = '\0';
  aa_queryseq[k2] = '\0';

  *strlen_g = k1;
  *strlen_c = k2;

  return;
}

static void
fill_aa_rev (int *strlen_g, int *strlen_c, int *netchars, char *aa_genomicseg, char *aa_queryseq,
	     char *nt_genomicseg, char *nt_queryseq, struct Pair_T *start, struct Pair_T *end) {
  int k1 = 0, k2 = 0, i1 = 0, i2 = 0;
  struct Pair_T *this;

  *netchars = 0;
  for (this = end; this >= start && i1 < MAXMUT && k1 < MAXMUT; --this) {
    if (this->gapp == false) {
      if (this->genome != ' ') {
	nt_genomicseg[i1++] = toupper(this->genome);
      } else {
	*netchars += 1;
      }
      if (this->aa_g != ' ') {
	aa_genomicseg[k1++] = this->aa_g;
      }
    }
  }

  for (this = end; this >= start && i2 < MAXMUT && k2 < MAXMUT; --this) {
    if (this->gapp == false) {
      if (this->cdna != ' ') {
	nt_queryseq[i2++] = toupper(this->cdna);
      } else {
	*netchars -= 1;
      }
      if (this->aa_e != ' ') {
	aa_queryseq[k2++] = this->aa_e;
      }
    }
  }

  if (i1 >= MAXMUT || k1 >= MAXMUT || i2 >= MAXMUT || k2 >= MAXMUT) {
    i1 = i2 = k1 = k2 = 0;
  }
  nt_genomicseg[i1] = '\0';
  aa_genomicseg[k1] = '\0';
  nt_queryseq[i2] = '\0';
  aa_queryseq[k2] = '\0';

  *strlen_g = k1;
  *strlen_c = k2;

  return;
}


static void
print_mutation (bool *printp, int aapos, int strlen_g, int strlen_c, int refquerypos,
		char *aa_genomicseg, char *aa_queryseq, char *nt_genomicseg, char *nt_queryseq) {
  bool print_nt_p = true, print_refquerypos_p = true;

  if (strlen_g > strlen_c) {
    if (*printp == true) printf(", "); else *printp = true;
    if (aa_genomicseg[0] == aa_queryseq[0]) {
      printf("del%s%d%s ",&(aa_genomicseg[1]),aapos+1,&(aa_queryseq[1]));
      refquerypos += 3;
    } else {
      printf("del%s%d%s ",aa_genomicseg,aapos,aa_queryseq);
    }
  } else if (strlen_g < strlen_c) {
    if (*printp == true) printf(", "); else *printp = true;
    if (strlen_c - strlen_g > 4) {
      printf("ins%d+%daa ",aapos,strlen_c-strlen_g);
      print_nt_p = false;
    } else if (aa_genomicseg[0] == aa_queryseq[0]) {
      printf("ins%s%d%s ",&(aa_genomicseg[1]),aapos,&(aa_queryseq[1]));
    } else {
      printf("ins%s%d%s ",aa_genomicseg,aapos,aa_queryseq);
    }
  } else if (aa_genomicseg[0] == 'X' || aa_queryseq[0] == 'X') {
    print_nt_p = false;
    print_refquerypos_p = false;
  } else {
    if (*printp == true) printf(", "); else *printp = true;
    printf("%s%d%s ",aa_genomicseg,aapos,aa_queryseq);
  }
#if 0
  if (print_nt_p == true) {
    printf("(%s>%s) ",nt_genomicseg,nt_queryseq);
  }
#endif
  if (print_refquerypos_p == true) {
#ifdef PMAP
    printf("[%d]",refquerypos+2);
#else    
    printf("[%d]",refquerypos);
#endif
  }
  
  return;
}

static void
print_large_deletion (bool *printp, int lastaapos, int nextaapos, int refquerypos) {

  if (*printp == true) printf(", "); else *printp = true;
  printf("del%d-%daa ",lastaapos+1,nextaapos-lastaapos-1);
  printf("[%d]",refquerypos+3);
  
  return;
}


void
Translation_compare (struct Pair_T *pairs, int npairs, struct Pair_T *refpairs, int nrefpairs,
		     int cdna_direction, int relaastart, int relaaend, int maxmutations) {
  int i, j;
  int aapos, strlen_g, strlen_c, netchars;
  bool printp = false;
  char nt_genomicseg[MAXMUT], aa_genomicseg[MAXMUT], nt_queryseq[MAXMUT], aa_queryseq[MAXMUT];

  printf("    Amino acid changes: ");

  if (relaastart < relaaend) {
    if ((aapos = pairs[i=0].aapos) == 0) {
      i = next_aapos_fwd(pairs,0,npairs,0);
    }
    while (i < npairs) {
      aapos = pairs[i].aapos;
      j = next_aapos_fwd(pairs,i,npairs,aapos);
      if (pairs[i].aa_g != ' ' && pairs[i].aa_e != ' ') {
	/* not a frameshift */
	fill_aa_fwd(&strlen_g,&strlen_c,&netchars,aa_genomicseg,aa_queryseq,nt_genomicseg,nt_queryseq,
		    &(pairs[i]),&(pairs[j-1]));
	if (strcmp(aa_genomicseg,aa_queryseq) && strcmp(nt_genomicseg,nt_queryseq)) {
	  if (netchars % 3 == 0 || netchars > 12 || netchars < -12) {
	    print_mutation(&printp,aapos,strlen_g,strlen_c,pairs[i].refquerypos,
			   aa_genomicseg,aa_queryseq,nt_genomicseg,nt_queryseq);
	  }
	} else if (j < npairs && pairs[j].aapos - aapos > 4) {
	  print_large_deletion(&printp,aapos,pairs[j].aapos,pairs[i].refquerypos);
	}
      }
      i = j;
    }

  } else {
    if ((aapos = pairs[i=npairs-1].aapos) == 0) {
      i = next_aapos_rev(pairs,0,npairs,0);
    }
    while (i >= 0) {
      aapos = pairs[i].aapos;
      j = next_aapos_rev(pairs,i,npairs,aapos);
      if (pairs[i].aa_g != ' ' && pairs[i].aa_e != ' ') {
	/* not a frameshift */
	fill_aa_rev(&strlen_g,&strlen_c,&netchars,aa_genomicseg,aa_queryseq,nt_genomicseg,nt_queryseq,
		    &(pairs[i]),&(pairs[j+1]));
	if (strcmp(aa_genomicseg,aa_queryseq) && strcmp(nt_genomicseg,nt_queryseq)) {
	  if (netchars % 3 == 0 || netchars > 12 || netchars < -12) {
	    print_mutation(&printp,aapos,strlen_g,strlen_c,pairs[i].refquerypos,
			   aa_genomicseg,aa_queryseq,nt_genomicseg,nt_queryseq);
	  }
	} else if (j >= 0 && pairs[j].aapos - aapos > 4) {
	  print_large_deletion(&printp,aapos,pairs[j].aapos,pairs[i].refquerypos);
	}
      }
      i = j;
    }
  }

  printf("\n");

  return;
}


