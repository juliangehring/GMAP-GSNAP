static char rcsid[] = "$Id: smooth.c,v 1.23 2006/04/04 23:21:53 twu Exp $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "smooth.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>		/* For pow() */
#include "mem.h"
#include "pair.h"
#include "pairdef.h"
#include "intlist.h"

/* Will mark internal exons smaller than this for dual genome gap */
#define SHORTEXONLEN_MID 20
#define ZERONETGAP 9

#ifdef DEBUG
#define debug(x) x
#else 
#define debug(x)
#endif


typedef enum {KEEP, DELETE, MARK} Exonstatus_T;


static int *
get_exonlengths (int *nexons, List_T pairs) {
  int *exonlengths;
  Intlist_T list = NULL;
  Pair_T firstpair, pair;
  int querypos, i;

  /* Initialize */
  firstpair = List_head(pairs);

  /* Find short first and middle exons */
  while (pairs != NULL) {
    pair = List_head(pairs);
    if (pair->gapp == true) {
      list = Intlist_push(list,querypos-firstpair->querypos+1);
      pairs = List_next(pairs);
      firstpair = List_head(pairs);
    } else {
      querypos = pair->querypos;
    }
    pairs = List_next(pairs);
  }

  list = Intlist_push(list,querypos-firstpair->querypos+1);
  list = Intlist_reverse(list);
  exonlengths = Intlist_to_array(&(*nexons),list);
  Intlist_free(&list);

  return exonlengths;
}

static int *
get_intronlengths (int *nintrons, List_T pairs) {
  int *intronlengths, i;
  Intlist_T list = NULL;
  Pair_T pair;

  while (pairs != NULL) {
    pair = List_head(pairs);
    if (pair->gapp == true) {
      list = Intlist_push(list,pair->genomejump - pair->queryjump);
    }
    pairs = List_next(pairs);
  }

  list = Intlist_reverse(list);
  intronlengths = Intlist_to_array(&(*nintrons),list);
  Intlist_free(&list);

  return intronlengths;
}



static bool
short_end_exon (int exonlen, int intronlen) {
  double prob;

  if (exonlen < 0) {
    /*
    fprintf(stderr,"Bug in short_end_exon.  exonlen = %d.  Please report to twu@gene.com\n",exonlen);
    abort();
    */
    return false;
  } else {
    prob = 1 - pow(1.0-pow(4.0,(double) -exonlen),(double) intronlen);
    debug(printf("Probability of exon of length %d with intron of length %d is %g\n",
		 exonlen,intronlen,prob));
    if (prob < 0.001) {
      return false;
    } else {
      return true;
    }
  }
}


static void
zero_net_gap (int *starti, int *startj, int i, int j, int *intronlengths) {
  int netgap, bestnetgap = 1000000;
  int k, l, adji;

  if (i == 0) {
    adji = 0;
  } else {
    adji = i-1;
  }

  for (k = adji; k < j; k++) {
    netgap = intronlengths[k];
    for (l = k+1; l < j; l++) {
      netgap += intronlengths[l];
      debug(printf("zero_net_gap: netgap from %d to %d is %d\n",k,l,netgap));
      if (abs(netgap) < bestnetgap) {
	bestnetgap = abs(netgap);
	*starti = k+1;
	*startj = l;
      }
    }
  }

  debug(printf("zero_net_gap: best result is %d from %d to %d\n",bestnetgap,*starti,*startj));

  if (bestnetgap > ZERONETGAP) {
    debug(printf("zero_net_gap: not recommending any deletions\n"));
    *starti = *startj = -1;
  }

  return;
}

static int
longest_middle_exon (int i, int j, int *exonlengths, int indexsize) {
  int length, bestlength = indexsize;
  int besti = -1, k;

  for (k = i; k < j; k++) {
    debug(printf("longest_middle_exon: looking at exon %d with length %d\n",k,exonlengths[k]));
    if (exonlengths[k] > bestlength) {
      bestlength = exonlengths[k];
      besti = k;
    }
  }

  debug(printf("longest_middle_exon: best result is exon %d\n",besti));

  return besti;
}



static int *
get_exonstatus (int *nexons, int *nshortexons, int *ndeleteexons, List_T pairs, int indexsize) {
  int *exonstatus, *exonlengths, *intronlengths, nintrons, i, j, k;
  int starti, startj, longi;

  exonlengths = get_exonlengths(&(*nexons),pairs);
  intronlengths = get_intronlengths(&nintrons,pairs);
  /* Debugging */
  debug(
	for (i = 0; i < *nexons-1; i++) {
	  printf("Exon %d of length %d\n",i,exonlengths[i]);
	  printf("Intron %d of length %d\n",i,intronlengths[i]);
	}
	printf("Exon %d of length %d\n",*nexons-1,exonlengths[*nexons-1]);
	);

  exonstatus = (int *) CALLOC(*nexons,sizeof(int));
  for (i = 0; i < *nexons; i++) {
    exonstatus[i] = KEEP;
  }
  
  /* Mark middle exons */
  for (i = 0; i < *nexons; i++) {
    if (exonlengths[i] < SHORTEXONLEN_MID) {
      exonstatus[i] = MARK;
    }
  }

  /* Trim 5' end */
  i = 0;
  j = 0;
  while (i < *nexons && j < nintrons && short_end_exon(exonlengths[i],intronlengths[j])) {
    exonstatus[i] = DELETE;
    i++;
    j++;
  }
  if (i < *nexons) {
    /* Keep next exon regardless of length */
    exonstatus[i] = KEEP;
  }

  /* Trim 3' end */
  i = (*nexons) - 1;
  j = nintrons - 1;
  while (i >= 0 && j >= 0 && short_end_exon(exonlengths[i],intronlengths[j])) {
    exonstatus[i] = DELETE;
    --i;
    --j;
  }
  if (i >= 0) {
    /* Keep next exon regardless of length */
    exonstatus[i] = KEEP;
  }

  /* Remove short exons */
  i = 0;
  while (i < *nexons) {
    if (exonstatus[i] == MARK) {
      j = i;
      while (j < *nexons && exonstatus[j] == MARK) {
	j++;
      }
      debug(printf("Calling zero_net_gap with %d exons\n",j-i));
      zero_net_gap(&starti,&startj,i,j,intronlengths);
      if (starti >= 0) {
	for (j = starti; j <= startj; j++) {
	  exonstatus[j] = DELETE;
	}
      } else if (j - i == 1) {
	exonstatus[i] = MARK;
      } else {
	for (k = i; k < j; k++) {
	  exonstatus[k] = DELETE;
	}
	if ((longi = longest_middle_exon(i,j,exonlengths,indexsize)) >= 0) {
	  exonstatus[longi] = MARK;
	}
      }
      i = j;
    } else {
      i++;
    }
  }

  /* Count */
  *nshortexons = 0;
  *ndeleteexons = 0;
  for (i = 0; i < *nexons; i++) {
    if (exonstatus[i] == MARK) {
      (*nshortexons)++;
    } else if (exonstatus[i] == DELETE) {
      (*ndeleteexons)++;
    }
  }

  /* Debugging */
  debug(
	for (i = 0; i < *nexons; i++) {
	  if (exonstatus[i] == KEEP) {
	    printf("Long exon %d of length %d => keep\n",i,exonlengths[i]);
	  } else if (exonstatus[i] == MARK) {
	    printf("Short exon %d of length %d => mark\n",i,exonlengths[i]);
	  } else if (exonstatus[i] == DELETE) {
	    printf("Short exon %d of length %d => delete\n",i,exonlengths[i]);
	  } else {
	    abort();
	  }
	}
	);

  if (nintrons > 0) {
    FREE(intronlengths);
  }
  FREE(exonlengths);

  return exonstatus;
}


/* Assumes pairs are from 1..querylength.  Reverses the pairs to be querylength..1 */
List_T
Smooth_pairs (int *nshortexons, int *ndeleteexons, List_T pairs, Pairpool_T pairpool,
	      int indexsize) {
  List_T newpairs = NULL, p;
  int *exonstatus, nexons, i;
  Pair_T pair;
  int prevstatus, currstatus;

  exonstatus = get_exonstatus(&nexons,&(*nshortexons),&(*ndeleteexons),pairs,indexsize);

  if (*nshortexons == 0 && *ndeleteexons == 0) {
    FREE(exonstatus);
    return pairs;
  }

  i = 0;
  currstatus = exonstatus[i];
  while (pairs != NULL) {
    pair = List_head(pairs);
    if (pair->gapp == true) {
      prevstatus = currstatus;
      currstatus = exonstatus[++i];
      debug(printf("Gap observed\n"));
      if (prevstatus != DELETE && currstatus != DELETE) {
	newpairs = Pairpool_push_existing(newpairs,pairpool,pair);
      }

    } else if (currstatus == KEEP) {
      /* debug(printf("Marking position %d as keep\n",pair->querypos)); */
      newpairs = Pairpool_push_existing(newpairs,pairpool,pair);
    } else if (currstatus == MARK) {
      debug(printf("Marking position %d as short\n",pair->querypos));
      pair->shortexonp = true;
      newpairs = Pairpool_push_existing(newpairs,pairpool,pair);
    } else {
      debug(printf("Marking position %d for deletion\n",pair->querypos));
    }
    pairs = List_next(pairs);
  }
	
  debug(printf("\n"));
  FREE(exonstatus);

  return List_reverse(newpairs);
}


void
Smooth_reset (List_T pairs) {
  List_T p;
  Pair_T pair;

  for (p = pairs; p != NULL; p = List_next(p)) {
    pair = List_head(p);
    pair->shortexonp = false;
  }
  return;
}


