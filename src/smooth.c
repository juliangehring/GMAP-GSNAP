static char rcsid[] = "$Id: smooth.c,v 1.18 2005/07/08 07:58:35 twu Exp $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "smooth.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>		/* For pow() */
#include "bool.h"
#include "mem.h"
#include "pair.h"
#include "pairdef.h"
#include "intlist.h"

/* Will mark internal exons smaller than this for dual genome gap */
#define SHORTEXONLEN_MID 80


#ifdef DEBUG
#define debug(x) x
#else 
#define debug(x)
#endif


typedef enum {KEEP, DELETE, MARK} Exonstatus_T;


static Intlist_T
get_exonlengths (List_T pairs) {
  Intlist_T exonlengths = NULL;
  List_T p;
  Pair_T firstpair, pair;
  int exonlen = 0;
  int lastquerypos, lastgenomepos, queryjump, genomejump;

  /* Initialize */
  firstpair = List_head(pairs);
  lastgenomepos = firstpair->genomepos - 1;
  lastquerypos = firstpair->querypos - 1;

  /* Find short first and middle exons */
  for (p = pairs; p != NULL; p = List_next(p)) {
    pair = List_head(p);
    queryjump = pair->querypos - lastquerypos - 1;
    genomejump = pair->genomepos - lastgenomepos - 1;

    if (queryjump <= 0 && genomejump <= 0) {
      exonlen += pair->querypos - lastquerypos;
    } else {
      exonlengths = Intlist_push(exonlengths,exonlen);
      exonlen = 1;
    }
    lastquerypos = pair->querypos;
    lastgenomepos = pair->genomepos;
  }

  exonlengths = Intlist_push(exonlengths,exonlen);
  exonlengths = Intlist_reverse(exonlengths);

  /* Debugging */
  debug(Pair_dump_list(pairs,true));
  debug(
	l = exonlengths;
	while (l != NULL) {
	  printf("Exon of length %d\n",Intlist_head(l));
	  l = Intlist_next(l);
	}
	);

  return exonlengths;
}

static Intlist_T
get_intronlengths (List_T pairs) {
  Intlist_T intronlengths = NULL;
  List_T p;
  Pair_T firstpair, pair;
  int lastquerypos, lastgenomepos, queryjump, genomejump;

  /* Initialize */
  firstpair = List_head(pairs);
  lastgenomepos = firstpair->genomepos - 1;
  lastquerypos = firstpair->querypos - 1;

  for (p = pairs; p != NULL; p = List_next(p)) {
    pair = List_head(p);
    queryjump = pair->querypos - lastquerypos - 1;
    genomejump = pair->genomepos - lastgenomepos - 1;

    if (queryjump <= 0 && genomejump <= 0) {
    } else {
      intronlengths = Intlist_push(intronlengths,pair->genomepos - lastgenomepos - 1);
    }
    lastquerypos = pair->querypos;
    lastgenomepos = pair->genomepos;
  }

  intronlengths = Intlist_reverse(intronlengths);

  /* Debugging */
  debug(
	l = intronlengths;
	while (l != NULL) {
	  printf("Intron of length %d\n",Intlist_head(l));
	  l = Intlist_next(l);
	}
	);

  return intronlengths;
}


static bool
short_end_exon (int exonlen, int intronlen) {
  double prob;

  if (exonlen < 0) {
    fprintf(stderr,"Bug in short_end_exon.  exonlen = %d.  Please report to twu@gene.com\n",exonlen);
    abort();
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


static Intlist_T
get_exonstatus (int *nshortexons, bool *deletep, List_T pairs) {
  Intlist_T exonstatus = NULL, exonlengths, intronlengths, q, l, k;

  *deletep = false;

  exonlengths = get_exonlengths(pairs);
  intronlengths = get_intronlengths(pairs);
  for (l = exonlengths; l != NULL; l = Intlist_next(l)) {
    exonstatus = Intlist_push(exonstatus,KEEP);
  }
  exonstatus = Intlist_reverse(exonstatus);
  
  /* Mark middle exons */
  q = exonstatus;
  l = exonlengths;
  while (l != NULL) {
    if (Intlist_head(l) < SHORTEXONLEN_MID) {
      Intlist_head_set(q,MARK);
    }
    q = Intlist_next(q);
    l = Intlist_next(l);
  }

  /* Trim 5' end */
  q = exonstatus;
  l = exonlengths;
  k = intronlengths;
  while (l != NULL && k != NULL && short_end_exon(Intlist_head(l),Intlist_head(k))) {
    Intlist_head_set(q,DELETE);
    *deletep = true;
    q = Intlist_next(q);
    l = Intlist_next(l);
    k = Intlist_next(k);
  }
  if (q != NULL) {
    /* Keep next exon regardless of length */
    Intlist_head_set(q,KEEP);
  }

  /* Trim 3' end */
  exonstatus = Intlist_reverse(exonstatus);
  exonlengths = Intlist_reverse(exonlengths);
  intronlengths = Intlist_reverse(intronlengths);
  
  q = exonstatus;
  l = exonlengths;
  k = intronlengths;
  while (l != NULL && k != NULL && short_end_exon(Intlist_head(l),Intlist_head(k))) {
    Intlist_head_set(q,DELETE);
    *deletep = true;
    q = Intlist_next(q);
    l = Intlist_next(l);
    k = Intlist_next(k);
  }
  if (q != NULL) {
    /* Keep next exon regardless of length */
    Intlist_head_set(q,KEEP);
  }

  *nshortexons = 0;
  for (q = exonstatus; q != NULL; q = Intlist_next(q)) {
    if (Intlist_head(q) == MARK) {
      (*nshortexons)++;
    }
  }

  exonstatus = Intlist_reverse(exonstatus);

  intronlengths = Intlist_reverse(intronlengths);
  exonlengths = Intlist_reverse(exonlengths);

  /* Debugging */
  debug(
	q = exonstatus;
	l = exonlengths;
	while (l != NULL) {
	  if (Intlist_head(q) == KEEP) {
	    printf("Long exon of length %d => keep\n",Intlist_head(l));
	  } else if (Intlist_head(q) == MARK) {
	    printf("Short middle exon of length %d => mark\n",Intlist_head(l));
	  } else if (Intlist_head(q) == DELETE) {
	    printf("Short end exon of length %d => delete\n",Intlist_head(l));
	  } else {
	    abort();
	  }
	  q = Intlist_next(q);
	  l = Intlist_next(l);
	}
	);

  Intlist_free(&intronlengths);
  Intlist_free(&exonlengths);
  return exonstatus;
}


/* Assumes pairs are from 1..querylength.  Reverses the pairs to be querylength..1 */
List_T
Smooth_pairs (int *nshortexons, List_T pairs, Pairpool_T pairpool) {
  List_T newpairs = NULL, p;
  Intlist_T exonstatus, q;
  Pair_T pair, firstpair;
  int currstatus;
  int lastquerypos, lastgenomepos, queryjump, genomejump;
  bool deletep;

  exonstatus = get_exonstatus(&(*nshortexons),&deletep,pairs);

  if (*nshortexons == 0 && deletep == false) {
    Intlist_free(&exonstatus);
    return pairs;
  }

  /* Initialize again */
  firstpair = List_head(pairs);
  lastgenomepos = firstpair->genomepos - 1;
  lastquerypos = firstpair->querypos - 1;

  /* Mark exons */
  q = exonstatus;
  currstatus = Intlist_head(q);
  for (p = pairs; p != NULL; p = List_next(p)) {
    pair = List_head(p);
    queryjump = pair->querypos - lastquerypos - 1;
    genomejump = pair->genomepos - lastgenomepos - 1;

    if (queryjump <= 0 && genomejump <= 0) {
      if (currstatus == KEEP) {
	newpairs = Pairpool_push_existing(newpairs,pairpool,pair);
      } else if (currstatus == MARK) {
	debug(printf("Marking position %d as short\n",pair->querypos));
	pair->shortexonp = true;
	newpairs = Pairpool_push_existing(newpairs,pairpool,pair);
      } else {
	debug(printf("Marking position %d for deletion\n",pair->querypos));
      }
    } else {
      q = Intlist_next(q);
      currstatus = Intlist_head(q);
      if (currstatus == KEEP) {
	newpairs = Pairpool_push_existing(newpairs,pairpool,pair);
      } else if (currstatus == MARK) {
	debug(printf("Marking position %d as short\n",pair->querypos));
	pair->shortexonp = true;
	newpairs = Pairpool_push_existing(newpairs,pairpool,pair);
      } else {
	debug(printf("Marking position %d for deletion\n",pair->querypos));
      }
    }
    lastquerypos = pair->querypos;
    lastgenomepos = pair->genomepos;
  }
  
  debug(printf("\n"));
  Intlist_free(&exonstatus);

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
