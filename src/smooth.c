static char rcsid[] = "$Id: smooth.c 90983 2013-04-01 19:42:40Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "smooth.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>		/* For pow(), log() */
#include "bool.h"
#include "mem.h"
#include "except.h"
#include "comp.h"
#include "pair.h"
#include "pairdef.h"
#include "intlist.h"


/* Will examine internal exons smaller than this prior to single gaps */
#define ZERONETGAP 9
#define SHORTEXONLEN_NETGAP 15

/* Will delete/mark internal exons smaller than this for dual genome gap */
#define DELETE_THRESHOLD 0.1
#define MARK_THRESHOLD 1e-7

/* Will delete internal exons smaller than this at ends */
#define SHORTEXONLEN_END 10

#ifdef GSNAP
/* Allow more intron predictions in ends of short reads */
#define SHORTEXONPROB_END 0.10
#else
#define SHORTEXONPROB_END 0.05
#endif

#define BIGGAP 9


#ifdef DEBUG
#define debug(x) x
#else 
#define debug(x)
#endif


typedef enum {KEEP, DELETE, MARK} Exonstatus_T;


static bool
big_gap_p (Pair_T pair, bool bysizep) {

  /* When bysizep is true, single gaps should have been solved already
     (in pass 2), so any gap remaining must be a poor one (must have
     been removed because of a negative score), and we should consider
     all gaps, regardless of size.  When bysizep is false, we must be
     doing initial smoothing, before any single gaps have been solved,
     so we want to consider only big gaps.  This change was motivated
     by alignment of mouse HER2 against human genome, where a 58/58
     section was deleted in stage 3, but subsequent smoothing failed
     to recognize the short exon between this section and the next
     intron. */

  if (bysizep == true) {
    return true;
  } else if (pair->genomejump > pair->queryjump) {
    return (pair->genomejump - pair->queryjump) > BIGGAP ? true : false;
  } else {
    return false;
  }
}


static int *
get_exonlengths (int **exonmatches, int *nexons, List_T pairs, bool bysizep) {
  int *exonlengths;
  Intlist_T list = NULL, matchlist = NULL;
  Pair_T firstpair, pair;
  int querypos, nmatches = 0;
#ifdef DEBUG
  int i;
#endif

  firstpair = List_head(pairs);
  querypos = firstpair->querypos;

  while (pairs != NULL) {
    pair = List_head(pairs);
    if (pair->gapp == true && big_gap_p(pair,bysizep) == true) {
      list = Intlist_push(list,querypos-firstpair->querypos+1);
      matchlist = Intlist_push(matchlist,nmatches);

      pairs = List_next(pairs);
      if (pairs != NULL) {
	firstpair = List_head(pairs);
	querypos = firstpair->querypos;
	nmatches = 0;
      }
    } else {
      querypos = pair->querypos;
      if (pair->comp == MATCH_COMP || pair->comp == DYNPROG_MATCH_COMP) {
	nmatches++;
      }
      pairs = List_next(pairs);
    }
  }

  list = Intlist_push(list,querypos-firstpair->querypos+1);
  list = Intlist_reverse(list);
  matchlist = Intlist_push(matchlist,nmatches);
  matchlist = Intlist_reverse(matchlist);

  exonlengths = Intlist_to_array(&(*nexons),list);
  *exonmatches = Intlist_to_array(&(*nexons),matchlist);
  debug(
	printf("%d exons: ",*nexons);
	for (i = 0; i < *nexons; i++) {
	  printf("%d (%d matches),",exonlengths[i],(*exonmatches)[i]);
	}
	printf("\n");
	);
  Intlist_free(&list);
  Intlist_free(&matchlist);

  return exonlengths;
}

static int *
get_intronlengths (int *nintrons, List_T pairs, bool bysizep) {
  int *intronlengths, length;
  Intlist_T list = NULL;
  Pair_T pair;
#ifdef DEBUG
  int i;
#endif

  while (pairs != NULL) {
    pair = List_head(pairs);
    if (pair->gapp == true && big_gap_p(pair,bysizep) == true) {
#if 0
      if (pair->genomejump > pair->queryjump) {
	list = Intlist_push(list,pair->genomejump);
      } else {
	list = Intlist_push(list,pair->queryjump);
      }
#else
      if ((length = pair->genomejump - pair->queryjump) < 0) {
	/* cDNA insertion */
	list = Intlist_push(list,-length);
      } else {
	list = Intlist_push(list,length);
      }
#endif
    }
    pairs = List_next(pairs);
  }

  list = Intlist_reverse(list);
  intronlengths = Intlist_to_array(&(*nintrons),list);
  debug(
	printf("%d introns: ",*nintrons);
	for (i = 0; i < *nintrons; i++) {
	  printf("%d,",intronlengths[i]);
	}
	printf("\n");
	);
  Intlist_free(&list);

  return intronlengths;
}


static const Except_T length_error = {"Negative exon or intron length"};

static double
compute_prob (int exonlen, int intronlen, int indexsize) {
  double prob;

  if (exonlen < indexsize) {
    prob = 1.0;
  } else {
    prob = 1 - pow(1.0-pow(4.0,(double) -exonlen),(double) intronlen);
  }
  debug(printf("Probability of exon of length %d (indexsize %d) with intron of length %d is %g\n",
	       exonlen,indexsize,intronlen,prob));
  return prob;
}

static bool
short_exon_byprob (int exonlen, int intronlen, int indexsize, double prob_threshold) {
  double prob;

  prob = compute_prob(exonlen,intronlen,indexsize);
  if (prob < prob_threshold) {
    return false;
  } else {
    return true;
  }
}

static bool
short_exon_bylength (int exonlen, int length_threshold) {

  if (exonlen >= length_threshold) {
    return false;
  } else {
    return true;
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


static int *
find_internal_shorts_by_netgap (bool *deletep, int *exonmatches, int nexons,
				int *intronlengths) {
  int *exonstatus;
  int starti, startj, i, j;
  int exonlen;

  *deletep = false;
  exonstatus = (int *) CALLOC(nexons,sizeof(int));
  for (i = 0; i < nexons; i++) {
    exonstatus[i] = KEEP;
  }
  
  /* Mark short middle exons */
  for (i = 1; i < nexons - 1; i++) {
    exonlen = exonmatches[i];
    if (exonlen < SHORTEXONLEN_NETGAP) {
      exonstatus[i] = MARK;
    }
  }

  /* Find internal shorts */
  i = 0;
  while (i < nexons) {
    if (exonstatus[i] == MARK) {
      j = i;
      while (j < nexons && exonstatus[j] == MARK) {
	j++;
      }
      debug(printf("Calling zero_net_gap with %d exons\n",j-i));
      zero_net_gap(&starti,&startj,i,j,intronlengths);
      if (starti >= 0) {
	for (j = starti; j <= startj; j++) {
	  *deletep = true;
	  exonstatus[j] = DELETE;
	}
      } else if (j - i == 1) {
	exonstatus[i] = MARK;
      }
      i = j;
    } else {
      i++;
    }
  }

  return exonstatus;
}

static int *
find_internal_shorts_by_size (bool *shortp, bool *deletep, int *exonmatches, int nexons,
			      int *intronlengths, int stage2_indexsize) {
  int *exonstatus;
  int i;
  int exonlen, intronlen;
  double prob;

  *shortp = *deletep = false;
  exonstatus = (int *) CALLOC(nexons,sizeof(int));
  for (i = 0; i < nexons; i++) {
    exonstatus[i] = KEEP;
  }
  
  /* Mark short middle exons */
  for (i = 1; i < nexons - 1; i++) {
    exonlen = exonmatches[i];
    intronlen = intronlengths[i-1]+intronlengths[i];
    prob = compute_prob(exonlen+4,intronlen,stage2_indexsize); /* Hack: add 4 for possible canonical dinucleotides */
    if (prob > DELETE_THRESHOLD) {
      *deletep = true;
      exonstatus[i] = DELETE;
    } else if (prob > MARK_THRESHOLD) {
      *shortp = true;
      exonstatus[i] = MARK;
    }
  }

  return exonstatus;
}


/* For ends, we turn off the indexsize parameter */
static int *
find_end_shorts (bool *deletep, int *exonmatches, int nexons, int *intronlengths) {
  int *exonstatus, i;
  bool shortp;

  *deletep = false;
  exonstatus = (int *) CALLOC(nexons,sizeof(int));
  for (i = 0; i < nexons; i++) {
    exonstatus[i] = KEEP;
  }

  shortp = true;
  i = 0;
  while (i < nexons - 1 && shortp == true) {
    if (short_exon_bylength(exonmatches[i],SHORTEXONLEN_END) == true &&
	short_exon_byprob(exonmatches[i],intronlengths[i],/*indexsize*/0,SHORTEXONPROB_END) == true) {
      *deletep = true;
      exonstatus[i] = DELETE;
    } else {
      shortp = false;
    }
    i++;
  }
    
  shortp = true;
  i = nexons - 1;
  while (i > 0 && shortp == true) {
    if (short_exon_bylength(exonmatches[i],SHORTEXONLEN_END) == true &&
	short_exon_byprob(exonmatches[i],intronlengths[i-1],/*indexsize*/0,SHORTEXONPROB_END) == true) {
      *deletep = true;
      exonstatus[i] = DELETE;
    } else {
      shortp = false;
    }
    --i;
  }
    
  return exonstatus;
}


static List_T
delete_and_mark_exons (List_T pairs,
#ifdef WASTE
		       Pairpool_T pairpool,
#endif
		       int *exonstatus, bool markp, bool bysizep) {
  List_T newpairs = NULL, pairptr;
  Pair_T pair;
  int currstatus, prevstatus;
  int i;

  debug(
	for (i = 0; i < nexons; i++) {
	  if (exonstatus[i] == KEEP) {
	    printf("Long exon %d of %d matches => keep\n",i,exonmatches[i]);
	  } else if (exonstatus[i] == MARK) {
	    printf("Short exon %d of %d matches => mark\n",i,exonmatches[i]);
	  } else if (exonstatus[i] == DELETE) {
	    printf("Exon %d of %d matches => delete\n",i,exonmatches[i]);
	  } else {
	    abort();
	  }
	}
	);

  i = 0;
  currstatus = exonstatus[i];
  while (pairs != NULL) {
    pairptr = pairs;
    pairs = Pairpool_pop(pairs,&pair);
    if (pair->gapp == true && big_gap_p(pair,bysizep) == true) {
      prevstatus = currstatus;
      currstatus = exonstatus[++i];
      debug(printf("Gap observed\n"));
      if (prevstatus != DELETE && currstatus != DELETE) {
#ifdef WASTE
	newpairs = Pairpool_push_existing(newpairs,pairpool,pair);
#else
	newpairs = List_push_existing(newpairs,pairptr);
#endif
      }

    } else if (currstatus == KEEP) {
      /* debug(printf("Marking position %d as keep\n",pair->querypos)); */
#ifdef WASTE
      newpairs = Pairpool_push_existing(newpairs,pairpool,pair);
#else
      newpairs = List_push_existing(newpairs,pairptr);
#endif
    } else if (currstatus == MARK) {
      debug(printf("Marking position %d as short in pair %p\n",pair->querypos,pair));
      if (markp == true) {
	pair->shortexonp = true;
      }
#ifdef WASTE
      newpairs = Pairpool_push_existing(newpairs,pairpool,pair);
#else
      newpairs = List_push_existing(newpairs,pairptr);
#endif
    } else {
      debug(printf("Marking position %d for deletion\n",pair->querypos));
    }
  }

  /* Remove gaps at end */
  if (newpairs != NULL) {
    pair = List_head(newpairs);
  }
  while (newpairs != NULL && pair->gapp == true) {
    debug(printf("Popping gap at end\n"));
    newpairs = Pairpool_pop(newpairs,&pair);
  }

  /* Remove gaps at beginning */
  newpairs = List_reverse(newpairs);
  if (newpairs != NULL) {
    pair = List_head(newpairs);
  }
  while (newpairs != NULL && pair->gapp == true) {
    debug(printf("Popping gap at beginning\n"));
    newpairs = Pairpool_pop(newpairs,&pair);
  }

  debug(printf("Result of delete_and_mark_exons:\n"));
  debug(Pair_dump_list(newpairs,/*zerobasedp*/true));
  debug(printf("\n"));

  return newpairs;
}


#if 0
static List_T
mark_exons (List_T pairs,
#ifdef WASTE
	    Pairpool_T pairpool,
#endif
	    int *exonstatus) {
  List_T newpairs = NULL, pairptr;
  Pair_T pair;
  int currstatus, prevstatus;
  int i;

  debug(
	for (i = 0; i < nexons; i++) {
	  if (exonstatus[i] == KEEP) {
	    printf("Long exon %d of %d matches => keep\n",i,exonmatches[i]);
	  } else if (exonstatus[i] == MARK) {
	    printf("Short exon %d of %d matches => mark\n",i,exonmatches[i]);
	  } else if (exonstatus[i] == DELETE) {
	    printf("Exon %d of %d matches => delete\n",i,exonmatches[i]);
	  } else {
	    abort();
	  }
	}
	);

  i = 0;
  currstatus = exonstatus[i];
  while (pairs != NULL) {
    pairptr = pairs;
    pairs = Pairpool_pop(pairs,&pair);
    if (pair->gapp == true && big_gap_p(pair,/*bysizep*/true) == true) {
      prevstatus = currstatus;
      currstatus = exonstatus[++i];
      debug(printf("Gap observed\n"));
      if (prevstatus != DELETE && currstatus != DELETE) {
#ifdef WASTE
	newpairs = Pairpool_push_existing(newpairs,pairpool,pair);
#else
	newpairs = List_push_existing(newpairs,pairptr);
#endif
      }

    } else if (currstatus == KEEP) {
      /* debug(printf("Marking position %d as keep\n",pair->querypos)); */
#ifdef WASTE
      newpairs = Pairpool_push_existing(newpairs,pairpool,pair);
#else
      newpairs = List_push_existing(newpairs,pairptr);
#endif
    } else if (currstatus == MARK) {
      debug(printf("Marking position %d as short in pair %p\n",pair->querypos,pair));
      pair->shortexonp = true;

#ifdef WASTE
      newpairs = Pairpool_push_existing(newpairs,pairpool,pair);
#else
      newpairs = List_push_existing(newpairs,pairptr);
#endif
    } else {
      /* Normally would delete */
#ifdef WASTE
      newpairs = Pairpool_push_existing(newpairs,pairpool,pair);
#else
      newpairs = List_push_existing(newpairs,pairptr);
#endif
    }
  }

  debug(printf("Result of mark_exons:\n"));
  debug(Pair_dump_list(newpairs,/*zerobasedp*/true));
  debug(printf("\n"));

  return List_reverse(newpairs);
}
#endif


static void
smooth_reset (List_T pairs) {
  List_T p;
  Pair_T pair;

  for (p = pairs; p != NULL; p = List_next(p)) {
    pair = List_head(p);
    pair->shortexonp = false;
  }
  return;
}



/* Assumes pairs are from 1..querylength.  Reverses the pairs to be querylength..1 */
List_T
Smooth_pairs_by_netgap (bool *deletep, List_T pairs, Pairpool_T pairpool) {
  int *exonstatus;
  int *exonlengths, *exonmatches, *intronlengths, nexons, nintrons;
#ifdef DEBUG
  int i;
#endif

  *deletep = false;
  smooth_reset(pairs);
  if (pairs != NULL) {
    /* Remove internal shorts */
    exonlengths = get_exonlengths(&exonmatches,&nexons,pairs,/*bysizep*/false);
    intronlengths = get_intronlengths(&nintrons,pairs,/*bysizep*/false);

    debug(
	  printf("Beginning of smoothing.  Initial structure:\n");
	  for (i = 0; i < nexons-1; i++) {
	    printf("Exon %d of length %d (%d matches)\n",i,exonlengths[i],exonmatches[i]);
	    printf("Intron %d of length %d\n",i,intronlengths[i]);
	  }
	  printf("Exon %d of length %d (%d matches)\n",nexons-1,exonlengths[nexons-1],exonmatches[nexons-1]);
	  );

    debug(printf("\nFind internal shorts\n"));
    exonstatus = find_internal_shorts_by_netgap(&(*deletep),exonmatches,nexons,intronlengths);
    debug(printf("\nRemove internal shorts\n"));
    if (*deletep == true) {
      pairs = delete_and_mark_exons(pairs,
#ifdef WASTE
				    pairpool,
#endif
				    exonstatus,/*markp*/false,/*bysizep*/false);
    }
    
#if 0
    /* This is not correct */
    debug(
	  printf("After removing internal shorts:\n");
	  for (i = 0; i < nexons-1; i++) {
	    printf("Exon %d of length %d (%d matches)\n",i,exonlengths[i],exonmatches[i]);
	    printf("Intron %d of length %d\n",i,intronlengths[i]);
	  }
	  printf("Exon %d of length %d\n",nexons-1,exonlengths[nexons-1]);
	  );
#endif

    FREE(exonstatus);
    FREE(intronlengths);
    FREE(exonmatches);
    FREE(exonlengths);

    debug(printf("Ending of smoothing\n\n"));
  }

  return pairs;
}


/* Assumes pairs are from 1..querylength.  Reverses the pairs to be querylength..1 */
List_T
Smooth_pairs_by_size (bool *shortp, bool *deletep, List_T pairs, Pairpool_T pairpool, int stage2_indexsize) {
  int *exonstatus;
  int *exonlengths, *exonmatches, *intronlengths, nexons, nintrons;
  bool delete1p, delete2p;
#ifdef DEBUG
  int i;
#endif

  *shortp = *deletep = false;
  smooth_reset(pairs);
  if (pairs != NULL) {
    /* Trim ends */
    exonlengths = get_exonlengths(&exonmatches,&nexons,pairs,/*bysizep*/true);
    intronlengths = get_intronlengths(&nintrons,pairs,/*bysizep*/true);

    debug(printf("\nFind end shorts\n"));
    exonstatus = find_end_shorts(&delete1p,exonmatches,nexons,intronlengths);
    if (delete1p == true) {
      *deletep = true;
      pairs = delete_and_mark_exons(pairs,
#ifdef WASTE
				    pairpool,
#endif
				    exonstatus,/*markp*/false,/*bysizep*/true);
    }

    FREE(exonstatus);

    FREE(intronlengths);
    FREE(exonmatches);
    FREE(exonlengths);
  }

  if (pairs != NULL) {
    /* Remove internal shorts */
    exonlengths = get_exonlengths(&exonmatches,&nexons,pairs,/*bysizep*/true);
    intronlengths = get_intronlengths(&nintrons,pairs,/*bysizep*/true);

    debug(
	  printf("Beginning of smoothing.  Initial structure:\n");
	  for (i = 0; i < nexons-1; i++) {
	    printf("Exon %d of length %d (%d matches)\n",i,exonlengths[i],exonmatches[i]);
	    printf("Intron %d of length %d\n",i,intronlengths[i]);
	  }
	  printf("Exon %d of length %d (%d matches)\n",nexons-1,exonlengths[nexons-1],exonmatches[nexons-1]);
	  );

    debug(printf("\nFind internal shorts\n"));
    exonstatus = find_internal_shorts_by_size(&(*shortp),&delete2p,exonmatches,nexons,intronlengths,stage2_indexsize);
    debug(printf("\nRemove internal shorts\n"));
    if (delete2p == true) {
      *deletep = true;
    }
    if (delete2p == true || *shortp == true) {
      pairs = delete_and_mark_exons(pairs,
#ifdef WASTE
				    pairpool,
#endif
				    exonstatus,/*markp*/true,/*bysizep*/true);
    }

    debug(
	  printf("After removing internal shorts:\n");
	  for (i = 0; i < nexons-1; i++) {
	    printf("Exon %d of length %d (%d matches)\n",i,exonlengths[i],exonmatches[i]);
	    printf("Intron %d of length %d\n",i,intronlengths[i]);
	  }
	  printf("Exon %d of length %d\n",nexons-1,exonlengths[nexons-1]);
	  );

    FREE(exonstatus);
    FREE(intronlengths);
    FREE(exonmatches);
    FREE(exonlengths);

  }

  debug(printf("Ending of smoothing\n\n"));

  return pairs;
}


#if 0
List_T
Smooth_mark_short_exons (List_T pairs, Pairpool_T pairpool, int stage2_indexsize) {
  int *exonstatus;
  int *exonlengths, *exonmatches, *intronlengths, nexons, nintrons;
  bool shortp;
  bool delete1p, delete2p;
#ifdef DEBUG
  int i;
#endif

  shortp = false;
  smooth_reset(pairs);
  if (pairs != NULL) {
    /* Trim ends */
    exonlengths = get_exonlengths(&exonmatches,&nexons,pairs,/*bysizep*/true);
    intronlengths = get_intronlengths(&nintrons,pairs,/*bysizep*/true);

    debug(printf("\nFind end shorts\n"));
    exonstatus = find_end_shorts(&delete1p,exonmatches,nexons,intronlengths);
    pairs = mark_exons(pairs,
#ifdef WASTE
		       pairpool,
#endif
		       exonstatus);

    FREE(exonstatus);

    FREE(intronlengths);
    FREE(exonmatches);
    FREE(exonlengths);
  }

  if (pairs != NULL) {
    /* Find internal shorts */
    exonlengths = get_exonlengths(&exonmatches,&nexons,pairs,/*bysizep*/true);
    intronlengths = get_intronlengths(&nintrons,pairs,/*bysizep*/true);

    debug(
	  printf("Beginning of smoothing.  Initial structure:\n");
	  for (i = 0; i < nexons-1; i++) {
	    printf("Exon %d of length %d (%d matches)\n",i,exonlengths[i],exonmatches[i]);
	    printf("Intron %d of length %d\n",i,intronlengths[i]);
	  }
	  printf("Exon %d of length %d (%d matches)\n",nexons-1,exonlengths[nexons-1],exonmatches[nexons-1]);
	  );

    debug(printf("\nFind internal shorts\n"));
    exonstatus = find_internal_shorts_by_size(&shortp,&delete2p,exonmatches,nexons,intronlengths,stage2_indexsize);
    debug(printf("\nMark internal shorts\n"));
    pairs = mark_exons(pairs,
#ifdef WASTE
		       pairpool,
#endif
		       exonstatus);

    debug(
	  printf("After marking internal shorts:\n");
	  for (i = 0; i < nexons-1; i++) {
	    printf("Exon %d of length %d (%d matches)\n",i,exonlengths[i],exonmatches[i]);
	    printf("Intron %d of length %d\n",i,intronlengths[i]);
	  }
	  printf("Exon %d of length %d\n",nexons-1,exonlengths[nexons-1]);
	  );

    FREE(exonstatus);
    FREE(intronlengths);
    FREE(exonmatches);
    FREE(exonlengths);

  }

  debug(printf("End of Smooth_mark_short_exons\n\n"));

  return pairs;
}
#endif


