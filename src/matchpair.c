static char rcsid[] = "$Id: matchpair.c,v 1.24 2005/03/01 20:22:07 twu Exp $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "matchpair.h"
#include <stdio.h>
#include <stdlib.h>
#include "mem.h"
#include "list.h"
#include "listdef.h"
#include "matchdef.h"

#define MAXCANDIDATES 10

#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif


/* Debugging of path finding */
#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif

/* Summary of clustering */
#ifdef DEBUG2
#define debug2(x) x
#else
#define debug2(x)
#endif


#define T Matchpair_T
struct T {
  int matchsize;
  Match_T bound5;
  Match_T bound3;
  unsigned int low;		/* low position */
  unsigned int high;		/* high position */
  bool watsonp;
  int clustersize;		/* for cluster algorithm */
  Matchpairend_T matchpairend;
};

Match_T
Matchpair_bound5 (T this) {
  return this->bound5;
}

Match_T
Matchpair_bound3 (T this) {
  return this->bound3;
}

bool
Matchpair_watsonp (T this) {
  return this->watsonp;
}

int
Matchpair_clustersize (T this) {
  return this->clustersize;
}

Matchpairend_T
Matchpair_matchpairend (T this) {
  return this->matchpairend;
}

T
Matchpair_new (Match_T bound5, Match_T bound3, int matchsize, int clustersize,
	       Matchpairend_T matchpairend) {
  T new = (T) MALLOC(sizeof(*new));
  unsigned int position_x1, position_x2, position_y1, position_y2, temp;

  new->matchsize = matchsize;
  new->bound5 = bound5;
  new->bound3 = bound3;
  new->watsonp = Match_forwardp(bound5);
  new->clustersize = clustersize;
  new->matchpairend = matchpairend;

  new->low = Match_position(bound5);
  new->high = Match_position(bound3);
  if (new->high < new->low) {
    temp = new->high;
    new->high = new->low;
    new->low = temp;
  }
  
  return new;
}

void
Matchpair_free (T *old) {
  if (*old) {
    FREE(*old);
  }
  return;
}


/* Not intended for qsort.  Returns 0 when not comparable. */
static int
Matchpair_dominate (const void *a, const void *b) {
  T x = * (T *) a;
  T y = * (T *) b;
  
  if (x->watsonp != y->watsonp) {
    return 0;			/* Different strands */

  } else if (y->low > x->high || x->low > y->high) {
    return 0;			/* No overlap */
  } else if (x->low < y->low) {
    if (x->high < y->high) {
      return 0;			/* Crossing */
    } else {
      return -1;
    }
  } else if (y->low < x->low) {
    if (y->high < x->high) {
      return 0;			/* Crossing */
    } else {
      return 1;
    }
  } else {
    if (x->high < y->high) {
      return 1;
    } else if (y->high < x->high) {
      return -1;
    } else {
      return 0;			/* Equal */
    }
  }
}

static bool
Matchpair_equal (const void *a, const void *b) {
  T x = * (T *) a;
  T y = * (T *) b;
  
  if (x->watsonp != y->watsonp) {
    return false;		/* Different strands */

  } else if (y->low != x->low) {
    return false;
  } else if (y->high != x->high) {
    return false;
  } else {
    return true;
  }
}

int
Matchpair_size_cmp (const void *a, const void *b) {
  T x = * (T *) a;
  T y = * (T *) b;

  if (x->clustersize > y->clustersize) {
    return -1;
  } else if (x->clustersize < y->clustersize) {
    return +1;
  } else {
    return 0;
  }
}


List_T
Matchpair_filter_duplicates (List_T matchpairlist) {
  List_T unique = NULL, p, q;
  T x, y, matchpair;
  int n, i, j, cmp;
  bool *eliminate;

  n = List_length(matchpairlist);
  if (n == 0) {
    return NULL;
  }

  debug(
	for (p = matchpairlist, i = 0; p != NULL; p = p->rest, i++) {
	  matchpair = p->first;
	  printf("  Initial %d: %u-%u (watsonp = %d)\n",
		 i,matchpair->low,matchpair->high,matchpair->watsonp);
	}
	);

  eliminate = (bool *) CALLOC(n,sizeof(bool));

  /* Check for duplicates */
  for (p = matchpairlist, i = 0; p != NULL; p = p->rest, i++) {
    x = p->first;
    for (q = p->rest, j = i+1; q != NULL; q = q->rest, j++) {
      y = q->first;
      if (Matchpair_equal(&x,&y) == true) {
	eliminate[j] = true;
      }
    }
  }

  for (p = matchpairlist, i = 0; p != NULL; p = p->rest, i++) {
    if (eliminate[i] == false) {
      debug(matchpair = p->first);
      debug(printf("  Keeping %u-%u (watsonp = %d)\n",matchpair->low,matchpair->high,matchpair->watsonp));
      unique = List_push(unique,p->first);
    } else {
      matchpair = p->first;
      debug(printf("  Eliminating %u-%u (watsonp = %d)\n",matchpair->low,matchpair->high,matchpair->watsonp));
      Matchpair_free(&matchpair);
    }
  }

  FREE(eliminate);
  List_free(&matchpairlist);

  debug(
	for (p = unique, i = 0; p != NULL; p = p->rest, i++) {
	  matchpair = p->first;
	  printf("  Final %d: %u-%u (watsonp = %d)\n",
		 i,matchpair->low,matchpair->high,matchpair->watsonp);
	}
	);

  return unique;
}



List_T
Matchpair_filter_unique (List_T matchpairlist) {
  List_T unique = NULL, p, q;
  T x, y, matchpair;
  int n, i, j, cmp;
  bool *eliminate;

  n = List_length(matchpairlist);
  if (n == 0) {
    return NULL;
  }

  debug(
	for (p = matchpairlist, i = 0; p != NULL; p = p->rest, i++) {
	  matchpair = p->first;
	  printf("  Initial %d: %u-%u (watsonp = %d)\n",
		 i,matchpair->low,matchpair->high,matchpair->watsonp);
	}
	);

  eliminate = (bool *) CALLOC(n,sizeof(bool));
  /* Check for subsumption */
  for (p = matchpairlist, i = 0; p != NULL; p = p->rest, i++) {
    x = p->first;
    for (q = p->rest, j = i+1; q != NULL; q = q->rest, j++) {
      y = q->first;
      if ((cmp = Matchpair_dominate(&x,&y)) == -1) {
	eliminate[j] = true;
      } else if (cmp == 1) {
	eliminate[i] = true;
      }
    }
  }
    
  for (p = matchpairlist, i = 0; p != NULL; p = p->rest, i++) {
    if (eliminate[i] == false) {
      debug(matchpair = p->first);
      debug(printf("  Keeping %u-%u (watsonp = %d)\n",matchpair->low,matchpair->high,matchpair->watsonp));
      unique = List_push(unique,p->first);
    } else {
      matchpair = p->first;
      debug(printf("  Eliminating %u-%u (watsonp = %d)\n",matchpair->low,matchpair->high,matchpair->watsonp));
      Matchpair_free(&matchpair);
    }
  }

  FREE(eliminate);
  List_free(&matchpairlist);

  debug(
	for (p = unique, i = 0; p != NULL; p = p->rest, i++) {
	  matchpair = p->first;
	  printf("  Final %d: %u-%u (watsonp = %d)\n",
		 i,matchpair->low,matchpair->high,matchpair->watsonp);
	}
	);

  return unique;
}

void
Matchpair_get_coords (Genomicpos_T *chrpos, Genomicpos_T *genomicstart, Genomicpos_T *genomiclength,
		      bool *watsonp, T this, Stage1_T stage1, Sequence_T queryseq, Genomicpos_T chrlength, 
		      int maxextension) {
  Match_T match1, match2;
  int stage1size;
  Genomicpos_T chrpos1, chrpos2, genomicpos1, genomicpos2, left, right;
  int querylength;

  querylength = Sequence_length(queryseq);

  if (Match_position(this->bound3) > Match_position(this->bound5)) {
    *watsonp = true;
    Stage1_find_extensions(&left,&right,stage1,this->bound5,this->bound3,maxextension);
    match1 = this->bound5;
    match2 = this->bound3;
  } else if (Match_position(this->bound3) < Match_position(this->bound5)) {
    *watsonp = false;
    Stage1_find_extensions(&right,&left,stage1,this->bound5,this->bound3,maxextension);
    match1 = this->bound3;
    match2 = this->bound5;
  } else {
    /* Special case of a single match */
    if (Match_forwardp(this->bound5) == true) {
      *watsonp = true;
      Stage1_find_extensions(&left,&right,stage1,this->bound5,this->bound3,maxextension);
    } else {
      *watsonp = false;
      Stage1_find_extensions(&right,&left,stage1,this->bound5,this->bound3,maxextension);
    }
    match1 = match2 = this->bound5;
  }

  if (Match_chrpos(match1) < left) {
    debug(printf("Proposed left is negative: %u < %u.  Problem!\n",
		 Match_chrpos(match1),left));
    chrpos1 = 0;		/* Match_chrpos(match1) - Match_chrpos(match1); */
    genomicpos1 = Match_position(match1) - Match_chrpos(match1);
  } else {
    debug(printf("Proposed left is non-negative: %u >= %u.  Okay.\n",
		 Match_chrpos(match1),left));
    chrpos1 = Match_chrpos(match1) - left;
    genomicpos1 = Match_position(match1) - left;
  }

  if (Match_chrpos(match2) + right >= chrlength) {
    debug(printf("Proposed right is too high: %u + %u >= %u.  Problem!\n",
		 Match_chrpos(match2),right,chrlength));
    chrpos2 = chrlength - 1;	/* Match_chrpos(match2) + (chrlength - 1 - Match_chrpos(match2)); */
    genomicpos2 = Match_position(match2) + (chrlength - 1 - Match_chrpos(match2));
  } else {
    debug(printf("Proposed right is not too high: %u + %u < %u.  Okay.\n",
		 Match_chrpos(match2),right,chrlength));
    chrpos2 = Match_chrpos(match2) + right;
    genomicpos2 = Match_position(match2) + right;
  }

  *genomicstart = genomicpos1;
  *genomiclength = genomicpos2 - genomicpos1;
  *chrpos = chrpos1;
  debug(printf("Coordinates are %u, length %u\n",*genomicstart,*genomiclength));

  return;
}


/************************************************************************
 *   Clustering for stage 1
 ************************************************************************/

static T
find_best_path (int *size, Match_T *matches, int n, int matchsize, int minsize, 
		Match_T prevstart, Match_T prevend, bool plusp) {
  T matchpair;
  int *prev, *score, j, i, bestj, besti, bestscore;
  bool fivep = false, threep = false;
  Matchpairend_T matchpairend;
  Genomicpos_T endpos;

  prev = (int *) CALLOC(n,sizeof(int));
  score = (int *) CALLOC(n,sizeof(int));
  prev[0] = -1;

  debug1(
	 for (j = 0; j < n; j++) {
	   if (plusp == true) {
	     printf("+");
	   } else {
	     printf("-");
	   }
	   printf("#%d:%u (%d,%d,%d)\n",Match_chrnum(matches[j]),
		  Match_chrpos(matches[j]),Match_querypos(matches[j]),
		  Match_forwardp(matches[j]),Match_fivep(matches[j]));
	 }
	 );

  if (plusp == true) {
    for (j = 1; j < n; j++) {
      endpos = matches[j]->querypos;
      bestscore = 0;
      besti = -1;
      for (i = 0; i < j; i++) {
	if (matches[i]->querypos < endpos) {
	  if (score[i] + 1 > bestscore) {
	    bestscore = score[i] + 1;
	    besti = i;
	  }
	}
      }
      score[j] = bestscore;
      prev[j] = besti;
    }
  } else {
    for (j = 1; j < n; j++) {
      endpos = matches[j]->querypos;
      bestscore = 0;
      besti = -1;
      for (i = 0; i < j; i++) {
	if (matches[i]->querypos > endpos) {
	  if (score[i] + 1 > bestscore) {
	    bestscore = score[i] + 1;
	    besti = i;
	  }
	}
      }
      score[j] = bestscore;
      prev[j] = besti;
    }
  }

  bestscore = 0;
  bestj = 0;
  for (j = 0; j < n; j++) {
    if (score[j] > bestscore) {
      bestscore = score[j];
      bestj = j;
    }
  }

  /* Traceback */
  besti = bestj;
  if (Match_fivep(matches[besti]) == true) {
    fivep = true;
  } else {
    threep = true;
  }
  while (prev[besti] >= 0) {
    besti = prev[besti];
    if (Match_fivep(matches[besti]) == true) {
      fivep = true;
    } else {
      threep = true;
    }
  }
  debug1(printf("Best path is %d to %d, size=%d\n",besti,bestj,bestscore+1));

  if ((*size = bestscore + 1) < minsize) {
    debug1(printf("Size is smaller than minsize of %d.  Discarding.\n",minsize));
    matchpair = NULL;
  } else if (matches[besti] == prevstart && matches[bestj] == prevend) {
    debug1(printf("Path equal to previous one.  Discarding.\n"));
    matchpair = NULL;
  } else {
    if (fivep == true && threep == true) {
      matchpairend = MIXED;
    } else if (fivep == true) {
      matchpairend = FIVEONLY;
    } else if (threep == true) {
      matchpairend = THREEONLY;
    } else {
      abort();
    }
    if (plusp == true) {
      matchpair = Matchpair_new(matches[besti],matches[bestj],matchsize,bestscore+1,matchpairend);
    } else {
      matchpair = Matchpair_new(matches[bestj],matches[besti],matchsize,bestscore+1,matchpairend);
    }
  }
  debug1(printf("\n"));

  FREE(score);
  FREE(prev);

  return matchpair;
}

static List_T
find_paths_bounded (int *bestsize, List_T clusterlist, Match_T *matches, int npositions, 
		    int maxintronlen, int matchsize, int minclustersize, double sizebound, bool plusp) {
  int starti, endi;
  Chrnum_T startchrnum;
  Match_T prevstart = NULL, prevend = NULL;
  T matchpair;
  Genomicpos_T startpos;
  int minsize, size;

  minsize = minclustersize;
  if ((int) ((*bestsize) * sizebound) > minsize) {
    minsize = (int) ((*bestsize) * sizebound);
  }
  debug2(printf("minsize = %d\n",minsize));

  endi = -1;
  for (starti = 0; starti < npositions - 1; starti++) {
    startchrnum = Match_chrnum(matches[starti]);
    startpos = Match_position(matches[starti]);
    while (endi+1 < npositions && Match_chrnum(matches[endi+1]) == startchrnum &&
	   Match_position(matches[endi+1]) < startpos+maxintronlen) {
      endi = endi+1;
    }
    if (endi - starti + 1 >= minsize) {
      if ((matchpair = find_best_path(&size,&(matches[starti]),endi-starti+1,matchsize,
				      minsize,prevstart,prevend,plusp)) == NULL) {
	prevstart = prevend = NULL;
      } else {
	clusterlist = List_push(clusterlist,(void *) matchpair);
	if (plusp == true) {
	  prevstart = matchpair->bound5;
	  prevend = matchpair->bound3;
	} else {
	  prevstart = matchpair->bound3;
	  prevend = matchpair->bound5;
	}
	if (size > *bestsize) {
	  *bestsize = size;
	  if (((int) (*bestsize) * sizebound) > minsize) {
	    minsize = (int) ((*bestsize) * sizebound);
	    debug2(printf("resetting minsize = %d\n",minsize));
	  }
	}
      }
    }
  }

  return clusterlist;
}


static void
separate_strands (Match_T **plus_matches, int *plus_npositions,
		  Match_T **minus_matches, int *minus_npositions,
		  List_T matches5, List_T matches3) {
  List_T p;
  Match_T match;
  int i = 0, j = 0;

  *plus_npositions = *minus_npositions = 0;

  /* Count */
  for (p = matches5; p != NULL; p = p->rest) {
    match = p->first;
    if (Match_forwardp(match) == true) {
      (*plus_npositions)++;
    } else {
      (*minus_npositions)++;
    }
  }

  for (p = matches3; p != NULL; p = p->rest) {
    match = p->first;
    if (Match_forwardp(match) == true) {
      (*plus_npositions)++;
    } else {
      (*minus_npositions)++;
    }
  }

  /* Allocate */
  if (*plus_npositions == 0) {
    *plus_matches = NULL;
  } else {
    *plus_matches = (Match_T *) CALLOC(*plus_npositions,sizeof(Match_T));
  }

  if (*minus_npositions == 0) {
    *minus_matches = NULL;
  } else {
    *minus_matches = (Match_T *) CALLOC(*minus_npositions,sizeof(Match_T));
  }

  /* Fill */
  for (p = matches5; p != NULL; p = p->rest) {
    match = p->first;
    if (Match_forwardp(match) == true) {
      (*plus_matches)[i++] = match;
    } else {
      (*minus_matches)[j++] = match;
    }
  }

  for (p = matches3; p != NULL; p = p->rest) {
    match = p->first;
    if (Match_forwardp(match) == true) {
      (*plus_matches)[i++] = match;
    } else {
      (*minus_matches)[j++] = match;
    }
  }

  /* Sort */
  if (*plus_matches) {
    qsort(*plus_matches,*plus_npositions,sizeof(Match_T),Match_cmp);
  }
  if (*minus_matches) {
    qsort(*minus_matches,*minus_npositions,sizeof(Match_T),Match_cmp);
  }

  return;
}


static List_T
bound_results_testing (List_T clusterlist, int bestsize, double sizebound) {
  List_T boundedlist = NULL, p;
  T matchpair;
  int minsize;

  minsize = (int) (bestsize * sizebound);
  debug2(printf("minsize = %d\n",minsize));
  for (p = clusterlist; p != NULL; p = p->rest) {
    matchpair = (T) List_head(p);
    if (matchpair->clustersize >= minsize) {
      boundedlist = List_push(boundedlist,(void *) matchpair);
    } else {
      Matchpair_free(&matchpair);
    }
  }
  List_free(&clusterlist);

  boundedlist = Matchpair_filter_unique(boundedlist);
  boundedlist = Matchpair_filter_duplicates(boundedlist);

  return boundedlist;
}


static List_T
bound_results_limited (List_T clusterlist, int bestsize, double sizebound) {
  List_T boundedlist = NULL, p;
  T matchpair;
  int *counts, minsize, total;

  clusterlist = Matchpair_filter_unique(clusterlist);
  clusterlist = Matchpair_filter_duplicates(clusterlist);

  counts = (int *) CALLOC(bestsize+1,sizeof(int));
  for (p = clusterlist; p != NULL; p = p->rest) {
    matchpair = (T) List_head(p);
    counts[matchpair->clustersize]++;
  }

  minsize = bestsize;
  total = counts[bestsize];
  debug2(printf("Setting minsize = %d, total = %d\n",minsize,total));
  while (minsize-1 >= 0 && minsize-1 >= (int) (bestsize * sizebound) &&
	 total+counts[minsize-1] <= MAXCANDIDATES) {
    total = total+counts[minsize-1];
    minsize = minsize-1;
    debug2(printf("Setting minsize = %d, total = %d\n",minsize,total));
  }
  FREE(counts);

  for (p = clusterlist; p != NULL; p = p->rest) {
    matchpair = (T) List_head(p);
    if (matchpair->clustersize >= minsize) {
      boundedlist = List_push(boundedlist,(void *) matchpair);
    } else {
      Matchpair_free(&matchpair);
    }
  }
  List_free(&clusterlist);
  
  return boundedlist;
}


static void
print_clusters (List_T list) {
  List_T p;
  T matchpair;

  for (p = list; p != NULL; p = p->rest) {
    matchpair = (T) List_head(p);
    if (matchpair->watsonp == true) {
      printf("+");
    } else {
      printf("-");
    }
    printf("#%d:%u - %u (%d entries)\n",
	   Match_chrnum(matchpair->bound5),
	   Match_chrpos(matchpair->bound5),Match_chrpos(matchpair->bound3),
	   matchpair->clustersize);
  }
  return;
}


List_T
Matchpair_find_clusters (List_T matches5, List_T matches3, int stage1size, 
			 int maxintronlen, int minclustersize, double sizebound, 
			 Boundmethod_T boundmethod) {
  List_T boundedlist = NULL, clusterlist = NULL, p;
  Match_T *plus_matches = NULL, *minus_matches = NULL, match, start;
  T matchpair, *matchpairarray;
  Genomicpos_T startpos;
  int plus_npositions, minus_npositions, i = 0, j = 0, 
    nmatchpairs, minbestsize, bestsize = 0;

  separate_strands(&plus_matches,&plus_npositions,&minus_matches,&minus_npositions,matches5,matches3);

  if (plus_npositions > 0) {
    clusterlist = find_paths_bounded(&bestsize,clusterlist,plus_matches,plus_npositions,maxintronlen,
				     stage1size,minclustersize,sizebound,true);
  }

  if (minus_npositions > 0) {
    clusterlist = find_paths_bounded(&bestsize,clusterlist,minus_matches,minus_npositions,maxintronlen,
				     stage1size,minclustersize,sizebound,false);
  }
	
  FREE(plus_matches);
  FREE(minus_matches);

  if (boundmethod == NO_BOUND) {
    return clusterlist;
  } else if (boundmethod == MOVING_THRESHOLD) {
    boundedlist = bound_results_limited(clusterlist,bestsize,sizebound);

    debug2(printf("Bounded clusters (%d):\n",List_length(boundedlist)));
    debug2(print_clusters(boundedlist));
    return boundedlist;

  } else if (boundmethod == BY_CANDIDATES) {
    boundedlist = bound_results_testing(clusterlist,bestsize,sizebound);

    if (List_length(boundedlist) <= MAXCANDIDATES) {
      debug2(printf("Enough clusters (%d)\n",List_length(boundedlist)));
      debug2(print_clusters(boundedlist));
      return boundedlist;

    } else {
      debug2(printf("Too many clusters (%d)\n",List_length(boundedlist)));
      debug2(print_clusters(boundedlist));

      for (p = boundedlist; p != NULL; p = p->rest) {
	matchpair = (T) List_head(p);
	Matchpair_free(&matchpair);
      }
      List_free(&boundedlist);
      
      return NULL;
    }
  } else {
    abort();
  }
}


