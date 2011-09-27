static char rcsid[] = "$Id: gregion.c,v 1.12 2010-07-10 01:28:53 twu Exp $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "gregion.h"
#include <stdlib.h>
#include "mem.h"
#include "indexdb.h"		/* For SUFFICIENT_SUPPORT */
#include "interval.h"
#include "chrnum.h"
#include "listdef.h"


#define EXTRA_SHORTEND  30000
#define EXTRA_LONGEND   100000


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/* Sufficient support */
#ifdef DEBUG5
#define debug5(x) x
#else
#define debug5(x)
#endif


#define T Gregion_T
struct T {
  int nexons;

  Genomicpos_T genomicstart;
  Genomicpos_T genomicend;
  Genomicpos_T genomiclength;
  bool plusp;

  Genomicpos_T extension5;
  Genomicpos_T extension3;
  bool extendedp;

  Chrnum_T chrnum;
  Genomicpos_T chrpos;
  Genomicpos_T chroffset;	/* This is for chr, not the segment */
  Genomicpos_T chrlength;	/* This is for chr, not the segment */
  int querystart;		/* Used only for maponly mode */
  int queryend;			/* Used only for maponly mode */
  int matchsize;		/* Used only for maponly mode */

  int trimstart;
  int trimend;
  bool sufficient_support_p;

  double weight;		/* Product of match weights */
  int support;

  int ncovered;
  int source;			/* Oligoindex source in stage 2 */
};



void
Gregion_print (T this) {

  printf(" %d..%d  %u..%u  #%d:%u..%u  length:%u  weight:%.2f  support:%d",
	 this->querystart,this->queryend,this->genomicstart,this->genomicend,this->chrnum,this->chrpos,
	 this->chrpos+this->genomiclength,this->genomiclength,this->weight,this->support);
  printf("\n");

  return;
}


void
Gregion_free (T *old) {
  FREE(*old);
  return;
}

Genomicpos_T
Gregion_genomicstart (T this) {
  return this->genomicstart;
}

Genomicpos_T
Gregion_genomicend (T this) {
  return this->genomicend;
}

Genomicpos_T
Gregion_genomiclength (T this) {
  return this->genomiclength;
}

bool
Gregion_plusp (T this) {
  return this->plusp;
}

bool
Gregion_revcompp (T this) {
  return !(this->plusp);
}

Chrnum_T
Gregion_chrnum (T this) {
  return this->chrnum;
}

/* Used only for debugging.  String is allocated and should be freed. */
char *
Gregion_chr (T this, IIT_T chromosome_iit) {
  return Chrnum_to_string(this->chrnum,chromosome_iit);
}

Genomicpos_T
Gregion_chrpos (T this) {
  return this->chrpos;
}

Genomicpos_T
Gregion_chroffset (T this) {
  return this->chroffset;
}

Genomicpos_T
Gregion_chrlength (T this) {
  return this->chrlength;
}

int
Gregion_querystart (T this) {
  return this->querystart;
}

int
Gregion_queryend (T this) {
  return this->queryend;
}

int
Gregion_matchsize (T this) {
  return this->matchsize;
}

double
Gregion_weight (T this) {
  return this->weight;
}

int
Gregion_support (T this) {
  return this->support;
}

bool 
Gregion_extendedp (T this) {
  return this->extendedp;
}

void
Gregion_set_ncovered (T this, int ncovered, int source) {
  this->ncovered = ncovered;
  this->source = source;
  return;
}

int
Gregion_ncovered (T this) {
  return this->ncovered;
}



T
Gregion_new (int nexons, Genomicpos_T genomicstart, Genomicpos_T genomicend,
	     bool plusp, IIT_T chromosome_iit, int querystart, int queryend, int querylength,
	     int matchsize, int trimstart, int trimend) {
  T new = (T) MALLOC(sizeof(*new));

  debug(printf("Creating gregion with genomicstart %u, genomicend %u\n",
	       genomicstart,genomicend));

  new->nexons = nexons;
  if (chromosome_iit == NULL) {
    new->chrnum = 0;
    new->chroffset = 0U;
    new->chrlength = 0U;
  } else {
    new->chrnum = IIT_get_one(chromosome_iit,/*divstring*/NULL,genomicstart,genomicstart);
    new->chroffset = Interval_low(IIT_interval(chromosome_iit,new->chrnum));
    new->chrlength = Interval_high(IIT_interval(chromosome_iit,new->chrnum)) - new->chroffset;
  }
  new->chrpos = genomicstart - new->chroffset;
  
  new->genomicstart = genomicstart;
  new->genomicend = genomicend;
  if (genomicend < genomicstart) {
    fprintf(stderr,"genomicend %u < genomicstart %u\n",genomicend,genomicstart);
    abort();
  } else {
    new->genomiclength = genomicend - genomicstart;
  }
  new->plusp = plusp;

  new->extension5 = 0U;
  new->extension3 = 0U;
  new->extendedp = false;

  new->querystart = querystart;
  new->queryend = queryend;
  new->matchsize = matchsize;

  new->trimstart = trimstart;
  new->trimend = trimend;

#ifdef PMAP
  debug5(printf("  Testing bound5+extension5 = %d - %u < %d + %d, bound3+extension3 = %d + %u > %d - %d\n",
		new->querystart,new->extension5,trimstart,SUFFICIENT_SUPPORT,
		new->queryend,new->extension3,trimend,SUFFICIENT_SUPPORT));
  if (new->querystart - new->extension5 < trimstart + SUFFICIENT_SUPPORT && 
      new->queryend + new->extension3 > trimend - SUFFICIENT_SUPPORT) {
    new->sufficient_support_p = true;
  } else {
    new->sufficient_support_p = false;
  }
#else
  debug5(printf("  Testing bound5+extension5 = %d - %u < %d + %d, bound3 = %d + %u > %d - %d\n",
		new->querystart,new->extension5,trimstart,SUFFICIENT_SUPPORT,
		new->queryend,new->extension3,trimend,SUFFICIENT_SUPPORT));
  if (new->querystart - new->extension5 < trimstart + SUFFICIENT_SUPPORT && 
      new->queryend + new->extension3 > trimend - SUFFICIENT_SUPPORT) {
    new->sufficient_support_p = true;
  } else {
    new->sufficient_support_p = false;
  }
#endif

  new->weight = 0.0;
  new->support = queryend - querystart + matchsize;

  return new;
}


T
Gregion_new_from_matches (Match_T match5, Match_T match3, IIT_T chromosome_iit, 
			  int querylength, int matchsize, int trimstart, int trimend) {
  T gregion;
  Genomicpos_T genomicstart, genomicend;
  int querystart, queryend;

  querystart = Match_querypos(match5);
  queryend = Match_querypos(match3);

  if (Match_forwardp(match5)) {
    genomicstart = Match_position(match5);
    genomicend = Match_position(match3) + 1U;
    /* chrnum = Match_chrnum(match5); */

  } else {
    genomicstart = Match_position(match3);
    genomicend = Match_position(match5) + 1U;
    /* chrnum = Match_chrnum(match3); */
  }

#if 0
  if (chromosome_iit == NULL) {
    chroffset = 0U;
    chrlength = 0U;
  } else {
    chroffset = Interval_low(IIT_interval(chromosome_iit,chrnum));
    chrlength = Interval_high(IIT_interval(chromosome_iit,chrnum)) - chroffset;
  }
#endif

  debug(printf("Coordinates are %u .. %u\n",genomicstart,genomicend));

  gregion = Gregion_new(/*nexons*/0,genomicstart,genomicend,Match_forwardp(match5),chromosome_iit,
			querystart,queryend,querylength,matchsize,trimstart,trimend);

  gregion->weight = Match_weight(match5) * Match_weight(match3);
  Match_incr_npairings(match5);
  Match_incr_npairings(match3);

  return gregion;
}



static int
weight_cmp (const void *x, const void *y) {
  T a = * (T *) x;
  T b = * (T *) y;

  if (a->weight > b->weight) {
    return -1;
  } else if (a->weight < b->weight) {
    return +1;
  } else if (a->support > b->support) {
    return -1;
  } else if (a->support < b->support) {
    return +1;
  } else if (a->genomiclength > b->genomiclength) {
    return -1;
  } else if (a->genomiclength < b->genomiclength) {
    return +1;
  } else {
    return 0;
  }
}


/* Not intended for qsort.  Returns 0 when not comparable. */
static bool
gregion_overlap_p (T x, T y) {

  if (x->plusp != y->plusp) {
    return false;			/* Different strands */

  } else if (y->genomicstart > x->genomicend || x->genomicstart > y->genomicend) {
    return false;		/* No overlap */

  } else {
    return true;
  }
}


List_T
Gregion_filter_unique (List_T gregionlist) {
  List_T unique = NULL;
  T x, y, gregion, *array;
  int n, i, j;
  bool *eliminate;
#ifdef DEBUG
  List_T p, q;
#endif

  n = List_length(gregionlist);
  if (n == 0) {
    return NULL;
  }

  debug(
	for (p = gregionlist, i = 0; p != NULL; p = p->rest, i++) {
	  gregion = (T) p->first;
	  printf("  Initial %d: %d..%d %u-%u (plusp = %d)\n",
		 i,gregion->querystart,gregion->queryend,gregion->genomicstart,gregion->genomicend,gregion->plusp);
	}
	);

  eliminate = (bool *) CALLOC(n,sizeof(bool));

  /* Not necessary if false is zero */
  /*
  for (i = 0; i < n; i++) {
    eliminate[i] = false;
  }
  */

  array = (T *) List_to_array(gregionlist,NULL);
  List_free(&gregionlist);
  qsort(array,n,sizeof(T),weight_cmp);

  for (i = 0; i < n; i++) {
    x = array[i];
    for (j = i+1; j < n; j++) {
      y = array[j];
      if (gregion_overlap_p(x,y) == true) {
	eliminate[j] = true;
      }
    }
  }

  for (i = n-1; i >= 0; i--) {
    gregion = array[i];
    if (eliminate[i] == false) {
      debug(printf("  Keeping %u-%u (plusp = %d)\n",
		   gregion->genomicstart,gregion->genomicend,gregion->plusp));
      unique = List_push(unique,(void *) gregion);
    } else {
      debug(printf("  Eliminating %u-%u (plusp = %d)\n",
		   gregion->genomicstart,gregion->genomicend,gregion->plusp));
      /*
      if (gregion->match5 != NULL) {
	Match_decr_npairings(gregion->match5);
	Match_decr_npairings(gregion->match3);
      }
      */
      Gregion_free(&gregion);
    }
  }

  FREE(eliminate);
  FREE(array);

  debug(
	for (p = unique, i = 0; p != NULL; p = p->rest, i++) {
	  gregion = (T) p->first;
	  printf("  Final %d: %u-%u (plusp = %d)\n",
		 i,gregion->genomicstart,gregion->genomicend,gregion->plusp);
	}
	);

  return unique;
}


List_T
Gregion_filter_support (List_T gregionlist, int boundary, double pct_max, int diff_max) {
  List_T good = NULL, p;
  int threshold, maxsupport = 0;
  T gregion;

  if (gregionlist == NULL) {
    return NULL;
  }

  for (p = gregionlist; p != NULL; p = List_next(p)) {
    gregion = (T) List_head(p);
    if (gregion->support > maxsupport) {
      maxsupport = gregion->support;
    }
  }

  if (maxsupport > boundary) {
    /* Use diff_max */
    for (p = gregionlist; p != NULL; p = List_next(p)) {
      gregion = (T) List_head(p);
      if (maxsupport - gregion->support < diff_max) {
	good = List_push(good,gregion);
      } else {
	Gregion_free(&gregion);
      }
    }
  } else {
    threshold = (int) ((double) maxsupport * pct_max);
    /* Use threshold only */
    for (p = gregionlist; p != NULL; p = List_next(p)) {
      gregion = (T) List_head(p);
      if (gregion->support >= threshold) {
	good = List_push(good,gregion);
      } else {
	Gregion_free(&gregion);
      }
    }
  }

  List_free(&gregionlist);
  return List_reverse(good);
}


double
Gregion_best_weight (List_T gregionlist) {
  double best_weight = 0.0;
  T gregion;
  List_T p;

  for (p = gregionlist; p != NULL; p = List_next(p)) {
    gregion = (T) List_head(p);
    if (gregion->weight > best_weight) {
      best_weight = gregion->weight;
    }
  }

  return best_weight;
}



#if 0

/* Not intended for qsort.  Returns 0 when not comparable. */
static int
gregion_dominate (const void *a, const void *b) {
  T x = * (T *) a;
  T y = * (T *) b;

  if (x->plusp != y->plusp) {
    return 0;			/* Different strands */

  } else if (y->genomicstart > x->genomicend || x->genomicstart > y->genomicend) {
    return 0;			/* No overlap */

  } else if (x->genomicstart < y->genomicstart) {
    if (x->genomicend < y->genomicend) {
      return 0;			/* Crossing */
    } else {
      return -1;
    }

  } else if (y->genomicstart < x->genomicstart) {
    if (y->genomicend < x->genomicend) {
      return 0;			/* Crossing */
    } else {
      return 1;
    }

  } else {
    /* equal genomicstart */
    if (x->genomicend < y->genomicend) {
      return 1;
    } else if (y->genomicend < x->genomicend) {
      return -1;

    } else if (y->querystart > x->queryend || x->querystart > y->queryend) {
      return 0;			/* Equal */
    } else if (x->querystart < y->querystart) {
      if (x->queryend < y->queryend) {
	return 0;
      } else {
	return -1;
      }
    } else if (y->querystart < x->querystart) {
      if (y->queryend < x->queryend) {
	return 0;
      } else {
	return 1;
      }
    } else {
      return 0;
    }
  }
}

static bool
gregion_equal (const void *a, const void *b) {
  T x = * (T *) a;
  T y = * (T *) b;
  
  if (x->plusp != y->plusp) {
    return false;		/* Different strands */

  } else if (y->genomicstart != x->genomicstart) {
    return false;
  } else if (y->genomicend != x->genomicend) {
    return false;
  } else {
    return true;
  }
}

List_T
Gregion_filter_unique_old (List_T gregionlist) {
  List_T unique = NULL, p, q;
  T x, y, gregion;
  int n, i, j, cmp;
  bool *eliminate;

  n = List_length(gregionlist);
  if (n == 0) {
    return NULL;
  }

  debug(
	for (p = gregionlist, i = 0; p != NULL; p = p->rest, i++) {
	  gregion = (T) p->first;
	  printf("  Initial %d: %d..%d %u-%u (plusp = %d), %d*%d points\n",
		 i,gregion->querystart,gregion->queryend,gregion->genomicstart,gregion->genomicend,gregion->plusp,
		 gregion->points5,gregion->points3);
	}
	);

  eliminate = (bool *) CALLOC(n,sizeof(bool));

  /* Not necessary if false is zero */
  /*
  for (i = 0; i < n; i++) {
    eliminate[i] = false;
  }
  */

  /* Check for subsumption */
  for (p = gregionlist, i = 0; p != NULL; p = p->rest, i++) {
    x = (T) p->first;
    for (q = p->rest, j = i+1; q != NULL; q = q->rest, j++) {
      y = q->first;
      if ((cmp = gregion_dominate(&x,&y)) == -1) {
	debug(printf("  Gregion %d..%d %u-%u (plusp = %d), %d*%d points eliminates %d..%d %u-%u (plusp = %d), %d*%d points\n",
		     x->querystart,x->queryend,x->genomicstart,x->genomicend,x->plusp,x->points5,x->points3,
		     y->querystart,y->queryend,y->genomicstart,y->genomicend,y->plusp,y->points5,y->points3));
	eliminate[j] = true;
	if (x->querypositions[y->querystart] == false) {
	  x->querypositions[y->querystart] = true;
	  x->points5 += 1;
	}
	if (x->querypositions[y->queryend] == false) {
	  x->querypositions[y->queryend] = true;
	  x->points3 += 1;
	}
      } else if (cmp == 1) {
	debug(printf("  Gregion %d..%d %u-%u (plusp = %d), %d*%d points eliminates %d..%d %u-%u (plusp = %d), %d*%d points\n",
		     y->querystart,y->queryend,y->genomicstart,y->genomicend,y->plusp,y->points5,y->points3,
		     x->querystart,x->queryend,x->genomicstart,x->genomicend,x->plusp,x->points5,x->points3));
	eliminate[i] = true;
	if (y->querypositions[x->querystart] == false) {
	  y->querypositions[x->querystart] = true;
	  y->points5 += 1;
	}
	if (y->querypositions[x->queryend] == false) {
	  y->querypositions[x->queryend] = true;
	  y->points3 += 1;
	}
      } else if (gregion_equal(&x,&y) == true) {
	/* Eliminate second occurrence */
	debug(printf("  Gregion %d..%d %u-%u (plusp = %d), %d*%d points eliminates %d..%d %u-%u (plusp = %d), %d*%d points\n",
		     x->querystart,x->queryend,x->genomicstart,x->genomicend,x->plusp,x->points5,x->points3,
		     y->querystart,y->queryend,y->genomicstart,y->genomicend,y->plusp,y->points5,y->points3));
	eliminate[j] = true;
	if (x->querypositions[y->querystart] == false) {
	  x->querypositions[y->querystart] = true;
	  x->points5 += 1;
	}
	if (x->querypositions[y->queryend] == false) {
	  x->querypositions[y->queryend] = true;
	  x->points3 += 1;
	}
      }
    }
  }
    
  for (p = gregionlist, i = 0; p != NULL; p = p->rest, i++) {
    if (eliminate[i] == false) {
      debug(gregion = p->first);
      debug(printf("  Keeping %u-%u (plusp = %d)\n",
		   gregion->genomicstart,gregion->genomicend,gregion->plusp));
      unique = List_push(unique,(void *) p->first);
    } else {
      gregion = (T) p->first;
      debug(printf("  Eliminating %u-%u (plusp = %d)\n",
		   gregion->genomicstart,gregion->genomicend,gregion->plusp));
      /*
      if (gregion->match5 != NULL) {
	Match_decr_npairings(gregion->match5);
	Match_decr_npairings(gregion->match3);
      }
      */
      Gregion_free(&gregion);
    }
  }

  FREE(eliminate);
  List_free(&gregionlist);

  debug(
	for (p = unique, i = 0; p != NULL; p = p->rest, i++) {
	  gregion = (T) p->first;
	  printf("  Final %d: %u-%u (plusp = %d)\n",
		 i,gregion->genomicstart,gregion->genomicend,gregion->plusp);
	}
	);

  return unique;
}



List_T
Gregion_filter_by_evidence (List_T gregionlist) {
  List_T goodlist = NULL, p;
  T gregion;
  int threshold_evidence, bestevidence = 0, evidence;

  for (p = gregionlist; p != NULL; p = List_next(p)) {
    gregion = (T) List_head(p);
    if (gregion->match5 != NULL) {
      evidence = Match_npairings(gregion->match5) * Match_npairings(gregion->match3);
      if (evidence > bestevidence) {
	bestevidence = evidence;
      }
    }
  }

  threshold_evidence = bestevidence / 2;

  for (p = gregionlist; p != NULL; p = p->rest) {
    gregion = (T) List_head(p);
    if (gregion->match5 == NULL) {
      goodlist = List_push(goodlist,(void *) p->first);
    } else {
      evidence = Match_npairings(gregion->match5) * Match_npairings(gregion->match3);
#if 0
      if (evidence > threshold_evidence) {
#else
	if (Match_npairings(gregion->match5) > 2 && Match_npairings(gregion->match3) > 2) {
#endif
	goodlist = List_push(goodlist,(void *) p->first);
      } else {
	fprintf(stderr,"  Eliminating %u-%u (plusp = %d)\n",
		gregion->genomicstart,gregion->genomicend,gregion->plusp);
	debug(printf("  Eliminating %u-%u (plusp = %d)\n",
		     gregion->genomicstart,gregion->genomicend,gregion->plusp));
	Gregion_free(&gregion);
      }
    }
  }

  List_free(&gregionlist);
  return goodlist;
}
#endif


bool
Gregion_sufficient_support (T this) {
  return this->sufficient_support_p;
}


void
Gregion_extend (T this, Genomicpos_T extension5, Genomicpos_T extension3, int querylength) {
  Genomicpos_T chrend, left, right;

  debug(printf("Entering Gregion_extend with extension5 %u and extension3 %u\n",extension5,extension3));
  debug(printf("  genomicstart %u, genomiclength %u\n",this->genomicstart,this->genomiclength));
  this->extension5 = extension5;
  this->extension3 = extension3;
  this->extendedp = true;

  if (this->nexons == 1 || Gregion_sufficient_support(this) == true || this->support < 100) {
    if (this->plusp == true) {
      left = extension5 + querylength + EXTRA_SHORTEND;
      right = extension3 + querylength + EXTRA_SHORTEND;
    } else {
      left = extension3 + querylength + EXTRA_SHORTEND;
      right = extension5 + querylength + EXTRA_SHORTEND;
    }
  } else {
    if (this->plusp == true) {
      left = extension5 + EXTRA_LONGEND;
      right = extension3 + EXTRA_LONGEND;
    } else {
      left = extension3 + EXTRA_LONGEND;
      right = extension5 + EXTRA_LONGEND;
    }
  }

  chrend = this->chrpos + this->genomiclength;

  if (this->chrpos < left) {
    /* At beginning of chromosome */
    this->genomicstart -= this->chrpos;
    this->chrpos = 0U;
  } else {
    this->genomicstart -= left;
    this->chrpos -= left;
  }

  if (chrend + right >= this->chrlength) {
    /* At end of chromosome */
    this->genomicend = this->chroffset + this->chrlength;
  } else {
    this->genomicend += right;
  }

  this->genomiclength = this->genomicend - this->genomicstart + 1U;

#ifdef PMAP
  debug5(printf("  Testing bound5+extension5 = %d - %u < %d + %d, bound3+extension3 = %d + %u > %d - %d\n",
		this->querystart,this->extension5,this->trimstart,SUFFICIENT_SUPPORT,
		this->queryend,this->extension3,this->trimend,SUFFICIENT_SUPPORT));
  if (this->querystart - this->extension5 < this->trimstart + SUFFICIENT_SUPPORT && 
      this->queryend + this->extension3 > this->trimend - SUFFICIENT_SUPPORT) {
    this->sufficient_support_p = true;
  } else {
    this->sufficient_support_p = false;
  }
#else
  debug5(printf("  Testing bound5+extension5 = %d - %u < %d + %d, bound3 = %d + %u > %d - %d\n",
		this->querystart,this->extension5,this->trimstart,SUFFICIENT_SUPPORT,
		this->queryend,this->extension3,this->trimend,SUFFICIENT_SUPPORT));
  if (this->querystart - this->extension5 < this->trimstart + SUFFICIENT_SUPPORT && 
      this->queryend + this->extension3 > this->trimend - SUFFICIENT_SUPPORT) {
    this->sufficient_support_p = true;
  } else {
    this->sufficient_support_p = false;
  }
#endif

  debug(printf("  genomicstart %u, genomiclength %u\n",this->genomicstart,this->genomiclength));
  return;
}


int
Gregion_cmp (const void *a, const void *b) {
  T x = * (T *) a;
  T y = * (T *) b;

  if (x->ncovered > y->ncovered) {
    return -1;
  } else if (y->ncovered > x->ncovered) {
    return +1;
  } else {
    return 0;
  }
}
