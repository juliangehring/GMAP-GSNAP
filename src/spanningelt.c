static char rcsid[] = "$Id: spanningelt.c 100404 2013-07-03 21:32:22Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "spanningelt.h"
#include <stdlib.h>
#include <math.h>		/* For qsort */
#include "mem.h"
#include "indexdbdef.h"

#ifdef WORDS_BIGENDIAN
#include "bigendian.h"
#endif

#define T Spanningelt_T

#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/* List all positions */
#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif

#ifdef CHECK
#define check(x) x
#else
#define check(x)
#endif

/* binary_search */
#ifdef DEBUG10
#define debug10(x) x
#else
#define debug10(x)
#endif

static bool free_positions_p;	/* Needs to be true if Indexdb positions are FILEIO */

static Width_T index1part;
static Width_T index1interval;
static Width_T spansize;

/* spansize is the interval at which spanningelts must be spaced, due
   to sampling (i.e., largest multiple of index1interval >=
   index1part.  We want spansize to exceed index1part, because we want
   a lower bound on the number of mismatches (to be checked by actual
   alignment against the genome).  If spansize is smaller than
   index1part, this lower bound won't hold. */

Width_T
Spanningelt_setup (Width_T index1part_in, Width_T index1interval_in) {
  index1part = index1part_in;
  index1interval = index1interval_in;

  spansize = 0;			/* Refer to global spansize */
  while (spansize + index1interval_in < index1part) {
    spansize += index1interval_in;
  }
  spansize += index1interval_in;

  return spansize;
}
  

void
Spanningelt_init_positions_free (bool positions_fileio_p) {
  if (positions_fileio_p == true) {
    free_positions_p = true;
  } else {
    free_positions_p = false;
  }
  return;
}

void
Spanningelt_free (T *old) {

  if ((*old)->intersection_diagonals_reset != NULL) {
    FREE((*old)->intersection_diagonals_reset);
  }
  if ((*old)->compoundpos != NULL) {
    Compoundpos_free(&((*old)->compoundpos));
  }
  FREE((*old)->positions_allocated);
  /* Should not free partner_positions */

  FREE(*old);
  return;
}


T
Spanningelt_reset (T this) {
  this->intersection_diagonals = this->intersection_diagonals_reset;
  this->intersection_ndiagonals = this->intersection_ndiagonals_reset;

  this->partner_positions = this->partner_positions_reset;
  this->partner_npositions = this->partner_npositions_reset;

  this->positions = this->positions_reset;
  this->npositions = this->npositions_reset;

  if (this->compoundpos != NULL) {
    Compoundpos_reset(this->compoundpos);
    Compoundpos_heap_init(this->compoundpos,this->querylength,this->compoundpos_diagterm);
  }

  return this;
}

void
Spanningelt_print (T this) {
#ifdef DEBUG1
  int i;
#endif

  if (this->intersection_diagonals != NULL) {
    printf("Intersection (%d diagonals), ",this->intersection_ndiagonals);
#ifdef DEBUG1
    for (i = 0; i < this->intersection_ndiagonals; i++) {
      printf("%u ",this->intersection_diagonals[i]);
    }
#endif
  }

  if (this->partnerp == true) {
    printf("Partner @ %d (%d positions, diagterm %d), ",
	   this->partner_querypos,this->partner_npositions,this->partner_diagterm);
#ifdef DEBUG1
    for (i = 0; i < this->partner_npositions; i++) {
      printf("%u ",this->partner_positions[i]);
    }
#endif
  }

  if (this->compoundpos != NULL) {
    printf("Compound @ %d (diagterm %d): ",this->querypos,this->compoundpos_diagterm);
    Compoundpos_print_sizes(this->compoundpos);

  } else {
    printf("Ordinary @ %d (%d positions, diagterm %d)",
	   this->querypos,this->npositions,this->diagterm);
#ifdef DEBUG1
    for (i = 0; i < this->npositions; i++) {
      printf(" %u",this->positions[i]);
    }
#endif
  }

  printf("\n");
  return;
}

void
Spanningelt_print_set (List_T spanningset) {
  List_T p;

  printf("Set of length %d\n",List_length(spanningset));
  for (p = spanningset; p != NULL; p = List_next(p)) {
    Spanningelt_print((T) List_head(p));
  }
  printf("\n");
  return;
}



static T
Spanningelt_new (Storedoligomer_T *stage1_oligos, bool **stage1_retrievedp,
		 Univcoord_T ***stage1_positions, int **stage1_npositions,
		 Indexdb_T indexdb, bool partnerp, int querypos1, int querypos2,
		 int query_lastpos, int querylength, bool plusp) {
  T new = (T) MALLOC(sizeof(*new));
  int partnerpos, querypos;

  debug(printf("Entered Spanningelt_new with querypos1 %d, querypos2 %d, query_lastpos %d, querylength %d\n",
	       querypos1,querypos2,query_lastpos,querylength));
  new->querylength = querylength;
  new->intersection_diagonals = (Univcoord_T *) NULL;
  new->intersection_ndiagonals = 0;
  new->intersection_diagonals_reset = (Univcoord_T *) NULL;
  new->intersection_ndiagonals_reset = 0;

  /* Handle partner */
  if ((new->partnerp = partnerp) == false) {
    debug(printf("partnerp is false, so querypos = querypos1 %d\n",querypos1));
    new->partner_positions = (Univcoord_T *) NULL;
    new->partner_npositions = 0;
    new->querypos = querypos = querypos1;

  } else {
    if (querypos1 < 0 || querypos1 > query_lastpos) {
      /* Treat querypos2 as partner */
      debug(printf("partnerp is true and querypos1 is out of bounds, so querypos = querypos1 %d\n",querypos1));
      new->querypos = querypos = querypos1;
      new->partner_querypos = partnerpos = querypos2;

    } else {
      /* Treat querypos1 as partner */
      debug(printf("partnerp is true and querypos1 is in bounds, so querypos = querypos2 %d\n",querypos2));
      new->querypos = querypos = querypos2;
      new->partner_querypos = partnerpos = querypos1;
    }
    new->partner_positions = 
      Indexdb_read_inplace(&new->partner_npositions,indexdb,stage1_oligos[partnerpos]);
    new->partner_diagterm = plusp ? querylength - partnerpos : partnerpos + index1part; /* FORMULA */

    (*stage1_retrievedp)[partnerpos] = true;
    (*stage1_positions)[partnerpos] = new->partner_positions;
    (*stage1_npositions)[partnerpos] = new->partner_npositions;
    
  }

  /* Handle main querypos */
  if (querypos >= 0 && querypos <= query_lastpos) {
    /* Ordinary position */
    new->compoundpos = NULL;
    new->positions = 
      Indexdb_read_inplace(&new->npositions,indexdb,stage1_oligos[querypos]);
    if (free_positions_p == true) {
      new->positions_allocated = new->positions;
    } else {
      new->positions_allocated = (Univcoord_T *) NULL;
    }
    new->diagterm = plusp ? querylength - querypos : querypos + index1part; /* FORMULA */

  } else if (plusp) {
    /* Plus compoundpos */
    new->positions = NULL;
    new->positions_allocated = NULL;
    new->npositions = 0;

    if (querypos == -2) {
      new->compoundpos = Indexdb_compoundpos_left_subst_2(indexdb,stage1_oligos[0]);
      new->compoundpos_diagterm = querylength+2; /* FORMULA */

    } else if (querypos == -1) {
      new->compoundpos = Indexdb_compoundpos_left_subst_1(indexdb,stage1_oligos[0]);
      new->compoundpos_diagterm = querylength+1; /* FORMULA */
      
    } else if (querypos == query_lastpos + 1) {
      new->compoundpos = Indexdb_compoundpos_right_subst_1(indexdb,stage1_oligos[query_lastpos]);
      new->compoundpos_diagterm = index1part-1; /* FORMULA */
      
    } else if (querypos == query_lastpos + 2) {
      new->compoundpos = Indexdb_compoundpos_right_subst_2(indexdb,stage1_oligos[query_lastpos]);
      new->compoundpos_diagterm = index1part-2; /* FORMULA */
      
    }
    Compoundpos_heap_init(new->compoundpos,querylength,new->compoundpos_diagterm);
    Compoundpos_set(new->compoundpos);

  } else {
    /* Minus compoundpos */
    new->positions = NULL;
    new->positions_allocated = NULL;
    new->npositions = 0;

    if (querypos == -2) {
      new->compoundpos = Indexdb_compoundpos_right_subst_2(indexdb,stage1_oligos[0]);
      new->compoundpos_diagterm = index1part-2; /* FORMULA */

    } else if (querypos == -1) {
      new->compoundpos = Indexdb_compoundpos_right_subst_1(indexdb,stage1_oligos[0]);
      new->compoundpos_diagterm = index1part-1; /* FORMULA */
      
    } else if (querypos == query_lastpos + 1) {
      new->compoundpos = Indexdb_compoundpos_left_subst_1(indexdb,stage1_oligos[query_lastpos]);
      new->compoundpos_diagterm = querylength+1; /* FORMULA */
      
    } else if (querypos == query_lastpos + 2) {
      new->compoundpos = Indexdb_compoundpos_left_subst_2(indexdb,stage1_oligos[query_lastpos]);
      new->compoundpos_diagterm = querylength+2; /* FORMULA */
      
    }
    Compoundpos_heap_init(new->compoundpos,querylength,new->compoundpos_diagterm);
    Compoundpos_set(new->compoundpos);
  }
      
  if (partnerp == true) {
    new->candidates_score = new->partner_npositions + (*stage1_npositions)[querypos];
    new->pruning_score = new->partner_npositions;
  } else {
    new->candidates_score = new->pruning_score = (*stage1_npositions)[querypos];
  }

  if (plusp) {
    if (partnerp == false) {
      new->miss_querypos5 = querypos;
      new->miss_querypos3 = querypos + index1part;
    } else if (partnerpos < querypos) {
      new->miss_querypos5 = partnerpos;
      new->miss_querypos3 = querypos + index1part;
    } else {
      new->miss_querypos5 = querypos;
      new->miss_querypos3 = partnerpos + index1part;
    }
  } else {
    if (partnerp == false) {
      new->miss_querypos5 = query_lastpos - querypos;
      new->miss_querypos3 = querylength - querypos;
    } else if (partnerpos < querypos) {
      new->miss_querypos5 = query_lastpos - querypos;
      new->miss_querypos3 = querylength - partnerpos;
    } else {
      new->miss_querypos5 = query_lastpos - partnerpos;
      new->miss_querypos3 = querylength - querypos;
    }
  }

  if (new->miss_querypos5 < 0) {
    new->miss_querypos5 = 0;
  }
  if (new->miss_querypos3 > querylength) {
    new->miss_querypos3 = querylength;
  }

  /* Set for later reset */
  new->partner_positions_reset = new->partner_positions;
  new->partner_npositions_reset = new->partner_npositions;

  new->positions_reset = new->positions;
  new->npositions_reset = new->npositions;

  return new;
}



static List_T
make_spanningset_aux (int *minscore, Storedoligomer_T *stage1_oligos, bool **stage1_retrievedp,
		      Univcoord_T ***stage1_positions, int **stage1_npositions,
		      Indexdb_T indexdb, int query_lastpos, int querylength,
		      int first, int last, bool plusp) {
  List_T spanningset;
  T elt;
  int worstpos, partnerpos, leftpos, rightpos, querypos, diff;
  int maxpositions = 0;
  bool first_anchored_p = true;


  worstpos = first;
  for (querypos = first; querypos < last; querypos += spansize) {
    if ((*stage1_npositions)[querypos] > maxpositions) {
      maxpositions = (*stage1_npositions)[querypos];
      worstpos = querypos;
    }
  }

  for (querypos = last; querypos > first; querypos -= spansize) {
    if ((*stage1_npositions)[querypos] > maxpositions) {
      maxpositions = (*stage1_npositions)[querypos];
      worstpos = querypos;
      first_anchored_p = false;
    }
  }

  /* Handle middle position first */
  diff = (last - first) % spansize;
  if (diff == 0) {
    elt = Spanningelt_new(stage1_oligos,&(*stage1_retrievedp),&(*stage1_positions),
			  &(*stage1_npositions),indexdb,/*partnerp*/false,/*querypos1*/worstpos,
			  /*querypos2*/0,query_lastpos,querylength,plusp);
    leftpos = rightpos = worstpos;
  } else {
    if (first_anchored_p == true) {
      partnerpos = worstpos + diff;
      leftpos = worstpos;
      rightpos = partnerpos;
    } else {
      partnerpos = worstpos - diff;
      leftpos = partnerpos;
      rightpos = worstpos;
    }
    elt = Spanningelt_new(stage1_oligos,&(*stage1_retrievedp),&(*stage1_positions),
			  &(*stage1_npositions),indexdb,/*partnerp*/true,/*querypos1*/partnerpos,
			  /*querypos2*/worstpos,query_lastpos,querylength,plusp);
  }
  spanningset = List_push(NULL,elt);
  *minscore = elt->candidates_score;

  /* Add left positions */
  for (querypos = first; querypos < leftpos; querypos += spansize) {
    elt = Spanningelt_new(stage1_oligos,&(*stage1_retrievedp),&(*stage1_positions),
			  &(*stage1_npositions),indexdb,/*partnerp*/false,/*querypos1*/querypos,
			  /*querypos2*/0,query_lastpos,querylength,plusp);
    spanningset = List_push(spanningset,elt);
    if (elt->candidates_score < *minscore) {
      *minscore = elt->candidates_score;
    }
  }

  /* Add right positions */
  for (querypos = last; querypos > rightpos; querypos -= spansize) {
    elt = Spanningelt_new(stage1_oligos,&(*stage1_retrievedp),&(*stage1_positions),
			  &(*stage1_npositions),indexdb,/*partnerp*/false,/*querypos1*/querypos,
			  /*querypos2*/0,query_lastpos,querylength,plusp);
    spanningset = List_push(spanningset,elt);
    if (elt->candidates_score < *minscore) {
      *minscore = elt->candidates_score;
    }
  }

  return spanningset;
}


List_T
Spanningelt_set (int *minscore, Storedoligomer_T *stage1_oligos, bool **stage1_retrievedp,
		 Univcoord_T ***stage1_positions, int **stage1_npositions,
		 Indexdb_T indexdb, int query_lastpos, int querylength,
		 int mod, bool plusp) {
  int first, last;

  if (mod == 0) {
    first = 0;
  } else {
    first = mod - index1interval;
  }

  if (query_lastpos % index1interval == mod) {
    last = query_lastpos;
  } else if ((query_lastpos + 1) % index1interval == mod) {
    last = query_lastpos + 1;
  } else {
    last = query_lastpos + 2;
  }

  return make_spanningset_aux(&(*minscore),stage1_oligos,&(*stage1_retrievedp),
			      &(*stage1_positions),&(*stage1_npositions),
			      indexdb,query_lastpos,querylength,
			      first,last,plusp);
}


int
Spanningelt_candidates_cmp (const void *a, const void *b) {
  T x = * (T *) a;
  T y = * (T *) b;

  if (x->candidates_score < y->candidates_score) {
    return -1;
  } else if (y->candidates_score < x->candidates_score) {
    return +1;
  } else {
    return 0;
  }
}

int
Spanningelt_pruning_cmp (const void *a, const void *b) {
  T x = * (T *) a;
  T y = * (T *) b;

  if (x->pruning_score < y->pruning_score) {
    return -1;
  } else if (y->pruning_score < x->pruning_score) {
    return +1;
  } else {
    return 0;
  }
}


/* Called only by exact/sub:1 procedures, so need to do Bigendian conversion */
#ifdef WORDS_BIGENDIAN
static int
binary_search (int lowi, int highi, Univcoord_T *positions, Univcoord_T goal) {
  int middlei;

  debug10(printf("entered binary search with lowi=%d, highi=%d, goal=%u\n",lowi,highi,goal));

  while (lowi < highi) {
    middlei = (lowi+highi)/2;
    debug10(printf("  binary: %d:%u %d:%u %d:%u   vs. %u\n",
		   lowi,Bigendian_convert_univcoord(positions[lowi]),middlei,Bigendian_convert_univcoord(positions[middlei]),
		   highi,Bigendian_convert_univcoord(positions[highi]),goal));
    if (goal < Bigendian_convert_univcoord(positions[middlei])) {
      highi = middlei;
    } else if (goal > Bigendian_convert_univcoord(positions[middlei])) {
      lowi = middlei + 1;
    } else {
      debug10(printf("binary search returns %d\n",middlei));
      return middlei;
    }
  }

  debug10(printf("binary search returns %d\n",highi));
  return highi;
}
#else
static int
binary_search (int lowi, int highi, Univcoord_T *positions, Univcoord_T goal) {
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
#endif


/* This procedure needs to eliminate duplicates, which can happen with a SNP-tolerant indexdb */
static Univcoord_T *
compute_intersection (int *ndiagonals, Univcoord_T *positionsa, int diagterma, int npositionsa, 
		      Univcoord_T *positionsb, int diagtermb, int npositionsb) {
  Univcoord_T *diagonals, *positions0, *positions1, local_goal, last_diagonal = 0U, this_diagonal;
  int npositions0, npositions1, delta, j, diagterm;

  if (npositionsa < npositionsb) {
    positions0 = positionsa;
    npositions0 = npositionsa;
    positions1 = positionsb;
    npositions1 = npositionsb;
    diagterm = diagtermb;	/* local_goal based on larger list */
    delta = diagterma - diagtermb; /* list0 + (diagterm0 - diagterm1) = list1 */
  } else {
    positions0 = positionsb;
    npositions0 = npositionsb;
    positions1 = positionsa;
    npositions1 = npositionsa;
    diagterm = diagterma;	/* local_goal based on larger list */
    delta = diagtermb - diagterma; /* list0 + (diagterm0 - diagterm1) = list1 */
  }

  debug(printf("compute_intersection with %d positions <= %d positions.  diagterm %d.\n",
		npositions0,npositions1,diagterm));

  *ndiagonals = 0;
  if (npositions0 == 0) {
    debug(printf("intersection is null\n"));
    return (Univcoord_T *) NULL;
  } else {
    /* Allocate maximum possible size */
    diagonals = (Univcoord_T *) CALLOC(npositions0,sizeof(Univcoord_T));
  }

  while (npositions0 > 0) {
#ifdef WORDS_BIGENDIAN
    local_goal = Bigendian_convert_univcoord(*positions0) + delta;
    debug(printf("intersection list 0: %d:%u => local_goal %u\n",
		 npositions0,Bigendian_convert_univcoord(*positions0),local_goal));
    if (npositions1 > 0 && Bigendian_convert_univcoord(*positions1) < local_goal) {
      j = 1;
      while (j < npositions1 && Bigendian_convert_univcoord(positions1[j]) < local_goal) {
	j <<= 1;		/* gallop by 2 */
      }
      if (j >= npositions1) {
	j = binary_search(j >> 1,npositions1,positions1,local_goal);
      } else {
	j = binary_search(j >> 1,j,positions1,local_goal);
      }
      positions1 += j;
      npositions1 -= j;
    }
#else
    local_goal = (*positions0) + delta;
    debug(printf("intersection list 0: %d:%u => local_goal %u\n",npositions0,*positions0,local_goal));
    if (npositions1 > 0 && *positions1 < local_goal) {
      j = 1;
      while (j < npositions1 && positions1[j] < local_goal) {
	j <<= 1;		/* gallop by 2 */
      }
      if (j >= npositions1) {
	j = binary_search(j >> 1,npositions1,positions1,local_goal);
      } else {
	j = binary_search(j >> 1,j,positions1,local_goal);
      }
      positions1 += j;
      npositions1 -= j;
    }
#endif

#ifdef WORDS_BIGENDIAN
    if (npositions1 <= 0) {
      return diagonals;
    } else if (Bigendian_convert_univcoord(*positions1) == local_goal) {
      /* Found local goal.  Save and advance */
      debug(printf("    intersection list 1: %d:%u  found\n",
		   npositions1,Bigendian_convert_univcoord(*positions1)));
      if ((this_diagonal = local_goal + diagterm) != last_diagonal) {
	diagonals[(*ndiagonals)++] = this_diagonal;
      }
      last_diagonal = this_diagonal;
      ++positions1;
      --npositions1;
    }

#else
    if (npositions1 <= 0) {
      return diagonals;
    } else if ((*positions1) == local_goal) {
      /* Found local goal.  Save and advance */
      debug(printf("    intersection list 1: %d:%u  found\n",npositions1,*positions1));
      if ((this_diagonal = local_goal + diagterm) != last_diagonal) {
	diagonals[(*ndiagonals)++] = this_diagonal;
      }
      last_diagonal = this_diagonal;
      ++positions1;
      --npositions1;
    }
#endif

    ++positions0;
    --npositions0;
  }
  debug(printf("\n"));

  return diagonals;
}


/* This procedure needs to eliminate duplicates, which can happen with a SNP-tolerant indexdb */
static Univcoord_T *
compoundpos_intersect (int *ndiagonals, Univcoord_T *positions0, int diagterm0, int npositions0, 
		       Compoundpos_T compoundpos, int diagterm1) {
  Univcoord_T *diagonals, local_goal, last_local_goal;
  int delta;
  bool emptyp;

  delta = diagterm0 - diagterm1; /* list0 + (diagterm0 - diagterm1) = list1 */

  *ndiagonals = 0;
  if (npositions0 == 0) {
    return (Univcoord_T *) NULL;
  } else {
  /* Could add up compoundpos->npositions to see if we could allocate less memory */
    diagonals = (Univcoord_T *) CALLOC(npositions0,sizeof(Univcoord_T));
  }
  
  last_local_goal = 0U;
  while (npositions0 > 0) {
#ifdef WORDS_BIGENDIAN
    local_goal = Bigendian_convert_univcoord(*positions0) + delta;
#else
    local_goal = (*positions0) + delta;
#endif
    if (local_goal != last_local_goal) {
      debug(printf("Local goal: %u\n",local_goal));
      if (Compoundpos_find(&emptyp,compoundpos,local_goal) == true) {
	/* Found local_goal.  Save and advance */
	debug(printf("Found.  Pushing %u onto positions\n",local_goal));
	diagonals[(*ndiagonals)++] = local_goal + diagterm1; /* = list0 + diagterm0 */
	/* Could potentially advance compoundpos, but let next iteration do it */

      } else if (emptyp == true) {
	/* Empty, so finished */
	debug(printf("Returning %d positions\n",*ndiagonals));
	return diagonals;
      }
    }
    last_local_goal = local_goal;
    
    ++positions0;
    --npositions0;
  }

  return diagonals;
}


#ifdef CHECK
static void
check_diagonals (Univcoord_T *diagonals, int ndiagonals) {
  Univcoord_T last_diagonal = 0U;
  int i;
  
  for (i = 0; i < ndiagonals; i++) {
    if (diagonals[i] == last_diagonal) {
      fprintf(stderr,"Saw repeat of %u\n",diagonals[i]);
      exit(9);
    } else {
      last_diagonal = diagonals[i];
    }
  }
  return;
}
#endif


/* Assumes that we want a single list of diagonals for use in
   generating candidates.  Hides all issues about partner and
   compoundpos. */
Univcoord_T *
Spanningelt_diagonals (int *ndiagonals, T this, int *miss_querypos5, int *miss_querypos3) {
  Univcoord_T *p, *q, last_diagonal;
  int diagterm, i;

  *miss_querypos5 = this->miss_querypos5;
  *miss_querypos3 = this->miss_querypos3;

  debug(printf("** Spanningelt_diagonals\n"));

  if (this->intersection_diagonals != NULL) {
    debug(printf("Intersection diagonals previously computed (%d diagonals)\n",this->intersection_ndiagonals));
    /* Previously computed a result */
    *ndiagonals = this->intersection_ndiagonals;
    check(check_diagonals(this->intersection_diagonals,this->intersection_ndiagonals));
    return this->intersection_diagonals;

  } else if (this->partnerp == false) {
    if (this->compoundpos == NULL) {
      debug(printf("Positions only.  Converting to diagonals\n"));
      /* Positions only.  Convert to diagonals */
      if (this->npositions == 0) {
	this->intersection_diagonals = (Univcoord_T *) NULL;
	*ndiagonals = this->intersection_ndiagonals = 0;
      } else {
	q = this->intersection_diagonals = (Univcoord_T *) CALLOC(this->npositions,sizeof(Univcoord_T));
	p = this->positions;
	diagterm = this->diagterm;

	last_diagonal = 0U;
#ifdef WORDS_BIGENDIAN
	for (i = 0; i < this->npositions; i++) {
	  if (Bigendian_convert_univcoord(*p) != last_diagonal) {
	    *q++ = Bigendian_convert_univcoord(*p) + diagterm;
	  }
	  last_diagonal = Bigendian_convert_univcoord(*p++);
	}
#else
	for (i = 0; i < this->npositions; i++) {
	  if (*p != last_diagonal) {
	    *q++ = *p + diagterm;
	  }
	  last_diagonal = *p++;
	}
#endif

	/* Should be no change in candidates_score or pruning_score */
	*ndiagonals = this->intersection_ndiagonals = (q - this->intersection_diagonals);
      }

      this->intersection_diagonals_reset = this->intersection_diagonals;
      this->intersection_ndiagonals_reset = this->intersection_ndiagonals;

      debug(printf("Returning %p (%d diagonals)\n",this->intersection_diagonals,this->intersection_ndiagonals));
      check(check_diagonals(this->intersection_diagonals,this->intersection_ndiagonals));
      return this->intersection_diagonals;

    } else {
      debug(printf("Compoundpos only.  Converting to diagonals\n"));
      /* Compoundpos */
      this->intersection_diagonals = 
	Indexdb_merge_compoundpos(&this->intersection_ndiagonals,this->compoundpos,this->compoundpos_diagterm);
      *ndiagonals = this->candidates_score = this->pruning_score = this->intersection_ndiagonals;

      this->intersection_diagonals_reset = this->intersection_diagonals;
      this->intersection_ndiagonals_reset = this->intersection_ndiagonals;

      check(check_diagonals(this->intersection_diagonals,this->intersection_ndiagonals));
      return this->intersection_diagonals;
    }

  } else {
    if (this->compoundpos == NULL) {
      debug(printf("Two positions.  Converting to diagonals\n"));
      /* Two positions */
      this->intersection_diagonals = 
	compute_intersection(&this->intersection_ndiagonals,
			     this->partner_positions,this->partner_diagterm,this->partner_npositions,
			     this->positions,this->diagterm,this->npositions);
      *ndiagonals = this->candidates_score = this->pruning_score = this->intersection_ndiagonals;

      this->intersection_diagonals_reset = this->intersection_diagonals;
      this->intersection_ndiagonals_reset = this->intersection_ndiagonals;

      check(check_diagonals(this->intersection_diagonals,this->intersection_ndiagonals));
      return this->intersection_diagonals;

    } else {
      debug(printf("Compoundpos and positions.  Converting to diagonals\n"));
      /* Compoundpos + positions */
      this->intersection_diagonals = 
	compoundpos_intersect(&this->intersection_ndiagonals,this->partner_positions,
			      this->partner_diagterm,this->partner_npositions,
			      this->compoundpos,this->compoundpos_diagterm);
      *ndiagonals = this->candidates_score = this->pruning_score = this->intersection_ndiagonals;

      this->intersection_diagonals_reset = this->intersection_diagonals;
      this->intersection_ndiagonals_reset = this->intersection_ndiagonals;

      check(check_diagonals(this->intersection_diagonals,this->intersection_ndiagonals));
      return this->intersection_diagonals;
    }      
  }
}


