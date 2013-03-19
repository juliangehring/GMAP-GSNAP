static char rcsid[] = "$Id: pairpool.c 82070 2012-12-19 21:42:59Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "pairpool.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>		/* For memcpy */
#include "assert.h"
#include "mem.h"
#include "comp.h"
#include "pairdef.h"
#include "listdef.h"
#include "intron.h"

#define CHUNKSIZE 20000

#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/* For mechanics of memory allocation and deallocation */
#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif

/* For popping */
#ifdef DEBUG2
#define debug2(x) x
#else
#define debug2(x)
#endif


#define T Pairpool_T
struct T {
  int npairs;
  int pairctr;
  struct Pair_T *pairptr;
  List_T pairchunks;

  int nlistcells;
  int listcellctr;
  struct List_T *listcellptr;
  List_T listcellchunks;
};

void
Pairpool_free (T *old) {
  List_T p;
  struct Pair_T *pairptr;
  struct List_T *listcellptr;

  if (*old) {
    for (p = (*old)->pairchunks; p != NULL; p = List_next(p)) {
      pairptr = (struct Pair_T *) List_head(p);
      FREE(pairptr);
    }
    List_free(&(*old)->pairchunks);
    for (p = (*old)->listcellchunks; p != NULL; p = List_next(p)) {
      listcellptr = (struct List_T *) List_head(p);
      FREE(listcellptr);
    }
    List_free(&(*old)->listcellchunks);
    FREE(*old);
  }
  return;
}

void
Pairpool_free_memory (T this) {
  List_T p;
  struct Pair_T *pairptr;
  struct List_T *listcellptr;

  for (p = this->pairchunks; p != NULL; p = List_next(p)) {
    pairptr = (struct Pair_T *) List_head(p);
    FREE_KEEP(pairptr);
  }
  List_free_keep(&this->pairchunks);
  for (p = this->listcellchunks; p != NULL; p = List_next(p)) {
    listcellptr = (struct List_T *) List_head(p);
    FREE_KEEP(listcellptr);
  }
  List_free_keep(&this->listcellchunks);

  this->npairs = 0;
  this->pairctr = 0;
  this->pairchunks = NULL;
  /* this->pairptr = add_new_pairchunk(this); */

  this->nlistcells = 0;
  this->listcellctr = 0;
  this->listcellchunks = NULL;
  /* this->listcellptr = add_new_listcellchunk(this); */

  return;
}


void
Pairpool_report_memory (T this) {
  printf("Pairpool has %d pairchunks and %d listcellchunks\n",
	 List_length(this->pairchunks),List_length(this->listcellchunks));
  return;
}


static struct Pair_T *
add_new_pairchunk (T this) {
  struct Pair_T *chunk;

  chunk = (struct Pair_T *) MALLOC_KEEP(CHUNKSIZE*sizeof(struct Pair_T));
  this->pairchunks = List_push_keep(this->pairchunks,(void *) chunk);
  debug1(printf("Adding a new chunk of pairs.  Ptr for pair %d is %p\n",
		this->npairs,chunk));

  this->npairs += CHUNKSIZE;

  return chunk;
}

static struct List_T *
add_new_listcellchunk (T this) {
  struct List_T *chunk;

  chunk = (struct List_T *) MALLOC_KEEP(CHUNKSIZE*sizeof(struct List_T));
  this->listcellchunks = List_push_keep(this->listcellchunks,(void *) chunk);
  debug1(printf("Adding a new chunk of listcells.  Ptr for listcell %d is %p\n",
	       this->nlistcells,chunk));

  this->nlistcells += CHUNKSIZE;

  return chunk;
}

T
Pairpool_new (void) {
  T new = (T) MALLOC(sizeof(*new));

  new->npairs = 0;
  new->pairctr = 0;
  new->pairchunks = NULL;
  /* new->pairptr = add_new_pairchunk(new); */

  new->nlistcells = 0;
  new->listcellctr = 0;
  new->listcellchunks = NULL;
  /* new->listcellptr = add_new_listcellchunk(new); */

  return new;
}

void
Pairpool_reset (T this) {
  this->pairctr = 0;
  this->listcellctr = 0;
  return;
}

/* gapp should be false for the following comps: MATCH_COMP,
   DYNPROG_MATCH_COMP, AMBIGUOUS_COMP, MISMATCH_COMP, INDEL_COMP,
   SHORTGAP_COMP */

List_T
Pairpool_push (List_T list, T this, int querypos, int genomepos, char cdna, char comp,
	       char genome, char genomealt, int dynprogindex) {
  List_T listcell;
  Pair_T pair;
  List_T p;
  int n;

  assert(querypos >= 0);

  if (this->pairctr >= this->npairs) {
    this->pairptr = add_new_pairchunk(this);
  } else if ((this->pairctr % CHUNKSIZE) == 0) {
    for (n = this->npairs - CHUNKSIZE, p = this->pairchunks;
	 n > this->pairctr; p = p->rest, n -= CHUNKSIZE) ;
    this->pairptr = (struct Pair_T *) p->first;
    debug1(printf("Located pair %d at %p\n",this->pairctr,this->pairptr));
  }    
  pair = this->pairptr++;
  this->pairctr++;

  pair->querypos = querypos;
  pair->genomepos = genomepos;
  pair->aapos = 0;
  pair->aaphase_g = -1;
  pair->aaphase_e = -1;
  pair->cdna = cdna;
  pair->comp = comp;
  pair->genome = genome;
  pair->genomealt = genomealt;
  pair->dynprogindex = dynprogindex;

  pair->aa_g = ' ';
  pair->aa_e = ' ';
  pair->shortexonp = false;
  pair->gapp = false;
  pair->knowngapp = false;
  pair->introntype = NONINTRON;
  if (comp == EXTRAEXON_COMP) {
    pair->extraexonp = true;
  } else {
    pair->extraexonp = false;
  }
  
  pair->queryjump = 0;
  pair->genomejump = 0;

  pair->state = GOOD;
  pair->protectedp = false;
  pair->disallowedp = false;
  pair->donor_prob = 0.0;
  pair->acceptor_prob = 0.0;
  pair->end_intron_p = false;

  debug(
	printf("Creating %p: %d %d %c %c %c\n",
	       pair,pair->querypos,pair->genomepos,pair->cdna,pair->comp,pair->genome);
	);

	if (this->listcellctr >= this->nlistcells) {
    this->listcellptr = add_new_listcellchunk(this);
  } else if ((this->listcellctr % CHUNKSIZE) == 0) {
    for (n = this->nlistcells - CHUNKSIZE, p = this->listcellchunks;
	 n > this->listcellctr; p = p->rest, n -= CHUNKSIZE) ;
    this->listcellptr = (struct List_T *) p->first;
    debug1(printf("Located listcell %d at %p\n",this->listcellctr,this->listcellptr));
  }
  listcell = this->listcellptr++;
  this->listcellctr++;

  listcell->first = (void *) pair;
  listcell->rest = list;

  return listcell;
}


List_T
Pairpool_push_copy (List_T list, T this, Pair_T orig) {
  List_T listcell;
  Pair_T pair;
  List_T p;
  int n;

  if (this->pairctr >= this->npairs) {
    this->pairptr = add_new_pairchunk(this);
  } else if ((this->pairctr % CHUNKSIZE) == 0) {
    for (n = this->npairs - CHUNKSIZE, p = this->pairchunks;
	 n > this->pairctr; p = p->rest, n -= CHUNKSIZE) ;
    this->pairptr = (struct Pair_T *) p->first;
    debug1(printf("Located pair %d at %p\n",this->pairctr,this->pairptr));
  }    
  pair = this->pairptr++;
  this->pairctr++;

  memcpy(pair,orig,sizeof(struct Pair_T));

  debug(
	printf("Copying %p: %d %d %c %c %c\n",
	       pair,pair->querypos,pair->genomepos,pair->cdna,pair->comp,pair->genome);
	);

	if (this->listcellctr >= this->nlistcells) {
    this->listcellptr = add_new_listcellchunk(this);
  } else if ((this->listcellctr % CHUNKSIZE) == 0) {
    for (n = this->nlistcells - CHUNKSIZE, p = this->listcellchunks;
	 n > this->listcellctr; p = p->rest, n -= CHUNKSIZE) ;
    this->listcellptr = (struct List_T *) p->first;
    debug1(printf("Located listcell %d at %p\n",this->listcellctr,this->listcellptr));
  }
  listcell = this->listcellptr++;
  this->listcellctr++;

  listcell->first = (void *) pair;
  listcell->rest = list;

  return listcell;
}


List_T
Pairpool_push_gapalign (List_T list, T this, int querypos, int genomepos, char cdna, char comp,
			int introntype, char genome, char genomealt, bool extraexonp) {
  List_T listcell;
  Pair_T pair;
  List_T p;
  int n;

  if (this->pairctr >= this->npairs) {
    this->pairptr = add_new_pairchunk(this);
  } else if ((this->pairctr % CHUNKSIZE) == 0) {
    for (n = this->npairs - CHUNKSIZE, p = this->pairchunks;
	 n > this->pairctr; p = p->rest, n -= CHUNKSIZE) ;
    this->pairptr = (struct Pair_T *) p->first;
    debug1(printf("Located pair %d at %p\n",this->pairctr,this->pairptr));
  }    
  pair = this->pairptr++;
  this->pairctr++;

  pair->querypos = querypos;
  pair->genomepos = genomepos;
  pair->aapos = 0;
  pair->aaphase_g = -1;
  pair->aaphase_e = -1;
  pair->cdna = cdna;
  pair->comp = comp;
  pair->genome = genome;
  pair->genomealt = genomealt;
  pair->dynprogindex = 0;

  pair->aa_g = ' ';
  pair->aa_e = ' ';
  pair->shortexonp = false;
  pair->gapp = true;
  pair->knowngapp = false;
  pair->introntype = introntype;
  pair->extraexonp = extraexonp;
  
  pair->queryjump = 0;
  pair->genomejump = 0;

  pair->state = GOOD;
  pair->protectedp = false;
  pair->disallowedp = false;
  pair->donor_prob = 0.0;
  pair->acceptor_prob = 0.0;
  pair->end_intron_p = false;

  debug(
	printf("Creating %p: %d %d %c %c %c introntype %d\n",
	       pair,pair->querypos,pair->genomepos,pair->cdna,pair->comp,pair->genome,pair->introntype);
	);

	if (this->listcellctr >= this->nlistcells) {
    this->listcellptr = add_new_listcellchunk(this);
  } else if ((this->listcellctr % CHUNKSIZE) == 0) {
    for (n = this->nlistcells - CHUNKSIZE, p = this->listcellchunks;
	 n > this->listcellctr; p = p->rest, n -= CHUNKSIZE) ;
    this->listcellptr = (struct List_T *) p->first;
    debug1(printf("Located listcell %d at %p\n",this->listcellctr,this->listcellptr));
  }
  listcell = this->listcellptr++;
  this->listcellctr++;

  listcell->first = (void *) pair;
  listcell->rest = list;

  return listcell;
}

List_T
Pairpool_push_gapholder (List_T list, T this, int queryjump, int genomejump, bool knownp) {
  List_T listcell;
  Pair_T pair;
  List_T p;
  int n;

  if (this->pairctr >= this->npairs) {
    this->pairptr = add_new_pairchunk(this);
  } else if ((this->pairctr % CHUNKSIZE) == 0) {
    for (n = this->npairs - CHUNKSIZE, p = this->pairchunks;
	 n > this->pairctr; p = p->rest, n -= CHUNKSIZE) ;
    this->pairptr = (struct Pair_T *) p->first;
    debug1(printf("Located pair %d at %p\n",this->pairctr,this->pairptr));
  }    
  pair = this->pairptr++;
  this->pairctr++;

  pair->querypos = -1;
  pair->genomepos = -1;

  pair->aapos = 0;
  pair->aaphase_g = -1;
  pair->aaphase_e = -1;
  pair->cdna = ' ';
  pair->comp = ' ';
  pair->genome = ' ';
  pair->genomealt = ' ';
  pair->dynprogindex = 0;

  pair->aa_g = ' ';
  pair->aa_e = ' ';
  pair->shortexonp = false;
  pair->gapp = true;
  if (knownp == true) {
    pair->knowngapp = true;
    pair->donor_prob = 2.0;
    pair->acceptor_prob = 2.0;
  } else {
    pair->knowngapp = false;
    pair->donor_prob = 0.0;
    pair->acceptor_prob = 0.0;
  }
  pair->introntype = NONINTRON;
  pair->extraexonp = false;

  pair->queryjump = queryjump;
  pair->genomejump = genomejump;

  pair->state = GOOD;
  pair->protectedp = false;
  pair->disallowedp = false;
  pair->end_intron_p = false;

  debug(printf("Creating gap %p, queryjump=%d, genomejump=%d\n",pair,queryjump,genomejump));

  if (this->listcellctr >= this->nlistcells) {
    this->listcellptr = add_new_listcellchunk(this);
  } else if ((this->listcellctr % CHUNKSIZE) == 0) {
    for (n = this->nlistcells - CHUNKSIZE, p = this->listcellchunks;
	 n > this->listcellctr; p = p->rest, n -= CHUNKSIZE) ;
    this->listcellptr = (struct List_T *) p->first;
    debug1(printf("Located listcell %d at %p\n",this->listcellctr,this->listcellptr));
  }
  listcell = this->listcellptr++;
  this->listcellctr++;

  listcell->first = (void *) pair;
  listcell->rest = list;

  return listcell;
}

List_T
Pairpool_push_existing (List_T list, T this, Pair_T pair) {
  List_T listcell;
  List_T p;
  int n;

  debug(
	Pair_T head;
	if (pair->gapp == true) {
	  printf("Pushing gap %p: queryjump=%d, genomejump=%d onto ",
		 pair,pair->queryjump,pair->genomejump);
	} else {
	  printf("Pushing %p: %d %d %c %c %c onto ",
		 pair,pair->querypos,pair->genomepos,pair->cdna,pair->comp,pair->genome);
	}
	if (list == NULL) {
	  printf("NULL\n");
	} else {
	  head = list->first;
	  if (head->gapp == true) {
	    printf("gap %p: queryjump=%d, genomejump=%d\n",
		   head,head->queryjump,head->genomejump);
	  } else {
	    printf("%p: %d %d %c %c %c\n",
		   head,head->querypos,head->genomepos,head->cdna,head->comp,head->genome);
	  }
	}
	);
  
  if (this->listcellctr >= this->nlistcells) {
    this->listcellptr = add_new_listcellchunk(this);
  } else if ((this->listcellctr % CHUNKSIZE) == 0) {
    for (n = this->nlistcells - CHUNKSIZE, p = this->listcellchunks;
	 n > this->listcellctr; p = p->rest, n -= CHUNKSIZE) ;
    this->listcellptr = (struct List_T *) p->first;
    debug1(printf("Located listcell %d at %p\n",this->listcellctr,this->listcellptr));
  }
  listcell = this->listcellptr++;
  this->listcellctr++;

  listcell->first = (void *) pair;
  listcell->rest = list;

  return listcell;
}


/* Note: this does not free the list cell */
List_T
Pairpool_pop (List_T list, Pair_T *x) {
  List_T head;

  if (list != NULL) {
    head = list->rest;
    *x = (Pair_T) list->first;
    debug2(
	   if ((*x)->gapp == true) {
	     printf("Popping gap: queryjump=%d, genomejump=%d\n",
		    (*x)->queryjump,(*x)->genomejump);
	   } else {
	     printf("Popping: %d %d %c %c %c\n",
		    (*x)->querypos,(*x)->genomepos,(*x)->cdna,(*x)->comp,(*x)->genome);
	   }
	   );
    return head;
  } else {
    return list;
  }
}


List_T
Pairpool_transfer (List_T dest, List_T source) {
  List_T p, next;
#ifdef DEBUG
  Pair_T pair;
#endif

  for (p = source; p != NULL; p = next) {
    debug(
	  pair = List_head(p);
	  if (pair->cdna == '\0' || pair->genome == '\0') {
	    abort();
	  }
	  if (pair->gapp) {
	    printf("Transferring gap %p: queryjump=%d, genomejump=%d\n",
		   pair,pair->queryjump,pair->genomejump);
	  } else {
	    printf("Transferring %p: %d %d %c %c %c\n",
		   pair,pair->querypos,pair->genomepos,pair->cdna,pair->comp,pair->genome);
	  }
	  );
    next = p->rest;
    p->rest = dest;
    dest = p;
  }
  return dest;
}

List_T
Pairpool_transfer_n (List_T dest, List_T source, int n) {
  List_T p, next;
#ifdef DEBUG
  Pair_T pair;
#endif

  for (p = source; p != NULL && --n >= 0; p = next) {
    debug(
	  pair = List_head(p);
	  if (pair->cdna == '\0' || pair->genome == '\0') {
	    abort();
	  }
	  if (pair->gapp) {
	    printf("Transferring gap %p: queryjump=%d, genomejump=%d\n",
		   pair,pair->queryjump,pair->genomejump);
	  } else {
	    printf("Transferring %p: %d %d %c %c %c\n",
		   pair,pair->querypos,pair->genomepos,pair->cdna,pair->comp,pair->genome);
	  }
	  );
    next = p->rest;
    p->rest = dest;
    dest = p;
  }
  return dest;
}

int
Pairpool_count_bounded (int *nstart, List_T source, int minpos, int maxpos) {
  int npairs = 0;
  List_T p, next;
  Pair_T pair;

  *nstart = 0;
  for (p = source; p != NULL; p = next) {
    pair = List_head(p);
    next = p->rest;
    if (pair->querypos < minpos) {
      *nstart += 1;
    } else if (pair->querypos < maxpos) {
      npairs++;
    } else {
      p = NULL;			/* Terminate transfer */
    }
  }
  return npairs;
}


/* Note: This code is designed to handle source, which may still have
   gaps with querypos undefined */
List_T
Pairpool_clip_bounded (List_T source, int minpos, int maxpos) {
  List_T dest, *prev, p;
  Pair_T pair;
  int starti = -1, endi = -1, i;

  for (p = source, i = 0; p != NULL; p = p->rest, i++) {
    pair = (Pair_T) List_head(p);
    if (pair->querypos == minpos) {
      starti = i;		/* Advances in case of ties */
    } else if (pair->querypos > minpos && starti < 0) {
      starti = i;		/* Handles case where minpos was skipped */
    }

    if (pair->querypos == maxpos && endi < 0) {
      endi = i + 1;		/* Does not advance in case of tie */
    } else if (pair->querypos > maxpos && endi < 0) {
      endi = i;	   /* Handles case where maxpos was skipped */
    }
  }

  if (starti < 0) {
    starti = 0;
  }
  if (endi < 0) {
    endi = i;
  }

  p = source;
  i = 0;
  while (i < starti) {
    p = p->rest;
    i++;
  }

  dest = p;
  prev = &p->rest;
  while (i < endi) {
    prev = &p->rest;
    p = p->rest;
    i++;
  }

  *prev = NULL;		/* Clip rest of list */
  return dest;
}


#if 0
List_T
Pairpool_transfer_bounded (List_T dest, List_T source, int minpos, int maxpos) {
  List_T p, next;
  Pair_T pair;

  for (p = source; p != NULL; p = next) {
    debug(
	  pair = (Pair_T) List_head(p);
	  if (pair->cdna == '\0' || pair->genome == '\0') {
	    abort();
	  }
	  printf("Transferring %p: %d %d %c %c %c\n",
		 pair,pair->querypos,pair->genomepos,pair->cdna,pair->comp,pair->genome);
	  );
    pair = (Pair_T) List_head(p);
    next = p->rest;
    if (pair->querypos == minpos) {
      if (dest != NULL) {
	/* Pop last querypos off the stack, because we want only one of them */
	dest = dest->rest;
      }
      p->rest = dest;
      dest = p;
    } else if (pair->querypos == maxpos) {
      p->rest = dest;
      dest = p;
      p = NULL;			/* Terminate transfer */
    } else if (pair->querypos > minpos && pair->querypos < maxpos) {
      p->rest = dest;
      dest = p;
    }
  }

  return dest;
}
#endif


/* Originally prohibited copying of gaps */
List_T
Pairpool_copy (List_T source, T this) {
  List_T dest = NULL;

  while (source != NULL) {
    dest = Pairpool_push_copy(dest,this,/*orig*/source->first);
    source = source->rest;
  }
  return List_reverse(dest);
}


struct Pair_T *
Pairpool_copy_array (struct Pair_T *source, int npairs) {
  struct Pair_T *dest;

  dest = (struct Pair_T *) CALLOC_OUT(npairs,sizeof(struct Pair_T));
  memcpy(dest,source,npairs*sizeof(struct Pair_T));
  return dest;
}


