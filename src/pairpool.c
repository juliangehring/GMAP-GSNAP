static char rcsid[] = "$Id: pairpool.c,v 1.28 2005/10/06 20:42:05 twu Exp $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "pairpool.h"
#include <stdio.h>
#include <stdlib.h>
#include "mem.h"
#include "comp.h"
#include "pairdef.h"
#include "listdef.h"

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


static struct Pair_T *
add_new_pairchunk (T this) {
  struct Pair_T *chunk;

  chunk = (struct Pair_T *) MALLOC(CHUNKSIZE*sizeof(struct Pair_T));
  this->pairchunks = List_push(this->pairchunks,(void *) chunk);
  debug1(printf("Adding a new chunk of pairs.  Ptr for pair %d is %p\n",
		this->npairs,chunk));

  this->npairs += CHUNKSIZE;
  return chunk;
}

static struct List_T *
add_new_listcellchunk (T this) {
  struct List_T *chunk;

  chunk = (struct List_T *) MALLOC(CHUNKSIZE*sizeof(struct List_T));
  this->listcellchunks = List_push(this->listcellchunks,(void *) chunk);
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
  new->pairptr = add_new_pairchunk(new);

  new->nlistcells = 0;
  new->listcellctr = 0;
  new->listcellchunks = NULL;
  new->listcellptr = add_new_listcellchunk(new);

  return new;
}

void
Pairpool_reset (T this) {
  this->pairctr = 0;
  this->listcellctr = 0;
  return;
}

List_T
Pairpool_push (List_T list, T this, int querypos, int genomepos, char cdna, char comp, char genome) {
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
  pair->aamarker = false;
  pair->cdna = cdna;
  pair->comp = comp;
  pair->genome = genome;
  pair->shortexonp = false;
  switch (comp) {
  case MATCH_COMP: case DYNPROG_MATCH_COMP: case AMBIGUOUS_COMP: case MISMATCH_COMP: case INDEL_COMP: case SHORTGAP_COMP:
    pair->gapp = false; break;
  default: pair->gapp = true;
  }

  debug(
	printf("Creating: %d %d %c %c %c\n",
	       pair->querypos,pair->genomepos,pair->cdna,pair->comp,pair->genome);
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
Pairpool_push_existing (List_T list, T this, Pair_T pair) {
  List_T listcell;
  List_T p;
  int n;

  debug(
	printf("Pushing: %d %d %c %c %c",
	       pair->querypos,pair->genomepos,pair->cdna,pair->comp,pair->genome);
	if (pair->gapp == true) {
	  printf(" (gap)");
	}
	printf("\n");
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
Pairpool_pop (List_T list, Pair_T *x) {
  List_T head;

  if (list != NULL) {
    head = list->rest;
    *x = (Pair_T) list->first;
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
	  printf("Transferring: %d %d %c %c %c\n",
		 pair->querypos,pair->genomepos,pair->cdna,pair->comp,pair->genome);
	  );
    next = p->rest;
    p->rest = dest;
    dest = p;
  }
  return dest;
}

List_T
Pairpool_transfer_copy (List_T dest, List_T source, T this) {
  Pair_T pair;

  while (source != NULL) {
    source = Pairpool_pop(source,&pair);
    dest = Pairpool_push(dest,this,pair->querypos,pair->genomepos,pair->cdna,pair->comp,pair->genome);
    debug(printf("Copying: %d %d %c %c %c\n",
		 pair->querypos,pair->genomepos,pair->cdna,pair->comp,pair->genome));
  }
  return dest;
}


List_T
Pairpool_transfer_bounded (List_T dest, List_T source, int minpos, int maxpos) {
  List_T p, next;
  Pair_T pair;

  for (p = source; p != NULL; p = next) {
    debug(
	  pair = List_head(p);
	  if (pair->cdna == '\0' || pair->genome == '\0') {
	    abort();
	  }
	  printf("Transferring: %d %d %c %c %c\n",
		 pair->querypos,pair->genomepos,pair->cdna,pair->comp,pair->genome);
	  );
    pair = List_head(p);
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

#if 0
List_T
Pairpool_transfer_copy_bounded (List_T dest, List_T source, T this, int minpos, int maxpos) {
  Pair_T pair;

  while (source != NULL) {
    source = Pairpool_pop(source,&pair);
    if (pair->querypos >= minpos && pair->querypos <= maxpos) {
      dest = Pairpool_push(dest,this,pair->querypos,pair->genomepos,pair->cdna,pair->comp,pair->genome);
      debug(printf("Copying: %d %d %c %c %c\n",
		   pair->querypos,pair->genomepos,pair->cdna,pair->comp,pair->genome));
    }
  }
  return dest;
}
#endif


