static char rcsid[] = "$Id: intpool.c,v 1.1 2007/02/05 07:45:55 twu Exp $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "intpool.h"
#include <stdio.h>
#include <stdlib.h>
#include "mem.h"
#include "comp.h"
#include "intlistdef.h"
#include "listdef.h"
#include "list.h"

#define CHUNKSIZE 16384

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


#define T Intpool_T
struct T {
  int nintlistcells;
  int intlistcellctr;
  struct Intlist_T *intlistcellptr;
  List_T chunks;
};

void
Intpool_free (T *old) {
  List_T p;
  struct Intlist_T *intlistcellptr;

  if (*old) {
    for (p = (*old)->chunks; p != NULL; p = List_next(p)) {
      intlistcellptr = (struct Intlist_T *) List_head(p);
      FREE(intlistcellptr);
    }
    List_free(&(*old)->chunks);
    FREE(*old);
  }
  return;
}


static struct Intlist_T *
add_new_chunk (T this) {
  struct Intlist_T *chunk;

  chunk = (struct Intlist_T *) MALLOC(CHUNKSIZE*sizeof(struct Intlist_T));
  this->chunks = List_push(this->chunks,(void *) chunk);
  debug1(printf("Adding a new chunk of listcells.  Ptr for listcell %d is %p\n",
	       this->nintlistcells,chunk));

  this->nintlistcells += CHUNKSIZE;
  return chunk;
}

T
Intpool_new (void) {
  T new = (T) MALLOC(sizeof(*new));

  new->nintlistcells = 0;
  new->intlistcellctr = 0;
  new->chunks = NULL;
  /* new->listcellptr = add_new_listcellchunk(new); */

  return new;
}

void
Intpool_reset (T this) {
  this->intlistcellctr = 0;
  return;
}

Intlist_T
Intpool_push (Intlist_T intlist, T this, int integer) {
  Intlist_T intlistcell;
  int *integerptr;
  List_T p;
  int n;

  if (this->intlistcellctr >= this->nintlistcells) {
    this->intlistcellptr = add_new_chunk(this);
  } else if ((this->intlistcellctr % CHUNKSIZE) == 0) {
    for (n = this->nintlistcells - CHUNKSIZE, p = this->chunks;
	 n > this->intlistcellctr; p = p->rest, n -= CHUNKSIZE) ;
    this->intlistcellptr = (struct Intlist_T *) p->first;
    debug1(printf("Located listcell %d at %p\n",this->intlistcellctr,this->intlistcellptr));
  }
  intlistcell = this->intlistcellptr++;
  this->intlistcellctr++;

  intlistcell->first = integer;
  intlistcell->rest = intlist;

  return intlistcell;
}

Intlist_T
Intpool_pop (Intlist_T intlist, int *integer) {
  Intlist_T head;

  if (intlist != NULL) {
    head = intlist->rest;
    *integer = intlist->first;
    return head;
  } else {
    return intlist;
  }
}


