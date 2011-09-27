static char rcsid[] = "$Id: gbuffer.c,v 1.4 2006/11/17 02:40:22 twu Exp $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "gbuffer.h"

#include <stddef.h>
#include <stdlib.h>
#include "mem.h"
#include "pair.h"
#include "sequence.h"


#define T Gbuffer_T
struct T {
  int default_gbufferlen;
  int gbufferlen;		/* Actual gbufferlen */
  int alignmentlen;

  char *chars1;			/* For retrieval of genomic sequence after stage 1 */
  char *chars2;
  char *chars3;

  int *lastGT;			/* For dinucleotide positions in stage 2 */
  int *lastAG;
  int *lastCT;
  int *lastAC;

  int *matchscores;		/* For trimming of alignment in stage 3 */
};


int
Gbuffer_gbufferlen (T this) {
  return this->gbufferlen;
}

int
Gbuffer_alignmentlen (T this) {
  return this->alignmentlen;
}

char *
Gbuffer_chars1 (T this) {
  return this->chars1;
}

char *
Gbuffer_chars2 (T this) {
  return this->chars2;
}

char *
Gbuffer_chars3 (T this) {
  return this->chars3;
}

int *
Gbuffer_lastGT (T this) {
  return this->lastGT;
}

int *
Gbuffer_lastAG (T this) {
  return this->lastAG;
}

int *
Gbuffer_lastCT (T this) {
  return this->lastCT;
}

int *
Gbuffer_lastAC (T this) {
  return this->lastAC;
}

int *
Gbuffer_matchscores (T this) {
  return this->matchscores;
}


static void
free_contents (T this) {
  FREE(this->lastAC);
  FREE(this->lastCT);
  FREE(this->lastAG);
  FREE(this->lastGT);
  FREE(this->chars3);
  FREE(this->chars2);
  FREE(this->chars1);
  return;
}

static void
alloc_contents (T this, int gbufferlen) {

  this->chars1 = (char *) CALLOC(gbufferlen+1,sizeof(char));
  this->chars2 = (char *) CALLOC(gbufferlen+1,sizeof(char));
  this->chars3 = (char *) CALLOC(gbufferlen+1,sizeof(char));
  
  this->lastGT = (int *) CALLOC(gbufferlen-2,sizeof(int));
  this->lastAG = (int *) CALLOC(gbufferlen+1,sizeof(int));
  this->lastCT = (int *) CALLOC(gbufferlen-2,sizeof(int));
  this->lastAC = (int *) CALLOC(gbufferlen+1,sizeof(int));

  this->gbufferlen = gbufferlen;

  return;
}


void
Gbuffer_free (T *old) {
  if (*old) {
    FREE((*old)->matchscores);
    free_contents(*old);
    FREE(*old);
  }
  return;
}

T
Gbuffer_new (int default_gbufferlen) {
  T new = (T) MALLOC(sizeof(*new));

  new->default_gbufferlen = default_gbufferlen;

  /* Maximum alignment length occurs if there is one gap after each
     oligomer in stage 2, which is a 6-mer */
#ifdef PMAP
  new->alignmentlen = 3*MAXSEQLEN + (MATCHESPERGAP)*3*MAXSEQLEN/6;
#else
  new->alignmentlen = MAXSEQLEN + (MATCHESPERGAP)*MAXSEQLEN/6;
#endif
  new->matchscores = (int *) CALLOC(new->alignmentlen,sizeof(int));;

  alloc_contents(new,default_gbufferlen);

  return new;
}

void
Gbuffer_check_alloc (T this, int gbufferlen) {

  if (gbufferlen > this->gbufferlen) {
    free_contents(this);
    alloc_contents(this,gbufferlen);

  } else if (this->gbufferlen == this->default_gbufferlen) {
    /* Do nothing; already at default */

  } else if (gbufferlen > this->default_gbufferlen) {
    /* Do nothing; keep currently enlarged buffer */
    
  } else {
    /* Reduce buffer to default */
    free_contents(this);
    alloc_contents(this,this->default_gbufferlen);
  }

  return;
}

