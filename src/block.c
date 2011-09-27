static char rcsid[] = "$Id: block.c,v 1.46 2005/07/19 17:54:53 twu Exp $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "block.h"

#include <stdio.h>
#include "mem.h"
#ifdef PMAP
#include "oligop.h"
#else
#include "oligo.h"
#endif

#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

#define T Block_T
struct T {
  Reader_T reader;

  cDNAEnd_T cdnaend;

  int last_querypos;
  Oligostate_T last_state;

#ifdef PMAP
  unsigned int aaindex;
#else
  Storedoligomer_T forward;
  Storedoligomer_T revcomp;
#endif
};

int
Block_querypos (T this) {
  return this->last_querypos;
}

#ifdef PMAP
unsigned int
Block_aaindex (T this) {
  return this->aaindex;
}
#else
Storedoligomer_T
Block_forward (T this) {
  return this->forward;
}

Storedoligomer_T
Block_revcomp (T this) {
  return this->revcomp;
}
#endif


extern void
Block_reset_ends (T this) {
  Reader_reset_ends(this->reader);

  if (this->cdnaend == FIVE) {
    this->last_querypos = Reader_querystart(this->reader);
  } else {
    this->last_querypos = Reader_queryend(this->reader);
  }
  this->last_state = INIT;

#ifdef PMAP  
  this->aaindex = 0U;
#else
  this->forward = 0U;
  this->revcomp = 0U;
#endif

  return;
}


T
Block_new (cDNAEnd_T cdnaend, Reader_T reader) {
  T new = (T) MALLOC(sizeof(*new));

  new->reader = reader;
  new->cdnaend = cdnaend;

  if (cdnaend == FIVE) {
    new->last_querypos = Reader_startpos(reader);
  } else if (cdnaend == THREE) {
    new->last_querypos = Reader_endpos(reader);
  }
  new->last_state = INIT;

#ifdef PMAP  
  new->aaindex = 0U;
#else
  new->forward = 0U;
  new->revcomp = 0U;
#endif

  return new;
}

void
Block_free (T *old) {
  if (*old) {
    FREE(*old);
  }
  return;
}



bool
Block_next (T this) {
  if (this->last_state == DONE) {
    return false;
  } else {

#ifdef PMAP
    this->last_state = Oligo_next(this->last_state,&this->last_querypos,
				  &this->aaindex,this->reader,this->cdnaend);
    debug(printf("Block has aaindex %u at querypos %d\n",
		 this->aaindex,this->last_querypos));
#else
    this->last_state = Oligo_next(this->last_state,&this->last_querypos,
				  &this->forward,&this->revcomp,
				  this->reader,this->cdnaend);
    debug(printf("Block has oligo %08X at querypos %d\n",
		 this->forward,this->last_querypos));
#endif

    if (this->last_state == DONE) {
      return false;
    } else {
      return true;
    }
  }
}

bool
Block_skip (T this, int nskip) {
  if (this->last_state == DONE) {
    return false;
  } else {

#ifdef PMAP
    this->last_state = Oligo_skip(this->last_state,&this->last_querypos,
				  &this->aaindex,this->reader,this->cdnaend,nskip);
    debug(printf("Block has aaindex %u at querypos %d\n",
		 this->aaindex,this->last_querypos));
#else
    this->last_state = Oligo_skip(this->last_state,&this->last_querypos,
				  &this->forward,&this->revcomp,
				  this->reader,this->cdnaend,nskip);
    debug(printf("Block has oligo %08X at querypos %d\n",
		 this->forward,this->last_querypos));
#endif

    if (this->last_state == DONE) {
      return false;
    } else {
      return true;
    }
  }
}


/* Returns querypos */
#ifdef PMAP
int
Block_process_oligo (Genomicpos_T **fwdpositions, int *nfwdhits, 
		     Genomicpos_T **revpositions, int *nrevhits,
		     T this, Indexdb_T indexdb_fwd, Indexdb_T indexdb_rev) {

  *nfwdhits = Oligo_lookup(&(*fwdpositions),indexdb_fwd,this->aaindex);
  *nrevhits = Oligo_lookup(&(*revpositions),indexdb_rev,this->aaindex);

  return this->last_querypos;
}

#else

int
Block_process_oligo (Genomicpos_T **fwdpositions, int *nfwdhits,
		     Genomicpos_T **revpositions, int *nrevhits,
		     T this, Indexdb_T indexdb) {

  *nfwdhits = Oligo_lookup(&(*fwdpositions),indexdb,this->forward);
  *nrevhits = Oligo_lookup(&(*revpositions),indexdb,this->revcomp);

  return this->last_querypos;
}

#endif
