static char rcsid[] = "$Id: params.c,v 1.62 2005/07/21 16:54:49 twu Exp $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "params.h"
#include "mem.h"


#define T Params_T
struct T {
  Genome_T genome;
  IIT_T altstrain_iit;
#ifdef PMAP
  Indexdb_T indexdb_fwd;
  Indexdb_T indexdb_rev;
#else
  Indexdb_T indexdb;
#endif
  IIT_T chromosome_iit;
  Chrsubset_T chrsubset;
  IIT_T contig_iit;
  IIT_T map_iit;
  int maxextension;

  int stuttercycles;
  int stutterhits;
  int indexsize;
  int maxpeelback;
  int sufflookback;
  int nsufflookback;
  int nullgap;
  int extramaterial_end;
  int extramaterial_paired;
  int extraband_single;
  int extraband_end;
  int extraband_paired;
};


Genome_T
Params_genome (T this) {
  return this->genome;
}

IIT_T
Params_altstrain_iit (T this) {
  return this->altstrain_iit;
}

char *
Params_refstrain (T this) {
  return IIT_typestring(this->altstrain_iit,0);
}

#ifdef PMAP
Indexdb_T
Params_indexdb_fwd (T this) {
  return this->indexdb_fwd;
}

Indexdb_T
Params_indexdb_rev (T this) {
  return this->indexdb_rev;
}
#else
Indexdb_T
Params_indexdb (T this) {
  return this->indexdb;
}
#endif

IIT_T
Params_chromosome_iit (T this) {
  return this->chromosome_iit;
}

Chrsubset_T
Params_chrsubset (T this) {
  return this->chrsubset;
}

IIT_T
Params_contig_iit (T this) {
  return this->contig_iit;
}

IIT_T
Params_map_iit (T this) {
  return this->map_iit;
}

int
Params_maxextension (T this) {
  return this->maxextension;
}

int
Params_stuttercycles (T this) {
  return this->stuttercycles;
}

int
Params_stutterhits (T this) {
  return this->stutterhits;
}

int
Params_indexsize (T this) {
  return this->indexsize;
}

int
Params_maxpeelback (T this) {
  return this->maxpeelback;
}

int
Params_sufflookback (T this) {
  return this->sufflookback;
}

int
Params_nsufflookback (T this) {
  return this->nsufflookback;
}

int
Params_nullgap (T this) {
  return this->nullgap;
}

int
Params_extramaterial_end (T this) {
  return this->extramaterial_end;
}

int
Params_extramaterial_paired (T this) {
  return this->extramaterial_paired;
}

int
Params_extraband_single (T this) {
  return this->extraband_single;
}

int
Params_extraband_end (T this) {
  return this->extraband_end;
}

int
Params_extraband_paired (T this) {
  return this->extraband_paired;
}


T
Params_new (Genome_T genome, IIT_T altstrain_iit, 
#ifdef PMAP
	    Indexdb_T indexdb_fwd,
	    Indexdb_T indexdb_rev,
#else
	    Indexdb_T indexdb, 
#endif
	    IIT_T chromosome_iit, Chrsubset_T chrsubset, IIT_T contig_iit, IIT_T map_iit, 
	    int maxextension, int stuttercycles, int stutterhits, int indexsize,
	    int maxpeelback, int sufflookback, int nsufflookback, int nullgap, 
	    int extramaterial_end, int extramaterial_paired,
	    int extraband_single, int extraband_end, int extraband_paired) {
  T new = (T) MALLOC(sizeof(*new));

  new->genome = genome;
  new->altstrain_iit = altstrain_iit;
#ifdef PMAP
  new->indexdb_fwd = indexdb_fwd;
  new->indexdb_rev = indexdb_rev;
#else
  new->indexdb = indexdb;
#endif
  new->chromosome_iit = chromosome_iit;
  new->contig_iit = contig_iit;
  new->chrsubset = chrsubset;
  new->map_iit = map_iit;
  new->maxextension = maxextension;
  new->stuttercycles = stuttercycles;
  new->stutterhits = stutterhits;
  new->indexsize = indexsize;
  new->maxpeelback = maxpeelback;
  new->sufflookback = sufflookback;
  new->nsufflookback = nsufflookback;
  new->nullgap = nullgap;
  new->extramaterial_end = extramaterial_end;
  new->extramaterial_paired = extramaterial_paired;
  new->extraband_single = extraband_single;
  new->extraband_end = extraband_end;
  new->extraband_paired = extraband_paired;
  return new;
}

void
Params_free (T *old) {
  if (*old) {
    FREE(*old);
  }
  return;
}

