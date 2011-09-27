static char rcsid[] = "$Id: params.c,v 1.67 2006/11/28 00:50:11 twu Exp $";
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
  bool map_iit_universal_p;
  int map_iit_forward_type;
  int map_iit_reverse_type;

  int stuttercycles;
  int stutterhits;
  int maxoligohits;
  int minindexsize;
  int maxindexsize;
  int maxpeelback;
  int sufflookback;
  int nsufflookback;
  int nullgap;
  int extramaterial_end;
  int extramaterial_paired;
  int extraband_single;
  int extraband_end;
  int extraband_paired;
  int maxmutations;
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

bool
Params_map_iit_universal_p (T this) {
  return this->map_iit_universal_p;
}

int
Params_map_iit_forward_type (T this) {
  return this->map_iit_forward_type;
}

int
Params_map_iit_reverse_type (T this) {
  return this->map_iit_reverse_type;
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
Params_maxoligohits (T this) {
  return this->maxoligohits;
}

int
Params_minindexsize (T this) {
  return this->minindexsize;
}

int
Params_maxindexsize (T this) {
  return this->maxindexsize;
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

int
Params_maxmutations (T this) {
  return this->maxmutations;
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
	    bool map_iit_universal_p, int map_iit_forward_type, int map_iit_reverse_type,
	    int stuttercycles, int stutterhits, int maxoligohits,
	    int minindexsize, int maxindexsize, int maxpeelback,
	    int sufflookback, int nsufflookback, int nullgap, 
	    int extramaterial_end, int extramaterial_paired,
	    int extraband_single, int extraband_end, int extraband_paired,
	    int maxmutations) {
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
  new->map_iit_universal_p = map_iit_universal_p;
  new->map_iit_forward_type = map_iit_forward_type;
  new->map_iit_reverse_type = map_iit_reverse_type;
  new->stuttercycles = stuttercycles;
  new->stutterhits = stutterhits;
  new->maxoligohits = maxoligohits;
  new->minindexsize = minindexsize;
  new->maxindexsize = maxindexsize;
  new->maxpeelback = maxpeelback;
  new->sufflookback = sufflookback;
  new->nsufflookback = nsufflookback;
  new->nullgap = nullgap;
  new->extramaterial_end = extramaterial_end;
  new->extramaterial_paired = extramaterial_paired;
  new->extraband_single = extraband_single;
  new->extraband_end = extraband_end;
  new->extraband_paired = extraband_paired;
  new->maxmutations = maxmutations;
  return new;
}

void
Params_free (T *old) {
  if (*old) {
    FREE(*old);
  }
  return;
}

