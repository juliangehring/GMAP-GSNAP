static char rcsid[] = "$Id: match.c,v 1.66 2005/07/08 14:39:55 twu Exp $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#ifndef HAVE_MEMCPY
# define memcpy(d,s,n) bcopy((s),(d),(n))
#endif
#ifndef HAVE_MEMMOVE
# define memmove(d,s,n) bcopy((s),(d),(n))
#endif

#include "match.h"
#include "matchdef.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>		/* For memcpy */
#include "mem.h"
#include "segmentpos.h"
#include "separator.h"

#define T Match_T

#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

int
Match_querypos (T this) {
  return this->querypos;
}

bool
Match_forwardp (T this) {
  return this->forwardp;
}

bool
Match_fivep (T this) {
  return this->fivep;
}

Genomicpos_T
Match_position (T this) {
  return this->position;
}

Chrnum_T
Match_chrnum (T this) {
  return this->chrnum;
}

Genomicpos_T
Match_chrpos (T this) {
  return this->chrpos;
}

void
Match_set_pairedp (T this) {
  this->pairedp = true;
}

bool
Match_pairedp (T this) {
  return this->pairedp;
}

int
Match_cmp (const void *a, const void *b) {
  T x = * (T *) a;
  T y = * (T *) b;

  if (x->position < y->position) {
    return -1;
  } else if (x->position > y->position) {
    return 1;
  } else {
    return 0;
  }
}

T
Match_new (int querypos, bool forwardp, bool fivep,
	   Genomicpos_T position, IIT_T chromosome_iit) {
  T new = (T) MALLOC(sizeof(*new));
  int index;

  new->querypos = querypos;
  new->position = position;
  new->forwardp = forwardp;
  new->fivep = fivep;
  new->pairedp = false;

  if (chromosome_iit == NULL) {
    new->chrnum = 0;
    new->chrpos = position;
  } else {
    if (forwardp == true) {
      index = IIT_get_one(chromosome_iit,position,position);
    } else {
      index = IIT_get_one(chromosome_iit,position-1,position-1);
    }
    new->chrpos = position - Interval_low(IIT_interval(chromosome_iit,index));
    new->chrnum = index;
  }

  return new;
}

T
Match_copy (T this) {
  T new = (T) MALLOC(sizeof(*new));

  memcpy((void *) new,this,sizeof(struct T));
  return new;
}


void
Match_free (T *old) {
  if (*old) {
    FREE(*old);
  }
  return;
}

/* Matches format of Pair_print_two */
void
Match_print_two (int pathnum, T start, T end, IIT_T chromosome_iit, IIT_T contig_iit, 
		 char *dbversion, bool zerobasedp) {
  unsigned int querypos1, querypos2;
  Genomicpos_T chrpos1, chrpos2, position1, position2;
  char *chrstring, *comma1, *comma2;

  querypos1 = start->querypos;
  querypos2 = end->querypos;

  printf("  Path %d: ",pathnum);

  printf("query %u%s%u (%u bp) => ",
	 querypos1 + !zerobasedp,SEPARATOR,querypos2 + !zerobasedp,querypos2-querypos1+1);
  
  chrstring = Chrnum_to_string(start->chrnum,chromosome_iit);
  if (start->chrpos <= end->chrpos) {
    chrpos1 = start->chrpos;
    chrpos2 = end->chrpos;
  } else {
    chrpos1 = end->chrpos;
    chrpos2 = start->chrpos;
  }    
    
  comma1 = Genomicpos_commafmt(chrpos1 + !zerobasedp);
  comma2 = Genomicpos_commafmt(chrpos2 + !zerobasedp);
  printf("chr %s:%s%s%s (%u bp)\n",
	 chrstring,comma1,SEPARATOR,comma2,chrpos2-chrpos1+1);
  FREE(comma2);
  FREE(comma1);

  position1 = start->position;
  position2 = end->position;
  comma1 = Genomicpos_commafmt(position1 + !zerobasedp);
  comma2 = Genomicpos_commafmt(position2 + !zerobasedp);
  printf("    Genomic pos: %s:%s%s%s",dbversion,comma1,SEPARATOR,comma2);
  FREE(comma2);
  FREE(comma1);

  if (start->chrpos <= end->chrpos) {
    printf(" (+ strand)\n");
  } else {
    printf(" (- strand)\n");
  }
  
  if (contig_iit != NULL) {
    if (position1 <= position2) {
      Segmentpos_print_accessions(contig_iit,position1,position2,/*referencealignp*/true,/*align_strain*/NULL,zerobasedp);
    } else {
      Segmentpos_print_accessions(contig_iit,position2,position1,/*referencealignp*/true,/*align_strain*/NULL,zerobasedp);
    }
  }
    
  printf("\n");
  FREE(chrstring);

  return;
}

