/* $Id: matchdef.h 27450 2010-08-05 19:02:48Z twu $ */
#ifndef MATCHDEF_INCLUDED
#define MATCHDEF_INCLUDED
#include "bool.h"
#include "chrnum.h"
#include "genomicpos.h"

#define T Match_T
struct T {
  Genomicpos_T position;
  Chrnum_T chrnum;
  Genomicpos_T chrpos;
  double weight;		/* equal to 1/nentries */
  bool has_weight_p;
  int querypos;
  int npairings;		/* number of matchpairs made with
				   other matches */
  bool forwardp;
  bool fivep;
};

#undef T
#endif

