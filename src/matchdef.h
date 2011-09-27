/* $Id: matchdef.h,v 1.6 2010-07-10 01:35:55 twu Exp $ */
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

