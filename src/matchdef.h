/* $Id: matchdef.h,v 1.5 2006/02/26 19:07:11 twu Exp $ */
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
  int querypos;
  int npairings;		/* number of matchpairs made with
				   other matches */
  bool forwardp;
  bool fivep;
};

#undef T
#endif

