/* $Id: matchdef.h,v 1.2 2005/01/25 06:21:59 twu Exp $ */
#ifndef MATCHDEF_INCLUDED
#define MATCHDEF_INCLUDED
#include "bool.h"
#include "chrnum.h"
#include "genomicpos.h"

#define T Match_T
struct T {
  int querypos;
  Genomicpos_T position;
  Chrnum_T chrnum;
  Genomicpos_T chrpos;
  bool forwardp;
  bool fivep;
  bool pairedp;
};

#undef T
#endif

