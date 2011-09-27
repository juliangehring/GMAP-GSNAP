/* $Id: diagdef.h,v 1.6 2008/09/15 17:43:57 twu Exp $ */
#ifndef DIAGDEF_INCLUDED
#define DIAGDEF_INCLUDED
#include "bool.h"

#define T Diag_T
struct T {
  int diagonal;
  int querystart;
  int queryend;
  int nconsecutive;
  bool dominatedp;
  double score;
};

#undef T
#endif

