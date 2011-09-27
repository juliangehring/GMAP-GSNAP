/* $Id: diagdef.h,v 1.2 2007/02/06 14:30:44 twu Exp $ */
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
};

#undef T
#endif

