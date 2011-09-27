/* $Id: diagdef.h 40271 2011-05-28 02:29:18Z twu $ */
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

