/* $Id: diagdef.h 99737 2013-06-27 19:33:03Z twu $ */
#ifndef DIAGDEF_INCLUDED
#define DIAGDEF_INCLUDED
#include "bool.h"

#define T Diag_T
struct T {
  Chrpos_T diagonal;
  int querystart;
  int queryend;
  int nconsecutive;
  bool dominatedp;
  double score;
};

#undef T
#endif

