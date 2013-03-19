/* $Id: diagdef.h 79302 2012-11-15 23:55:58Z twu $ */
#ifndef DIAGDEF_INCLUDED
#define DIAGDEF_INCLUDED
#include "bool.h"

#define T Diag_T
struct T {
  Genomicpos_T diagonal;
  int querystart;
  int queryend;
  int nconsecutive;
  bool dominatedp;
  double score;
};

#undef T
#endif

