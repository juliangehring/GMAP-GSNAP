/* $Id: pairdef.h,v 1.9 2005/01/22 14:51:06 twu Exp $ */
#ifndef PAIRDEF_INCLUDED
#define PAIRDEF_INCLUDED
#include "bool.h"

#define T Pair_T
struct T {
  int querypos;
  int genomepos;
  int refquerypos;
  int aapos;
  char cdna;
  char comp;
  char genome;
  char aa_g;			/* Genomic aa */
  char aa_e;			/* EST aa */
  bool aamarker;
  bool gapp;			/* True if comp is in a big gap:
                                   >])([<#= (but not '-' or '~'). */
  bool shortexonp;
};

#undef T
#endif
