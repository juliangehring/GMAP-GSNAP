/* $Id: pairdef.h,v 1.11 2006/02/20 03:11:20 twu Exp $ */
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
  bool aamarker_g;
  bool aamarker_e;
  bool gapp;			/* True if comp is in a big gap:
                                   >])([<#= (but not '-' or '~'). */
  int queryjump;		/* Used only for gaps */
  int genomejump;		/* Used only for gaps */

  bool shortexonp;
};

#undef T
#endif
