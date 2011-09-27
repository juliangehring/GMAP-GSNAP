/* $Id: genomicpos.h 40271 2011-05-28 02:29:18Z twu $ */
#ifndef GENOMICPOS_INCLUDED
#define GENOMICPOS_INCLUDED
#include "types.h"

/* A genomic position, typically 3 billion or less, requiring 32 bits
   or 4 bytes */
typedef UINT4 Genomicpos_T;

#define T Genomicpos_T

extern char *
Genomicpos_commafmt (T N);
extern int
Genomicpos_compare (const void *a, const void *b);

#undef T
#endif
