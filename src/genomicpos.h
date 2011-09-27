/* $Id: genomicpos.h 27450 2010-08-05 19:02:48Z twu $ */
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
