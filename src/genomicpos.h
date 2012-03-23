/* $Id: genomicpos.h 56964 2012-02-02 17:57:52Z twu $ */
#ifndef GENOMICPOS_INCLUDED
#define GENOMICPOS_INCLUDED
#include <stdlib.h>
#include "types.h"

/* A genomic position, typically 3 billion or less, requiring 32 bits
   or 4 bytes */
typedef UINT4 Genomicpos_T;

#define T Genomicpos_T

extern char *
Genomicpos_commafmt (size_t N);
extern int
Genomicpos_compare (const void *a, const void *b);

#undef T
#endif
