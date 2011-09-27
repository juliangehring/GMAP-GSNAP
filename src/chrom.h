/* $Id: chrom.h,v 1.5 2008-03-31 23:49:01 twu Exp $ */
#ifndef CHROM_INCLUDED
#define CHROM_INCLUDED
#include "bool.h"
#include "genomicpos.h"

#define T Chrom_T
typedef struct T *T;

extern void
Chrom_free (T *old);
extern char *
Chrom_string (T this);
extern T
Chrom_from_string (char *string);
extern int
Chrom_cmp (T a, T b);
extern int
Chrom_compare (const void *x, const void *y);
extern int
Chrom_compare_table (const void *x, const void *y);
extern unsigned int
Chrom_hash_table (const void *key);

#undef T
#endif
