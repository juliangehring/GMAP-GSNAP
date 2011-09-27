/* $Id: chrom.h 44063 2011-08-01 18:04:15Z twu $ */
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
Chrom_from_string (char *string, char *mitochondrial_string, unsigned int order);

extern int
Chrom_cmp_alpha (T a, T b);
extern int
Chrom_cmp_numeric_alpha (T a, T b);
extern int
Chrom_cmp_chrom (T a, T b);

extern int
Chrom_compare_order (const void *x, const void *y);
extern int
Chrom_compare_alpha (const void *x, const void *y);
extern int
Chrom_compare_numeric_alpha (const void *x, const void *y);
extern int
Chrom_compare_chrom (const void *x, const void *y);

extern int
Chrom_compare_table (const void *x, const void *y);
extern unsigned int
Chrom_hash_table (const void *key);

#undef T
#endif
