/* $Id: chrom.h,v 1.3 2005/02/15 01:26:18 twu Exp $ */
#ifndef CHROM_INCLUDED
#define CHROM_INCLUDED
#include "bool.h"
#include "genomicpos.h"
#include "iit-read.h"

#define T Chrom_T
typedef struct T *T;

extern void
Chrom_free (T *old);
extern char *
Chrom_to_string (T this);
extern char *
Chrom_to_string_signed (T this, bool watsonp);
extern T
Chrom_from_string (char *string);
extern char *
Chrom_string_from_position (Genomicpos_T *chrpos, Genomicpos_T position, 
			    IIT_T chromosome_iit);
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
