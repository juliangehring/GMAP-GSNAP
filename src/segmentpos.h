/* $Id: segmentpos.h,v 1.27 2005/07/08 07:58:34 twu Exp $ */
#ifndef SEGMENTPOS_INCLUDED
#define SEGMENTPOS_INCLUDED
#include <stdio.h>
#include "bool.h"
#include "genomicpos.h"
#include "chrom.h"
#include "iit-read.h"

#define T Segmentpos_T
typedef struct T *T;

extern Chrom_T
Segmentpos_chrom (T this);
extern Genomicpos_T
Segmentpos_chrpos1 (T this);
extern Genomicpos_T
Segmentpos_chrpos2 (T this);
extern Genomicpos_T
Segmentpos_length (T this);
extern int
Segmentpos_type (T this);
extern bool
Segmentpos_revcompp (T this);
extern T
Segmentpos_new (Chrom_T chrom, Genomicpos_T chrpos1, Genomicpos_T chrpos2, 
		bool revcompp, Genomicpos_T length, int type);
extern void
Segmentpos_free (T *old);
extern void
Segmentpos_print (FILE *fp, T this, char *acc, Genomicpos_T offset);
extern int
Segmentpos_compare (const void *x, const void *y);

extern void
Segmentpos_print_accessions (IIT_T contig_iit, Genomicpos_T position1,
			     Genomicpos_T position2, bool referencealignp, 
			     char *align_strain, bool zerobasedp);

#undef T
#endif


