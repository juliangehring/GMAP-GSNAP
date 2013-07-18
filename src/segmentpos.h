/* $Id: segmentpos.h 99737 2013-06-27 19:33:03Z twu $ */
#ifndef SEGMENTPOS_INCLUDED
#define SEGMENTPOS_INCLUDED
#include <stdio.h>
#include "bool.h"
#include "genomicpos.h"
#include "types.h"
#include "chrom.h"
#include "iit-read-univ.h"

#define T Segmentpos_T
typedef struct T *T;

extern Chrom_T
Segmentpos_chrom (T this);
extern Chrpos_T
Segmentpos_chrpos1 (T this);
extern Chrpos_T
Segmentpos_chrpos2 (T this);
extern Chrpos_T
Segmentpos_length (T this);
extern int
Segmentpos_type (T this);
extern bool
Segmentpos_revcompp (T this);
extern T
Segmentpos_new (Chrom_T chrom, Chrpos_T chrpos1, Chrpos_T chrpos2, 
		bool revcompp, Chrpos_T length, int type);
extern void
Segmentpos_free (T *old);
extern void
Segmentpos_print (FILE *fp, T this, char *acc, Univcoord_T chroffset);
extern int
Segmentpos_compare_alpha (const void *x, const void *y);
extern int
Segmentpos_compare_numeric_alpha (const void *x, const void *y);
extern int
Segmentpos_compare_chrom (const void *x, const void *y);

extern void
Segmentpos_print_accessions (FILE *fp, Univ_IIT_T contig_iit, Univcoord_T position1,
			     Univcoord_T position2, bool referencealignp, 
			     char *align_strain);

#undef T
#endif


