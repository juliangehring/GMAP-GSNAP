/* $Id: genomicpos.h 101488 2013-07-15 16:52:36Z twu $ */
#ifndef GENOMICPOS_INCLUDED
#define GENOMICPOS_INCLUDED
#include <stdlib.h>
#include "types.h"

/* A genomic position */
#ifdef LARGE_GENOMES
#include "uint8list.h"
typedef Uint8list_T Genomicposlist_T;
#define Genomicposlist_length(x) Uint8list_length(x)
#define Genomicposlist_to_array(x,y) Uint8list_to_array(x,y)
#define Genomicposlist_free(x) Uint8list_free(x)
#else
#include "uintlist.h"
typedef Uintlist_T Genomicposlist_T;
#define Genomicposlist_length(x) Uintlist_length(x)
#define Genomicposlist_to_array(x,y) Uintlist_to_array(x,y)
#define Genomicposlist_free(x) Uintlist_free(x)
#endif

/* A chromosomal position */
typedef UINT4 Chrpos_T;

extern char *
Genomicpos_commafmt (
#ifdef HAVE_64_BIT
		     UINT8 N
#else
		     UINT4 N
#endif
		     );
extern int
UINT8_compare (const void *a, const void *b);
extern int
UINT4_compare (const void *a, const void *b);
extern int
Univcoord_compare (const void *a, const void *b);
extern int
Chrpos_compare (const void *a, const void *b);

#endif
