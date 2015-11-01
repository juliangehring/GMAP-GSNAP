#ifndef BITPACK64_WRITE_INCLUDED
#define BITPACK64_WRITE_INCLUDED
#include <stdio.h>
#include "types.h"

/* Stores values 0..n */
extern void
Bitpack64_write_differential (char *ptrsfile, char *compfile, UINT4 *ascending, UINT4 n);
extern void
Bitpack64_write_differential_paired (char *ptrsfile, char *compfile, UINT4 *ascending, UINT4 n);
extern void
Bitpack64_write_fixed10 (char *ptrsfile, char *compfile, UINT4 *ascending, UINT4 n);
extern void
Bitpack64_write_differential_huge (char *pagesfile, char *ptrsfile, char *compfile,
				   UINT8 *ascending, UINT4 n);
extern void
Bitpack64_write_fixed10_huge (char *pagesfile, char *ptrsfile, char *compfile,
			      UINT8 *ascending, UINT4 n);

/* Stores values 0..(n-1) */
extern void
Bitpack64_write_direct (char *ptrsfile, char *compfile, UINT4 *direct, UINT4 n);

#endif
