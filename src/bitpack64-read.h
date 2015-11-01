#ifndef BITPACK64_READ_INCLUDED
#define BITPACK64_READ_INCLUDED
#include "types.h"

/* For reading differential-coded bitstreams */

extern void
Bitpack64_read_setup ();

extern UINT4
Bitpack64_offsetptr (UINT4 *end0, Storedoligomer_T oligo,
		     UINT4 *bitpackptrs, UINT4 *bitpackcomp);

#ifdef LARGE_GENOMES
UINT8
Bitpack64_offsetptr_huge (UINT8 *end0, Storedoligomer_T oligo,
			  UINT4 *bitpackpages, UINT4 *bitpackptrs, UINT4 *bitpackcomp);
#endif


extern UINT4
Bitpack64_offsetptr_only (Storedoligomer_T oligo,
			  UINT4 *bitpackptrs, UINT4 *bitpackcomp);

#ifdef LARGE_GENOMES
extern UINT8
Bitpack64_offsetptr_only_huge (Storedoligomer_T oligo, UINT4 *bitpackpages,
			       UINT4 *bitpackptrs, UINT4 *bitpackcomp);
#endif

#ifndef PMAP
extern void
Bitpack64_block_offsets (UINT4 *offsets, Storedoligomer_T oligo,
			 UINT4 *bitpackptrs, Offsetscomp_T *bitpackcomp);
#endif


#ifndef PMAP
#if defined(HAVE_64_BIT) && (defined(UTILITYP) || defined(LARGE_GENOMES))
void
Bitpack64_block_offsets_huge (UINT8 *offsets, Storedoligomer_T oligo,
			      UINT4 *bitpackpages, UINT4 *bitpackptrs, Offsetscomp_T *bitpackcomp);
#endif
#endif

#endif
