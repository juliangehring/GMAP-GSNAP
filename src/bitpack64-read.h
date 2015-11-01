#ifndef BITPACK64_READ_INCLUDED
#define BITPACK64_READ_INCLUDED
#include "types.h"

extern void
Bitpack64_read_setup (UINT4 *bitpackptrs_in, UINT4 *offsetscomp_in,
		      Blocksize_T offsetscomp_blocksize_in);

extern Positionsptr_T
Bitpack64_offsetptr (Positionsptr_T *end0, Storedoligomer_T oligo);

extern Positionsptr_T
Bitpack64_offsetptr_only (Storedoligomer_T oligo);

extern void
Bitpack64_block_offsets (Positionsptr_T *offsets, UINT4 *bitpackptrs, Offsetscomp_T *offsetscomp,
			 Blocksize_T offsetscomp_blocksize, Storedoligomer_T oligo);
#endif
