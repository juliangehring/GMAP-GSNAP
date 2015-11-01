#ifndef BITPACK64_WRITE_INCLUDED
#define BITPACK64_WRITE_INCLUDED
#include <stdio.h>
#include "types.h"

extern void
Bitpack64_write_setup ();

extern int
Bitpack64_write_vert (FILE *offsetscomp_fp, Positionsptr_T *offsets_buffer, int offsets_buffer_size, int offsets_buffer_i,
		      const UINT4 *horizontal, int packsize);

extern int
Bitpack64_write_horiz (FILE *offsetscomp_fp, Positionsptr_T *offsets_buffer, int offsets_buffer_size, int offsets_buffer_i,
		       const UINT4 *horizontal, int packsize);

#endif
