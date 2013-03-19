#ifndef OLIGOINDEX_HR_INCLUDED
#define OLIGOINDEX_HR_INCLUDED
#include "bool.h"
#include "types.h"
#include "mode.h"
#include "genomicpos.h"
#include "oligoindex.h"

#define T Oligoindex_T

extern void
Oligoindex_hr_setup (UINT4 *ref_blocks_in, Mode_T mode_in);

extern void
Oligoindex_hr_tally (T this, Genomicpos_T mappingstart, Genomicpos_T mappingend, bool plusp,
		     char *queryuc_ptr, int querylength, Genomicpos_T chrpos, int genestrand);

#undef T
#endif

