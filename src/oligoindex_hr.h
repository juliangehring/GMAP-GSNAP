#ifndef OLIGOINDEX_HR_INCLUDED
#define OLIGOINDEX_HR_INCLUDED
#include "bool.h"
#include "types.h"
#include "genomicpos.h"
#include "oligoindex.h"

#define T Oligoindex_T

extern void
Oligoindex_hr_setup (UINT4 *ref_blocks_in);

extern void
Oligoindex_hr_tally (T this, Genomicpos_T genomicstart, Genomicpos_T genomicend,
		     Genomicpos_T mappingstart, Genomicpos_T mappingend, bool plusp,
		     char *queryuc_ptr, int querylength);

#undef T
#endif

