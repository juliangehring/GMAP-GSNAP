/* $Id: block.h,v 1.35 2005/05/04 18:03:48 twu Exp $ */
#ifndef BLOCK_INCLUDED
#define BLOCK_INCLUDED
#include "bool.h"
#include "genomicpos.h"
#include "indexdb.h"
#include "reader.h"

#define T Block_T
typedef struct T *T;

extern int
Block_querypos (T this);
extern Storedoligomer_T
Block_forward (T this);
extern Storedoligomer_T
Block_revcomp (T this);
extern void
Block_reset_ends (T this);

extern T
Block_new (cDNAEnd_T cdnaend, Reader_T reader);
extern void
Block_free (T *old);
extern bool
Block_next (T this);
extern bool
Block_skip (T this, int nskip);

extern int
Block_process_oligo (Genomicpos_T **fwdpositions, int *nfwdhits, Genomicpos_T **revpositions, int *nrevhits,
		     T this, Indexdb_T indexdb);

#undef T
#endif
