/* $Id: block.h 27450 2010-08-05 19:02:48Z twu $ */
#ifndef BLOCK_INCLUDED
#define BLOCK_INCLUDED
#include "bool.h"
#include "genomicpos.h"
#include "indexdb.h"
#include "reader.h"

#define T Block_T
typedef struct T *T;

extern Reader_T
Block_reader (T this);

extern int
Block_querypos (T this);

#ifdef PMAP
unsigned int
Block_aaindex (T this);
#else
extern Storedoligomer_T
Block_forward (T this);
extern Storedoligomer_T
Block_revcomp (T this);
#endif

extern bool
Block_donep (T this);

extern void
Block_save (T this);
extern void
Block_restore (T this);
extern void
Block_reset_ends (T this);

extern T
Block_new (cDNAEnd_T cdnaend,
#ifndef PMAP
	   int oligosize, int leftreadshift,
#endif
	   Reader_T reader, int querylength);
extern void
Block_free (T *old);
extern bool
Block_next (T this);
extern bool
Block_skip (T this, int nskip);
extern bool
Block_skipto (T this, int querypos);

#ifdef PMAP
extern int
Block_process_oligo (Genomicpos_T **fwdpositions, int *nfwdhits, 
		     Genomicpos_T **revpositions, int *nrevhits,
		     T this, Indexdb_T indexdb_fwd, Indexdb_T indexdb_rev);
#else
extern int
Block_process_oligo (Genomicpos_T **fwdpositions, int *nfwdhits,
		     Genomicpos_T **revpositions, int *nrevhits,
		     T this, Indexdb_T indexdb);
#endif

#undef T
#endif
