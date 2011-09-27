/* $Id: oligo.h,v 1.25 2005/07/08 14:37:40 twu Exp $ */
#ifndef OLIGO_INCLUDED
#define OLIGO_INCLUDED
#include "bool.h"
#include "genomicpos.h"
#include "indexdb.h"
#include "reader.h"

typedef enum {INIT, DONE, INVALID, VALID} Oligostate_T;

extern int
Oligo_lookup (Genomicpos_T **positions, Indexdb_T indexdb, Storedoligomer_T oligo);

extern Oligostate_T
Oligo_next (Oligostate_T last_state, int *querypos, Storedoligomer_T *forward, 
	    Storedoligomer_T *revcomp, Reader_T reader, cDNAEnd_T cdnaend);
extern Oligostate_T
Oligo_skip (Oligostate_T last_state, int *querypos, Storedoligomer_T *forward,
	    Storedoligomer_T *revcomp, Reader_T reader, cDNAEnd_T cdnaend, int nskip);

extern char *
Oligo_nt (Storedoligomer_T oligo1, Storedoligomer_T oligo2, int oligosize);

#endif
