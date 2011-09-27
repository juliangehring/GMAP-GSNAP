/* $Id: oligo.h,v 1.34 2009/05/16 14:06:40 twu Exp $ */
#ifndef OLIGO_INCLUDED
#define OLIGO_INCLUDED
#include "bool.h"
#include "genomicpos.h"
#include "indexdb.h"
#include "reader.h"

typedef enum {INIT, DONE, INVALID, VALID} Oligostate_T;

extern char *
Oligo_one_nt (Storedoligomer_T oligo, int oligosize);

extern int
Oligo_lookup (Genomicpos_T **positions, Indexdb_T indexdb, Storedoligomer_T storedoligo);

extern Oligostate_T
Oligo_next (Oligostate_T last_state, int *querypos, Storedoligomer_T *forward, 
	    Storedoligomer_T *revcomp, int oligosize, Reader_T reader, cDNAEnd_T cdnaend);
extern Oligostate_T
Oligo_skip (Oligostate_T last_state, int *querypos, Storedoligomer_T *forward,
	    Storedoligomer_T *revcomp, int oligosize, Reader_T reader, cDNAEnd_T cdnaend, int nskip);

extern char *
Oligo_nt (Storedoligomer_T oligo1, Storedoligomer_T oligo2, int oligosize);

extern bool
Oligo_repetitive_p (Storedoligomer_T oligo);

extern bool
Oligo_mark_repetitive (bool **repetitivep, Storedoligomer_T *oligos,
		       int first_querypos, int last_querypos);

#endif
