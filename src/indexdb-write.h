/* $Id: indexdb-write.h 99737 2013-06-27 19:33:03Z twu $ */
#ifndef INDEXDB_WRITE_INCLUDED
#define INDEXDB_WRITE_INCLUDED
#include <stdio.h>
#include "access.h"
#include "types.h"
#include "mode.h"
#include "genomicpos.h"
#include "bool.h"
#include "iit-read-univ.h"

#ifdef PMAP
#include "alphabet.h"
#endif


#ifdef PMAP

#define FWD_FILESUFFIX "pf"
#define REV_FILESUFFIX "pr"

#else
#define IDX_FILESUFFIX "ref"
#endif

#define OFFSETS_FILESUFFIX "offsets"
#define POSITIONS_FILESUFFIX "positions"



extern void
Indexdb_write_gammaptrs (char *gammaptrsfile, char *offsetscompfile, Positionsptr_T *offsets,
			 Oligospace_T oligospace, int blocksize);

extern void
Indexdb_write_offsets (char *gammaptrsfile, char *offsetscompfile, FILE *sequence_fp, Univ_IIT_T chromosome_iit,
		       int offsetscomp_basesize,
#ifdef PMAP
		       int alphabet_size, int index1part_aa, bool watsonp,
#else
		       int index1part,
#endif
		       int index1interval, bool genome_lc_p, char *fileroot, bool mask_lowercase_p);

extern void
Indexdb_write_positions (char *positionsfile, char *gammaptrsfile, char *offsetscompfile,
			 FILE *sequence_fp, Univ_IIT_T chromosome_iit, int offsetscomp_basesize,
#ifdef PMAP
			 int alphabet_size, int index1part_aa, bool watsonp,
#else
			 int index1part,
#endif
			 int index1interval, bool genome_lc_p, bool writefilep,
			 char *fileroot, bool mask_lowercase_p, bool coord_values_8p);

#endif

