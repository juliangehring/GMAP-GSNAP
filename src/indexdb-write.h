/* $Id: indexdb-write.h 100260 2013-07-02 23:48:24Z twu $ */
#ifndef INDEXDB_WRITE_INCLUDED
#define INDEXDB_WRITE_INCLUDED

#include "bool.h"
#include "types.h"		/* For Oligospace_T */
#include "iit-read-univ.h"

#ifdef PMAP

#define FWD_FILESUFFIX "pf"
#define REV_FILESUFFIX "pr"

#else
#define IDX_FILESUFFIX "ref"
#endif

#define OFFSETS_FILESUFFIX "offsets"
#define POSITIONS_FILESUFFIX "positions"

extern void
Indexdb_write_bitpackptrs (char *bitpackptrsfile, char *offsetsfile, Positionsptr_T *offsets,
			   Oligospace_T oligospace, Blocksize_T blocksize);

extern void
Indexdb_write_gammaptrs (char *gammaptrsfile, char *offsetscompfile, Positionsptr_T *offsets,
			 Oligospace_T oligospace, Blocksize_T blocksize);

extern void
Indexdb_write_offsets (char *destdir, char interval_char, FILE *sequence_fp, Univ_IIT_T chromosome_iit,
		       Width_T offsetscomp_basesize,
#ifdef PMAP
		       int alphabet_size, Width_T index1part_aa, bool watsonp,
#else
		       Width_T index1part,
#endif
		       Width_T index1interval, bool genome_lc_p, char *fileroot, bool mask_lowercase_p,
		       int compression_types);

extern void
Indexdb_write_positions (char *positionsfile, char *pointersfile, char *offsetsfile,
			 FILE *sequence_fp, Univ_IIT_T chromosome_iit, Width_T offsetscomp_basesize,
#ifdef PMAP
			 int alphabet_size, Width_T index1part_aa, bool watsonp,
#else
			 Width_T index1part,
#endif
			 Width_T index1interval, bool genome_lc_p, bool writefilep,
			 char *fileroot, bool mask_lowercase_p, int compression_type,
			 bool coord_values_8p);

#endif

