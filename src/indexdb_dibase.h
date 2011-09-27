/* $Id: indexdb_dibase.h 31011 2010-10-28 16:42:01Z twu $ */
#ifndef INDEXDB_DIBASE_INCLUDED
#define INDEXDB_DIBASE_INCLUDED
#include <stdio.h>
#include "types.h"
#include "genomicpos.h"
#include "bool.h"
#include "iit-read.h"


#define SUFFICIENT_SUPPORT 18

#define DIBASE_FILESUFFIX "dibase"
#define INDEX1PART 12		/* In dibases */

#define OFFSETS_FILESUFFIX "offsets"
#define POSITIONS_FILESUFFIX "positions"


/* Typically 12 nt or 24 bits, requiring 3 bytes */
typedef UINT4 Storedoligomer_T;

#define T Indexdb_T
typedef struct T *T;

extern void
Indexdb_free (T *old);
extern int
Indexdb_interval (T this);
extern double
Indexdb_mean_size (T this, bool cmetp);
extern T
Indexdb_new_genome (char *genomesubdir, char *fileroot, char *idx_filesuffix, char *snps_root,
		    int required_interval, bool batch_offsets_p, bool batch_positions_p);
extern T
Indexdb_new_segment (char *genomicseg, int index1interval);

extern Genomicpos_T *
Indexdb_read (int *nentries, T this, Storedoligomer_T oligo);
extern Genomicpos_T *
Indexdb_read_inplace (int *nentries, T this, Storedoligomer_T oligo);

extern Genomicpos_T *
Indexdb_read_with_diagterm (int *nentries, T this, Storedoligomer_T oligo, int diagterm);
extern Genomicpos_T *
Indexdb_read_with_diagterm_sizelimit (int *nentries, T this, Storedoligomer_T oligo, int diagterm,
				      int size_threshold);

extern void
Indexdb_write_offsets (FILE *offsets_fp, FILE *sequence_fp, IIT_T chromosome_iit,
		       int index1interval, bool genome_lc_p, char *fileroot, bool mask_lowercase_p);

extern void
Indexdb_write_positions (char *positionsfile, FILE *offsets_fp, FILE *sequence_fp, 
			 IIT_T chromosome_iit, int index1interval,
			 bool genome_lc_p, bool writefilep, char *fileroot, bool mask_lowercase_p);

extern int
Storedoligomer_compare (const void *a, const void *b);

#undef T
#endif

