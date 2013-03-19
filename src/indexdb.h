/* $Id: indexdb.h 78839 2012-11-09 21:34:05Z twu $ */
#ifndef INDEXDB_INCLUDED
#define INDEXDB_INCLUDED
#include <stdio.h>
#include "access.h"
#include "types.h"
#include "mode.h"
#include "genomicpos.h"
#include "bool.h"
#include "iit-read.h"
#include "indexdbdef.h"

#ifdef PMAP
#include "alphabet.h"
#endif



#ifdef PMAP
#define SUFFICIENT_SUPPORT 9
#else
#define SUFFICIENT_SUPPORT 18
#endif

#ifdef PMAP

#define FWD_FILESUFFIX "pf"
#define REV_FILESUFFIX "pr"

#else
#define IDX_FILESUFFIX "ref"
#endif

#define OFFSETS_FILESUFFIX "offsets"
#define POSITIONS_FILESUFFIX "positions"


#define T Indexdb_T
typedef struct T *T;

#ifdef PMAP
extern void
Indexdb_setup (int index1part_aa_in);
#else
extern void
Indexdb_setup (int index1part_in);
#endif

extern void
Indexdb_free (T *old);
#ifndef PMAP
extern int
Indexdb_interval (T this);
#endif
extern bool
Indexdb_positions_fileio_p (T this);
extern double
Indexdb_mean_size (T this, Mode_T mode, int index1part);

extern bool
Indexdb_get_filenames (char **gammaptrs_filename, char **offsetscomp_filename, char **positions_filename,
		       char **gammaptrs_basename_ptr, char **offsetscomp_basename_ptr, char **positions_basename_ptr,
		       char **gammaptrs_index1info_ptr, char **offsetscomp_index1info_ptr, char **positions_index1info_ptr,
#ifdef PMAP
		       Alphabet_T *alphabet, Alphabet_T required_alphabet,
#endif
		       int *basesize, int *index1part, int *index1interval, char *genomesubdir,
		       char *fileroot, char *idx_filesuffix, char *snps_root,
		       int required_basesize, int required_index1part, int required_interval);

extern Genomicpos_T *
Indexdb_point_one_shift (int *nentries, T this, Storedoligomer_T subst);
extern int
Indexdb_count_one_shift (T this, Storedoligomer_T subst, int nadjacent);


extern Positionsptr_T *
Indexdb_offsets_from_gammas (char *gammaptrsfile, char *offsetscompfile, int offsetscomp_basesize
#ifdef PMAP
			     , int alphabet_size, int index1part_aa
#else
			     , int index1part
#endif
			     );

extern T
Indexdb_new_genome (int *basesize, int *index1part, int *index1interval,
		    char *genomesubdir, char *fileroot, char *idx_filesuffix, char *snps_root,
#ifdef PMAP
		    Alphabet_T *alphabet, int *alphabet_size, Alphabet_T required_alphabet,
#endif
		    int required_basesize, int required_index1part, int required_interval, bool expand_offsets_p,
		    Access_mode_T offsetscomp_access, Access_mode_T positions_access);
extern T
Indexdb_new_segment (char *genomicseg,
#ifdef PMAP
		     int alphabet_size, int index1part_aa, bool watsonp,
#else
		     int index1part,
#endif
		     int index1interval);

#ifdef PMAP
extern Genomicpos_T *
Indexdb_read (int *nentries, T this, unsigned int aaindex);
#else
extern Genomicpos_T *
Indexdb_read (int *nentries, T this, Storedoligomer_T oligo);
extern Genomicpos_T *
Indexdb_read_inplace (int *nentries, T this, Storedoligomer_T oligo);
#endif

extern Genomicpos_T *
Indexdb_read_with_diagterm (int *nentries, T this, Storedoligomer_T oligo, int diagterm);
extern Genomicpos_T *
Indexdb_read_with_diagterm_sizelimit (int *nentries, T this, Storedoligomer_T oligo, int diagterm,
				      int size_threshold);

extern void
Indexdb_write_gammaptrs (char *gammaptrsfile, char *offsetscompfile, Positionsptr_T *offsets,
			 Oligospace_T oligospace, int blocksize);

extern void
Indexdb_write_offsets (char *gammaptrsfile, char *offsetscompfile, FILE *sequence_fp, IIT_T chromosome_iit,
		       int offsetscomp_basesize,
#ifdef PMAP
		       int alphabet_size, int index1part_aa, bool watsonp,
#else
		       int index1part,
#endif
		       int index1interval, bool genome_lc_p, char *fileroot, bool mask_lowercase_p);

extern void
Indexdb_write_positions (char *positionsfile, char *gammaptrsfile, char *offsetscompfile,
			 FILE *sequence_fp, IIT_T chromosome_iit, int offsetscomp_basesize,
#ifdef PMAP
			 int alphabet_size, int index1part_aa, bool watsonp,
#else
			 int index1part,
#endif
			 int index1interval, bool genome_lc_p, bool writefilep,
			 char *fileroot, bool mask_lowercase_p);

extern int
Storedoligomer_compare (const void *a, const void *b);

#undef T
#endif

