/* $Id: indexdb.h 99753 2013-06-27 21:13:11Z twu $ */
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
Indexdb_setup (Width_T index1part_aa_in);
#else
extern void
Indexdb_setup (Width_T index1part_in);
#endif

extern void
Indexdb_free (T *old);
#ifndef PMAP
extern Width_T
Indexdb_interval (T this);
#endif
extern bool
Indexdb_positions_fileio_p (T this);
extern double
Indexdb_mean_size (T this, Mode_T mode, Width_T index1part);

extern bool
Indexdb_get_filenames (char **gammaptrs_filename, char **offsetscomp_filename, char **positions_filename,
		       char **gammaptrs_basename_ptr, char **offsetscomp_basename_ptr, char **positions_basename_ptr,
		       char **gammaptrs_index1info_ptr, char **offsetscomp_index1info_ptr, char **positions_index1info_ptr,
#ifdef PMAP
		       Alphabet_T *alphabet, Alphabet_T required_alphabet,
#endif
		       Width_T *basesize, Width_T *index1part, Width_T *index1interval, char *genomesubdir,
		       char *fileroot, char *idx_filesuffix, char *snps_root,
		       Width_T required_basesize, Width_T required_index1part, Width_T required_interval);

extern Univcoord_T *
Indexdb_point_one_shift (int *nentries, T this, Storedoligomer_T subst);
extern int
Indexdb_count_one_shift (T this, Storedoligomer_T subst, int nadjacent);


extern Positionsptr_T *
Indexdb_offsets_from_gammas (char *gammaptrsfile, char *offsetscompfile, Width_T offsetscomp_basesize
#ifdef PMAP
			     , int alphabet_size, Width_T index1part_aa
#else
			     , Width_T index1part
#endif
			     );

extern T
Indexdb_new_genome (Width_T *basesize, Width_T *index1part, Width_T *index1interval,
		    char *genomesubdir, char *fileroot, char *idx_filesuffix, char *snps_root,
#ifdef PMAP
		    Alphabet_T *alphabet, int *alphabet_size, Alphabet_T required_alphabet,
#endif
		    Width_T required_basesize, Width_T required_index1part, Width_T required_interval, bool expand_offsets_p,
		    Access_mode_T offsetscomp_access, Access_mode_T positions_access);
extern T
Indexdb_new_segment (char *genomicseg,
#ifdef PMAP
		     int alphabet_size, Width_T index1part_aa, bool watsonp,
#else
		     Width_T index1part,
#endif
		     Width_T index1interval);

#ifdef PMAP
extern Univcoord_T *
Indexdb_read (int *nentries, T this, Storedoligomer_T aaindex);
#else
extern Univcoord_T *
Indexdb_read (int *nentries, T this, Storedoligomer_T oligo);
extern Univcoord_T *
Indexdb_read_inplace (int *nentries, T this, Storedoligomer_T oligo);
#endif

extern Univcoord_T *
Indexdb_read_with_diagterm (int *nentries, T this, Storedoligomer_T oligo, int diagterm);
extern Univcoord_T *
Indexdb_read_with_diagterm_sizelimit (int *nentries, T this, Storedoligomer_T oligo, int diagterm,
				      int size_threshold);

extern int
Storedoligomer_compare (const void *a, const void *b);

#undef T
#endif

