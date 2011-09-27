/* $Id: indexdb.h,v 1.52 2009/08/14 15:05:19 twu Exp $ */
#ifndef INDEXDB_INCLUDED
#define INDEXDB_INCLUDED
#include <stdio.h>
#include "types.h"
#include "genomicpos.h"
#include "bool.h"
#include "iit-read.h"


#ifdef PMAP
#define SUFFICIENT_SUPPORT 9
#else
#define SUFFICIENT_SUPPORT 18
#endif

#ifdef PMAP

#define FWD_FILESUFFIX "pf"
#define REV_FILESUFFIX "pr"

#define INDEX1PART_AA 6		/* In amino acids */
#define INDEX1PART_NT 18	/* In nucleotides */

#if INDEX1PART_AA == 5
#define NAMINOACIDS 20
#define MSB 160000		/* 20^(5-1) */

#elif INDEX1PART_AA == 6
#define NAMINOACIDS 20
#define MSB 3200000		/* 20^(6-1) */

#elif INDEX1PART_AA == 7
#ifdef NOCOLLAPSE
#define NAMINOACIDS 20
#define MSB 64000000		/* 20^(7-1) */
#else
#define NAMINOACIDS 15
#define MSB 11390625		/* 15^(7-1) */
#endif

#elif INDEX1PART_AA == 8
#ifdef NOCOLLAPSE
#define NAMINOACIDS 20
#define MSB 1280000000		/* 20^(8-1) */
#else
#define NAMINOACIDS 12
#define MSB 35831808		/* 12^(8-1) */
#endif

#else
#error The given value of INDEX1PART_AA is not supported by indexdb.h
#endif	/* INDEX1PART */

#else  /* PMAP */
#define IDX_FILESUFFIX "ref"
#define INDEX1PART 12		/* In nucleotides */
#endif	/* PMAP */

#define OFFSETS_FILESUFFIX "offsets"
#define POSITIONS_FILESUFFIX "positions"


/* Typically 12 nt or 24 bits, requiring 3 bytes */
typedef UINT4 Storedoligomer_T;

#define T Indexdb_T
typedef struct T *T;

extern void
Indexdb_free (T *old);
#ifndef PMAP
extern int
Indexdb_interval (T this);
#endif
extern double
Indexdb_mean_size (T this, bool cmetp);
extern T
Indexdb_new_genome (char *genomesubdir, char *fileroot, char *idx_filesuffix, char *snps_root,
		    int required_interval, bool batch_offsets_p, bool batch_positions_p);
extern T
Indexdb_new_segment (char *genomicseg, int index1interval
#ifdef PMAP
		     , bool watsonp
#endif
		     );

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
Indexdb_write_offsets (FILE *offsets_fp, FILE *sequence_fp, IIT_T chromosome_iit,
		       int index1interval,
#ifdef PMAP
		       bool watsonp,
#endif
		       bool genome_lc_p, char *fileroot, bool mask_lowercase_p);

extern void
Indexdb_write_positions (char *positionsfile, FILE *offsets_fp, FILE *sequence_fp, 
			 IIT_T chromosome_iit, int index1interval,
#ifdef PMAP
			 bool watsonp,
#endif	
			 bool genome_lc_p, bool writefilep, char *fileroot, bool mask_lowercase_p);

extern int
Storedoligomer_compare (const void *a, const void *b);

#undef T
#endif

