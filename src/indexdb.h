/* $Id: indexdb.h,v 1.26 2005/10/12 18:00:39 twu Exp $ */
#ifndef INDEXDB_INCLUDED
#define INDEXDB_INCLUDED
#include <stdio.h>
#include "types.h"
#include "genomicpos.h"
#include "bool.h"
#include "iit-read.h"


#ifdef PMAP

#define INDEX1PART_AA 7		/* In amino acids */

#if INDEX1PART_AA == 6
#define INDEX1PART_NT 18	/* In nucleotides */
#define INDEX1INTERVAL 6	/* In amino acids */
#elif INDEX1PART_AA == 7
#define INDEX1PART_NT 21	/* In nucleotides */
#define INDEX1INTERVAL 7	/* In amino acids */
#else
#error The given value of INDEX1PART_AA is not supported by indexdb.h
#endif

#ifdef NOCOLLAPSE
#define NAMINOACIDS 20
#if INDEX1PART_AA == 6
#define MSB 3200000		/* 20^(6-1) */
#elif INDEX1PART_AA == 7
#define MSB 64000000		/* 20^(7-1) */
#else
#error The given value of INDEX1PART_AA is not supported by indexdb.h
#endif
#else
#define NAMINOACIDS 15
#if INDEX1PART_AA == 6
#define MSB 759375		/* 15^(6-1) */
#elif INDEX1PART_AA == 7
#define MSB 11390625		/* 16^(7-1) */
#else
#error The given value of INDEX1PART_AA is not supported by indexdb.h
#endif
#endif

#else
#define INDEX1PART 12		/* In nucleotides */
#define INDEX1INTERVAL 6
#endif



/* Typically 12 nt or 24 bits, requiring 3 bytes */
typedef UINT4 Storedoligomer_T;

#define T Indexdb_T
typedef struct T *T;

extern void
Indexdb_free (T *old);
extern T
Indexdb_new_genome (char *genomesubdir, char *fileroot,
#ifdef PMAP
		    bool watsonp,
#endif
		    bool batch_offsets_p, bool batch_positions_p);
extern T
Indexdb_new_segment (char *genomicseg, int index1interval
#ifdef PMAP
		     , bool watsonp
#endif
		     );

#ifdef PMAP
extern Genomicpos_T *
Indexdb_read (int *nentries, T this, Storedoligomer_T oligo);
#else
extern Genomicpos_T *
Indexdb_read (int *nentries, T this, unsigned int aaindex);
#endif

extern void
Indexdb_write_offsets (FILE *offsets_fp, FILE *sequence_fp, IIT_T altstrain_iit,
		       int index1interval,
#ifdef PMAP
		       bool watsonp,
#endif
		       bool uncompressedp);

extern void
Indexdb_write_positions (char *positionsfile, FILE *offsets_fp, FILE *sequence_fp, 
			 IIT_T altstrain_iit, int index1interval, bool uncompressedp, 
#ifdef PMAP
			 bool watsonp,
#endif
			 bool writefilep);

#undef T
#endif

