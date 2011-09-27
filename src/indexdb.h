/* $Id: indexdb.h,v 1.33 2007/02/20 17:01:59 twu Exp $ */
#ifndef INDEXDB_INCLUDED
#define INDEXDB_INCLUDED
#include <stdio.h>
#include "types.h"
#include "genomicpos.h"
#include "bool.h"
#include "iit-read.h"


#ifdef PMAP

#define INDEX1PART_AA 6		/* In amino acids */

#if INDEX1PART_AA == 5
#define INDEX1PART_NT 15	/* In nucleotides */
#define INDEX1INTERVAL 5	/* In amino acids */
#elif INDEX1PART_AA == 6
#define INDEX1PART_NT 18	/* In nucleotides */
#define INDEX1INTERVAL 6	/* In amino acids */
#elif INDEX1PART_AA == 7
#define INDEX1PART_NT 21	/* In nucleotides */
#define INDEX1INTERVAL 7	/* In amino acids */
#elif INDEX1PART_AA == 8
#define INDEX1PART_NT 24	/* In nucleotides */
#define INDEX1INTERVAL 8	/* In amino acids */
#else
#error The given value of INDEX1PART_AA is not supported by indexdb.h
#endif

#if INDEX1PART_AA == 5
#define PFXOFFSETS "pf5offsets"
#define PFXPOSITIONS "pf5positions"
#define PRXOFFSETS "pr5offsets"
#define PRXPOSITIONS "pr5positions"
#define NAMINOACIDS 20
#define MSB 160000		/* 20^(5-1) */

#elif INDEX1PART_AA == 6
#define PFXOFFSETS "pf6offsets"
#define PFXPOSITIONS "pf6positions"
#define PRXOFFSETS "pr6offsets"
#define PRXPOSITIONS "pr6positions"
#define NAMINOACIDS 20
#define MSB 3200000		/* 20^(6-1) */

#elif INDEX1PART_AA == 7
#define PFXOFFSETS "pfxoffsets"
#define PFXPOSITIONS "pfxpositions"
#define PRXOFFSETS "prxoffsets"
#define PRXPOSITIONS "prxpositions"
#ifdef NOCOLLAPSE
#define NAMINOACIDS 20
#define MSB 64000000		/* 20^(7-1) */
#else
#define NAMINOACIDS 15
#define MSB 11390625		/* 15^(7-1) */
#endif

#elif INDEX1PART_AA == 8
#define PFXOFFSETS "pf8offsets"
#define PFXPOSITIONS "pf8positions"
#define PRXOFFSETS "pr8offsets"
#define PRXPOSITIONS "pr8positions"
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
#define IDXOFFSETS "idxoffsets"
#define IDXPOSITIONS "idxpositions"
#define INDEX1PART 12		/* In nucleotides */
#define INDEX1INTERVAL 6
#endif	/* PMAP */



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
extern Genomicpos_T *
Indexdb_read_inplace (int *nentries, T this, Storedoligomer_T oligo);
#endif

extern void
Indexdb_write_offsets (FILE *offsets_fp, FILE *sequence_fp, IIT_T altstrain_iit,
		       int index1interval,
#ifdef PMAP
		       bool watsonp,
#endif
		       bool uncompressedp, char *fileroot);

extern void
Indexdb_write_positions (char *positionsfile, FILE *offsets_fp, FILE *sequence_fp, 
			 IIT_T altstrain_iit, int index1interval, bool uncompressedp, 
#ifdef PMAP
			 bool watsonp,
#endif
			 bool writefilep, char *fileroot);

#undef T
#endif

