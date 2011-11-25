#ifndef SPLICETRIE_BUILD_INCLUDED
#define SPLICETRIE_BUILD_INCLUDED

typedef enum {DONOR, ANTIDONOR, ACCEPTOR, ANTIACCEPTOR} Splicetype_T;

#include "bool.h"
#include "types.h"
#include "genomicpos.h"
#include "iit-read.h"
#include "genome.h"
#include "list.h"

/* For offsets */
/* #define USE_2BYTE_RELOFFSETS 1 */
#ifdef USE_2BYTE_RELOFFSETS
#define NULL_POINTER 65535
#else
#define NULL_POINTER -1U    /* Note: 0 does not work */
#endif

#define MAX_DUPLICATES 1000
#define DUPLICATE_NODE -1000U	/* Needs to be -MAX_DUPLICATES */

#define INTERNAL_NODE -1U

#define single_leaf_p(x) (x) < DUPLICATE_NODE
#define multiple_leaf_p(x) (x) < INTERNAL_NODE


extern char *
Splicetype_string (Splicetype_T splicetype);

extern void
Splicestring_gc (List_T *splicestrings, int nsplicesites);

extern Genomicpos_T *
Splicetrie_retrieve_via_splicesites (bool *distances_observed_p,
#ifdef GSNAP
				     UINT4 **splicecomp,
#endif
				     Splicetype_T **splicetypes, Genomicpos_T **splicedists,
				     List_T **splicestrings, UINT4 **splicefrags_ref, UINT4 **splicefrags_alt,
				     int *nsplicesites, IIT_T splicing_iit, int *splicing_divint_crosstable,
				     int donor_typeint, int acceptor_typeint, IIT_T chromosome_iit,
				     Genome_T genome, Genome_T genomealt, Genomicpos_T shortsplicedist);
extern Genomicpos_T *
Splicetrie_retrieve_via_introns (
#ifdef GSNAP
				 UINT4 **splicecomp,
#endif
				 Splicetype_T **splicetypes, Genomicpos_T **splicedists,
				 List_T **splicestrings, UINT4 **splicefrags_ref, UINT4 **splicefrags_alt,
				 int *nsplicesites, IIT_T splicing_iit, int *splicing_divint_crosstable,
				 IIT_T chromosome_iit, Genome_T genome, Genome_T genomealt);
extern void
Splicetrie_npartners (int **nsplicepartners_skip, int **nsplicepartners_obs, int **nsplicepartners_max,
		      Genomicpos_T *splicesites, Splicetype_T *splicetypes,
		      Genomicpos_T *splicedists, List_T *splicestrings, int nsplicesites,
		      IIT_T chromosome_iit, Genomicpos_T max_distance,
		      bool distances_observed_p);

extern void
Splicetrie_build_via_splicesites (unsigned int **triecontents_obs, unsigned int **trieoffsets_obs,
				  unsigned int **triecontents_max, unsigned int **trieoffsets_max,
				  int *nsplicepartners_skip, int *nsplicepartners_obs, int *nsplicepartners_max,
				  Splicetype_T *splicetypes, List_T *splicestrings, int nsplicesites);

extern void
Splicetrie_build_via_introns (unsigned int **triecontents_obs, unsigned int **trieoffsets_obs,
			      Genomicpos_T *splicesites, Splicetype_T *splicetypes,
			      List_T *splicestrings, int nsplicesites,
			      IIT_T chromosome_iit, IIT_T splicing_iit, int *splicing_divint_crosstable);

#endif

