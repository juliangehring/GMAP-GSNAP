#ifndef SPLICETRIE_INCLUDED
#define SPLICETRIE_INCLUDED

#include "bool.h"
#include "types.h"
#include "genomicpos.h"
#include "sequence.h"
#include "iit-read.h"
#include "genome.h"
#include "compress.h"
#include "intlist.h"

typedef enum {DONOR, ANTIDONOR, ACCEPTOR, ANTIACCEPTOR} Splicetype_T;

extern void
Splicestring_gc (List_T *splicestrings, int nsplicesites);

extern Genomicpos_T *
Splicetrie_retrieve_splicesites (bool *distances_observed_p, Splicetype_T **splicetypes, Genomicpos_T **splicedists,
				 List_T **splicestrings, UINT4 **splicefrags_ref, UINT4 **splicefrags_alt,
				 int *nsplicesites, IIT_T splicesites_iit, int *splicesites_divint_crosstable,
				 int donor_typeint, int acceptor_typeint, IIT_T chromosome_iit,
				 Genome_T genome, Genome_T genomealt, Genomicpos_T shortsplicedist);

extern void
Splicetrie_npartners (int **nsplicepartners_skip, int **nsplicepartners_obs, int **nsplicepartners_max,
		      Genomicpos_T *splicesites, Splicetype_T *splicetypes,
		      Genomicpos_T *splicedists, List_T *splicestrings, int nsplicesites,
		      IIT_T chromosome_iit, Genomicpos_T max_distance,
		      bool distances_observed_p);

extern void
Splicetrie_build (unsigned int **triecontents_obs, unsigned int **trieoffsets_obs,
		  unsigned int **triecontents_max, unsigned int **trieoffsets_max,
		  int *nsplicepartners_skip, int *nsplicepartners_obs, int *nsplicepartners_max,
		  Splicetype_T *splicetypes, List_T *splicestrings, int nsplicesites);

extern Intlist_T
Splicetrie_find_left (int *best_nmismatches, Intlist_T *nmismatches_list, int i,
		      Genomicpos_T origleft, int pos5, int pos3,
		      Genomicpos_T *splicesites, UINT4 *splicefrags_ref, UINT4 *splicefrags_alt, int *nsplicepartners_skip,
		      unsigned int *trieoffsets_obs, unsigned int *triecontents_obs, int *nsplicepartners_obs,
		      unsigned int *trieoffsets_max, unsigned int *triecontents_max, int *nsplicepartners_max,
		      Splicetype_T *splicetypes, List_T *splicestrings,
		      Compress_T query_compress, UINT4 *genome_blocks, UINT4 *snp_blocks,
		      char *query, char *queryptr, int max_mismatches_allowed,
		      bool dibasep, bool cmetp, bool plusp);

extern Intlist_T
Splicetrie_find_right (int *best_nmismatches, Intlist_T *nmismatches_list, int i,
		       Genomicpos_T origleft, int pos5, int pos3,
		       Genomicpos_T *splicesites, UINT4 *splicefrags_ref, UINT4 *splicefrags_alt, int *nsplicepartners_skip,
		       unsigned int *trieoffsets_obs, unsigned int *triecontents_obs, int *nsplicepartners_obs,
		       unsigned int *trieoffsets_max, unsigned int *triecontents_max, int *nsplicepartners_max,
		       Splicetype_T *splicetypes, List_T *splicestrings,
		       Compress_T query_compress, UINT4 *genome_blocks, UINT4 *snp_blocks,
		       char *query, char *queryptr, int max_mismatches_allowed,
		       bool dibasep, bool cmetp, bool plusp);

#endif
