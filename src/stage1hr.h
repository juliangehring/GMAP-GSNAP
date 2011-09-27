/* $Id: stage1hr.h 36302 2011-03-09 05:29:24Z twu $ */
#ifndef STAGE1HR_INCLUDED
#define STAGE1HR_INCLUDED
#include "bool.h"
#include "types.h"
#include "genomicpos.h"
#include "indexdb.h"
#include "shortread.h"
#include "list.h"
#include "iit-read.h"
#include "genome.h"
#include "splicetrie.h"
#include "stage3hr.h"


#define MAX_QUERYLENGTH 500

typedef enum {MASK_NONE, MASK_FREQUENT, MASK_REPETITIVE, MASK_GREEDY_FREQUENT, MASK_GREEDY_REPETITIVE} Masktype_T;


typedef struct Floors_T *Floors_T;

extern void
Floors_free (Floors_T *old);


#define T Stage1_T
typedef struct T *T;

extern void
Stage1_init_positions_free (bool positions_fileio_p);
extern void
Stage1_free (T *old, int querylength);


extern Stage3_T *
Stage1_single_read (int *npaths, Shortread_T queryseq, Indexdb_T indexdb, Indexdb_T indexdb2,
		    int indexdb_size_threshold, IIT_T chromosome_iit, Genome_T genome,
		    Genome_T genomealt, Floors_T *floors_array,
		    bool knownsplicingp, bool novelsplicingp, bool canonicalp,
		    int maxpaths, int maxchimerapaths, double usermax_level_float, int subopt_levels,
		    Masktype_T masktype, int terminal_penalty, int max_terminal_length,
		    int indel_penalty, int max_middle_insertions, int max_middle_deletions,
		    bool allow_end_indels_p, int max_end_insertions, int max_end_deletions, int min_indel_end_matches,
		    Genomicpos_T shortsplicedist,
		    int localsplicing_penalty, int distantsplicing_penalty, int min_localsplicing_end_matches,
		    int min_distantsplicing_end_matches, double min_distantsplicing_identity, int min_shortend,
		    bool find_novel_doublesplices_p, Genomicpos_T *splicesites, Splicetype_T *splicetypes,
		    Genomicpos_T *splicedists, UINT4 *splicefrags_ref, UINT4 *splicefrags_alt, int nsplicesites,
		    int *nsplicepartners_skip,
		    unsigned int *trieoffsets_obs, unsigned int *triecontents_obs, int *nsplicepartners_obs,
		    unsigned int *trieoffsets_max, unsigned int *triecontents_max, int *nsplicepartners_max,
		    List_T *splicestrings, bool dibasep, bool cmetp);

extern Stage3pair_T *
Stage1_paired_read (int *npaths, bool *concordantp, Stage3_T **stage3array5, int *nhits5, Stage3_T **stage3array3, int *nhits3,
		    Shortread_T queryseq5, Shortread_T queryseq3,
		    Indexdb_T indexdb, Indexdb_T indexdb2, int indexdb_size_threshold,
		    IIT_T chromosome_iit, Genome_T genome, Genome_T genomealt, Floors_T *floors_array,
		    bool knownsplicingp, bool novelsplicingp, bool canonicalp,
		    int maxpaths, int maxchimerapaths, double usermax_level_float, int subopt_levels, 
		    Masktype_T masktype, int terminal_penalty, int max_terminal_length,
		    int indel_penalty, int max_middle_insertions, int max_middle_deletions,
		    bool allow_end_indels_p, int max_end_insertions, int max_end_deletions, int min_indel_end_matches,
		    Genomicpos_T shortsplicedist,
		    int localsplicing_penalty, int distantsplicing_penalty, int min_localsplicing_end_matches,
		    int min_distantsplicing_end_matches, double min_distantsplicing_identity, int min_shortend,
		    bool find_novel_doublesplices_p, Genomicpos_T *splicesites, Splicetype_T *splicetypes,
		    Genomicpos_T *splicedists, UINT4 *splicefrags_ref, UINT4 *splicefrags_alt, int nsplicesites,
		    int *nsplicepartners_skip,
		    unsigned int *trieoffsets_obs, unsigned int *triecontents_obs, int *nsplicepartners_obs,
		    unsigned int *trieoffsets_max, unsigned int *triecontents_max, int *nsplicepartners_max,
		    List_T *splicestrings, Genomicpos_T pairmax, Genomicpos_T expected_pairlength,
		    Genomicpos_T pairlength_deviation, bool dibasep, bool cmetp);

#undef T
#endif

