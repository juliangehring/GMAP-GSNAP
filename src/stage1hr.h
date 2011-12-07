/* $Id: stage1hr.h 52862 2011-11-20 23:53:04Z twu $ */
#ifndef STAGE1HR_INCLUDED
#define STAGE1HR_INCLUDED
#include "bool.h"
#include "types.h"
#include "mode.h"
#include "genomicpos.h"
#include "indexdb.h"
#include "shortread.h"
#include "iit-read.h"
#include "genome.h"
#include "splicetrie.h"
#include "resulthr.h"		/* For Pairtype_T */
#include "stage3hr.h"

#include "oligoindex.h"
#include "pairpool.h"
#include "diagpool.h"
#include "dynprog.h"


#if 0
/* Obsolete */
typedef enum {MASK_NONE, MASK_FREQUENT, MASK_REPETITIVE, MASK_GREEDY_FREQUENT, MASK_GREEDY_REPETITIVE} Masktype_T;
#endif

#define GMAP_IMPROVEMENT 1
#define GMAP_TERMINAL 2
#define GMAP_PAIRSEARCH 4


typedef struct Floors_T *Floors_T;

extern void
Floors_free (Floors_T *old);
extern void
Floors_free_keep (Floors_T *old);


#define T Stage1_T
typedef struct T *T;

extern void
Stage1_init_positions_free (bool positions_fileio_p);
extern void
Stage1_free (T *old, int querylength);


extern Stage3end_T *
Stage1_single_read (int *npaths, int *second_mapq,
		    Shortread_T queryseq, Indexdb_T indexdb, Indexdb_T indexdb2,
		    int indexdb_size_threshold, Genome_T genome, Floors_T *floors_array,
		    int maxpaths, double usermax_level_float, int subopt_levels,
		    int indel_penalty_middle, int indel_penalty_end, int max_middle_insertions, int max_middle_deletions,
		    bool allow_end_indels_p, int max_end_insertions, int max_end_deletions, int min_indel_end_matches,
		    Genomicpos_T shortsplicedist,
		    int localsplicing_penalty, int distantsplicing_penalty,
		    int min_distantsplicing_end_matches, double min_distantsplicing_identity, int min_shortend,
		    Oligoindex_T *oligoindices_major, int noligoindices_major,
		    Oligoindex_T *oligoindices_minor, int noligoindices_minor,
		    Pairpool_T pairpool, Diagpool_T diagpool,
		    Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,
		    bool keep_floors_p);

extern Stage3pair_T *
Stage1_paired_read (int *npaths, int *second_mapq, Pairtype_T *final_pairtype,
		    Stage3end_T **stage3array5, int *nhits5, int *second_mapq5,
		    Stage3end_T **stage3array3, int *nhits3, int *second_mapq3,
		    Shortread_T queryseq5, Shortread_T queryseq3,
		    Indexdb_T indexdb, Indexdb_T indexdb2, int indexdb_size_threshold,
		    Genome_T genome, Floors_T *floors_array,
		    int maxpaths, double usermax_level_float, int subopt_levels,
		    int indel_penalty_middle, int indel_penalty_end, int max_middle_insertions, int max_middle_deletions,
		    bool allow_end_indels_p, int max_end_insertions, int max_end_deletions, int min_indel_end_matches,
		    Genomicpos_T shortsplicedist, int localsplicing_penalty, int distantsplicing_penalty,
		    int min_distantsplicing_end_matches, double min_distantsplicing_identity, int min_shortend,
		    Oligoindex_T *oligoindices_major, int noligoindices_major,
		    Oligoindex_T *oligoindices_minor, int noligoindices_minor,
		    Pairpool_T pairpool, Diagpool_T diagpool,
		    Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,
		    Genomicpos_T pairmax, bool keep_floors_p);

extern void
Stage1hr_setup (int index1part_in, IIT_T chromosome_iit_in, int nchromosomes_in,
		Genome_T genomealt, Mode_T mode_in,
		int terminal_threshold_in,

		Genomicpos_T *splicesites_in, Splicetype_T *splicetypes_in,
		Genomicpos_T *splicedists_in, int nsplicesites_in,

		bool novelsplicingp_in, bool knownsplicingp_in,	int shortsplicedist_known_in,

		int nullgap_in, int maxpeelback_in, int maxpeelback_distalmedial_in,
		int extramaterial_end_in, int extramaterial_paired_in,
		int gmap_mode, int trigger_score_for_gmap_in,
		int max_gmap_pairsearch_in, int max_gmap_terminal_in,
		int max_gmap_improvement_in, int antistranded_penalty_in);


#undef T
#endif

