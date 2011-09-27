/* $Id: stage1hr.h,v 1.62 2010-07-27 00:01:18 twu Exp $ */
#ifndef STAGE1HR_INCLUDED
#define STAGE1HR_INCLUDED
#include "bool.h"
#include "genomicpos.h"
#include "indexdb.h"
#include "sequence.h"
#include "list.h"
#include "iit-read.h"
#include "genome.h"


#define MAX_QUERYLENGTH 200

typedef enum {MASK_NONE, MASK_FREQUENT, MASK_REPETITIVE, MASK_GREEDY_FREQUENT, MASK_GREEDY_REPETITIVE} Masktype_T;
typedef enum {DONOR, ANTIDONOR, ACCEPTOR, ANTIACCEPTOR} Splicetype_T;


typedef struct Floors_T *Floors_T;

extern void
Floors_free (Floors_T *old);


#define T Stage1_T
typedef struct T *T;

extern void
Stage1_init_positions_free (bool positions_fileio_p);
extern void
Stage1_free (T *old, int querylength);


extern Genomicpos_T *
Stage1_retrieve_splicesites (Splicetype_T **splicetypes, int *nsplicesites,
			     IIT_T splicesites_iit, int *splicesites_divint_crosstable,
			     int donor_typeint, int acceptor_typeint, IIT_T chromosome_iit);
extern List_T
Stage1_single_read (Sequence_T queryseq, Indexdb_T indexdb, Indexdb_T indexdb2,
		    int indexdb_size_threshold, IIT_T geneprob_iit, IIT_T chromosome_iit, Genome_T genome,
		    Genome_T genomealt, Floors_T *floors_array,
		    bool knownsplicingp, bool novelsplicingp, bool canonicalp, int trim_maxlength,
		    int maxpaths, int maxchimerapaths, double usermax_level_float, int subopt_levels,
		    Masktype_T masktype, int indel_penalty, int max_middle_insertions, int max_middle_deletions,
		    bool allow_end_indels_p, int max_end_insertions, int max_end_deletions, int min_indel_end_matches,
		    Genomicpos_T shortsplicedist,
		    int localsplicing_penalty, int distantsplicing_penalty, int min_localsplicing_end_matches,
		    int min_distantsplicing_end_matches, double min_distantsplicing_identity,
		    Genomicpos_T *splicesites, Splicetype_T *splicetypes, int nsplicesites,
		    IIT_T splicesites_iit, int *splicesites_divint_crosstable,
		    int donor_typeint, int acceptor_typeint, bool dibasep, bool cmetp);
extern List_T
Stage1_paired_read (List_T *singlehits5, List_T *singlehits3,
		    Sequence_T queryseq5, Sequence_T queryseq3,
		    Indexdb_T indexdb, Indexdb_T indexdb2, int indexdb_size_threshold,
		    IIT_T geneprob_iit, IIT_T chromosome_iit, Genome_T genome, Genome_T genomealt, Floors_T *floors_array,
		    bool knownsplicingp, bool novelsplicingp, bool canonicalp, int trim_maxlength,
		    int maxpaths, int maxchimerapaths, double usermax_level_float, int subopt_levels, 
		    Masktype_T masktype, int indel_penalty, int max_middle_insertions, int max_middle_deletions,
		    bool allow_end_indels_p, int max_end_insertions, int max_end_deletions, int min_indel_end_matches,
		    Genomicpos_T shortsplicedist,
		    int localsplicing_penalty, int distantsplicing_penalty, int min_localsplicing_end_matches,
		    int min_distantsplicing_end_matches, double min_distantsplicing_identity,
		    Genomicpos_T *splicesites, Splicetype_T *splicetypes, int nsplicesites,
		    IIT_T splicesites_iit, int *splicesites_divint_crosstable,
		    int donor_typeint, int acceptor_typeint,
		    bool dibasep, bool cmetp, int pairmax, int expected_pairlength);

#undef T
#endif

