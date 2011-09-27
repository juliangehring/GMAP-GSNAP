/* $Id: params.h,v 1.68 2006/11/28 00:50:11 twu Exp $ */
#ifndef PARAMS_INCLUDED
#define PARAMS_INCLUDED
#include "bool.h"
#include "genome.h"
#include "indexdb.h"
#include "iit-read.h"
#include "chrsubset.h"

#define T Params_T
typedef struct T *T;

extern Genome_T
Params_genome (T this);
extern IIT_T
Params_altstrain_iit (T this);
extern char *
Params_refstrain (T this);
#ifdef PMAP
extern Indexdb_T
Params_indexdb_fwd (T this);
extern Indexdb_T
Params_indexdb_rev (T this);
#else
extern Indexdb_T
Params_indexdb (T this);
#endif
extern IIT_T
Params_chromosome_iit (T this);
extern Chrsubset_T
Params_chrsubset (T this);
extern IIT_T
Params_contig_iit (T this);
extern IIT_T
Params_map_iit (T this);
extern bool
Params_map_iit_universal_p (T this);
extern int
Params_map_iit_forward_type (T this);
extern int
Params_map_iit_reverse_type (T this);
extern int
Params_stuttercycles (T this);
extern int
Params_stutterhits (T this);
extern int
Params_maxoligohits (T this);
extern int
Params_minindexsize (T this);
extern int
Params_maxindexsize (T this);
extern int
Params_maxpeelback (T this);
extern int
Params_nullgap (T this);
extern int
Params_sufflookback (T this);
extern int
Params_nsufflookback (T this);
extern int
Params_extramaterial_end (T this);
extern int
Params_extramaterial_paired (T this);
extern int
Params_extraband_single (T this);
extern int
Params_extraband_end (T this);
extern int
Params_extraband_paired (T this);
extern int
Params_maxmutations (T this);
extern T
Params_new (Genome_T genome, IIT_T altstrain_iit, 
#ifdef PMAP
	    Indexdb_T indexdb_fwd,
	    Indexdb_T indexdb_rev,
#else
	    Indexdb_T indexdb, 
#endif
	    IIT_T chromosome_iit, Chrsubset_T chrsubset, IIT_T contig_iit, IIT_T map_iit, 
	    bool map_iit_universal_p, int map_iit_forward_type, int map_iit_reverse_type,
	    int stuttercycles, int stutterhits, int maxoligohits,
	    int minindexsize, int maxindexsize, int maxpeelback,
	    int sufflookback, int nsufflookback, int nullgap, 
	    int extramaterial_end, int extramaterial_paired,
	    int extraband_single, int extraband_end, int extraband_paired,
	    int maxmutations);
extern void
Params_free (T *old);

#undef T
#endif

