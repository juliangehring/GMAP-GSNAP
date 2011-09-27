/* $Id: params.h,v 1.62 2005/07/08 07:58:34 twu Exp $ */
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
extern Indexdb_T
Params_indexdb (T this);
extern IIT_T
Params_chromosome_iit (T this);
extern Chrsubset_T
Params_chrsubset (T this);
extern IIT_T
Params_contig_iit (T this);
extern IIT_T
Params_map_iit (T this);
extern int
Params_maxextension (T this);
extern int
Params_stuttercycles (T this);
extern int
Params_stutterhits (T this);
extern int
Params_indexsize (T this);
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
extern T
Params_new (Genome_T genome, IIT_T altstrain_iit, Indexdb_T indexdb, 
	    IIT_T chromosome_iit, Chrsubset_T chrsubset, IIT_T contig_iit, IIT_T map_iit,
	    int maxextension, int stuttercycles, int stutterhits, 
	    int indexsize, int maxpeelback, int sufflookback, int nsufflookback, int nullgap, 
	    int extramaterial_end, int extramaterial_paired,
	    int extraband_single, int extraband_end, int extraband_paired);
extern void
Params_free (T *old);

#undef T
#endif

