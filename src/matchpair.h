/* $Id: matchpair.h,v 1.30 2006/11/28 00:52:21 twu Exp $ */
#ifndef MATCHPAIR_INCLUDED
#define MATCHPAIR_INCLUDED
#include "bool.h"
#include "match.h"
#include "list.h"
#include "sequence.h"
#include "genomicpos.h"
#include "stage1.h"
#include "pairpool.h"
#include "genome.h"
#include "gbuffer.h"

typedef enum {FIVEONLY, THREEONLY, MIXED} Matchpairend_T;
typedef enum {NO_BOUND, BY_CANDIDATES, MOVING_THRESHOLD} Boundmethod_T;

#define T Matchpair_T
typedef struct T *T;

extern Match_T
Matchpair_bound5 (T this);
extern Match_T
Matchpair_bound3 (T this);
extern unsigned int
Matchpair_extension5 (T this);
extern unsigned int
Matchpair_extension3 (T this);
extern bool
Matchpair_usep (T this);
extern void
Matchpair_set_usep (T this);
extern void
Matchpair_clear_usep (T this);
extern bool
Matchpair_watsonp (T this);
extern int
Matchpair_clustersize (T this);
extern Matchpairend_T
Matchpair_matchpairend (T this);
extern int
Matchpair_support (T this);
extern double
Matchpair_fsupport (T this);
extern double
Matchpair_compute_stretch (Match_T bound5, Match_T bound3);

extern T
Matchpair_new (Match_T match5, Match_T match3, int stage1size, int clustersize, 
	       int trimstart, int trimend, int trimlength, Matchpairend_T matchpairend);
extern void
Matchpair_free (T *old);

extern int
Matchpair_size_cmp (const void *a, const void *b);
extern int
Matchpair_support_cmp (const void *a, const void *b);
extern List_T
Matchpair_filter_unique (List_T matchpairlist, bool removedupsp);
extern bool
Matchpair_sufficient_support (T this);
extern void
Matchpair_continue (T this, Stage1_T stage1, int maxextension);
extern void
Matchpair_get_coords (Genomicpos_T *chrpos, Genomicpos_T *genomicstart, Genomicpos_T *genomiclength,
		      bool *watsonp, T this, Stage1_T stage1, Genomicpos_T chrlength, int maxextension,
		      bool maponlyp);

extern List_T
Matchpair_find_clusters (List_T matches5, List_T matches3, int stage1size, 
			 int maxintronlen, int minclustersize, double sizebound, 
			 int trimstart, int trimend, int trimlength, Boundmethod_T boundmethod);
extern List_T
Matchpair_salvage_hits (List_T matches5, List_T matches3, int stage1size, int trimstart, int trimend, int trimlength);

extern List_T
Matchpair_make_path (T this, char *queryseq_ptr, Pairpool_T pairpool, Genome_T genome, 
		     Genomicpos_T chroffset, Genomicpos_T chrpos, Gbuffer_T gbuffer);

#undef T
#endif
