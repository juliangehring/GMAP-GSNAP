/* $Id: match.h 40271 2011-05-28 02:29:18Z twu $ */
#ifndef MATCH_INCLUDED
#define MATCH_INCLUDED
#include "bool.h"
#include "genomicpos.h"
#include "iit-read.h"
#include "list.h"
#include "genome.h"
#include "chrnum.h"

#define T Match_T
typedef struct T *T;

extern int
Match_querypos (T this);
extern bool
Match_forwardp (T this);
extern bool
Match_fivep (T this);
extern Genomicpos_T
Match_position (T this);
extern Chrnum_T
Match_chrnum (T this);
extern char *
Match_chr (T this, IIT_T chromosome_iit);
extern Genomicpos_T
Match_chrpos (T this);
extern int
Match_incr_npairings (T this);
extern int
Match_decr_npairings (T this);
extern int
Match_npairings (T this);
extern void
Match_set_weight (T this, double weight);
extern double
Match_weight (T this);
extern bool
Match_has_weight_p (T this);

extern int
Match_cmp (const void *a, const void *b);


#ifndef USE_MATCHPOOL
extern T
Match_new (int querypos, bool forwardp, bool fivep,
	   Genomicpos_T position, IIT_T chromosome_iit);
extern void
Match_free (T *old);
#endif

extern void
Match_print_mer (T this, char *queryseq_ptr, Genome_T genome, IIT_T chromosome_iit, int stage1size);
extern void
Match_print (T this, IIT_T chromosome_iit);

extern bool
Match_acceptable_pair (T match5, T match3, int trimlength, int stage1size);


#undef T
#endif

