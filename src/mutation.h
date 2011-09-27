/* $Id: mutation.h,v 1.7 2005/02/07 23:56:56 twu Exp $ */
#ifndef MUTATION_INCLUDED
#define MUTATION_INCLUDED
#include "bool.h"

typedef enum {INVALID, SUBSTITUTION, INSERTION, DELETION} Muttype_T;

#define T Mutation_T
typedef struct T *T;

extern Muttype_T
Mutation_type (T this);
extern int
Mutation_refquerypos (T this);
extern int
Mutation_aapos (T this);
extern int
Mutation_naa (T this);
extern bool
Mutation_at_ends_p (T segmental, T single);

extern void
Mutation_print (T this);

extern int
Mutation_cmp (const void *a, const void *b);

extern T
Mutation_substitution_new (int refquerypos, int aapos, char aa_g, char aa_e);
extern T
Mutation_insertion_new (int refquerypos, int aapos, char aa_e, int ninsertions);
extern T
Mutation_deletion_new (int refquerypos, int deletionstart, char aa_g1, int deletionend, char aa_g2);


extern void
Mutation_free (T *old);

extern T
Mutation_merge_insertions (T mutation1, T mutation2);
extern T
Mutation_merge_deletions (T mutation1, T mutation2);

extern void
Mutation_invalidate (T this);
extern bool
Mutation_validp (T this);

#undef T
#endif

