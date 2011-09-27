/* $Id: result.h,v 1.33 2005/05/06 18:44:03 twu Exp $ */
#ifndef RESULT_INCLUDED
#define RESULT_INCLUDED
#include "stage1.h"
#include "stage3.h"

#define T Result_T
typedef struct T *T;

extern int
Result_id (T this);
extern int
Result_chimerapos (T this);
extern int
Result_nonchimera_matches (T this);
extern int
Result_nonchimera_mismatches (T this);
extern int
Result_nonchimera_indels (T this);
extern Stage1_T
Result_stage1 (T this);
extern Stage3_T *
Result_array (int *npaths, T this);

extern T
Result_new (int id, int chimerapos, int nonchimera_matches, int nonchimera_mismatches,
	    int nonchimera_indels, Stage1_T stage1, Stage3_T *array, int npaths);
extern void
Result_free (T *old);

#undef T
#endif

