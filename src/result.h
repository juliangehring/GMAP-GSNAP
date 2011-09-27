/* $Id: result.h,v 1.32 2005/02/15 01:58:50 twu Exp $ */
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
extern Stage1_T
Result_stage1 (T this);
extern Stage3_T *
Result_array (int *npaths, T this);

extern T
Result_new (int id, int chimerapos, Stage1_T stage1, Stage3_T *array, int npaths);
extern void
Result_free (T *old);

#undef T
#endif

