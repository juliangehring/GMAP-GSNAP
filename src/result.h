/* $Id: result.h,v 1.36 2005/06/23 22:47:42 twu Exp $ */
#ifndef RESULT_INCLUDED
#define RESULT_INCLUDED

#ifndef BETATEST
typedef struct Chimera_T *Chimera_T;
#else
#include "chimera.h"
#endif

#include "stage1.h"
#include "stage3.h"

typedef enum {NO_FAILURE, EMPTY_SEQUENCE, REPETITIVE} Failure_T;

#define T Result_T
typedef struct T *T;

extern int
Result_id (T this);
extern Chimera_T
Result_chimera (T this);
extern Stage1_T
Result_stage1 (T this);
extern Stage3_T *
Result_array (int *npaths, T this);
extern Failure_T
Result_failuretype (T this);

extern T
Result_new (int id, Chimera_T chimera, Stage1_T stage1, 
	    Stage3_T *array, int npaths, Failure_T failuretype);
extern void
Result_free (T *old);

#undef T
#endif

