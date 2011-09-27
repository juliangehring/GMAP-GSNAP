/* $Id: result.h,v 1.45 2008/12/13 02:06:45 twu Exp $ */
#ifndef RESULT_INCLUDED
#define RESULT_INCLUDED

#ifndef BETATEST
typedef struct Chimera_T *Chimera_T;
#else
#include "chimera.h"
#endif

#include "diagnostic.h"
#include "stage3.h"

typedef enum {NO_FAILURE, EMPTY_SEQUENCE, SHORT_SEQUENCE, POOR_SEQUENCE, REPETITIVE} Failure_T;

#define T Result_T
typedef struct T *T;

extern int
Result_id (T this);
extern int
Result_worker_id (T this);
extern Chimera_T
Result_chimera (T this);
extern Stage3_T *
Result_array (int *npaths, T this);
extern Diagnostic_T
Result_diagnostic (T this);
extern Failure_T
Result_failuretype (T this);

extern T
Result_new (int id, int worker_id, Chimera_T chimera, Stage3_T *array,
	    int npaths, Diagnostic_T diagnostic, Failure_T failuretype);
extern T
Result_blank (int id, int worker_id);
extern void
Result_free (T *old);

#undef T
#endif

