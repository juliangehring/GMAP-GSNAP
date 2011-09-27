/* $Id: result.h 40271 2011-05-28 02:29:18Z twu $ */
#ifndef RESULT_INCLUDED
#define RESULT_INCLUDED

#include "chimera.h"
#include "diagnostic.h"
#include "stage3.h"

typedef enum {NO_FAILURE, EMPTY_SEQUENCE, SHORT_SEQUENCE, POOR_SEQUENCE, REPETITIVE} Failure_T;

#define T Result_T
typedef struct T *T;

extern int
Result_id (T this);
extern Chimera_T
Result_chimera (T this);
extern Stage3_T *
Result_array (int *npaths, T this);
extern List_T
Result_gregionlist (T this);
extern List_T
Result_diagonals (T this);
extern Diagnostic_T
Result_diagnostic (T this);
extern Failure_T
Result_failuretype (T this);

extern T
Result_new (int id, Chimera_T chimera, Stage3_T *array,
	    int npaths, Diagnostic_T diagnostic, Failure_T failuretype);
extern T
Result_new_stage1debug (int id, List_T gregionlist,
			Diagnostic_T diagnostic, Failure_T failuretype);
extern T
Result_new_diag_debug (int id, List_T diagonals,
		       Diagnostic_T diagnostic, Failure_T failuretype);
extern void
Result_free (T *old);

#undef T
#endif

