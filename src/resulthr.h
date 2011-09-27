/* $Id: resulthr.h,v 1.12 2010/03/05 05:52:28 twu Exp $ */
#ifndef RESULTHR_INCLUDED
#define RESULTHR_INCLUDED
#include "bool.h"

typedef enum {SINGLEEND_READ, PAIREDEND_CONCORDANT, PAIREDEND_SAMECHR, 
	      PAIREDEND_AS_SINGLES, PAIREDEND_AS_SINGLES_UNIQUE} Resulttype_T;

#define T Result_T
typedef struct T *T;

extern Resulttype_T
Result_resulttype (T this);
extern int
Result_id (T this);
extern int
Result_worker_id (T this);
extern void **
Result_array (int *npaths, T this);
extern void **
Result_array2 (int *npaths, T this);
extern T
Result_single_read_new (int id, int worker_id, void **resultarray, int npaths);
extern T
Result_paired_read_new (int id, int worker_id, void **resultarray, int npaths);
extern T
Result_paired_as_singles_new (int id, int worker_id, void **hits5, int npaths5, void **hits3, int npaths3);
extern T
Result_new (int id, int worker_id, void **resultarray, int npaths, bool pairedp);
extern void
Result_free (T *old);

#undef T
#endif

