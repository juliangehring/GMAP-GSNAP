/* $Id: resulthr.h 36798 2011-03-18 20:33:23Z twu $ */
#ifndef RESULTHR_INCLUDED
#define RESULTHR_INCLUDED
#include "bool.h"

typedef enum {SINGLEEND_NOMAPPING, PAIREDEND_NOMAPPING,
	      SINGLEEND_UNIQ, SINGLEEND_MULT,
	      PAIRED_UNIQ, PAIRED_MULT,
	      CONCORDANT_UNIQ, CONCORDANT_MULT,
	      HALFMAPPING_UNIQ, HALFMAPPING_MULT,
	      UNPAIRED_UNIQ, UNPAIRED_MULT} Resulttype_T;

#define T Result_T
typedef struct T *T;

extern Resulttype_T
Result_resulttype (T this);
extern bool
Result_translocationp (T this);
extern char *
Resulttype_string (Resulttype_T resulttype);
extern int
Result_id (T this);
extern int
Result_worker_id (T this);
extern void **
Result_array (int *npaths, T this);
extern void **
Result_array2 (int *npaths, T this);
extern T
Result_single_read_new (int id, void **resultarray, int npaths);
extern T
Result_paired_read_new (int id, void **resultarray, int npaths, bool concordantp);
extern T
Result_paired_as_singles_new (int id, void **hits5, int npaths5, void **hits3, int npaths3);
extern void
Result_free (T *old);

#undef T
#endif

