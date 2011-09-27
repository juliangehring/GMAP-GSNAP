static char rcsid[] = "$Id: result.c,v 1.56 2008/07/18 01:00:58 twu Exp $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "result.h"
#include <stdlib.h>
#include "mem.h"


#define T Result_T
struct T {
  int id;
  int worker_id;

  Chimera_T chimera;		/* NULL indicates not a chimera */
  Stage3_T *array;
  int npaths;
  Diagnostic_T diagnostic;
  Failure_T failuretype;
};


int
Result_id (T this) {
  return this->id;
}


int
Result_worker_id (T this) {
  return this->worker_id;
}


Chimera_T
Result_chimera (T this) {
  return this->chimera;
}


Stage3_T *
Result_array (int *npaths, T this) {
  *npaths = this->npaths;
  return this->array;
}


Diagnostic_T
Result_diagnostic (T this) {
  return this->diagnostic;
}


Failure_T
Result_failuretype (T this) {
  return this->failuretype;
}


T
Result_new (int id, int worker_id, Chimera_T chimera, Stage3_T *array,
	    int npaths, Diagnostic_T diagnostic, Failure_T failuretype) {
  T new = (T) MALLOC(sizeof(*new));

  new->id = id;
  new->worker_id = worker_id;
  new->chimera = chimera;
  new->array = array;
  new->npaths = npaths;
  new->diagnostic = diagnostic;
  new->failuretype = failuretype;

  return new;
}

T
Result_blank (int id, int worker_id) {
  T new = (T) MALLOC(sizeof(*new));

  new->id = id;
  new->worker_id = worker_id;
  new->chimera = NULL;
  new->array = NULL;

  return new;
}

void
Result_free (T *old) {
#ifdef BETATEST
  Chimera_T chimera;
#endif
  Stage3_T stage3;
  int i;

  if (*old) {
#ifdef BETATEST    
    if ((chimera = (*old)->chimera) != NULL) {
      Chimera_free(&chimera);
    }
#endif
    if ((*old)->diagnostic != NULL) {
      Diagnostic_free(&(*old)->diagnostic);
    }
    if ((*old)->array) {
      for (i = 0; i < (*old)->npaths; i++) {
	stage3 = (*old)->array[i];
	Stage3_free(&stage3);
      }
      FREE((*old)->array);
    }

    FREE(*old);
  }
  return;
}

