static char rcsid[] = "$Id: result.c,v 1.57 2010-07-10 14:51:40 twu Exp $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "result.h"
#include <stdlib.h>
#include "mem.h"
#include "diag.h"


#define T Result_T
struct T {
  int id;
  int worker_id;

  Chimera_T chimera;		/* NULL indicates not a chimera */
  List_T gregionlist;		/* For debugging of stage 1 */
  List_T diagonals;		/* For debugging of diag */
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


List_T
Result_gregionlist (T this) {
  return this->gregionlist;
}

List_T
Result_diagonals (T this) {
  return this->diagonals;
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
  new->gregionlist = (List_T) NULL;
  new->diagonals = (List_T) NULL;
  new->array = array;
  new->npaths = npaths;
  new->diagnostic = diagnostic;
  new->failuretype = failuretype;

  return new;
}

T
Result_new_stage1debug (int id, int worker_id, List_T gregionlist,
			Diagnostic_T diagnostic, Failure_T failuretype) {
  T new = (T) MALLOC(sizeof(*new));

  new->id = id;
  new->worker_id = worker_id;
  new->chimera = (Chimera_T) NULL;
  new->gregionlist = gregionlist;
  new->diagonals = (List_T) NULL;
  new->array = (Stage3_T *) NULL;
  new->npaths = List_length(gregionlist);
  new->diagnostic = diagnostic;
  new->failuretype = failuretype;

  return new;
}

T
Result_new_diag_debug (int id, int worker_id, List_T diagonals,
		       Diagnostic_T diagnostic, Failure_T failuretype) {
  T new = (T) MALLOC(sizeof(*new));

  new->id = id;
  new->worker_id = worker_id;
  new->chimera = (Chimera_T) NULL;
  new->gregionlist = (List_T) NULL;
  new->diagonals = diagonals;
  new->array = (Stage3_T *) NULL;
  new->npaths = List_length(diagonals);
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
  List_T p;
  Gregion_T gregion;
#ifndef USE_DIAGPOOL
  Diag_T diag;
#endif

  if (*old) {
#ifdef BETATEST    
    if ((chimera = (*old)->chimera) != NULL) {
      Chimera_free(&chimera);
    }
#endif
    if ((*old)->diagnostic != NULL) {
      Diagnostic_free(&(*old)->diagnostic);
    }
    if ((*old)->gregionlist) {
      for (p = (*old)->gregionlist; p != NULL; p = List_next(p)) {
	gregion = (Gregion_T) List_head(p);
	Gregion_free(&gregion);
      }
      List_free(&((*old)->gregionlist));
    }

#ifndef USE_DIAGPOOL
    /* No need to free since memory is allocated separately */
    if ((*old)->diagonals) {
      for (p = (*old)->diagonals; p != NULL; p = List_next(p)) {
	diag = (Diag_T) List_head(p);
	Diag_free(&diag);
      }
      List_free(&((*old)->diagonals));
    }
#endif

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

