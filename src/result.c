static char rcsid[] = "$Id: result.c,v 1.47 2005/10/14 19:35:58 twu Exp $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "result.h"
#include <stdlib.h>
#include "mem.h"


#define T Result_T
struct T {
  int id;
  Chimera_T chimera;		/* NULL indicates not a chimera */
  Stage1_T stage1;
  Stage3_T *array;
  int npaths;
  Failure_T failuretype;
};


int
Result_id (T this) {
  return this->id;
}


Chimera_T
Result_chimera (T this) {
  return this->chimera;
}


Stage1_T
Result_stage1 (T this) {
  return this->stage1;
}


Stage3_T *
Result_array (int *npaths, T this) {
  *npaths = this->npaths;
  return this->array;
}


Failure_T
Result_failuretype (T this) {
  return this->failuretype;
}


T
Result_new (int id, Chimera_T chimera, Stage1_T stage1, 
	    Stage3_T *array, int npaths, Failure_T failuretype) {
  T new = (T) MALLOC(sizeof(*new));

  new->id = id;
  new->chimera = chimera;
  new->stage1 = stage1;
  new->array = array;
  new->npaths = npaths;
  new->failuretype = failuretype;
  return new;
}

void
Result_free (T *old) {
#ifdef BETATEST
  Chimera_T chimera;
#endif
  Stage1_T stage1;
  Stage3_T stage3;
  int i;

  if (*old) {
#ifdef BETATEST    
    if ((chimera = (*old)->chimera) != NULL) {
      Chimera_free(&chimera);
    }
#endif
    if ((stage1 = (*old)->stage1) != NULL) {
      Stage1_free(&stage1);
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

