static char rcsid[] = "$Id: result.c,v 1.42 2005/02/07 23:56:57 twu Exp $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "result.h"
#include <stdlib.h>
#include "mem.h"


#define T Result_T
struct T {
  int id;
  int chimerapos;		/* -1 indicates not a chimera */
  Stage1_T stage1;
  Stage3_T *array;
  int npaths;
};


int
Result_id (T this) {
  return this->id;
}


int
Result_chimerapos (T this) {
  return this->chimerapos;
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


T
Result_new (int id, int chimerapos, Stage1_T stage1, Stage3_T *array, int npaths) {
  T new = (T) MALLOC(sizeof(*new));

  new->id = id;
  new->chimerapos = chimerapos;
  new->stage1 = stage1;
  new->array = array;
  new->npaths = npaths;
  return new;
}

void
Result_free (T *old) {
  Stage1_T stage1;
  Stage3_T stage3;
  int i;

  if (*old) {
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

