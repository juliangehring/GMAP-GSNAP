static char rcsid[] = "$Id: diag.c,v 1.5 2007/02/07 08:35:27 twu Exp $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "diag.h"
#include "diagdef.h"
#include <stdlib.h>


#define T Diag_T

int
Diag_diagonal (T this) {
  return this->diagonal;
}

int
Diag_querystart (T this) {
  return this->querystart;
}

int
Diag_queryend (T this) {
  return this->queryend;
}

int
Diag_nconsecutive (T this) {
  return this->nconsecutive;
}

bool
Diag_dominatedp (T this) {
  return this->dominatedp;
}

void
Diag_set_dominatedp (T this) {
  this->dominatedp = true;
}


int
Diag_compare_nconsecutive (const void *x, const void *y) {
  T a = * (T *) x;
  T b = * (T *) y;

  if (a->nconsecutive > b->nconsecutive) {
    return -1;
  } else if (b->nconsecutive > a->nconsecutive) {
    return 1;
  } else {
    return 0;
  }
}

int
Diag_compare_diagonal (const void *x, const void *y) {
  T a = * (T *) x;
  T b = * (T *) y;

  if (a->diagonal < b->diagonal) {
    return -1;
  } else if (b->diagonal < a->diagonal) {
    return +1;
  } else {
    return 0;
  }
}

