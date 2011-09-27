static char rcsid[] = "$Id: datum.c,v 1.1 2005/06/15 22:38:49 twu Exp $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "datum.h"
#include <stdlib.h>
#include "mem.h"


#define T Datum_T
struct T {
  unsigned int chrpos;
  double value;
};


unsigned int
Datum_chrpos (T this) {
  return this->chrpos;
}

double
Datum_value (T this) {
  return this->value;
}


T
Datum_new (unsigned int chrpos, double value) {
  T new = (T) MALLOC(sizeof(*new));

  new->chrpos = chrpos;
  new->value = value;
  return new;
}

void
Datum_free (T *old) {
  FREE(*old);
  return;
}


int
Datum_cmp (const void *a, const void *b) {
  T x = * (T *) a;
  T y = * (T *) b;

  if (x->chrpos < y->chrpos) {
    return -1;
  } else if (x->chrpos > y->chrpos) {
    return +1;
  } else {
    return 0;
  }
}


