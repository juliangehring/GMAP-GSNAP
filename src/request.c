static char rcsid[] = "$Id: request.c,v 1.18 2005/02/15 01:58:50 twu Exp $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "request.h"
#include "mem.h"

#define T Request_T
struct T {
  int id;
  Sequence_T queryseq;
  Sequence_T usersegment;
};


int
Request_id (T this) {
  return this->id;
}

Sequence_T
Request_queryseq (T this) {
  return this->queryseq;
}

Sequence_T
Request_usersegment (T this) {
  return this->usersegment;
}

T
Request_new (int id, Sequence_T queryseq, Sequence_T usersegment) {
  T new = (T) MALLOC(sizeof(*new));

  new->id = id;
  new->queryseq = queryseq;
  new->usersegment = usersegment;
  return new;
}

void
Request_free (T *old) {
  if (*old) {
    Sequence_free(&(*old)->queryseq);
    FREE(*old);
  }
  return;
}

