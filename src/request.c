static char rcsid[] = "$Id: request.c,v 1.19 2008/01/08 01:43:38 twu Exp $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "request.h"
#include "mem.h"

#define T Request_T
struct T {
  int id;
#ifdef GSNAP
  Sequence_T queryseq1;
  Sequence_T queryseq2;
#else
  Sequence_T queryseq;
  Sequence_T usersegment;
#endif  
};


int
Request_id (T this) {
  return this->id;
}

#ifdef GSNAP

Sequence_T
Request_queryseq1 (T this) {
  return this->queryseq1;
}

Sequence_T
Request_queryseq2 (T this) {
  return this->queryseq2;
}

T
Request_new (int id, Sequence_T queryseq1, Sequence_T queryseq2) {
  T new = (T) MALLOC(sizeof(*new));

  new->id = id;
  new->queryseq1 = queryseq1;
  new->queryseq2 = queryseq2;
  return new;
}

void
Request_free (T *old) {
  if (*old) {
    Sequence_free(&(*old)->queryseq1);
    if ((*old)->queryseq2) {
      Sequence_free(&(*old)->queryseq2);
    }
    FREE(*old);
  }
  return;
}

#else

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

#endif

