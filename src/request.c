static char rcsid[] = "$Id: request.c 40330 2011-05-30 17:40:46Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "request.h"
#include "mem.h"

#define T Request_T
struct T {
  int id;
#ifdef GSNAP
  Shortread_T queryseq1;
  Shortread_T queryseq2;
#else
  Sequence_T queryseq;
#endif  
};


int
Request_id (T this) {
  return this->id;
}

#ifdef GSNAP

Shortread_T
Request_queryseq1 (T this) {
  return this->queryseq1;
}

Shortread_T
Request_queryseq2 (T this) {
  return this->queryseq2;
}


T
Request_new (int id, Shortread_T queryseq1, Shortread_T queryseq2) {
  T new = (T) MALLOC_IN(sizeof(*new));

  new->id = id;
  new->queryseq1 = queryseq1;
  new->queryseq2 = queryseq2;
  return new;
}

void
Request_free (T *old) {
  if (*old) {
    Shortread_free(&(*old)->queryseq1);
    if ((*old)->queryseq2) {
      Shortread_free(&(*old)->queryseq2);
    }
    FREE_IN(*old);
  }
  return;
}

#else

Sequence_T
Request_queryseq (T this) {
  return this->queryseq;
}

T
Request_new (int id, Sequence_T queryseq) {
  T new = (T) MALLOC_IN(sizeof(*new));

  new->id = id;
  new->queryseq = queryseq;
  return new;
}

void
Request_free (T *old) {
  if (*old) {
    Sequence_free(&(*old)->queryseq);
    FREE_IN(*old);
  }
  return;
}

#endif

