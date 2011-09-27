/* $Id: request.h,v 1.18 2005/02/07 23:56:57 twu Exp $ */
#ifndef REQUEST_INCLUDED
#define REQUEST_INCLUDED
#include "sequence.h"

#define T Request_T
typedef struct T *T;

extern int
Request_id (T this);
extern Sequence_T
Request_queryseq (T this);
extern Sequence_T
Request_usersegment (T this);
extern T
Request_new (int id, Sequence_T queryseq, Sequence_T usersegment);
extern void
Request_free (T *old);

#undef T
#endif
