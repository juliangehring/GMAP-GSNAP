/* $Id: request.h,v 1.19 2008/01/08 01:43:38 twu Exp $ */
#ifndef REQUEST_INCLUDED
#define REQUEST_INCLUDED
#include "sequence.h"

#define T Request_T
typedef struct T *T;

extern int
Request_id (T this);

#ifdef GSNAP
extern Sequence_T
Request_queryseq1 (T this);
extern Sequence_T
Request_queryseq2 (T this);
extern T
Request_new (int id, Sequence_T queryseq1, Sequence_T queryseq2);
#else
extern Sequence_T
Request_queryseq (T this);
extern Sequence_T
Request_usersegment (T this);
extern T
Request_new (int id, Sequence_T queryseq, Sequence_T usersegment);
#endif

extern void
Request_free (T *old);

#undef T
#endif
