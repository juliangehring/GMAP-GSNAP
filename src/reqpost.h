/* $Id: reqpost.h,v 1.15 2008-02-28 18:12:06 twu Exp $ */
#ifndef REQPOST_INCLUDED
#define REQPOST_INCLUDED
#include "blackboard.h"
#include "request.h"
#ifdef GSNAP
#include "resulthr.h"
#else
#include "result.h"
#endif

#define T Reqpost_T
typedef struct T *T;

extern T
Reqpost_new (Blackboard_T blackboard, int id);
extern void
Reqpost_free (T *old);
extern int
Reqpost_id (T this);

extern void
Reqpost_put_request (T this, Request_T request);
extern Request_T
Reqpost_get_request (T this);
extern void
Reqpost_put_result (T this, Result_T result);
Result_T
Reqpost_get_result (Request_T *request, T this);
extern void
Reqpost_reset (T this);

#undef T
#endif
