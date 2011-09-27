/* $Id: reqpost.h,v 1.13 2005/12/06 02:17:18 twu Exp $ */
#ifndef REQPOST_INCLUDED
#define REQPOST_INCLUDED
#include "blackboard.h"
#include "request.h"
#include "result.h"
#include "params.h"

#define T Reqpost_T
typedef struct T *T;

extern T
Reqpost_new (Blackboard_T blackboard, int id, Params_T params);
extern void
Reqpost_free (T *old);
extern int
Reqpost_id (T this);
extern Params_T
Reqpost_params (T this);

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
