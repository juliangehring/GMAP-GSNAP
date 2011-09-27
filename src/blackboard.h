/* $Id: blackboard.h,v 1.15 2008/02/28 18:12:06 twu Exp $ */
#ifndef BLACKBOARD_INCLUDED
#define BLACKBOARD_INCLUDED
#include <stdio.h>

/* Avoid circularity because reqpost.h includes blackboard.h
#include "reqpost.h"
*/

#include "bool.h"
#include "request.h"
#ifdef GSNAP
#include "resulthr.h"
#else
#include "result.h"
#endif
#include "sequence.h"

#define T Blackboard_T
typedef struct T *T;

extern T
Blackboard_new (FILE *input, char **files, int nfiles, int nextchar,
		Sequence_T usersegment, int nworkers);
extern void
Blackboard_free (T *old);

extern bool
Blackboard_donep (T this);
extern FILE *
Blackboard_input (T this);
extern char **
Blackboard_files (T this);
extern int
Blackboard_nfiles (T this);
extern int
Blackboard_nextchar (T this);
extern Sequence_T
Blackboard_usersegment (T this);

/* extern Reqpost_T -- Avoid circularity because reqpost.h includes blackboard.h */
extern struct Reqpost_T *
Blackboard_get_reqpost (T this, int i);

extern void
Blackboard_put_request (T this, Request_T request);
extern void
Blackboard_set_inputdone (T this);
extern void
Blackboard_notify_output_ready (T this, int i);
extern Result_T
Blackboard_get_result (Request_T *request, T this);

#undef T
#endif

