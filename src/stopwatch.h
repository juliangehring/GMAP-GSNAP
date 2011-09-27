/* $Id: stopwatch.h,v 1.6 2006-10-04 19:24:58 twu Exp $ */
#ifndef STOPWATCH_INCLUDED
#define STOPWATCH_INCLUDED

#define T Stopwatch_T
typedef struct T *T;

extern T
Stopwatch_new ();
extern void
Stopwatch_free (T *old);
extern void
Stopwatch_start (T this);
extern double 
Stopwatch_stop (T this);

#undef T
#endif


