/* $Id: assert.h 27450 2010-08-05 19:02:48Z twu $ */
#ifndef ASSERT_INCLUDED
#define ASSERT_INCLUDED

#include "except.h"

#undef assert

/* #define CHECK_ASSERTIONS 1 */
#ifdef CHECK_ASSERTIONS
extern void assert (int e);
#define assert(e) ((void) ((e) || (RAISE(Assert_Failed),0)))
#else
#define assert(e) ((void) 0)
#endif

#endif
