/* $Id: assert.h,v 1.6 2010-07-27 16:21:32 twu Exp $ */
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
