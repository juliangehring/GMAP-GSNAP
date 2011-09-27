/* $Id: assert.h,v 1.4 2005/02/07 23:56:55 twu Exp $ */
#ifndef ASSERT_INCLUDED
#define ASSERT_INCLUDED

#include "except.h"

#undef assert

/* #define NDEBUG */
#ifdef NDEBUG
#define assert(e) ((void) 0)
#else
extern void assert (int e);
#define assert(e) ((void) ((e) || (RAISE(Assert_Failed),0)))
#endif

#endif
