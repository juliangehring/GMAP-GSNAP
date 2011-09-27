/* $Id: types.h 27450 2010-08-05 19:02:48Z twu $ */
#ifndef TYPES_INCLUDED
#define TYPES_INCLUDED
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/* A 2-byte word */
typedef unsigned short UINT2;

/* A 4-byte word */
typedef unsigned int UINT4;
typedef int INT4;

/* An 8-byte word */
#if (SIZEOF_UNSIGNED_LONG == 8)
#define HAVE_64_BIT
typedef unsigned long UINT8;
#elif (SIZEOF_UNSIGNED_LONG_LONG == 8)
#define HAVE_64_BIT
typedef unsigned long long UINT8;
#endif

#endif

