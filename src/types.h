/* $Id: types.h,v 1.7 2009-02-20 22:19:25 twu Exp $ */
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

