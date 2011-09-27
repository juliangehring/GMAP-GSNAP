/* $Id: types.h 44854 2011-08-13 06:42:25Z twu $ */
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

/* An offset into the positions file of an IndexDB.  Typically, 3
   billion divided by sampling interval, requiring a maximum of 32
   bits or 4 bytes */
typedef UINT4 Positionsptr_T;

#endif

