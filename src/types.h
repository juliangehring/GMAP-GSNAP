/* $Id: types.h 89124 2013-03-13 22:23:20Z twu $ */
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
/* Oligospace_T needs to hold 1 more than maximum Storedoligomer_T.
   If 8-byte words are not available, then maximum k-mer is 15 */
#if (SIZEOF_UNSIGNED_LONG == 8)
#define HAVE_64_BIT
#define MAXIMUM_KMER 16
typedef unsigned long UINT8;
typedef unsigned long Oligospace_T;
#elif (SIZEOF_UNSIGNED_LONG_LONG == 8)
#define HAVE_64_BIT
#define MAXIMUM_KMER 16
typedef unsigned long long UINT8;
typedef unsigned long long Oligospace_T;
#else
#define MAXIMUM_KMER 15
#define OLIGOSPACE_NOT_LONG
typedef unsigned int Oligospace_T;
#endif

/* An offset into the positions file of an IndexDB.  Typically, 3
   billion divided by sampling interval, requiring a maximum of 32
   bits or 4 bytes */
typedef UINT4 Positionsptr_T;

/* Typically 12 nt or 24 bits, requiring 3 bytes */
typedef UINT4 Storedoligomer_T;

typedef UINT4 Genomecomp_T;

#endif

