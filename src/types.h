/* $Id: types.h 102177 2013-07-20 00:51:23Z twu $ */
#ifndef TYPES_INCLUDED
#define TYPES_INCLUDED
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/* Number of bits, such as index1part or basesize.  Need to allow for negative values. */
typedef int Width_T;

/* Number of entries, such as offsetscomp_blocksize */
typedef unsigned int Blocksize_T;


/* A 2-byte word */
typedef unsigned short UINT2;

/* A 4-byte word */
typedef unsigned int UINT4;
typedef int INT4;


/* Compressed representation of genome (high, low, flags).  Always
   UINT4.  Can think of as a genome block unit.  */
typedef UINT4 Genomecomp_T;


/* An 8-byte word */
/* Oligospace_T needs to hold 1 more than maximum Storedoligomer_T.
   If 8-byte words are not available, then maximum k-mer is 15 */
/* Prefer to use unsigned long long, whic should be 8 bytes on all systems */
#if (SIZEOF_UNSIGNED_LONG_LONG == 8)
#define HAVE_64_BIT 1
#define MAXIMUM_KMER 16
typedef unsigned long long UINT8;
typedef unsigned long long Oligospace_T;
#elif (SIZEOF_UNSIGNED_LONG == 8)
#define HAVE_64_BIT 1
#define MAXIMUM_KMER 16
typedef unsigned long UINT8;
typedef unsigned long Oligospace_T;
#else
#define MAXIMUM_KMER 15
#define OLIGOSPACE_NOT_LONG
typedef unsigned int Oligospace_T;
#endif

/* Pointer into compressed offsets file.  Can be UINT4 as long as file
   has fewer than 4 billion words, which should hold for k <= 16. */
typedef UINT4 Gammaptr_T;

/* Contents of compressed offsets file.  Storing as UINT4, even for
   large genomes, to reduce zero-padding of bitstreams.  For large
   genomes, need to store 64-bit Positionsptr_T quantity in 2 UINT4
   words. */
typedef UINT4 Offsetscomp_T;

/* Holds a k-mer.  Can be UINT4 as long as k <= 16. */
/* Some procedures use Shortoligomer_T, which should be the same */
typedef UINT4 Storedoligomer_T;

/* An offset into the positions file of an IndexDB.  For small genomes
   < 2^32 bp such as human, need 3 billion divided by sampling
   interval (default 3), requiring a maximum of 32 bits or 4 bytes */
typedef UINT4 Positionsptr_T;


/* For definition of Univcoord_T and Chrpos_T, see genomicpos.h */

/* For intervals and IIT files */
#ifdef HAVE_64_BIT

#ifdef UTILITYP
typedef UINT8 Univcoord_T;
#elif defined LARGE_GENOMES
typedef UINT8 Univcoord_T;
#else
typedef UINT4 Univcoord_T;
#endif

#else
typedef UINT4 Univcoord_T;
#endif

/* For splicetrie */
typedef UINT4 Trieoffset_T;
typedef UINT4 Triecontent_T;

/* For suffix array */
typedef UINT4 Sarrayptr_T;

#endif

