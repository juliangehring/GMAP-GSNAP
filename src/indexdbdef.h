/* $Id: indexdbdef.h 44878 2011-08-13 23:18:46Z twu $ */
#ifndef INDEXDBDEF_INCLUDED
#define INDEXDBDEF_INCLUDED

#include "genomicpos.h"
#include "access.h"
#include "types.h"


/* Typically 12 nt or 24 bits, requiring 3 bytes */
typedef UINT4 Storedoligomer_T;

#define BADVAL (Genomicpos_T) -1

#define T Indexdb_T
struct T {
  int index1part;
  int index1interval;
  int offsetscomp_basesize;		/* e.g., 12 */
  int offsetscomp_blocksize;		/* e.g., 64 = 4^(15-12) */

  /* Access_T gammaptrs_access; -- Always ALLOCATED */ 
  int gammaptrs_fd;
  size_t gammaptrs_len;
  Positionsptr_T *gammaptrs;

  Access_T offsetscomp_access;
  int offsetscomp_fd;
  size_t offsetscomp_len;
  Positionsptr_T *offsetscomp;


  /* With gamma encoding, we can expand into offsets */
  Access_T offsets_access;
  int offsets_fd;
  size_t offsets_len;
  Positionsptr_T *offsets;
#ifdef HAVE_PTHREAD
#ifdef PMAP
  pthread_mutex_t offsets_read_mutex;
#endif
#endif


  Access_T positions_access;
  int positions_fd;
  size_t positions_len;
  Genomicpos_T *positions;
#ifdef HAVE_PTHREAD
  pthread_mutex_t positions_read_mutex;
#endif
};

#undef T
#endif

