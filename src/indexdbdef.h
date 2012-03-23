/* $Id: indexdbdef.h 56653 2012-01-27 21:20:20Z twu $ */
#ifndef INDEXDBDEF_INCLUDED
#define INDEXDBDEF_INCLUDED

#include "genomicpos.h"
#include "access.h"
#include "types.h"

#ifdef PMAP
#include "alphabet.h"
#endif


#define BADVAL (Genomicpos_T) -1

#define T Indexdb_T
struct T {
#ifdef PMAP
  Alphabet_T alphabet;
  int alphabet_size;
#endif

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

