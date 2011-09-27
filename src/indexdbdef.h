/* $Id: indexdbdef.h 36299 2011-03-09 05:26:28Z twu $ */
#ifndef INDEXDBDEF_INCLUDED
#define INDEXDBDEF_INCLUDED

#include "indexdb.h"		/* For Storedoligomer_T */
#include "genomicpos.h"
#include "access.h"
#include "types.h"

/* An offset into the positions file of an IndexDB.  Typically, 3
   billion divided by sampling interval, requiring a maximum of 32
   bits or 4 bytes */
typedef UINT4 Positionsptr_T;

#define BADVAL (Genomicpos_T) -1

#define T Indexdb_T
struct T {
  int index1interval;

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

