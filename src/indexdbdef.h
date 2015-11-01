/* $Id: indexdbdef.h 132144 2014-04-02 16:02:28Z twu $ */
#ifndef INDEXDBDEF_INCLUDED
#define INDEXDBDEF_INCLUDED

#include "genomicpos.h"
#include "access.h"
#include "types.h"

#ifdef PMAP
#include "alphabet.h"
#endif


#define BADVAL (Univcoord_T) -1

/* Compression types */
#define NO_COMPRESSION 0
#define BITPACK64_COMPRESSION 1


#define T Indexdb_T
struct T {
#ifdef PMAP
  Alphabet_T alphabet;
  int alphabet_size;
#endif

  int compression_type;
  Width_T index1part;
  Width_T index1interval;
  Blocksize_T blocksize;	/* e.g., 64 = 4^(15-12) */

#ifdef LARGE_GENOMES
  UINT4 *offsetspages;
#endif

  int offsetsmeta_fd;
  size_t offsetsmeta_len;
  UINT4 *offsetsmeta;

  Access_T offsetsstrm_access;
  int offsetsstrm_fd;
  size_t offsetsstrm_len;
  UINT4 *offsetsstrm;

  Access_T positions_access;
#ifdef LARGE_GENOMES
  int positions_high_fd;
  size_t positions_high_len;
  int positions_low_fd;
  size_t positions_low_len;
  unsigned char *positions_high;
  UINT4 *positions_low;
#else
  int positions_fd;
  size_t positions_len;
  UINT4 *positions;		/* For small genomes, same as Univcoord_T */
#endif

#ifdef HAVE_PTHREAD
  pthread_mutex_t positions_read_mutex;
#endif
};

#undef T
#endif

