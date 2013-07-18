/* $Id: iitdef.h 99737 2013-06-27 19:33:03Z twu $ */
#ifndef IITDEF_INCLUDED
#define IITDEF_INCLUDED
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef HAVE_SYS_TYPES_H
# include <sys/types.h>		/* Needed to define pthread_t on Solaris */
#endif
#ifdef HAVE_PTHREAD
#include <pthread.h>
#endif
#include "bool.h"
#include "access.h"
#include "interval.h"
#include "types.h"


#define IIT_LATEST_VERSION 5


/* version 1 starts with nintervals (now handled separately as a Univ_IIT_T) */
/* version 2 starts with 0, then version number.  Also adds sign to each interval.  */
/* version 3 allows for multiple divs */
/* version 4 has label and annot pointers being 8-byte long unsigned ints */

/* version 5 has two extra fields indicating whether label pointers are
   4- or 8-bytes and whether annot pointers are 4- or 8-bytes.  Also
   stores rest of header line with annotation, so NULL => print '\n',
   otherwise print annotation. */


typedef enum {NO_SORT, ALPHA_SORT, NUMERIC_ALPHA_SORT, CHROM_SORT} Sorttype_T;

typedef struct FNode_T *FNode_T;
struct FNode_T {
  Chrpos_T value;
  int a;
  int b;
  int leftindex;
  int rightindex;
};

#define T IIT_T
struct T {
  char *name;			/* Name of IIT (optional) */
  int version;			
  bool label_pointers_8p;
  bool annot_pointers_8p;

  int fd;
  Access_T access;		/* access type */

#ifdef HAVE_PTHREAD
  pthread_mutex_t read_mutex;
#endif

  int ntypes;			/* Always >= 1 */
  int nfields;			/* Can be zero */

  int divsort;			/* Really Sorttype_T */
  int ndivs;
  UINT4 *divpointers;
  char *divstrings;

  int total_nintervals;
  int *nintervals;		/* Per div */
  int *cum_nintervals;
  int *nnodes;			/* Per div */
  int *cum_nnodes;

  int **alphas;			/* Strict ordering of Interval_low */
  int **betas;			/* Strict ordering of Interval_high */
  int **sigmas;			/* Ordering for IIT */
  int **omegas;			/* Ordering for IIT */

  struct FNode_T **nodes;	/* Per div */
  struct Interval_T **intervals; /* Per div */

  UINT4 *typepointers;
  char *typestrings;

  UINT4 *fieldpointers;
  char *fieldstrings;

  off_t labelorder_offset;
  size_t labelorder_length; /* mmap length (mmap uses size_t, not off_t) */
  char *labelorder_mmap;

  off_t labelpointers_offset;
  size_t labelpointers_length; /* mmap length (mmap uses size_t, not off_t) */
  char *labelpointers_mmap;

  off_t label_offset;
  size_t label_length;		/* mmap length (mmap uses size_t, not off_t) */
  char *label_mmap;

  off_t annotpointers_offset;
  size_t annotpointers_length; /* mmap length (mmap uses size_t, not off_t) */
  char *annotpointers_mmap;

  off_t annot_offset;
  size_t annot_length;		/* mmap length (mmap uses size_t, not off_t) */
  char *annot_mmap;

  int *labelorder;
  UINT4 *labelpointers;
#ifdef HAVE_64_BIT
  UINT8 *labelpointers8;
#endif
  char *labels;

  UINT4 *annotpointers;
#ifdef HAVE_64_BIT
  UINT8 *annotpointers8;
#endif
  char *annotations;

  void **datapointers;
};


#undef T
#endif
