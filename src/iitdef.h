/* $Id: iitdef.h,v 1.20 2009/08/29 00:36:12 twu Exp $ */
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


#ifdef HAVE_64_BIT
#define IIT_LATEST_VERSION 4
#define IIT_8BYTE_VERSION 4
#define IIT_4BYTE_VERSION 3
#else
#define IIT_LATEST_VERSION 3
#endif

/* version 1 starts with nintervals */
/* version 2 starts with 0, then version number.  Also adds sign to each interval.  */
/* version 3 allows for multiple divs */
/* version 4 has label and annot pointers being 8-byte long unsigned ints */


typedef enum {NO_SORT, ALPHA_SORT, CHROM_SORT} Sorttype_T;

typedef struct FNode_T *FNode_T;
struct FNode_T {
  unsigned int value;
  int a;
  int b;
  int leftindex;
  int rightindex;
};

#define T IIT_T
struct T {
  char *name;			/* Name of IIT (optional) */
  int version;			
  int fd;
  Access_T access;		/* access type */

#ifdef HAVE_PTHREAD
  pthread_mutex_t read_mutex;
#endif

  int ntypes;			/* Always >= 1 */
  int nfields;			/* Can be zero */

  int divsort;			/* Really Sorttype_T */
  int ndivs;
  unsigned int *divpointers;
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

  unsigned int *typepointers;
  char *typestrings;

  unsigned int *fieldpointers;
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
  unsigned int *labelpointers;
#ifdef HAVE_64_BIT
  UINT8 *labelpointers8;
#endif
  char *labels;

  unsigned int *annotpointers;
#ifdef HAVE_64_BIT
  UINT8 *annotpointers8;
#endif
  char *annotations;
};


#undef T
#endif
