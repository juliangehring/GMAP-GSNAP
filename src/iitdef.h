/* $Id: iitdef.h,v 1.14 2007/07/16 17:21:55 twu Exp $ */
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
#include "access.h"
#include "interval.h"


#define IIT_LATEST_VERSION 2
/* version 1 starts with nintervals */
/* version 2 starts with 0, then version number.  Also adds sign to each interval.  */


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

  Access_T access;		/* access type */

  int fd;			/* file descriptor */
  char *finfo;			/* result of mmap command */
  size_t flength;		/* mmap length (mmap uses size_t, not off_t) */
  off_t offset;			/* used only for fileio */
#ifdef HAVE_PTHREAD
  pthread_mutex_t read_mutex;
#endif

  int nintervals;
  int ntypes;			/* Always >= 1 */
  int nfields;			/* Can be zero */
  int nnodes;
  int *alphas;			/* Strict ordering of Interval_low */
  int *betas;			/* Strict ordering of Interval_high */
  int *sigmas;			/* Ordering for IIT */
  int *omegas;			/* Ordering for IIT */
  struct FNode_T *nodes;
  struct Interval_T *intervals;

  unsigned int *typepointers;
  char *typestrings;

  unsigned int *fieldpointers;
  char *fieldstrings;

  int *labelorder;
  unsigned int *labelpointers;
  char *labels;
  unsigned int *annotpointers;
  char *annotations;
};


#undef T
#endif
