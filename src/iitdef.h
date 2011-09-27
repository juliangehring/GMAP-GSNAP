/* $Id: iitdef.h,v 1.11 2005/10/19 03:48:58 twu Exp $ */
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

  Access_T access;		/* access type */

  int fd;			/* file descriptor */
  char *finfo;			/* result of mmap command */
  size_t flength;		/* mmap length (mmap uses size_t, not off_t) */
  off_t offset;			/* used only for fileio */
#ifdef HAVE_PTHREAD
  pthread_mutex_t read_mutex;
#endif

  int nintervals;
  int ntypes;
  int nnodes;
  int *sigmas;
  int *omegas;
  struct FNode_T *nodes;
  struct Interval_T *intervals;
  unsigned int *typepointers;
  char *typestrings;
  int *labelorder;
  unsigned int *labelpointers;
  char *labels;
  unsigned int *annotpointers;
  char *annotations;
};


#undef T
#endif
