/* $Id: iitdef.h,v 1.10 2005/07/08 19:16:53 twu Exp $ */
#ifndef IITDEF_INCLUDED
#define IITDEF_INCLUDED
#include <sys/types.h>
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

  int fd;			/* file descriptor */
  char *finfo;			/* result of mmap command */
  off_t flength;		/* mmap length */

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
