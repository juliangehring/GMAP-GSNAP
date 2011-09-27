/* $Id: iitdef.h,v 1.8 2005/02/15 01:55:34 twu Exp $ */
#ifndef IITDEF_INCLUDED
#define IITDEF_INCLUDED
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
  unsigned int flength;		/* mmap length */

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
