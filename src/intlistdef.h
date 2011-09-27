/* $Id: intlistdef.h,v 1.1 2007-02-05 07:13:21 twu Exp $ */
#ifndef INTLISTDEF_INCLUDED
#define INTLISTDEF_INCLUDED

#define T Intlist_T
struct T {
  int first;
  struct T *rest;
};

#undef T
#endif
