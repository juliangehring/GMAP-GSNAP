/* $Id: intlistdef.h 27450 2010-08-05 19:02:48Z twu $ */
#ifndef INTLISTDEF_INCLUDED
#define INTLISTDEF_INCLUDED

#define T Intlist_T
struct T {
  int first;
  struct T *rest;
};

#undef T
#endif
