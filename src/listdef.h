/* $Id: listdef.h 27450 2010-08-05 19:02:48Z twu $ */
#ifndef LISTDEF_INCLUDED
#define LISTDEF_INCLUDED

#define T List_T
struct T {
  void *first;
  struct T *rest;
};

#undef T
#endif

