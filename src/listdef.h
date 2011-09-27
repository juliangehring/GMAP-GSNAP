/* $Id: listdef.h,v 1.3 2005-02-15 01:55:58 twu Exp $ */
#ifndef LISTDEF_INCLUDED
#define LISTDEF_INCLUDED

#define T List_T
struct T {
  void *first;
  struct T *rest;
};

#undef T
#endif

