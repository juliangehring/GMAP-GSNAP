/* $Id: iit-write.h,v 1.12 2005/02/16 00:44:00 twu Exp $ */
#ifndef IIT_WRITE_INCLUDED
#define IIT_WRITE_INCLUDED
#include "list.h"
#include "uintlist.h"

#define T IIT_T
#ifndef IIT_TYPEDEF
#define IIT_TYPEDEF
typedef struct T *T;
#endif

extern void
IIT_free (T *old);
extern void
IIT_write (char *filename, List_T intervallist, List_T typelist, List_T labellist, List_T annotlist,
	   Uintlist_T annot_strlen_list);
extern void
IIT_backfill_sequence (T this, int index, int offset, char *Buffer);

#undef T
#endif

