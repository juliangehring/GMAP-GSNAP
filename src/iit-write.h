/* $Id: iit-write.h,v 1.18 2008-03-27 20:43:11 twu Exp $ */
#ifndef IIT_WRITE_INCLUDED
#define IIT_WRITE_INCLUDED
#include "list.h"
#include "uintlist.h"
#include "table.h"
#include "iitdef.h"		/* For Sorttype_T */

#define T IIT_T
#ifndef IIT_TYPEDEF
#define IIT_TYPEDEF
typedef struct T *T;
#endif

extern void
IIT_output_direct (char *iitfile, T this, int version);
extern void
IIT_write (char *iitfile, List_T divlist, List_T typelist, List_T fieldlist, Table_T intervaltable,
	   Table_T labeltable, Table_T annottable, Sorttype_T divsort, int version);
extern T
IIT_new (List_T intervallist);
extern void
IIT_backfill_sequence (T this, int index, int offset, char *Buffer);

#undef T
#endif

