/* $Id: iit-write-univ.h 99737 2013-06-27 19:33:03Z twu $ */
#ifndef IIT_WRITE_UNIV_INCLUDED
#define IIT_WRITE_UNIV_INCLUDED
#include "bool.h"
#include "list.h"
#include "uintlist.h"
#include "table.h"
#include "iitdef.h"

#define T IIT_T
#ifndef IIT_TYPEDEF
#define IIT_TYPEDEF
typedef struct T *T;
#endif

extern void
IIT_write_univ (char *iitfile, List_T divlist, List_T typelist, Table_T intervaltable,
		Table_T labeltable, Table_T annottable,
		bool coord_values_8p, bool label_pointers_8p, bool annot_pointers_8p);

#undef T
#endif

