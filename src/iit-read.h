/* $Id: iit-read.h,v 1.27 2005/03/11 17:55:20 twu Exp $ */
#ifndef IIT_READ_INCLUDED
#define IIT_READ_INCLUDED
#include <stdio.h>
#include "bool.h"
#include "intlist.h"
#include "uintlist.h"
#include "list.h"
#include "interval.h"

#define T IIT_T
#ifndef IIT_TYPEDEF
#define IIT_TYPEDEF
typedef struct T *T;
#endif

extern char *
IIT_name (T this);
extern int
IIT_nintervals (T this);
extern int
IIT_ntypes (T this);
extern unsigned int
IIT_totallength (T this);
extern Interval_T
IIT_interval (T this, int index);
extern char *
IIT_typestring (T this, int type);
extern int
IIT_typeint (T this, char *typestring);
extern char *
IIT_label (T this, int index);
extern char *
IIT_annotation (T this, int index);
extern unsigned int
IIT_annotation_strlen (T this, int index);

extern void
IIT_debug (char *filename);
extern void
IIT_dump_typestrings (FILE *fp, T this);
extern void
IIT_dump (T this);
extern void
IIT_dump_formatted (T this, bool directionalp);

extern void
IIT_free_mmapped (T *old);
extern T
IIT_read (char *filename, char *name, bool readonlyp);

extern Intlist_T
IIT_find (T this, char *label);
extern int
IIT_find_linear (T this, char *label);
extern int
IIT_find_one (T this, char *label);

extern Intlist_T
IIT_get (T this, unsigned int x, unsigned int y);
extern int
IIT_get_one (T this, unsigned int x, unsigned int y);
extern Intlist_T
IIT_get_typed (T this, unsigned int x, unsigned int y, int type);
extern int
IIT_get_exact (T this, unsigned int x, unsigned int y, int type);
extern void
IIT_print (T this, Intlist_T matches, bool map_bothstrands_p);

extern List_T
IIT_intervallist_typed (List_T *labellist, Uintlist_T *seglength_list, T this);
extern List_T
IIT_typelist (T this);


#undef T
#endif
