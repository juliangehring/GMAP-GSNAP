/* $Id: uinttable.h 40271 2011-05-28 02:29:18Z twu $ */
#ifndef UINTTABLE_INCLUDED
#define UINTTABLE_INCLUDED
#include "bool.h"

#define T Uinttable_T
typedef struct T *T;

extern T
Uinttable_new (int hint);
extern void 
Uinttable_free (T *table);
extern int   
Uinttable_length (T table);
extern void *
Uinttable_put (T table, const unsigned int key, void *value);
extern void *
Uinttable_get (T table, const unsigned int key);
extern void *
Uinttable_remove (T table, const unsigned int key);
extern void   
Uinttable_map (T table,
	       void (*apply)(const unsigned int key, void **value, void *cl),
	       void *cl);
extern unsigned int *
Uinttable_keys (T table, bool sortp);
extern unsigned int *
Uinttable_keys_by_timeindex (T table);
extern void **
Uinttable_values (T table);

#undef T
#endif
