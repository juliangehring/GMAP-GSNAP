/* $Id: table.h,v 1.12 2008/04/15 20:01:19 twu Exp $ */
#ifndef TABLE_INCLUDED
#define TABLE_INCLUDED


#define T Table_T
typedef struct T *T;

extern int
Table_string_compare (const void *x, const void *y);
extern unsigned int
Table_string_hash (const void *x);

extern T
Table_new (int hint,
	   int (*cmp)(const void *x, const void *y),
	   unsigned int hash(const void *key));
extern void 
Table_free (T *table);
extern int   
Table_length (T table);
extern void *
Table_put (T table, const void *key,
	   void *value);
extern void *
Table_get (T table, const void *key);
extern void *
Table_remove (T table, const void *key);
extern void   
Table_map (T table,
	   void (*apply)(const void *key, void **value, void *cl),
	   void *cl);
extern void **
Table_keys (T table, void *end);
extern void **
Table_keys_by_timeindex (T table, void *end);
extern void **
Table_values (T table, void *end);

#undef T
#endif
