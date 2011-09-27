/* $Id: tableint.h,v 1.4 2005/02/15 01:57:34 twu Exp $ */
#ifndef TABLEINT_INCLUDED
#define TABLEINT_INCLUDED

#define T Tableint_T
typedef struct T *T;

extern T
Tableint_new (int hint,
	      int (*cmp)(const void *x, const void *y),
	      unsigned int hash(const void *key));
extern void 
Tableint_free (T *table);
extern int   
Tableint_length (T table);
extern int
Tableint_put (T table, const void *key, int value);
extern int
Tableint_get (T table, const void *key);
extern int
Tableint_remove (T table, const void *key);
extern void   
Tableint_map (T table,
	      void (*apply)(const void *key, int *value, void *cl),
	      void *cl);
extern void **
Tableint_keys (T table, void *end);
extern int *
Tableint_values (T table, int end);

#undef T
#endif
