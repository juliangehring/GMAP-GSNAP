/* $Id: uint8list.h 148721 2014-09-24 00:45:45Z twu $ */
#ifndef UINT8LIST_INCLUDED
#define UINT8LIST_INCLUDED

typedef struct Uint8list_T *Uint8list_T;

#include "types.h"
#include "bool.h"

#define T Uint8list_T

extern T 
Uint8list_push (T list, UINT8 x);
extern T 
Uint8list_pop (T list, UINT8 *x);
extern UINT8 
Uint8list_head (T list);
extern T 
Uint8list_next (T list);
extern void 
Uint8list_head_set (T list, UINT8 x);
extern void 
Uint8list_free (T *list);
extern T 
Uint8list_reverse (T list);
extern int 
Uint8list_length (T list);
extern UINT8 *
Uint8list_to_array (int *n, T list);
extern UINT8 *
Uint8list_to_array_out (int *n, T list);
extern void
Uint8list_fill_array (UINT8 *array, T list);
extern void
Uint8list_fill_array_and_free (UINT8 *array, T *list);
extern T
Uint8list_from_array (UINT8 *array, int n);
extern T 
Uint8list_copy (T list);
extern T 
Uint8list_append (T list, T tail);
extern UINT8 
Uint8list_last_value (T this);
extern UINT8 
Uint8list_index (T this, int index);
extern bool
Uint8list_find (T this, UINT8 value);
extern char *
Uint8list_to_string (T this);
#undef T
#endif
