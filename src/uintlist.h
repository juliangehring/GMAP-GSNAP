/* $Id: uintlist.h,v 1.4 2005/02/07 23:56:57 twu Exp $ */
#ifndef UINTLIST_INCLUDED
#define UINTLIST_INCLUDED

#define T Uintlist_T
typedef struct T *T;

extern T 
Uintlist_push (T list, unsigned int x);
extern T 
Uintlist_pop (T list, unsigned int *x);
extern unsigned int 
Uintlist_head (T list);
extern T 
Uintlist_next (T list);
extern void 
Uintlist_head_set (T list, unsigned int x);
extern void 
Uintlist_free (T *list);
extern T 
Uintlist_reverse (T list);
extern int 
Uintlist_length (T list);
extern unsigned int *
Uintlist_to_array (int *n, T list);
extern T 
Uintlist_copy (T list);
extern T 
Uintlist_append (T list, T tail);
extern unsigned int 
Uintlist_last_value (T this);
extern unsigned int 
Uintlist_index (T this, int index);

#undef T
#endif
