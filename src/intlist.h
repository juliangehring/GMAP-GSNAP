/* $Id: intlist.h,v 1.8 2005/06/21 18:37:57 twu Exp $ */
#ifndef INTLIST_INCLUDED
#define INTLIST_INCLUDED

#define T Intlist_T
typedef struct T *T;

extern T 
Intlist_push (T list, int x);
extern T 
Intlist_pop (T list, int *x);
extern void
Intlist_delete (T prev, T this);
extern int 
Intlist_head (T list);
extern T 
Intlist_next (T list);
extern void 
Intlist_head_set (T list, int x);
extern void 
Intlist_free (T *list);
extern T 
Intlist_reverse (T list);
extern int 
Intlist_length (T list);
extern int *
Intlist_to_array (int *n, T list);
extern char *
Intlist_to_char_array (int *n, T list);
extern T 
Intlist_copy (T list);
extern T 
Intlist_append (T list, T tail);
extern int 
Intlist_last_value (T this);
extern int 
Intlist_index (T this, int index);

#undef T
#endif
