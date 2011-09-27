/* $Id: list.h,v 1.8 2005/02/07 23:56:56 twu Exp $ */
#ifndef LIST_INCLUDED
#define LIST_INCLUDED

#define T List_T
typedef struct T *T;

extern T List_push (T list, void *x);
extern T List_pop (T list, void **x);
extern void *List_head (T list);
extern T List_next (T list);
extern void List_head_set (T list, void *x);
extern void List_free (T *list);
extern T List_reverse (T list);
extern int List_length (T list);
extern void **List_to_array (T list, void *end);
extern T List_copy (T list);
extern T List_append (T list, T tail);
extern void *
List_last_value (T this);
extern void *
List_index (T this, int index);

#undef T
#endif
