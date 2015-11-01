/* $Id: doublelist.h 145990 2014-08-25 21:47:32Z twu $ */
#ifndef DOUBLELIST_INCLUDED
#define DOUBLELIST_INCLUDED

#define T Doublelist_T
typedef struct T *T;

extern T Doublelist_push (T list, double index);
extern T Doublelist_pop (T list, double *index);
extern double Doublelist_head (T list);
extern T Doublelist_next (T list);
extern void Doublelist_free (T *list);
extern T Doublelist_reverse (T list);
extern int Doublelist_length (T list);
extern double *
Doublelist_to_array (int *n, T list);
extern void
Doublelist_fill_array (double *array, T list);
extern T Doublelist_from_string (char *string);
extern double
Doublelist_max (T this);
extern double
Doublelist_min (T this);

#undef T
#endif
