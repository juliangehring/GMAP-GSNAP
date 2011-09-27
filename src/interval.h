/* $Id: interval.h 40271 2011-05-28 02:29:18Z twu $ */
#ifndef INTERVAL_INCLUDED
#define INTERVAL_INCLUDED
#include "bool.h"

#define T Interval_T
typedef struct T *T;
struct T {
  unsigned int low;		/* low <= high */
  unsigned int high;
  int sign;
  int type;
};

extern T
Interval_new (unsigned int low, unsigned int high, int type);
extern T
Interval_copy (T old);
extern void
Interval_free (T *old);
extern void
Interval_print (T this);

extern unsigned int
Interval_low (T this);
extern unsigned int
Interval_high (T this);
extern void
Interval_store_length (T this, unsigned int length);
extern int
Interval_sign (T this);
extern unsigned int
Interval_length (T this);
extern int
Interval_type (T this);

extern unsigned int
Interval_array_low (struct T *intervals, int index);
extern unsigned int
Interval_array_high (struct T *intervals, int index);

extern bool
Interval_is_contained (unsigned int x, struct T *intervals, int index);
extern bool
Interval_overlap_p (unsigned int x, unsigned int y, struct T *intervals, int index);

extern void
Interval_qsort_by_sigma (int *table, int i, int j, struct T *intervals);
extern void
Interval_qsort_by_omega (int *table, int i, int j, struct T *intervals);

extern int
Interval_cmp (const void *a, const void *b);
extern int
Interval_cmp_low (const void *a, const void *b);
extern int
Interval_cmp_high (const void *a, const void *b);


struct Interval_windex_T {
  int index;
  T interval;
};

extern int
Interval_windex_cmp (const void *a, const void *b);

#undef T
#endif


