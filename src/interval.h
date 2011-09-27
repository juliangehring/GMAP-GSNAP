/* $Id: interval.h,v 1.14 2005/05/09 22:33:33 twu Exp $ */
#ifndef INTERVAL_INCLUDED
#define INTERVAL_INCLUDED
#include "bool.h"

#define T Interval_T
typedef struct T *T;
struct T {
  unsigned int low;
  unsigned int high;
  int type;
};


extern T
Interval_new (unsigned int low, unsigned int high, int type);
extern T
Interval_copy (T old);
extern void
Interval_free (T *old);

extern unsigned int
Interval_low (T this);
extern unsigned int
Interval_high (T this);
extern unsigned int
Interval_length (T this);

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

#undef T
#endif


