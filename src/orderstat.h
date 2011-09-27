/* $Id: orderstat.h,v 1.2 2007/08/15 22:30:25 twu Exp $ */
#ifndef ORDERSTAT_INCLUDED
#define ORDERSTAT_INCLUDED

extern double
Orderstat_double_pct (double *set, int length, double pct);
extern double
Orderstat_double_pct_inplace (double *set, int length, double pct);
extern int
Orderstat_int_pct (int *set, int length, double pct);
extern int
Orderstat_int_pct_inplace (int *set, int length, double pct);

#endif
