/* $Id: orderstat.h,v 1.3 2010/02/03 18:12:25 twu Exp $ */
#ifndef ORDERSTAT_INCLUDED
#define ORDERSTAT_INCLUDED

extern double
Orderstat_double_pct (double *set, int length, double pct);
extern double
Orderstat_double_pct_inplace (double *set, int length, double pct);
extern int
Orderstat_int_pct (int *set, int length, double pct);
extern long int
Orderstat_long_int_pct (long int *set, int length, double pct);
extern int
Orderstat_int_pct_inplace (int *set, int length, double pct);

#endif
