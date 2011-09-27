#ifndef NR_X_INCLUDED
#define NR_X_INCLUDED

extern double
NR_pbinomc (double k, double size, double prob);
extern double
NR_pbinom (double k, double size, double prob);
extern double
NR_dbinom (double k, double size, double prob);
extern double
NR_pchisqc (double chisq, int df);
extern double
NR_pchisq (double chisq, int df);
extern double
NR_pnormc (double x);
extern double
NR_pnorm (double x);
extern double
NR_pfc (double F, double df1, double df2);
extern double
NR_pf (double F, double df1, double df2);
extern double
NR_ptc (double t, double df);
extern double
NR_pt (double t, double df);

#endif

