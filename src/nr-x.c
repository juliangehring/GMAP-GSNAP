#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "nr-x.h"

#define ITMAX 100
#define MAXIT 100
#define EPS 3.0e-7
#define FPMIN 1.0e-30

double
gammln (double xx) {
  double x, y, tmp, ser;
  static double cof[6] = {76.18009172947146,
			  -86.50532032941677,
			  24.01409824083091,
			  -1.231739572450155,
			  0.1208650973866179e-2,
			  -0.5395239384953e-5};
  int j;
  
  y = x = xx;
  tmp = x + 5.5;
  tmp -= (x+0.5)*log(tmp);
  ser = 1.000000000190015;
  for (j = 0; j <= 5; ++j) {
    ser += cof[j]/++y;
  }
  return -tmp + log(2.5066282746310005*ser/x);
}


static double
betacf (double a, double b, double x) {
  int m, m2;
  double aa, c, d, del, h, qab, qam, qap;

  qab = a+b;
  qap = a + 1.0;
  qam = a - 1.0;
  c = 1.0;
  d = 1.0-qab*x/qap;
  if (fabs(d) < FPMIN) {
    d = FPMIN;
  }
  d = 1.0/d;
  h = d;
  for (m = 1; m <= MAXIT; ++m) {
    m2 = 2*m;
    aa = m*(b-m)*x/((qam+m2)*(a+m2));
    d = 1.0 + aa*d;
    if (fabs(d) < FPMIN) {
      d = FPMIN;
    }
    c = 1.0 + aa/c;
    if (fabs(c) < FPMIN) {
      c = FPMIN;
    }
    d = 1.0/d;
    h *= d*c;
    aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
    d = 1.0 + aa*d;
    if (fabs(d) < FPMIN) {
      d = FPMIN;
    }
    c = 1.0 + aa/c;
    if (fabs(c) < FPMIN) {
      c = FPMIN;
    }
    d = 1.0/d;
    del = d*c;
    h *= del;
    if (fabs(del-1.0) < EPS) {
      break;
    }
  }
  if (m > MAXIT) {
    printf("betacf: a (%f) or b (%f) too big, or MAXIT too small\n",
	   a,b);
    abort();
  }
  return h;
}


static double
betai (double a, double b, double x) {
  double bt;

  if (x < 0.0 || x > 1.0) {
    printf("betai called with %f %f %f\n",a,b,x);
    abort();
  }
  if (x == 0.0 || x == 1.0) {
    bt = 0.0;
  } else {
    bt = exp(gammln(a+b)-gammln(a)-gammln(b)+a*log(x)+b*log(1.0-x));
  }
  if (x < (a+1.0)/(a+b+2.0)) {
    return bt*betacf(a,b,x)/a;
  } else {
    return 1.0-bt*betacf(b,a,1.0-x)/b;
  }
}

static void
gser (double *gamser, double a, double x, double *gln) {
  int n;
  double sum, del, ap;

  *gln = gammln(a);
  if (x <= 0.0) {
    if (x < 0.0) {
      printf("x less than 0 in routine gser\n");
    }
    *gamser=0.0;
    return;
  } else {
    ap = a;
    del = sum = 1.0/a;
    for (n = 1; n <= ITMAX; ++n) {
      ++ap;
      del *= x/ap;
      sum += del;
      if (fabs(del) < fabs(sum)*EPS) {
	*gamser = sum * exp(-x + a*log(x) - (*gln));
	return;
      }
    }
    printf("a too large, ITMAX too small in routine gser\n");
    return;
  }
}

static void
gcf (double *gammcf, double a, double x, double *gln) {
  int i;
  double an, b, c, d, del, h;

  *gln = gammln(a);
  b = x + 1.0 - a;
  c = 1.0/FPMIN;
  d = 1.0/b;
  h = d;
  for (i = 1; i <= ITMAX; ++i) {
    an = -i * (i-a);
    b += 2.0;
    d = an*d + b;
    if (fabs(d) < FPMIN) {
      d = FPMIN;
    }
    c = b + an/c;
    if (fabs(c) < FPMIN) {
      c = FPMIN;
    }
    d = 1.0/d;
    del = d*c;
    h *= del;
    if (fabs(del-1.0) < EPS) {
      break;
    }
  }
  if (i > ITMAX) {
    printf("gcf: a (%f) too large, ITMAX too small\n",a);
    printf("gcf: x = %f\n",x);
    abort();
  }
  *gammcf = exp(-x + a*log(x) - (*gln)) * h;
}

static double
gammp (double a, double x) {
  double gamser, gammcf, gln;

  if (x < 0.0 || a <= 0.0) {
    printf("Invalid arguments in routine gammp\n");
    exit(9);
    return 0.0;
  } else if (x < a+1.0) {
    gser(&gamser,a,x,&gln);
    return gamser;
  } else {
    gcf(&gammcf,a,x,&gln);
    return 1.0-gammcf;
  }
}

static double
gammq (double a, double x) {
  double gamser, gammcf, gln;

  if (x < 0.0 || a <= 0.0) {
    printf("Invalid arguments in routine gammq: x = %f, a = %f\n",x,a);
    abort();
    return 0.0;
  } else if (x < a+1.0) {
    gser(&gamser,a,x,&gln);
    return 1.0-gamser;
  } else {
    gcf(&gammcf,a,x,&gln);
    return gammcf;
  }
}


/************************************************************************
 *   Binomial
 ************************************************************************/

/* Probability of seeing k or more */
double
NR_pbinomc (double k, double size, double prob) {
  if (k == 0.0) {
    return 1.0;
  } else {
    return betai(k,size-k+1.0,prob);
  }
}

/* Same as Splus pbinom */
double
NR_pbinom (double k, double size, double prob) {
  if (k == size) {
    return 0.0;
  } else {
    return betai(size-k,k+1.0,1.0-prob);
  }
}

/* Same as Splus dbinom */
double
NR_dbinom (double k, double size, double prob) {
  double result;

  if (k == 0.0) {
    result = NR_pbinom(k,size,prob);
  } else {
    result = NR_pbinom(k,size,prob) - NR_pbinom(k-1,size,prob);
  }
  if (result == 0.0) {
    result = FPMIN;
  }
  return result;
}


/************************************************************************
 *  Chisq
 ************************************************************************/

double
NR_pchisqc (double chisq, int df) {
  if (chisq < 0.0) {
    return 1.0;
  }
  return gammq((double) df/2.0, chisq/2.0);
}

/* Does the same thing as Splus pchisq */
double
NR_pchisq (double chisq, int df) {
  if (chisq < 0.0) {
    return 0.0;
  }
  return gammp((double) df/2.0, chisq/2.0);
}


/************************************************************************
 *   Normal
 ************************************************************************/

/* This is 1.0 - pnorm */
double
NR_pnormc (double x) {
  if (x < -20.0) {
    return 1.0;
  } else if (x < 0.0) {
    return 0.5 + 0.5*gammp(0.5,x*x/2.0);
  } else if (x < 20.0) {
    return 0.5*gammq(0.5,x*x/2.0);
  } else {
    return 0.0;
  }
}

/* Does the same thing as Splus pnorm */
double
NR_pnorm (double x) {
  if (x < -20.0) {
    return 0.0;
  } else if (x < 0.0) {
    return 0.5*gammq(0.5,x*x/2.0);
  } else if (x < 20.0) {
    return 0.5 + 0.5*gammp(0.5,x*x/2.0);
  } else {
    return 1.0;
  }
}


/************************************************************************
 *   F
 ************************************************************************/

double
NR_pfc (double F, double df1, double df2) {
  return betai(df2/2.0,df1/2.0,df2/(df2+df1*F));
}


/* Same as Splus pf */
double
NR_pf (double F, double df1, double df2) {
  return betai(df1/2.0,df2/2.0,1.0-(df2/(df2+df1*F)));
}


/************************************************************************
 *   Student's t
 ************************************************************************/

double
NR_ptc (double t, double df) {
  double val;

  if (df > 4e5) {
    val = 1.0/(4.0*df);
    return NR_pnormc(t*(1.0-val)/sqrt(1.0+t*t*2.0*val));
  } else {
    return 0.5*betai(df/2.0,0.5,df/(df+t*t));
  }
}


/* Same as Splus pt */
double
NR_pt (double t, double df) {
  double val;

  if (df > 4e5) {
    val = 1.0/(4.0*df);
    return NR_pnorm(t*(1.0-val)/sqrt(1.0+t*t*2.0*val));
  } else {
    return 0.5 + 0.5*betai(0.5,df/2.0,1.0-(df/(df+t*t)));
  }
}

