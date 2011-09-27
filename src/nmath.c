static char rcsid[] = "$Id: nmath.c,v 1.1 2005/05/04 18:02:50 twu Exp $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "nmath.h"
#include <math.h>

#define SIXTEN 16
#define M_SQRT_32 5.656854249492380195206754896838 /* sqrt(32) */
/* #define DBL_MIN 2.247e-308 */
#define DBL_MIN 1e-30


double
Nmath_pnormc (double x) {
  double ccum, xden, xnum, temp, del, eps = 1e-16, xsq, y;
  int i, lower, upper;
  
  const double a[5] = {
    2.2352520354606839287,
    161.02823106855587881,
    1067.6894854603709582,
    18154.981253343561249,
    0.065682337918207449113
  };
  const double b[4] = {
    47.20258190468824187,
    976.09855173777669322,
    10260.932208618978205,
    45507.789335026729956
  };
  const double c[9] = {
    0.39894151208813466764,
    8.8831497943883759412,
    93.506656132177855979,
    597.27027639480026226,
    2494.5375852903726711,
    6848.1904505362823326,
    11602.651437647350124,
    9842.7148383839780218,
    1.0765576773720192317e-8
  };
  const double d[8] = {
    22.266688044328115691,
    235.38790178262499861,
    1519.377599407554805,
    6485.558298266760755,
    18615.571640885098091,
    34900.952721145977266,
    38912.003286093271411,
    19685.429676859990727
  };
  const double p[6] = {
    0.21589853405795699,
    0.1274011611602473639,
    0.022235277870649807,
    0.001421619193227893466,
    2.9112874951168792e-5,
    0.02307344176494017303
  };
  const double q[5] = {
    1.28426009614491121,
    0.468238212480865118,
    0.0659881378689285515,
    0.00378239633202758244,
    7.29751555083966205e-5
  };
  
  y = fabs(x);
  if (y <= 0.67448975) { /* qnorm(3/4) = .6744.... -- earlier had 0.66291 */
    if (y > eps) {
      xsq = x * x;
      xnum = a[4] * xsq;
      xden = xsq;
      for (i = 0; i < 3; ++i) {
	xnum = (xnum + a[i]) * xsq;
	xden = (xden + b[i]) * xsq;
      }
    } else {
      xnum = xden = 0.0;
    }
    
    temp = x * (xnum + a[3]) / (xden + b[3]);
    ccum = 0.5 - temp;
    
  } else if (y <= 8.0) {
    /* Threshold is sqrt(32) in R code */
    /* Evaluate pnorm for 0.674.. = qnorm(3/4) < |x| <= sqrt(32) ~= 5.657 */
    xnum = c[8] * y;
    xden = y;
    for (i = 0; i < 7; i++) {
      xnum = (xnum + c[i]) * y;
      xden = (xden + d[i]) * y;
    }
    temp = (xnum + c[7]) / (xden + d[7]);
    
    xsq = trunc(y * SIXTEN) / SIXTEN;
    del = (y - xsq) * (y + xsq);
    
    if (x > 0) {
      ccum = exp(-xsq * xsq * 0.5) * exp(-del * 0.5) * temp;
    } else {
      ccum = 1.0 - exp(-xsq * xsq * 0.5) * exp(-del * 0.5) * temp;
    }
    
  } else if (x > 0) {
    ccum = 0.0;
  } else {
    ccum = 1.0;
  }
  
  return ccum;
}


