#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "bool.h"
#include "mem.h"
#include "color.h"

#define HALFPI 1.50707963
#define PI 3.1415926
#define TWOPI 6.2831852


static double
g (double x) {
  if (x >= 6.0/29.0) {
    return x*x*x;
  } else {
    return 108.0/841.0 * (x - 4.0/29.0);
  }
}


/* Taken from the Adobe Postscript Reference Book */
/* D6500 white point is (0.9505 1.0 1.0890) */
void
Lab_to_XYZ (double *x, double *y, double *z, double L, double a, double b) {
  double yterm;

  *x = 0.9505 * g((L + 16.0)/116.0 + a/500.0);
  *y = 1.0000 * g((L + 16.0)/116.0);
  *z = 1.0890 * g((L + 16.0)/116.0 - b/200.0);
  return;
}

static void
XYZ_to_Lab (double *L, double *a, double *b, double x, double y, double z) {
  double yterm;

  if (y > 0.008856) {
    *L = 116.0*pow(y,1.0/3.0) - 16.0;
    yterm = pow(y,1.0/3.0);
  } else {
    *L = 903.3*y;
    yterm = 7.787*y + 16.0/116.0;
  }

  if (x > 0.008856) {
    *a = 500.0*(pow(x,1.0/3.0) - yterm);
  } else {
    *a = 500.0*(7.787*x + 16.0/116.0 - yterm);
  }

  if (z > 0.008856) {
    *b = 200.0*(yterm - pow(z,1.0/3.0));
  } else {
    *b = 200.0*(yterm - (7.787*z + 16.0/116.0));
  }
  return;
}
	    
static double
colorgamma_monitor (double x) {
  if (x < 0.018) {
    return 4.5*x;
  } else {
    return 1.099*pow(x,0.45) - 0.099;
  }
}

static double
colorgamma (double x) {
  if (x >= 0.0) {
    return pow(x,1.0/1.1);
  } else {
    return -pow(fabs(x),1.0/1.1);
  }
}

static double
colorgamma_alt (double x) {
  if (x >= 0.0) {
    return pow(x,0.45);
  } else {
    return -pow(fabs(x),0.45);
  }
}


/* Want det of 1 px py
               1 ax ay
               1 bx by */
static int
side (double ax, double ay, double bx, double by, double px, double py) {
  double det;

  det = (ax*by - ay*bx) + (bx*py - by*px) + (px*ay - py*ax);
  if (det > 0) {
    return 1;
  } else if (det == 0) {
    return 0;
  } else {
    return -1;
  }
}
  
static bool
intersectp (double ax, double ay, double bx, double by,
	    double px, double py, double qx, double qy) {
  int side1, side2, side3, side4;

  side1 = side(ax,ay,bx,by,px,py);
  side2 = side(ax,ay,bx,by,qx,qy);
  side3 = side(px,py,qx,qy,ax,ay);
  side4 = side(px,py,qx,qy,bx,by);
  if ((side1 != side2) && (side3 != side4)) {
    return true;
  } else {
    return false;
  }
}


/* (slope1, intercept1) => [intercept1; slope1, -1] */
/* (slope2, intercept2) => [intercept2; slope2, -1] */
static void
line_intersection (double *x, double *y, 
		   double slope1, double intercept1, 
		   double slope2, double intercept2) {
  double homo1, homo2, homo3;

  homo1 = -slope1 + slope2;
  homo2 = -intercept2 + intercept1;
  homo3 = intercept1*slope2 - intercept2*slope1;
  
  *x = homo2/homo1;
  *y = homo3/homo1;
  return;
}

static int
correct_color (double i) {
  if (i > 255) {
    return 255;
  } else if (i < 0) {
    return 0;
  } else {
    return i;
  }
}


/* red, green, blue coordinates found by taking solve(A) %*% unit vectors */
void
XYZ_to_rgb (int *r, int *g, int *b, double X, double Y, double Z) {
  double redx = 0.64, redy = 0.33;
  double greenx = 0.29, greeny = 0.60;
  double bluex = 0.15, bluey = 0.06;
  double whitex = 0.312713, whitey = 0.329016;

  double rgslope = (greeny-redy)/(greenx-redx);
  double rgintercept = (redy*greenx - greeny*redx)/(greenx-redx);
  double gbslope = (bluey-greeny)/(bluex-greenx);
  double gbintercept = (greeny*bluex - bluey*greenx)/(bluex-greenx);
  double brslope = (redy-bluey)/(redx-bluex);
  double brintercept = (bluey*redx - redy*bluex)/(redx-bluex);
  
  double pointx, pointy, wpslope, wpintercept;
  int rgcrossing, gbcrossing, brcrossing;
  double gamutx, gamuty, gamutz;
  double gamutX, gamutY, gamutZ;
  double R, G, B;

  pointx = X/(X+Y+Z);
  pointy = Y/(X+Y+Z);
  wpslope = (pointy - whitey)/(pointx - whitex);
  wpintercept = (whitey*pointx - pointy*whitex)/(pointx-whitex);

  if (intersectp(redx,redy,greenx,greeny,whitex,whitey,pointx,pointy)) {
    line_intersection(&gamutx,&gamuty,rgslope,rgintercept,wpslope,wpintercept);
    /*
    fprintf(stderr,"Adjusting %f %f for rg crossing to %f %f\n",
	    pointx,pointy,gamutx,gamuty);
    */
  } else if (intersectp(greenx,greeny,bluex,bluey,whitex,whitey,pointx,pointy)) {
    line_intersection(&gamutx,&gamuty,gbslope,gbintercept,wpslope,wpintercept);
    /*
    fprintf(stderr,"Adjusting %f %f for gb crossing to %f %f\n",
	    pointx,pointy,gamutx,gamuty);
    */
  } else if (intersectp(bluex,bluey,redx,redy,whitex,whitey,pointx,pointy)) {
    line_intersection(&gamutx,&gamuty,brslope,brintercept,wpslope,wpintercept);
    /*
    fprintf(stderr,"Adjusting %f %f for br crossing to %f %f\n",
	    pointx,pointy,gamutx,gamuty);
    */
  } else {
    gamutx = pointx;
    gamuty = pointy;
  }

  gamutz = 1.0 - gamutx - gamuty;
  gamutX = gamutx * Y/pointy;
  gamutY = gamuty * Y/pointy;
  gamutZ = gamutz * Y/pointy;
  
  /*
  R = +1.910*gamutX  -0.532*gamutY  -0.288*gamutZ;
  G = -0.985*gamutX  +1.999*gamutY  -0.028*gamutZ;
  B = +0.058*gamutX  -0.118*gamutY  +0.898*gamutZ;
  */

  /* These are the ITU standards: */
  R =  3.063*gamutX - 1.393*gamutY - 0.476*gamutZ;
  G = -0.969*gamutX + 1.876*gamutY + 0.042*gamutZ;
  B =  0.068*gamutX - 0.229*gamutY + 1.069*gamutZ;

  *r = correct_color(rint(colorgamma_monitor(R)*255));
  *g = correct_color(rint(colorgamma_monitor(G)*255));
  *b = correct_color(rint(colorgamma_monitor(B)*255));
  return;
}

static void
polar_to_rect (double *a, double *b, double r, double theta) {
  *a = r*cos(theta);
  *b = r*sin(theta);
  return;
}

static void
integertohex (char *hex1, char *hex0, int level) {
  char buffer[2];
  int bit1, bit0;

  bit1 = level/16;
  bit0 = level - bit1*16;
  sprintf(buffer,"%x",bit1);
  *hex1 = buffer[0];
  sprintf(buffer,"%x",bit0);
  *hex0 = buffer[0];
  return;
}

void
Color_rgb (double *red, double *green, double *blue,
	   double intensity, int group, int ngroups) {
  double r = 60.0, theta;
  double L_pure = 50.0, a_pure, b_pure, L, a, b, x, y, z;
  int ired, igreen, iblue;

  if (intensity > 1.0) {
    intensity = 1.0;
  } else if (intensity < 0.0) {
    intensity = 0.0;
  }

  theta = TWOPI * (double) group/(double) ngroups - PI;
  polar_to_rect(&a_pure,&b_pure,r,theta);

  L = 100.0 + (L_pure - 100)*intensity;
  if (L > 100.0) {
    L = 100.0;
  } else if (L < 50) {
    L = 50.0;
  }

  a = a_pure*intensity;
  b = b_pure*intensity;

  Lab_to_XYZ(&x,&y,&z,L,a,b);
  XYZ_to_rgb(&ired,&igreen,&iblue,x,y,z);
  *red = (double) ired/256.0;
  *green = (double) igreen/256.0;
  *blue = (double) iblue/256.0;
  
  return;
}

/* intensity expressed on 0 to 1 scale */
/* take point on radius at L=50 and interpolate towards (100,0,0) */
char *
Color_string (double intensity, int group, int ngroups) {
  char *rgb, hex1, hex0;
  double r = 40.0, theta;
  double L_pure = 50.0, a_pure, b_pure, L, a, b, x, y, z;
  int red, green, blue;

  /* intensity = 1.0 - cos(HALFPI*intensity); */

  if (intensity > 1.0) {
    intensity = 1.0;
  }

  rgb = (char *) CALLOC(7,sizeof(char));
  theta = TWOPI * (double) group/(double) ngroups - PI;
  polar_to_rect(&a_pure,&b_pure,r,theta);

  L = 100.0 + (L_pure - 100)*intensity;
  if (L > 100.0) {
    L = 100.0;
  } else if (L < 50) {
    L = 50.0;
  }

  a = a_pure*intensity;
  b = b_pure*intensity;

  Lab_to_XYZ(&x,&y,&z,L,a,b);
  XYZ_to_rgb(&red,&green,&blue,x,y,z);
  integertohex(&rgb[0],&rgb[1],red);
  integertohex(&rgb[2],&rgb[3],green);
  integertohex(&rgb[4],&rgb[5],blue);

  return rgb;
}


