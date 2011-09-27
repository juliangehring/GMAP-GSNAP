#ifndef COLOR_INCLUDED
#define COLOR_INCLUDED

extern void
Lab_to_XYZ (double *x, double *y, double *z, double L, double a, double b);
extern void
XYZ_to_rgb (int *r, int *g, int *b, double x, double y, double z);
extern void
Color_rgb (double *red, double *green, double *blue,
	   double intensity, int group, int ngroups);
extern char *
Color_string (double intensity, int group, int ngroups);

#endif
