#ifndef PLOTDATA_INCLUDED
#define PLOTDATA_INCLUDED
#include <stdio.h>
#include "bool.h"
#include "iit-read.h"
#include "chrsubset.h"
#include "uintlist.h"
#include "doublelist.h"

#define T Plotdata_T
typedef struct T *T;

extern int
Plotdata_nsamples (T this);
extern int
Plotdata_ngenes (T this);
extern int
Plotdata_nvalues (T this, int c);
extern unsigned int *
Plotdata_chrpositions (T this, int c);
extern double *
Plotdata_values (T this, int c);
extern double
Plotdata_get (unsigned int *chrpos, T this, int c, int g, int s);

extern T
Plotdata_read (FILE *fp, IIT_T chromosome_iit, Chrsubset_T chrsubset);
extern void
Plotdata_one_signature (T this, int s, Uintlist_T *segment_startpositions,
			Uintlist_T *segment_endpositions, Doublelist_T *segment_means,
			char *title, bool logp,
			IIT_T chromosome_iit, Chrsubset_T chrsubset,
			unsigned int mincoord, unsigned int maxcoord,
			int nincluded, double valuefactor, 
			double top, double bottom, double left, double right, 
			double annotheight, double annotwidth, 
			bool segmentp, bool skipvaluesp);

#undef T
#endif

