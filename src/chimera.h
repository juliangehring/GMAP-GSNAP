/* $Id: chimera.h,v 1.4 2005/05/04 23:50:22 twu Exp $ */
#ifndef CHIMERA_INCLUDED
#define CHIMERA_INCLUDED
#include "sequence.h"
#include "stage3.h"

extern double
Chimera_detect (int *breakpoint, int *margin, Stage3_T stage3, Sequence_T queryseq);
extern void
Chimera_bestpath (int *bestfrom, int *bestto, int *bestpos,
		  int *fromscore_3, int *toscore_5,
		  Stage3_T *stage3array, int npaths, int querylength);

#endif
