/* $Id: chimera.h,v 1.5 2005/05/09 22:26:24 twu Exp $ */
#ifndef CHIMERA_INCLUDED
#define CHIMERA_INCLUDED
#include "sequence.h"
#include "stage3.h"

extern int
Chimera_detect (int *margin, Stage3_T stage3, Sequence_T queryseq, double fthreshold);
extern void
Chimera_bestpath (int *bestfrom, int *bestto, int *bestpos, int *equivpos,
		  int *fromscore_3, int *toscore_5,
		  Stage3_T *stage3array, int npaths, int querylength);

#endif
