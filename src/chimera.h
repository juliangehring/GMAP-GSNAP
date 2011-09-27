/* $Id: chimera.h,v 1.2 2005/02/08 00:00:20 twu Exp $ */
#ifndef CHIMERA_INCLUDED
#define CHIMERA_INCLUDED
#include "stage3.h"

extern void
Chimera_bestpath (int *bestfrom, int *bestto, int *bestpos,
		  int *fromscore_3, int *toscore_5,
		  Stage3_T *stage3array, int npaths, int querylength);

#endif
