/* $Id: sarray-write.h 101477 2013-07-15 15:33:07Z twu $ */
#ifndef SARRAY_WRITE_INCLUDED
#define SARRAY_WRITE_INCLUDED
#include "types.h"
#include "genome.h"

extern void
Sarray_write_lcp (char *lcpfile, char *lcpptrsfile, char *lcpcompfile, char *saindexfile,
		  char *sarrayfile, Genome_T genome, UINT4 genomelength);
extern void
Sarray_write_all (char *genomesubdir, char *fileroot, Genome_T genome, UINT4 genomiclength);

#endif

