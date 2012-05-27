/* $Id: genome-write.h 64017 2012-05-14 22:35:15Z twu $ */
#ifndef GENOME_WRITE_INCLUDED
#define GENOME_WRITE_INCLUDED
#include <stdio.h>
#include "bool.h"
#include "iit-read.h"
#include "types.h"

extern void
Genome_write (char *genomesubdir, char *fileroot, FILE *input, 
	      IIT_T contig_iit, IIT_T altstrain_iit, bool uncompressedp, bool rawp,
	      bool writefilep, unsigned int genomelength, int index1part);

extern UINT4 *
Genome_create_blocks (char *genomicseg, unsigned int genomelength);

#endif
