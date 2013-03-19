/* $Id: genome-write.h 89125 2013-03-13 22:23:36Z twu $ */
#ifndef GENOME_WRITE_INCLUDED
#define GENOME_WRITE_INCLUDED
#include <stdio.h>
#include "bool.h"
#include "iit-read.h"
#include "types.h"

extern void
Genome_write (char *genomesubdir, char *fileroot, FILE *input, 
	      IIT_T contig_iit, IIT_T altstrain_iit, IIT_T chromosome_iit,
	      bool uncompressedp, bool rawp, bool writefilep,
	      unsigned int genomelength, int index1part, int nmessages);

extern Genomecomp_T *
Genome_create_blocks (char *genomicseg, unsigned int genomelength);

#endif
