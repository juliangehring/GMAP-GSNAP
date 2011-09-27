/* $Id: genome-write.h 42000 2011-06-30 17:12:04Z twu $ */
#ifndef GENOME_WRITE_INCLUDED
#define GENOME_WRITE_INCLUDED
#include <stdio.h>
#include "bool.h"
#include "iit-read.h"

extern void
Genome_write (char *genomesubdir, char *fileroot, FILE *input, 
	      IIT_T contig_iit, IIT_T altstrain_iit, bool uncompressedp, bool rawp,
	      bool writefilep, unsigned int genomelength, int index1part);

#endif
