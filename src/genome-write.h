/* $Id: genome-write.h,v 1.2 2005/11/09 22:04:23 twu Exp $ */
#ifndef GENOME_WRITE_INCLUDED
#define GENOME_WRITE_INCLUDED
#include <stdio.h>
#include "bool.h"
#include "iit-read.h"

extern void
Genome_write (char *genomesubdir, char *fileroot, FILE *input, 
	      IIT_T contig_iit, IIT_T altstrain_iit, bool uncompressedp, bool rawp,
	      bool writefilep, unsigned int genomelength);

#endif
