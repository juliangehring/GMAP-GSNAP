/* $Id: genome-write.h,v 1.1 2005/02/14 12:50:41 twu Exp $ */
#ifndef GENOME_WRITE_INCLUDED
#define GENOME_WRITE_INCLUDED
#include <stdio.h>
#include "bool.h"
#include "iit-read.h"

extern void
Genome_write (char *genomesubdir, char *fileroot, FILE *input, 
	      IIT_T contig_iit, IIT_T altstrain_iit, bool uncompressedp,
	      bool writefilep, unsigned int genomelength);

#endif
