/* $Id: genome-write.h 99737 2013-06-27 19:33:03Z twu $ */
#ifndef GENOME_WRITE_INCLUDED
#define GENOME_WRITE_INCLUDED
#include <stdio.h>
#include "bool.h"
#include "iit-read-univ.h"
#include "iit-read.h"
#include "types.h"

extern void
Genome_write (char *genomesubdir, char *fileroot, FILE *input, 
	      Univ_IIT_T contig_iit, IIT_T altstrain_iit, Univ_IIT_T chromosome_iit,
	      bool uncompressedp, bool rawp, bool writefilep,
	      Univcoord_T genomelength, int index1part, int nmessages);

extern Genomecomp_T *
Genome_create_blocks (char *genomicseg, Univcoord_T genomelength);

#endif
