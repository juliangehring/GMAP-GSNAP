/* $Id: sarray-write.h 122381 2013-12-24 01:22:55Z twu $ */
#ifndef SARRAY_WRITE_INCLUDED
#define SARRAY_WRITE_INCLUDED
#include "types.h"
#include "genome.h"

extern void
Sarray_write_array (char *sarrayfile, Genome_T genomecomp, UINT4 genomelength);
extern void
Sarray_write_plcp (char *plcpptrsfile, char *plcpcompfile, char *sarrayfile,
		   Genome_T genomecomp, UINT4 genomelength);
extern void
Sarray_write_child (
#ifdef USE_CHILD_BP
		    char *childbpfile, char *childfcfile, char *childs_pagesfile, char *childs_ptrsfile, char *childs_compfile,
		    char *childr_ptrsfile,  char *childr_compfile, char *childx_ptrsfile, char *childx_compfile,
		    char *pioneerbpfile, char *pior_ptrsfile, char *pior_compfile,
		    char *piom_ptrsfile, char *piom_compfile,
#else
		    char *childptrsfile, char *childcompfile, char *sanextpfile,
#endif
		    char *sarrayfile, char *plcpptrsfile, char *plcpcompfile, UINT4 genomelength);
extern void
Sarray_write_index (char *indexiptrsfile, char *indexicompfile, char *indexjptrsfile,char *indexjcompfile,
		    char *sarrayfile, Genome_T genomecomp, UINT4 genomelength);
extern void
Sarray_array_uncompress (Genome_T genomecomp, char *sarrayfile, char *plcpptrsfile, char *plcpcompfile,
			 UINT4 genomelength, UINT4 start, UINT4 end);

#ifdef USE_CHILD_BP
void
Sarray_child_uncompress (char *childbpfile, char *childfcfile, char *childs_pagesfile, char *childs_ptrsfile, char *childs_compfile,
			 char *childr_ptrsfile, char *childr_compfile, char *childx_ptrsfile, char *childx_compfile,
			 char *pioneerbpfile, char *pior_ptrsfile, char *pior_compfile,
			 char *piom_ptrsfile, char *piom_compfile, UINT4 genomelength, BP_size_t startl, BP_size_t endl);
#else
extern void
Sarray_child_uncompress (Genome_T genomecomp, char *plcpptrsfile, char *plcpcompfile, char *childptrsfile, char *childcompfile,
			 char *sanextpfile, char *sarrayfile, UINT4 genomelength, UINT4 start, UINT4 end);
#endif

#endif

