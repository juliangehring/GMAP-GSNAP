#ifndef PARSERANGE_INCLUDED
#define PARSERANGE_INCLUDED
#include "bool.h"
#include "genomicpos.h"
#include "iit-read.h"

extern bool
Parserange_query (char **divstring, unsigned int *coordstart, unsigned int *coordend, bool *revcomp,
		  char *query, char *filename);


/* genomicstart is 0-based, chrstart is 1-based */
extern bool
Parserange_universal (char **div, bool *revcomp,
		      Genomicpos_T *genomicstart, Genomicpos_T *genomiclength,
		      Genomicpos_T *chrstart, Genomicpos_T *chrend,
		      Genomicpos_T *chroffset, Genomicpos_T *chrlength,
		      char *query, char *genomesubdir, char *fileroot);

#if 0
  /* Old style */
extern char *
Parserange_universal (bool *revcomp, Genomicpos_T *genomicstart, Genomicpos_T *genomiclength,
		      Genomicpos_T *chroffset, Genomicpos_T *chrlength,
		      char *query, char *genomesubdir, char *fileroot);
#endif

extern bool
Parserange_universal_iit (char **div, bool *revcomp,
			  Genomicpos_T *genomicstart, Genomicpos_T *genomiclength,
			  Genomicpos_T *chrstart, Genomicpos_T *chrend,
			  Genomicpos_T *chroffset, Genomicpos_T *chrlength,
			  char *query, IIT_T chromosome_iit, IIT_T contig_iit);

extern bool
Parserange_simple (char **div, bool *revcomp, Genomicpos_T *chrstart, Genomicpos_T *chrend,
		   char *query);

#endif

