#ifndef GSNAPREAD_INCLUDED
#define GSNAPREAD_INCLUDED
#include "genomicpos.h"
#include "intlist.h"
#include "uintlist.h"

extern bool
Gsnapread_concordantp (char *header);

extern char *
Gsnapread_quality_string (char *header);

extern int
Gsnapread_nhits (char *header);

extern int
Gsnapread_readlength_preadapter (char *header);

extern int
Gsnapread_readlength (char *line);

extern char *
Gsnapread_accession (char *header);

extern void
Gsnapread_parse_line (char *line, int *query5, int *query3, char *firstend, char *secondend,
		      int *support, int *nmismatches, char *strand, char **chr,
		      Genomicpos_T *firstpos, Genomicpos_T *secondpos,
		      Genomicpos_T *donorpos, Genomicpos_T *acceptorpos, char *truestrand, bool *sensep);

#endif

