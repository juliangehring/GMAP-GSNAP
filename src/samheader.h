/* $Id: samheader.h 149320 2014-09-30 02:16:01Z twu $ */
#ifndef SAMHEADER_INCLUDED
#define SAMHEADER_INCLUDED

#include <stdio.h>
#include "bool.h"

extern void
SAM_header_change_HD_tosorted (FILE *fp, int headerlen);
extern void
SAM_header_print_HD (FILE *fp, int nworkers, bool orderedp);
extern void
SAM_header_print_PG (FILE *fp, int argc, char **argv, int optind);
extern int
SAM_header_length (int *lastchar, FILE *fp);

#endif

