/* $Id: samheader.h 154452 2014-12-02 19:28:04Z twu $ */
#ifndef SAMHEADER_INCLUDED
#define SAMHEADER_INCLUDED

#include <stdio.h>
#include "bool.h"

extern void
SAM_header_change_HD_tosorted_stdout (FILE *fp, int headerlen);
extern void
SAM_header_change_HD_tosorted_split (FILE *fp, int headerlen, FILE **outputs, int noutputs);
extern void
SAM_header_print_HD (FILE *fp, int nworkers, bool orderedp);
extern void
SAM_header_print_PG (FILE *fp, int argc, char **argv, int optind);
extern int
SAM_header_length (int *lastchar, FILE *fp);

#endif

