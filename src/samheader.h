/* $Id: samheader.h 103701 2013-08-02 23:04:54Z twu $ */
#ifndef SAMHEADER_INCLUDED
#define SAMHEADER_INCLUDED

#include <stdio.h>
#include "bool.h"

extern void
SAM_header_print_HD (FILE *fp, int nworkers, bool orderedp);
extern void
SAM_header_print_PG (FILE *fp, int argc, char **argv, int optind);

#endif

