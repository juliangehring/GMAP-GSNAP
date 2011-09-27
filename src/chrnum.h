/* $Id: chrnum.h 40271 2011-05-28 02:29:18Z twu $ */
#ifndef CHRNUM_INCLUDED
#define CHRNUM_INCLUDED
#include "bool.h"
#include "iit-read.h"
#include "genomicpos.h"

typedef int Chrnum_T;

extern char *
Chrnum_to_string (Chrnum_T chrnum, IIT_T chromosome_iit);
extern char *
Chrnum_to_string_signed (Chrnum_T chrnum, IIT_T chromosome_iit, bool watsonp);
extern unsigned int
Chrnum_length (Chrnum_T chrnum, IIT_T chromosome_iit);
extern unsigned int
Chrnum_offset (Chrnum_T chrnum, IIT_T chromosome_iit);
extern void
Chrnum_print_position (Genomicpos_T position, IIT_T chromosome_iit);

#endif
