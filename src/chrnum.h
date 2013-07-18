/* $Id: chrnum.h 99737 2013-06-27 19:33:03Z twu $ */
#ifndef CHRNUM_INCLUDED
#define CHRNUM_INCLUDED
#include "bool.h"
#include "iit-read-univ.h"
#include "types.h"
#include "genomicpos.h"

typedef int Chrnum_T;

extern char *
Chrnum_to_string (Chrnum_T chrnum, Univ_IIT_T chromosome_iit);
extern char *
Chrnum_to_string_signed (Chrnum_T chrnum, Univ_IIT_T chromosome_iit, bool watsonp);
extern Chrpos_T
Chrnum_length (Chrnum_T chrnum, Univ_IIT_T chromosome_iit);
extern Univcoord_T
Chrnum_offset (Chrnum_T chrnum, Univ_IIT_T chromosome_iit);
extern void
Chrnum_print_position (Univcoord_T position, Univ_IIT_T chromosome_iit);

#endif
