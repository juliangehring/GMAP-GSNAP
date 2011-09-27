/* $Id: chrnum.h,v 1.13 2006/05/11 17:15:00 twu Exp $ */
#ifndef CHRNUM_INCLUDED
#define CHRNUM_INCLUDED
#include "bool.h"
#include "iit-read.h"

typedef int Chrnum_T;

extern char *
Chrnum_to_string (Chrnum_T chrnum, IIT_T chromosome_iit, bool allocp);
extern char *
Chrnum_to_string_signed (Chrnum_T chrnum, IIT_T chromosome_iit, bool watsonp);
extern unsigned int
Chrnum_length (Chrnum_T chrnum, IIT_T chromosome_iit);
extern unsigned int
Chrnum_offset (Chrnum_T chrnum, IIT_T chromosome_iit);

#endif
