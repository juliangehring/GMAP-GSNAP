/* $Id: chrnum.h,v 1.11 2005/07/07 23:57:41 twu Exp $ */
#ifndef CHRNUM_INCLUDED
#define CHRNUM_INCLUDED
#include "bool.h"
#include "iit-read.h"

typedef int Chrnum_T;

extern char *
Chrnum_to_string (Chrnum_T chrnum, IIT_T chromosome_iit);
extern char *
Chrnum_to_string_signed (Chrnum_T chrnum, IIT_T chromosome_iit, bool watsonp);
extern unsigned int
Chrnum_length (Chrnum_T chrnum, IIT_T chromosome_iit);

#endif
