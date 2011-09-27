/* $Id: chrnum.h,v 1.10 2005/02/07 23:56:55 twu Exp $ */
#ifndef CHRNUM_INCLUDED
#define CHRNUM_INCLUDED
#include "bool.h"
#include "iit-read.h"

typedef int Chrnum_T;

extern char *
Chrnum_to_string (Chrnum_T chrnum, IIT_T chromosome_iit);
extern char *
Chrnum_to_string_signed (Chrnum_T chrnum, IIT_T chromosome_iit, bool watsonp);

#endif
