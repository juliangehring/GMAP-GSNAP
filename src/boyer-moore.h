/* $Id: boyer-moore.h,v 1.2 2004/05/27 00:38:56 twu Exp $ */
#ifndef BOYER_MOORE_INCLUDED
#define BOYER_MOORE_INCLUDED
#include "intlist.h"

extern Intlist_T
BoyerMoore (char *query, int querylen, char *text, int textlen);

#endif

