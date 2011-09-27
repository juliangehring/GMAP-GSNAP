/* $Id: boyer-moore.h,v 1.3 2010-07-21 17:13:39 twu Exp $ */
#ifndef BOYER_MOORE_INCLUDED
#define BOYER_MOORE_INCLUDED
#include "intlist.h"

extern Intlist_T
BoyerMoore (char *query, int querylen, char *text, int textlen);
extern int *
BoyerMoore_bad_char_shift (char *query, int querylen);
extern int
BoyerMoore_maxprefix (char *query, int querylen, char *text, int textlen,
		      int *bad_char_shift);

#endif

