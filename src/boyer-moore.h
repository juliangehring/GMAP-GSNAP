/* $Id: boyer-moore.h 79302 2012-11-15 23:55:58Z twu $ */
#ifndef BOYER_MOORE_INCLUDED
#define BOYER_MOORE_INCLUDED
#include "intlist.h"
#include "genomicpos.h"

extern Intlist_T
BoyerMoore (char *query, int querylen, char *text, int textlen);
extern Intlist_T
BoyerMoore_nt (char *query, int querylen, int textoffset, int textlen,
	       Genomicpos_T chroffset, Genomicpos_T chrhigh, bool watsonp);
extern int *
BoyerMoore_bad_char_shift (char *query, int querylen);
extern int
BoyerMoore_maxprefix (char *query, int querylen, char *text, int textlen,
		      int *bad_char_shift);

#endif

