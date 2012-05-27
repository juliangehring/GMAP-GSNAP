/* $Id: boyer-moore.h 64017 2012-05-14 22:35:15Z twu $ */
#ifndef BOYER_MOORE_INCLUDED
#define BOYER_MOORE_INCLUDED
#include "intlist.h"
#include "genomicpos.h"
#include "genome.h"

extern Intlist_T
BoyerMoore (char *query, int querylen, char *text, int textlen);
extern Intlist_T
BoyerMoore_nt (char *query, int querylen, int textoffset, int textlen, Genome_T genome,
	       Genomicpos_T chroffset, Genomicpos_T chrpos, Genomicpos_T genomiclength,
	       bool watsonp);
extern int *
BoyerMoore_bad_char_shift (char *query, int querylen);
extern int
BoyerMoore_maxprefix (char *query, int querylen, char *text, int textlen,
		      int *bad_char_shift);

#endif

