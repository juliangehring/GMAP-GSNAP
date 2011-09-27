/* $Id: sequence.h,v 1.23 2005/02/07 23:56:57 twu Exp $ */
#ifndef SEQUENCE_INCLUDED
#define SEQUENCE_INCLUDED
#include <stdio.h>
#include "bool.h"

#define MAXSEQLEN 1000000

#define T Sequence_T
typedef struct T *T;

extern char *
Sequence_pointer (T this);
extern char *
Sequence_pointer_full (T this);
extern int
Sequence_length (T this);
extern int
Sequence_length_full (T this);
extern int
Sequence_ntrimmed (T this);
extern char
Sequence_char (T this, int querypos);
extern void
Sequence_trim_polya (T this);

extern void
Sequence_free (T *old);
extern T
Sequence_genomic_new (char *contents, int length);
extern T
Sequence_read (FILE *input, bool polya_trim);
extern T
Sequence_read_unlimited (FILE *input);
extern int
Sequence_count_bad (T this, int pos, int max, int direction);

extern T
Sequence_subsequence (T this, int start, int end, bool fivep);
extern T
Sequence_revcomp (T this);
extern void
Sequence_endstream ();
extern void
Sequence_print_digest (T this);
extern void
Sequence_print_header (T this, bool checksump);
extern void
Sequence_print (T this, bool uppercasep, int wraplength);
extern char *
Sequence_accession (T this);

extern T
Sequence_substring (T usersegment, unsigned int left, unsigned int length, 
		    bool revcomp, char *gbuffer1, char *gbuffer2, int gbufferlen);

#undef T
#endif
