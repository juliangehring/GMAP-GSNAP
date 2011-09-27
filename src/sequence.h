/* $Id: sequence.h,v 1.27 2005/05/06 18:44:53 twu Exp $ */
#ifndef SEQUENCE_INCLUDED
#define SEQUENCE_INCLUDED
#include <stdio.h>
#include "bool.h"

#define MAXSEQLEN 1000000

#define T Sequence_T
typedef struct T *T;

extern char *
Sequence_fullpointer (T this);
extern char *
Sequence_trimpointer (T this);
extern int
Sequence_fulllength (T this);
extern int
Sequence_trimlength (T this);
extern void
Sequence_trim (T this, int trim_start, int trim_end);
extern int
Sequence_trim_start (T this);
extern int
Sequence_trim_end (T this);

extern void
Sequence_free (T *old);
extern T
Sequence_genomic_new (char *contents, int length);
extern T
Sequence_read (int *nextchar, FILE *input);
extern T
Sequence_read_unlimited (FILE *input);
extern int
Sequence_count_bad (T this, int pos, int max, int direction);

extern T
Sequence_subsequence (T this, int start, int end);
extern T
Sequence_revcomp (T this);
extern void
Sequence_endstream ();
extern void
Sequence_print_digest (T this);
extern void
Sequence_print_header (T this, bool checksump);
extern void
Sequence_print (T this, bool uppercasep, int wraplength, bool trimmedp);
extern char *
Sequence_accession (T this);

extern T
Sequence_substring (T usersegment, unsigned int left, unsigned int length, 
		    bool revcomp, char *gbuffer1, char *gbuffer2, int gbufferlen);

#undef T
#endif
