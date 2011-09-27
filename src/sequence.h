/* $Id: sequence.h,v 1.49 2009/08/21 19:26:11 twu Exp $ */
#ifndef SEQUENCE_INCLUDED
#define SEQUENCE_INCLUDED
#include <stdio.h>
#include "bool.h"

#ifdef PMAP
#define MAXSEQLEN 300000
#else
#define MAXSEQLEN 1000000
#endif
#define HALFLEN MAXSEQLEN/2

#define T Sequence_T
typedef struct T *T;

extern char *
Sequence_fullpointer (T this);
extern char *
Sequence_fullpointer_uc (T this);
extern char *
Sequence_trimpointer (T this);
extern int
Sequence_ntlength (T this);
extern int
Sequence_fulllength (T this);
extern int
Sequence_trimlength (T this);
extern int
Sequence_fulllength_given (T this);
extern void
Sequence_trim (T this, int trim_start, int trim_end);
extern int
Sequence_trim_start (T this);
extern int
Sequence_trim_end (T this);
extern int
Sequence_subseq_offset (T this);
extern int
Sequence_skiplength (T this);

extern void
Sequence_free (T *old);
extern T
Sequence_genomic_new (char *contents, int length);
extern T
Sequence_read (int *nextchar, FILE *input, bool maponlyp);
extern T
Sequence_read_multifile (int *nextchar, FILE **input, char ***files, int *nfiles, bool maponlyp);
extern T
Sequence_read_multifile_shortreads (int *nextchar, T *queryseq2, FILE **input, char ***files, int *nfiles,
				    bool circularp);
extern T
Sequence_read_unlimited (FILE *input);
#ifdef PMAP
extern T
Sequence_convert_to_nucleotides (T this);
#endif
extern int
Sequence_count_bad (T this, int pos, int max, int direction);

extern T
Sequence_subsequence (T this, int start, int end);
extern T
Sequence_revcomp (T this);
extern T
Sequence_uppercase (T this);
extern T
Sequence_alias (T this);


extern void
Sequence_print_digest (T this);
extern void
Sequence_print_header (T this, bool checksump);
extern void
Sequence_print_header_revcomp (T this);
extern void
Sequence_print (T this, bool uppercasep, int wraplength, bool trimmedp);
extern void
Sequence_print_two (T this, T alt, bool uppercasep, int wraplength);
extern void
Sequence_print_oneline (FILE *fp, T this);
extern void
Sequence_print_oneline_uc (FILE *fp, T this);
extern void
Sequence_print_oneline_revcomp (FILE *fp, T this);
extern void
Sequence_print_oneline_revcomp_uc (FILE *fp, T this);
extern void
Sequence_print_raw (T this);
extern char *
Sequence_accession (T this);

extern T
Sequence_substring (T usersegment, unsigned int left, unsigned int length, 
		    bool revcomp, char *gbuffer1, char *gbuffer2, int gbufferlen);

#undef T
#endif
