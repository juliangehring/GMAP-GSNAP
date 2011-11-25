/* $Id: inbuffer.h 49442 2011-10-08 01:32:38Z twu $ */
#ifndef INBUFFER_INCLUDED
#define INBUFFER_INCLUDED
#include <stdio.h>
#include "bool.h"
#include "outbuffer.h"
#include "sequence.h"
#include "request.h"

#ifdef HAVE_ZLIB
#include <zlib.h>
#endif


#define T Inbuffer_T
typedef struct T *T;

#ifndef GSNAP
extern T
Inbuffer_cmdline (char *contents, int length);
#endif

extern T
Inbuffer_new (int nextchar, FILE *input,
#ifdef GSNAP
	      FILE *input2,
#ifdef HAVE_ZLIB
	      gzFile gzipped, gzFile gzipped2,
#endif
#ifdef HAVE_GOBY
	      Gobyreader_T gobyreader,
#endif
#endif
	      char **files, int nfiles,
#ifdef GSNAP
	      bool fastq_format_p, bool creads_format_p,
	      int barcode_length, bool invert_first_p, bool invert_second_p,
	      bool chop_primers_p,
#else
	      bool maponlyp,
#endif
	      unsigned int nspaces, unsigned int maxchars, int part_interval, int part_modulus,
	      bool filter_if_both_p);

extern void
Inbuffer_set_outbuffer (T this, Outbuffer_T outbuffer);

extern void
Inbuffer_free (T *old);

extern unsigned int
Inbuffer_fill_init (T this);

extern Request_T
Inbuffer_get_request (T this);

#undef T
#endif

