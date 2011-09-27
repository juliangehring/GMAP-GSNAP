/* $Id: reader.h 40271 2011-05-28 02:29:18Z twu $ */
#ifndef READER_INCLUDED
#define READER_INCLUDED
#include <stdio.h>
#include "bool.h"

typedef enum {FIVE, THREE, MIDDLE} cDNAEnd_T;

#define T Reader_T
typedef struct T *T;

extern int
Reader_querystart (T this);
extern int
Reader_queryend (T this);
extern int
Reader_startpos (T this);
extern int
Reader_endpos (T this);
extern void
Reader_reset_start (T this, int querypos);
extern void
Reader_reset_end (T this, int querypos);
extern void
Reader_reset_ends (T this);

extern T
Reader_new (char *sequence, int querystart, int queryend, int blocksize);
extern void
Reader_free (T *old);
extern char
Reader_getc (T this, cDNAEnd_T cdnaend);

#undef T
#endif
