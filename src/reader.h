/* $Id: reader.h,v 1.13 2006/05/19 17:14:14 twu Exp $ */
#ifndef READER_INCLUDED
#define READER_INCLUDED
#include <stdio.h>

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
