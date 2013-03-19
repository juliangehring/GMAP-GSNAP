/* $Id: bzip2.h 83593 2013-01-16 22:59:40Z twu $ */
#ifndef BZIP2_INCLUDED
#define BZIP2_INCLUDED
#include "bool.h"

#define T Bzip2_T
typedef struct T *T;

extern T
Bzip2_new (char *filename);

extern void
Bzip2_free (T *old);

extern int
bzgetc (T this);

extern bool
bzeof (T this);

extern char *
bzgets (T this, char *buffer, int maxlength);

#undef T
#endif
