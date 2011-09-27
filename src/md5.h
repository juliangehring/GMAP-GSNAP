/* $Id: md5.h 30358 2010-10-14 15:59:25Z twu $ */
#ifndef MD5_INCLUDED
#define MD5_INCLUDED

#include <stdio.h>

extern unsigned char *
MD5_compute (unsigned char *input, int input_len);
extern void
MD5_print (FILE *fp, unsigned char *digest);

#endif

