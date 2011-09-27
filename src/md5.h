/* $Id: md5.h 40271 2011-05-28 02:29:18Z twu $ */
#ifndef MD5_INCLUDED
#define MD5_INCLUDED

#include <stdio.h>

extern unsigned char *
MD5_compute (unsigned char *input, int input_len);
extern void
MD5_print (FILE *fp, unsigned char *digest);

#endif

