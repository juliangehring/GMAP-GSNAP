/* $Id: md5.h,v 1.5 2005/02/15 01:58:50 twu Exp $ */
#ifndef MD5_INCLUDED
#define MD5_INCLUDED

extern unsigned char *
MD5_compute (unsigned char *input, int input_len);
extern void
MD5_print (unsigned char *digest);

#endif

