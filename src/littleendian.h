/* $Id: littleendian.h,v 1.1 2005/02/14 13:00:14 twu Exp $ */
#ifndef LITTLEENDIAN_INCLUDED
#define LITTLEENDIAN_INCLUDED
#include <stdio.h>

#define FREAD_INT(p,fp) fread(p,sizeof(int),1,fp)
#define FREAD_UINT(p,fp) fread(p,sizeof(unsigned int),1,fp)
#define FREAD_INTS(a,n,fp) fread(a,sizeof(int),n,fp)
#define FREAD_UINTS(a,n,fp) fread(a,sizeof(unsigned int),n,fp)
#define FWRITE_INT(x,fp) fwrite(&(x),sizeof(int),1,fp)
#define FWRITE_UINT(x,fp) fwrite(&(x),sizeof(unsigned int),1,fp)
#define FWRITE_INTS(a,n,fp) fwrite(a,sizeof(int),n,fp)
#define FWRITE_UINTS(a,n,fp) fwrite(a,sizeof(unsigned int),n,fp)

#define FREAD_CHARS(a,n,fp) fread(a,sizeof(char),n,fp)
#define FWRITE_CHARS(a,n,fp) fwrite(a,sizeof(char),n,fp)

#endif

