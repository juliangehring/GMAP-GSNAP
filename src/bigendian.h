/* $Id: bigendian.h,v 1.6 2005/02/14 13:02:12 twu Exp $ */
#ifndef BIGENDIAN_INCLUDED
#define BIGENDIAN_INCLUDED
#include <stdio.h>
#include <stddef.h>

extern int
Bigendian_convert_int (int littleendian);
extern size_t
Bigendian_fwrite_int (int value, FILE *fp);
extern size_t
Bigendian_fwrite_ints (int *array, int n, FILE *fp);
extern size_t
Bigendian_fread_int (int *value, FILE *fp);
extern size_t
Bigendian_fread_ints (int *array, int n, FILE *fp);

extern unsigned int
Bigendian_convert_uint (unsigned int littleendian);
extern size_t
Bigendian_fwrite_uint (unsigned int value, FILE *fp);
extern size_t
Bigendian_fwrite_uints (unsigned int *array, int n, FILE *fp);
extern size_t
Bigendian_fread_uint (unsigned int *value, FILE *fp);
extern size_t
Bigendian_fread_uints (unsigned int *array, int n, FILE *fp);

#define FREAD_INT(p,fp) Bigendian_fread_int(p,fp)
#define FREAD_UINT(p,fp) Bigendian_fread_uint(p,fp)
#define FREAD_INTS(a,n,fp) Bigendian_fread_ints(a,n,fp)
#define FREAD_UINTS(a,n,fp) Bigendian_fread_uints(a,n,fp)
#define FWRITE_INT(x,fp) Bigendian_fwrite_int(x,fp)
#define FWRITE_UINT(x,fp) Bigendian_fwrite_uint(x,fp)
#define FWRITE_INTS(a,n,fp) Bigendian_fwrite_ints(a,n,fp)
#define FWRITE_UINTS(a,n,fp) Bigendian_fwrite_uints(a,n,fp)

#define FREAD_CHARS(a,n,fp) fread(a,sizeof(char),n,fp)
#define FWRITE_CHARS(a,n,fp) fwrite(a,sizeof(char),n,fp)

#endif
