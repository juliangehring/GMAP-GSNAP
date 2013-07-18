/* $Id: bigendian.h 99737 2013-06-27 19:33:03Z twu $ */
#ifndef BIGENDIAN_INCLUDED
#define BIGENDIAN_INCLUDED
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stddef.h>
#include "types.h"

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

extern UINT4
Bigendian_convert_uint (UINT4 littleendian);
extern size_t
Bigendian_fwrite_uint (UINT4 value, FILE *fp);
extern void
Bigendian_write_uint (UINT4 value, int fd);
extern size_t
Bigendian_fwrite_uints (UINT4 *array, int n, FILE *fp);
extern size_t
Bigendian_fread_uint (UINT4 *value, FILE *fp);
extern size_t
Bigendian_fread_uints (UINT4 *array, int n, FILE *fp);
extern UINT4
Bigendian_fileio_read_uint (int fd);


#ifdef HAVE_64_BIT
extern UINT8
Bigendian_convert_uint8 (UINT8 littleendian);
extern size_t
Bigendian_fwrite_uint8 (UINT8 value, FILE *fp);
extern size_t
Bigendian_fwrite_uint8s (UINT8 *array, int n, FILE *fp);
extern size_t
Bigendian_fread_uint8 (UINT8 *value, FILE *fp);
extern size_t
Bigendian_fread_uint8s (UINT8 *array, int n, FILE *fp);
extern UINT8
Bigendian_fileio_read_uint8 (int fd);
#endif


#define FREAD_INT(p,fp) Bigendian_fread_int(p,fp)
#define FREAD_UINT(p,fp) Bigendian_fread_uint(p,fp)
#define FREAD_INTS(a,n,fp) Bigendian_fread_ints(a,n,fp)
#define FREAD_UINTS(a,n,fp) Bigendian_fread_uints(a,n,fp)
#ifdef HAVE_64_BIT
#define FREAD_UINT8(p,fp) Bigendian_fread_uint8(p,fp)
#define FREAD_UINT8S(a,n,fp) Bigendian_fread_uint8s(a,n,fp)
#endif

#define FWRITE_INT(x,fp) Bigendian_fwrite_int(x,fp)
#define FWRITE_UINT(x,fp) Bigendian_fwrite_uint(x,fp)
#define WRITE_UINT(x,fd) Bigendian_write_uint(x,fd)
#define WRITE_UINT8(x,fd) Bigendian_write_uin8t(x,fd)
#define FWRITE_INTS(a,n,fp) Bigendian_fwrite_ints(a,n,fp)
#define FWRITE_UINTS(a,n,fp) Bigendian_fwrite_uints(a,n,fp)
#ifdef HAVE_64_BIT
#define FWRITE_UINT8(x,fp) Bigendian_fwrite_uint8(x,fp)
#define FWRITE_UINT8S(a,n,fp) Bigendian_fwrite_uint8s(a,n,fp)
#endif

#define FREAD_CHARS(a,n,fp) fread(a,sizeof(char),n,fp)
#define FWRITE_CHARS(a,n,fp) fwrite(a,sizeof(char),n,fp)

#endif
