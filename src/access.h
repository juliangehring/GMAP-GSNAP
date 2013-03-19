/* $Id: access.h 77636 2012-10-26 00:14:01Z twu $ */
#ifndef ACCESS_INCLUDED
#define ACCESS_INCLUDED

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef HAVE_UNISTD_H
#include <unistd.h>		/* For size_t, and for mmap on Linux, lseek, and getpagesize */
#endif
#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>		/* For size_t, and for mmap and off_t */
#endif

#include "bool.h"

/* ALLOCATED implies bigendian conversion already done */
typedef enum {USE_ALLOCATE, USE_MMAP_ONLY, USE_MMAP_PRELOAD, USE_FILEIO} Access_mode_T;
typedef enum {ALLOCATED, MMAPPED, FILEIO} Access_T;
#define MAX32BIT 4294967295U	/* 2^32 - 1 */

extern bool
Access_file_exists_p (char *filename);

extern off_t
Access_filesize (char *filename);

extern size_t
Access_file_copy (char *dest_file, char *source_file);

extern bool
Access_file_equal (char *file1, char *file2);

extern int
Access_fileio (char *filename);

extern int
Access_fileio_rw (char *filename);

extern void *
Access_allocated (size_t *len, double *seconds, char *filename, size_t eltsize);

#ifdef HAVE_CADDR_T
extern caddr_t
#else
extern void *
#endif
Access_mmap (int *fd, size_t *len, char *filename, size_t eltsize, bool randomp);

#ifdef HAVE_CADDR_T
extern caddr_t
#else
extern void *
#endif
Access_mmap_offset (int *remainder, int fd, off_t offset, size_t length, size_t eltsize, bool randomp);

#ifdef HAVE_CADDR_T
extern caddr_t
#else
extern void *
#endif
Access_mmap_rw (int *fd, size_t *len, char *filename, size_t eltsize, bool randomp);

#ifdef HAVE_CADDR_T
extern caddr_t
#else
extern void *
#endif
Access_mmap_offset_rw (int *remainder, int fd, off_t offset, size_t length, size_t eltsize, bool randomp);

#ifdef HAVE_CADDR_T
extern caddr_t
#else
extern void *
#endif
Access_mmap_and_preload (int *fd, size_t *len, int *npages, double *seconds,
			 char *filename, size_t eltsize);

#endif
