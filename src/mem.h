/* $Id: mem.h,v 1.10 2010-07-16 20:10:47 twu Exp $ */
#ifndef MEM_INCLUDED
#define MEM_INCLUDED

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stddef.h>
#include "except.h"

extern const Except_T Mem_Failed;

extern void
Mem_trap_start (void **location, const char *file, int line);
extern void
Mem_trap_check (const char *file, int line);

/* #define MEMUSAGE 1 */
#ifdef MEMUSAGE
extern void
Mem_usage_init ();
extern long
Mem_usage_report ();
#endif


/* #define LEAKCHECK 1 */
/* For LEAKCHECK, may also want to turn on DEBUG */
#ifdef LEAKCHECK
extern void
Mem_leak_check_activate ();
extern void
Mem_leak_check_deactivate ();
extern void
Mem_leak_check_start (const char *file, int line);
extern void
Mem_leak_check_end (const char *file, int line);
#endif

extern void *Mem_alloc (size_t nbytes,
			const char *file, int line);
extern void *Mem_alloc_no_exception (size_t nbytes,
				     const char *file, int line);
extern void *Mem_calloc (size_t count, size_t nbytes,
			 const char *file, int line);
extern void *Mem_calloc_no_exception (size_t count, size_t nbytes,
				      const char *file, int line);
extern void Mem_free (void *ptr,
		      const char *file, int line);
extern void *Mem_resize (void *ptr, size_t nbytes,
			 const char *file, int line);

#define MTRAP(location) Mem_trap_start((location), __FILE__, __LINE__)
#define MCHECK() Mem_trap_check(__FILE__, __LINE__)

#define MALLOC(nbytes) \
	Mem_alloc((nbytes), __FILE__, __LINE__)
#define MALLOC_NO_EXCEPTION(nbytes) \
	Mem_alloc_no_exception((nbytes), __FILE__, __LINE__)
#define CALLOC(count, nbytes) \
	Mem_calloc((count), (nbytes), __FILE__, __LINE__)
#define CALLOC_NO_EXCEPTION(count, nbytes) \
	Mem_calloc((count), (nbytes), __FILE__, __LINE__)
#define  NEW(p) ((p) = MALLOC(sizeof *(p)))
#define NEW0(p) ((p) = CALLOC(1, sizeof *(p)))
#define FREE(ptr) ((void)(Mem_free((ptr), \
	__FILE__, __LINE__), (ptr) = 0))
#define RESIZE(ptr, nbytes) 	((ptr) = Mem_resize((ptr), \
	(nbytes), __FILE__, __LINE__))
#endif
