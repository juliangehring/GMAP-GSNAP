static char rcsid[] = "$Id: mem.c,v 1.10 2005/02/07 23:56:56 twu Exp $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "mem.h"
#include <stdio.h>
#include <stdlib.h>
#include "assert.h"
#include "except.h"

/*
#define TRAP
*/
static void *trap_contents;
static void **trap_location;

void
Mem_trap_start (void **location) {
  trap_location = location;
  trap_contents = * (void **) location;
  return;
}


const Except_T Mem_Failed = { "Allocation Failed" };
void *
Mem_alloc (size_t nbytes, const char *file, int line) {
  void *ptr;
  assert(nbytes > 0);
  ptr = malloc(nbytes);
#ifdef TRAP
  if (*trap_location != trap_contents) {
    printf("Trap violation in malloc at %p (%d bytes).  Location %p contains %p\n",
	   ptr,nbytes,*trap_location,trap_contents);
    fflush(stdout);
    abort();
  }
#endif

  if (ptr == NULL) {
    if (file == NULL) {
      RAISE(Mem_Failed);
    } else {
      Except_raise(&Mem_Failed, file, line);
    }
  }
  return ptr;
}

void *
Mem_alloc_no_exception (size_t nbytes, const char *file, int line) {
  void *ptr;
  assert(nbytes > 0);
  ptr = malloc(nbytes);
  return ptr;
}

void *
Mem_calloc (size_t count, size_t nbytes, const char *file, int line) {
  void *ptr;
  if (count <= 0) {
    if (file == NULL) {
      RAISE(Mem_Failed);
    } else {
      Except_raise(&Mem_Failed, file, line);
    }
  }
  assert(nbytes > 0);
  ptr = calloc(count, nbytes);
#ifdef TRAP
  if (*trap_location != trap_contents) {
    printf("Trap violation in malloc at %p (%d bytes).  Location %p contains %p\n",
	   ptr,nbytes,*trap_location,trap_contents);
    fflush(stdout);
    abort();
  }
#endif
  if (ptr == NULL) {
    if (file == NULL) {
      RAISE(Mem_Failed);
    } else {
      Except_raise(&Mem_Failed, file, line);
    }
  }
  return ptr;
}

void *
Mem_calloc_no_exception (size_t count, size_t nbytes, const char *file, int line) {
  void *ptr;

  if (count <= 0) {
    if (file == NULL) {
      RAISE(Mem_Failed);
    } else {
      Except_raise(&Mem_Failed, file, line);
    }
  }
  assert(nbytes > 0);
  ptr = calloc(count, nbytes);
  return ptr;
}

void 
Mem_free (void *ptr, const char *file, int line) {
#ifdef TRAP
  if (*trap_location != trap_contents) {
    printf("Trap violation in free at %p.  Location %p contains %p\n",
	   ptr,*trap_location,trap_contents);
    fflush(stdout);
    abort();
  }
#endif

  if (ptr) {
    free(ptr);
  }
}

void *
Mem_resize (void *ptr, size_t nbytes, const char *file, int line) {
  assert(ptr);
  assert(nbytes > 0);
  ptr = realloc(ptr, nbytes);
  if (ptr == NULL) {
    if (file == NULL) {
      RAISE(Mem_Failed);
    } else {
      Except_raise(&Mem_Failed, file, line);
    }
  }
  return ptr;
}
