static char rcsid[] = "$Id: mem.c,v 1.13 2005/07/13 19:25:15 twu Exp $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "mem.h"
#include <stdio.h>
#include <stdlib.h>
#include "assert.h"
#include "except.h"

#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/*
#define TRAP
*/

#ifdef TRAP
static void *trap_contents;
static void **trap_location;
static int startp = 0;

void
Mem_trap_start (void **location, const char *file, int line) {
  if (startp == 0) {
    trap_location = location;
    trap_contents = * (void **) location;
    startp = 1;
    printf("Initial value at location %p is %p from %s:%d\n",
	   trap_location,trap_contents,file,line);
    fflush(stdout);
  }
  return;
}

void
Mem_trap_check (const char *file, int line) {
  if (startp > 0 && *trap_location != trap_contents) {
      printf("Value changed at location %p.  Old value was %p.  New value is %p.  Observed during check at %s:%d\n",
	     trap_location,trap_contents,*trap_location,file,line);
      fflush(stdout);
      trap_contents = * (void **) trap_location;
  }
  return;
}
#endif


const Except_T Mem_Failed = { "Allocation Failed" };
void *
Mem_alloc (size_t nbytes, const char *file, int line) {
  void *ptr;

  assert(nbytes > 0);
  ptr = malloc(nbytes);
  debug(printf("Alloc of %d bytes requested from %s:%d => %p\n",nbytes,file,line,ptr));

#ifdef TRAP
  if (startp > 0 && *trap_location != trap_contents) {
      printf("Value changed at location %p.  Old value was %p.  New value is %p.  Observed during malloc at %s:%d\n",
	     trap_location,trap_contents,*trap_location,file,line);
      fflush(stdout);
      trap_contents = * (void **) trap_location;
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

  ptr = calloc(count,nbytes);

#ifdef TRAP
  if (startp > 0 && *trap_location != trap_contents) {
      printf("Value changed at location %p.  Old value is %p.  New value is %p.  Observed during calloc at %s:%d\n",
	     trap_location,trap_contents,*trap_location,file,line);
      fflush(stdout);
      trap_contents = * (void **) trap_location;
  }
#endif

  debug(printf("Calloc of %d x %d bytes requested from %s:%d => %p\n",count,nbytes,file,line,ptr));

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

  debug(printf("Location %p freed at %s:%d\n",ptr,file,line));

  if (ptr) {
    free(ptr);
  }

#ifdef TRAP
  if (startp > 0 && *trap_location != trap_contents) {
      printf("Value changed at location %p.  Old value was %p.  New value is %p.  Observed during free at %s:%d\n",
	     trap_location,trap_contents,*trap_location,file,line);
      fflush(stdout);
      trap_contents = * (void **) trap_location;
  }
#endif

  return;
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
