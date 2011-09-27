static char rcsid[] = "$Id: mem.c 36360 2011-03-10 17:02:50Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "mem.h"
#include <stdio.h>
#include <stdlib.h>
#include "assert.h"
#include "except.h"
#include "bool.h"

/* #define TRAP 1 */
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


#ifdef MEMUSAGE
static long memusage_nalloc = 0;

void
Mem_usage_init () {
  memusage_nalloc = 0;
}

long
Mem_usage_report () {
  return memusage_nalloc;
}

#define hash(p, t) (((unsigned long)(p)>>3) & (sizeof (t)/sizeof ((t)[0])-1))
struct descriptor {
  struct descriptor *link;
  const void *ptr;
  long size;
};
static struct descriptor *htab[2048];

/* Also removes element from linked list */
static struct descriptor *
find (const void *ptr) {
  struct descriptor *bp, **pp;

  pp = &(htab[hash(ptr, htab)]);
  while (*pp && (*pp)->ptr != ptr) {
    pp = &(*pp)->link;
  }
  if (*pp) {
    bp = *pp;
    *pp = bp->link;
    return bp;
  } else {
    return NULL;
  }
}
#endif



/* Prints out memory usage */
/* #define DEBUG 1 */
#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif


/* Stops program when leak is detected */
/* LEAKCHECK needs to be defined in mem.h */
#ifdef LEAKCHECK
static bool leak_check_p = false;
static int nalloc = 0;
static unsigned int total_alloc = 0U;

const Except_T Mem_Leak = { "Memory Leak" };

void
Mem_leak_check_start (const char *file, int line) {
  debug(printf("Starting leak check at %s:%d\n",file,line));
  leak_check_p = true;
  nalloc = 0;
  total_alloc = 0U;
  return;
}

void
Mem_leak_check_end (const char *file, int line) {
  if (nalloc != 0) {
    fprintf(stderr,"Leak check at %s:%d gives %d\n",file,line,nalloc);
    Except_raise(&Mem_Leak, file, line);
  } else {
    debug(printf("Ending leak check at %s:%d.  Total nalloc = %u\n",file,line,total_alloc));
    printf("Ending leak check at %s:%d.  Total nalloc = %u\n",file,line,total_alloc);
  }
  leak_check_p = false;
  return;
}

void
Mem_leak_check_activate () {
  leak_check_p = true;
  return;
}

void
Mem_leak_check_deactivate () {
  leak_check_p = false;
  return;
}

void
Mem_leak_check_add (unsigned int nbytes) {
  if (leak_check_p == true) {
    nalloc++;
    total_alloc += nbytes;
  }
  return;
}

void
Mem_leak_check_subtract () {
  if (leak_check_p == true) {
    nalloc--;
  }
  return;
}

#endif


const Except_T Mem_Failed = { "Allocation Failed" };


void *
Mem_alloc (size_t nbytes, const char *file, int line) {
  void *ptr;
#ifdef MEMUSAGE
  static struct descriptor *bp;
  unsigned h;
#endif

  assert(nbytes > 0);
  ptr = malloc(nbytes);

#ifdef MEMUSAGE
  memusage_nalloc += nbytes;
  h = hash(ptr,htab);
  bp = malloc(sizeof(*bp));
  bp->link = htab[h];
  bp->ptr = ptr;
  bp->size = nbytes;
  htab[h] = bp;
#endif

#ifdef LEAKCHECK
  Mem_leak_check_add(nbytes);
  debug(if (leak_check_p == true) {
	  printf("%d: Allocating %p to %p -- Malloc of %lu bytes requested from %s:%d\n",
		 nalloc,ptr,(char *) ptr + nbytes-1,nbytes,file,line);
	});
#else
  debug(printf("Allocating %p to %p -- Malloc of %lu bytes requested from %s:%d\n",
	       ptr,(char *) ptr + nbytes-1,nbytes,file,line));
#endif


#ifdef TRAP
  if (ptr == trap_location) {
    printf("Trap: Alloc of location %p by %s:%d\n",ptr,file,line);
  }
  if (startp > 0 && *trap_location != trap_contents) {
      printf("Value changed at location %p.  Old value was %p.  New value is %p.  Observed during malloc at %s:%d\n",
	     trap_location,trap_contents,*trap_location,file,line);
      fflush(stdout);
      trap_contents = * (void **) trap_location;
  }
#endif

  if (ptr == NULL) {
    fprintf(stderr,"Failed attempt to alloc %lu bytes\n",nbytes);
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
#ifdef MEMUSAGE
  static struct descriptor *bp;
  unsigned h;
#endif

  if (count <= 0) {
    fprintf(stderr,"Failed attempt to calloc %lu x %lu bytes\n",count,nbytes);
    if (file == NULL) {
      RAISE(Mem_Failed);
    } else {
      Except_raise(&Mem_Failed, file, line);
    }
  }
  assert(nbytes > 0);

  ptr = calloc(count,nbytes);

#ifdef TRAP
  if (ptr == trap_location) {
    printf("Trap: Calloc of location %p by %s:%d\n",ptr,file,line);
  }

  if (startp > 0 && *trap_location != trap_contents) {
      printf("Value changed at location %p.  Old value is %p.  New value is %p.  Observed during calloc at %s:%d\n",
	     trap_location,trap_contents,*trap_location,file,line);
      fflush(stdout);
      trap_contents = * (void **) trap_location;
  }
#endif

#ifdef MEMUSAGE
  memusage_nalloc += count*nbytes;
  h = hash(ptr,htab);
  bp = malloc(sizeof(*bp));
  bp->link = htab[h];
  bp->ptr = ptr;
  bp->size = count*nbytes;
  htab[h] = bp;
#endif

#ifdef LEAKCHECK
  Mem_leak_check_add(count*nbytes);
  debug(if (leak_check_p == true) {
	  printf("%d: Allocating %p to %p -- Calloc of %lu x %lu = %lu bytes requested from %s:%d\n",
		 nalloc,ptr,(char *) ptr + count*nbytes-1,count,nbytes,count*nbytes,file,line);
	});
#else
  debug(printf("Allocating %p to %p -- Calloc of %lu x %lu = %lu bytes requested from %s:%d\n",
	       ptr,(char *) ptr + count*nbytes-1,count,nbytes,count*nbytes,file,line));
#endif

  if (ptr == NULL) {
    fprintf(stderr,"Failed attempt to calloc %lu x %lu bytes\n",count,nbytes);
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
#ifdef MEMUSAGE
  static struct descriptor *bp;
  unsigned h;
#endif

  if (count <= 0) {
    fprintf(stderr,"Failed attempt to allocate %lu x %lu bytes\n",count,nbytes);
    if (file == NULL) {
      RAISE(Mem_Failed);
    } else {
      Except_raise(&Mem_Failed, file, line);
    }
  }
  assert(nbytes > 0);

  ptr = calloc(count, nbytes);

#ifdef MEMUSAGE
  memusage_nalloc += count*nbytes;
  h = hash(ptr,htab);
  bp = malloc(sizeof(*bp));
  bp->link = htab[h];
  bp->ptr = ptr;
  bp->size = count*nbytes;
  htab[h] = bp;
#endif

#ifdef LEAKCHECK
  Mem_leak_check_add(count*nbytes);
#endif

  return ptr;
}

void 
Mem_free (void *ptr, const char *file, int line) {
#ifdef MEMUSAGE
  struct descriptor *bp;
#endif

#ifdef TRAP
  if (ptr == trap_location) {
    printf("Trap: Location %p freed at %s:%d\n",ptr,file,line);
  }
#endif

  if (ptr) {
#ifdef MEMUSAGE
    if ((bp = find(ptr)) == NULL) {
      Except_raise(&Mem_Failed, file, line);
    } else {
      memusage_nalloc -= bp->size;
      free(bp);
    }
#endif

#ifdef LEAKCHECK
    Mem_leak_check_subtract();
    debug(if (leak_check_p == true) {
	printf("%d: Freeing %p at %s:%d\n",nalloc,ptr,file,line);
      });
#else
    debug(printf("Freeing %p at %s:%d\n",ptr,file,line));
#endif
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
    fprintf(stderr,"Failed attempt to realloc %lu bytes\n",nbytes);
    if (file == NULL) {
      RAISE(Mem_Failed);
    } else {
      Except_raise(&Mem_Failed, file, line);
    }
  }
  return ptr;
}
