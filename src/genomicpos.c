static char rcsid[] = "$Id: genomicpos.c 40271 2011-05-28 02:29:18Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#ifndef HAVE_MEMCPY
# define memcpy(d,s,n) bcopy((s),(d),(n))
#endif
#ifndef HAVE_MEMMOVE
# define memmove(d,s,n) bcopy((s),(d),(n))
#endif

#include "genomicpos.h"
#include <string.h>
#include "mem.h"

#define T Genomicpos_T

/* Adapted from public domain code by Bob Stout and Mark Kamradt */

#define BUFSIZE 100

char *
Genomicpos_commafmt (T N) {
  char *string, *buffer;
  int len, posn = 1;
  char *ptr, *start;

  buffer = (char *) CALLOC(BUFSIZE+1,sizeof(char));
  start = ptr = &(buffer[BUFSIZE]);
  buffer[BUFSIZE] = '\0';

  if (N == 0U) {
    *--ptr = '0';
  } else {
    while (N > 0U) {
      *--ptr = (char)((N % 10U) + '0');
      N /= 10U;
      if (N > 0U) {
	if ((posn % 3) == 0) {
	  *--ptr = ',';
	}
      }
      posn++;
    }
  }

  len = start - ptr;		/* Not including terminal '\0'. */
  string = (char *) CALLOC(len+1,sizeof(char));
  memcpy(string,ptr,len+1);
  FREE(buffer);
  return string;
}

int
Genomicpos_compare (const void *a, const void *b) {
  Genomicpos_T x = * (Genomicpos_T *) a;
  Genomicpos_T y = * (Genomicpos_T *) b;

  if (x < y) {
    return -1;
  } else if (y < x) {
    return 1;
  } else {
    return 0;
  }
}

