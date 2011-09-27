static char rcsid[] = "$Id: chrom.c,v 1.7 2008/03/31 23:49:01 twu Exp $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "chrom.h"
#include <stdio.h>
#include <stdlib.h>		/* For atoi */
#include <string.h>
#include "mem.h"
#include "interval.h"

#define T Chrom_T
struct T {
  bool numericp;
  char *string;			/* The original string */
  unsigned int num;		/* The initial numeric part; valid only if numericp == true */
  char *alpha;			/* The alphabetic part, possibly after the number; valid only if numericp == true */
};

void
Chrom_free (T *old) {
  if ((*old)->numericp == true) {
    FREE((*old)->alpha);
  }
  FREE((*old)->string);
  FREE(*old);
  return;
}


char *
Chrom_string (T this) {
  return this->string;
}

/* Largest number for an unsigned int is 4294967295, which is 10
   digits.  However, the organism with the most chromosomes is the
   Indian fern, with 1260.  Hence, more than 4 digits would suggest
   a non-chromosomal string.  Also, if first digit is '0', treat as
   a string. */

T
Chrom_from_string (char *string) {
  T new = (T) MALLOC(sizeof(*new));
  int ndigits = 0;
  char *p;

  new->string = (char *) CALLOC(strlen(string)+1,sizeof(char));
  strcpy(new->string,string);

  p = string;
  while (p != '\0' && *p >= '0' && *p <= '9') {
    ndigits++;
    p++;
  }

  if (ndigits > 0 && ndigits <= 4 && string[0] != '0') {
    new->numericp = true;
    new->num = atoi(string);
    new->alpha = (char *) CALLOC(strlen(p)+1,sizeof(char));
    strcpy(new->alpha,p);
  } else {
    new->numericp = false;
    new->num = 0;
    new->alpha = (char *) NULL;
  }
  return new;
}

int
Chrom_cmp (T a, T b) {

  if (a->numericp == true && b->numericp == false) {
    /* 1 and X */
    return -1;
  } else if (a->numericp == false && b->numericp == true) {
    /* X and 1 */
    return +1;
  } else if (a->numericp == true && b->numericp == true) {
    if (a->num < b->num) {
      /* 1 and 2U */
      return -1;
    } else if (a->num > b->num) {
      /* 2U and 1 */
      return +1;
    } else {
      return strcmp(a->alpha,b->alpha);
    }
  } else {
    return strcmp(a->string,b->string);
  }
}

/* For use by qsorting an array */
int
Chrom_compare (const void *x, const void *y) {
  T a = * (T *) x;
  T b = * (T *) y;

  return Chrom_cmp(a,b);
}

/* For use in table comparisons */
int
Chrom_compare_table (const void *x, const void *y) {
  T a = (T) x;
  T b = (T) y;

  return Chrom_cmp(a,b);
}

/* This is the X31 hash function */
static unsigned int
string_hash (char *a) {
  unsigned int h = 0U;
  char *p;
  
  for (p = a; *p != '\0'; p++) {
    h = (h << 5) - h + *p;
  }
  return h;
}

unsigned int
Chrom_hash_table (const void *key) {
  T this = (T) key;

  return string_hash(this->string);
}
