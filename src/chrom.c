static char rcsid[] = "$Id: chrom.c,v 1.5 2005/07/08 07:58:27 twu Exp $";
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
  unsigned int num;		/* Valid only if numericp == true */
  char *alpha;			/* The alphabetic part, possibly after the number */
};

void
Chrom_free (T *old) {
  FREE((*old)->alpha);
  FREE(*old);
  return;
}


char *
Chrom_to_string (T this) {
  char *string;
  unsigned int num, numplaces = 1, power = 10, i = 0, j, div;
  
  if (this->numericp == true) {
    num = this->num;
    while (power <= num) {
      power *= 10;
      numplaces += 1;
    }

    string = (char *) CALLOC(numplaces+strlen(this->alpha)+1,sizeof(char));
    power /= 10;
    while (i < numplaces) {
      div = num/power;
      string[i] = '0' + div;
      num -= div*power;
      power /= 10;
      i++;
    }
  } else {
    string = (char *) CALLOC(strlen(this->alpha)+1,sizeof(char));
  }

  for (j = 0; j < strlen(this->alpha); j++) {
    string[i++] = this->alpha[j];
  }
  string[i] = '\0';
  return string;
}


char *
Chrom_to_string_signed (T this, bool watsonp) {
  char *string;
  int num, numplaces = 2, power = 10, i = 1, j, div;
  
  num = this->num;
  while (power <= num) {
    power *= 10;
    numplaces += 1;
  }

  string = (char *) CALLOC(numplaces+strlen(this->alpha)+1,sizeof(char));
  if (watsonp) {
    string[0] = '+';
  } else {
    string[0] = '-';
  }
  power /= 10;
  while (i < numplaces) {
    div = num/power;
    string[i] = '0' + div;
    num -= div*power;
    power /= 10;
    i++;
  }

  for (j = 0; j < strlen(this->alpha); j++) {
    string[i++] = this->alpha[j];
  }
  string[i] = '\0';
  return string;
}

T
Chrom_from_string (char *string) {
  T new = (T) MALLOC(sizeof(*new));
  int i = 0;
  char *p;

  if (string[0] >= '0' && string[0] <= '9') {
    new->numericp = true;
    new->num = atoi(string);
    while (string[i] != '\0' && string[i] >= '0' && string[i] <= '9') {
      i++;
    }
    p = &(string[i]);
    new->alpha = (char *) CALLOC(strlen(p)+1,sizeof(char));
    strcpy(new->alpha,p);
  } else {
    new->numericp = false;
    new->num = 0;
    new->alpha = (char *) CALLOC(strlen(string)+1,sizeof(char));
    strcpy(new->alpha,string);
  }
  return new;
}

char *
Chrom_string_from_position (Genomicpos_T *chrpos, Genomicpos_T position, 
			    IIT_T chromosome_iit) {
  int index;

  index = IIT_get_one(chromosome_iit,position,position);
  *chrpos = position - Interval_low(IIT_interval(chromosome_iit,index));
  return IIT_label(chromosome_iit,index);
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
    return strcmp(a->alpha,b->alpha);
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
  const char *p;
  
  for (p = a; *p != '\0'; p++) {
    h = (h << 5) - h + *p;
  }
  return h;
}

unsigned int
Chrom_hash_table (const void *key) {
  T this = (T) key;

  if (this->numericp == true) {
    return this->num;
  } else {
    return string_hash(this->alpha);
  }
}
