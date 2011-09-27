static char rcsid[] = "$Id: boyer-moore.c,v 1.3 2005/02/07 23:56:55 twu Exp $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "boyer-moore.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "mem.h"
#include "bool.h"

#define ASIZE 5			/* A, C, G, T, other */
#define MAX(a,b) ((a) > (b) ? (a) : (b))

#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/* Successes */
#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif


static int
na_index (char c) {
  switch(toupper(c)) {
  case 'A': return 0;
  case 'C': return 1;
  case 'G': return 2;
  case 'T': return 3;
  default: return 4;
  }
}

static int *
precompute_bad_char_shift (char *query, int querylen, char *text, int textlen) {
  int *bad_char_shift;
  int i;
  char *p;

  bad_char_shift = (int *) CALLOC(ASIZE,sizeof(int));
  for (i = 0; i < ASIZE; i++) {
    bad_char_shift[i] = querylen;
  }
  for (i = 0, p = query; i < querylen - 1; i++, p++) {
    bad_char_shift[na_index(*p)] = querylen - i - 1;
  }
  return bad_char_shift;
}

static int *
precompute_suffix (char *query, int querylen) {
  int *suffix;
  int f, g, i;

  suffix = (int *) CALLOC(querylen,sizeof(int));
  suffix[querylen - 1] = querylen;
  g = querylen - 1;
  for (i = querylen - 2; i >= 0; i--) {
    if (i > g && suffix[i + querylen - 1 -f] < i - g) {
      suffix[i] = suffix[i + querylen - 1 - f];
    } else {
      if (i < g) {
	g = i;
      }
      f = i;
      while (g >= 0 && query[g] == query[g + querylen - 1 - f]) {
	g--;
      }
      suffix[i] = f - g;
    }
  }
  return suffix;
}
  
static int *
precompute_good_suffix_shift (char *query, int querylen, char *text, int textlen) {
  int *good_suffix_shift;
  int i, j, *suffix;

  good_suffix_shift = (int *) CALLOC(querylen,sizeof(int));
  suffix = precompute_suffix(query,querylen);
  
  for (i = 0; i < querylen; i++) {
    good_suffix_shift[i] = querylen;
  }
  j = 0;
  for (i = querylen - 1; i >= -1; i--) {
    if (i == -1 || suffix[i] == i + 1) {
      for ( ; j < querylen - 1 - i; j++) {
	if (good_suffix_shift[j] == querylen) {
	  good_suffix_shift[j] = querylen - 1 - i;
	}
      }
    }
  }
  for (i = 0; i <= querylen - 2; i++) {
    good_suffix_shift[querylen - 1 - suffix[i]] = querylen - 1 - i;
  }

  FREE(suffix);
  return good_suffix_shift;
}


static bool
query_okay (char *query, int querylen) {
  int i;
  char *p, c;

  for (i = 0, p = query; i < querylen; i++, p++) {
    c = toupper(*p);
    if (c != 'A' && c != 'C' && c != 'G' && c != 'T') {
      return false;
    }
  }
  return true;
}

Intlist_T
BoyerMoore (char *query, int querylen, char *text, int textlen) {
  Intlist_T hits = NULL;
  int i, j, *good_suffix_shift, *bad_char_shift;

  if (query_okay(query,querylen)) {
    good_suffix_shift = precompute_good_suffix_shift(query,querylen,text,textlen);
    bad_char_shift = precompute_bad_char_shift(query,querylen,text,textlen);

    debug(
	  for (i = 0; i < ASIZE; i++) {
	    printf("%d %d\n",i,bad_char_shift[i]);
	  }
	  printf("\n");
	  for (i = 0; i < querylen; i++) {
	    printf("%d %d\n",i,good_suffix_shift[i]);
	  }
	  );

    j = 0;
    while (j <= textlen - querylen) {
      for (i = querylen - 1; i >= 0 && toupper(query[i]) == toupper(text[i+j]); i--) ;
      if (i < 0) {
	hits = Intlist_push(hits,j);
	
	debug1(printf("Success at %d\n",j));
	debug(printf("Shift by %d (Gs[0])\n",good_suffix_shift[0]));
	j += good_suffix_shift[0];
      } else {
	debug(
	      if (good_suffix_shift[i] > 
		  bad_char_shift[na_index(text[i+j])] - querylen + 1 + i) {
		printf("Shift by %d (Gs[%d])\n",
		       good_suffix_shift[i],i);
	      } else {
		printf("Shift by %d (Gs[%d] == Bc[%c] - %d + %d)\n",
		       bad_char_shift[na_index(text[i+j])] - querylen + 1 + i,
		       i,text[i+j],querylen,i+1);
	      }
	      );
	j += MAX(good_suffix_shift[i],
		 bad_char_shift[na_index(text[i+j])] - querylen + 1 + i);
      }
    }
    FREE(bad_char_shift);
    FREE(good_suffix_shift);
  }

  return hits;
}


/*
int
main (int argc, char *argv[]) {
  char *query, *text;
  int querylen, textlen;

  text = argv[1];
  query = argv[2];
  querylen = strlen(query);
  textlen = strlen(text);
  BoyerMoore(query,querylen,text,textlen);
  return 0;
}
*/

