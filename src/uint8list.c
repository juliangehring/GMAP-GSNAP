static char rcsid[] = "$Id: uint8list.c 99737 2013-06-27 19:33:03Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "uint8list.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mem.h"

#define T Uint8list_T
struct T {
  UINT8 first;
  T rest;
};

T
Uint8list_push (T list, UINT8 x) {
  T new = (T) MALLOC(sizeof(*new));
  
  new->first = x;
  new->rest = list;
  return new;
}

T
Uint8list_pop (T list, UINT8 *x) {
  T head;

  if (list) {
    head = list->rest;
    *x = list->first;
    FREE(list);
    return head;
  } else {
    return list;
  }
}
  
UINT8
Uint8list_head (T list) {
  return list->first;
}

T
Uint8list_next (T list) {
  if (list) {
    return list->rest;
  } else {
    return NULL;
  }
}

void
Uint8list_head_set (T this, UINT8 x) {
  this->first = x;
  return;
}

void
Uint8list_free (T *list) {
  T prev;

  while ((prev = *list) != NULL) {
    *list = (*list)->rest;
    FREE(prev);
  }
}

T
Uint8list_reverse (T list) {
  T head = NULL, next;

  for ( ; list; list = next) {
    next = list->rest;
    list->rest = head;
    head = list;
  }
  return head;
}

int
Uint8list_length (T list) {
  int n;
  
  for (n = 0; list; list = list->rest) {
    n++;
  }
  return n;
}

UINT8 *
Uint8list_to_array (int *n, T list) {
  UINT8 *array;
  int i;

  *n = Uint8list_length(list);
  array = (UINT8 *) CALLOC(*n,sizeof(UINT8));
  for (i = 0; i < *n; i++) {
    array[i] = list->first;
    list = list->rest;
  }
  return array;
}

T
Uint8list_copy (T list) {
  T head, *p = &head;

  for ( ; list; list = list->rest) {
    *p = (T) MALLOC(sizeof(**p));
    (*p)->first = list->first;
    p = &(*p)->rest;
  }
  *p = NULL;
  return head;
}

T
Uint8list_append (T list, T tail) {
  T *p = &list;

  while (*p) {
    p = &(*p)->rest;
  }
  *p = tail;
  return list;
}

UINT8
Uint8list_last_value (T this) {
  T last = NULL, r;

  for (r = this; r != NULL; r = r->rest) {
    last = r;
  }
  return last->first;
}

UINT8
Uint8list_index (T this, int index) {
  while (index-- > 0) {
    this = this->rest;
  }
  return this->first;
}


bool
Uint8list_find (T this, UINT8 value) {
  T r;

  for (r = this; r != NULL; r = r->rest) {
    if (r->first == value) {
      return true;
    }
  }
  return false;
}

char *
Uint8list_to_string (T this) {
  char *string, Buffer[256];
  T p;
  int n, i, strlength;

  if ((n = Uint8list_length(this)) == 0) {
    string = (char *) CALLOC(1,sizeof(char));
    string[0] = '\0';
  } else {
    strlength = 0;
    for (i = 0, p = this; i < n-1; i++, p = Uint8list_next(p)) {
      sprintf(Buffer,"%u,",Uint8list_head(p));
      strlength += strlen(Buffer);
    }
    sprintf(Buffer,"%u",Uint8list_head(p));
    strlength += strlen(Buffer);

    string = (char *) CALLOC(strlength + 1,sizeof(char));
    string[0] = '\0';
    for (i = 0, p = this; i < n-1; i++, p = Uint8list_next(p)) {
      sprintf(Buffer,"%u,",Uint8list_head(p));
      strcat(string,Buffer);
    }
    sprintf(Buffer,"%u",Uint8list_head(p));
    strcat(string,Buffer);
  }

  return string;
}
