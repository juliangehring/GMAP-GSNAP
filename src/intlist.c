static char rcsid[] = "$Id: intlist.c,v 1.14 2007/02/05 07:13:21 twu Exp $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "intlist.h"
#include "intlistdef.h"
#include <stdio.h>		/* For sprintf */
#include <stdlib.h>
#include <string.h>		/* For strlen */
#include "mem.h"

#define T Intlist_T

T
Intlist_push (T list, int x) {
  T new = (T) MALLOC(sizeof(*new));
  
  new->first = x;
  new->rest = list;
  return new;
}

T
Intlist_pop (T list, int *x) {
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
  
void
Intlist_delete (T prev, T this) {
  prev->rest = this->rest;
  FREE(this);
  return;
}


int
Intlist_head (T list) {
  return list->first;
}

T
Intlist_next (T list) {
  if (list) {
    return list->rest;
  } else {
    return NULL;
  }
}

void
Intlist_head_set (T this, int x) {
  this->first = x;
  return;
}

void
Intlist_free (T *list) {
  T prev;

  while ((prev = *list) != NULL) {
    *list = (*list)->rest;
    FREE(prev);
  }
}

T
Intlist_reverse (T list) {
  T head = NULL, next;

  for ( ; list; list = next) {
    next = list->rest;
    list->rest = head;
    head = list;
  }
  return head;
}

int
Intlist_length (T list) {
  int n;
  
  for (n = 0; list; list = list->rest) {
    n++;
  }
  return n;
}

int *
Intlist_to_array (int *n, T list) {
  int *array;
  int i;

  *n = Intlist_length(list);
  if (*n == 0) {
    return NULL;
  } else {
    array = (int *) CALLOC(*n,sizeof(int));
    for (i = 0; i < *n; i++) {
      array[i] = list->first;
      list = list->rest;
    }
    return array;
  }
}

char *
Intlist_to_char_array (int *n, T list) {
  char *array;
  int i;

  *n = Intlist_length(list);
  if (*n == 0) {
    return NULL;
  } else {
    array = (char *) CALLOC(*n + 1,sizeof(char));
    for (i = 0; i < *n; i++) {
      array[i] = (char) list->first;
      list = list->rest;
    }
    array[*n] = '\0';
    return array;
  }
}

T
Intlist_copy (T list) {
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
Intlist_append (T list, T tail) {
  T *p = &list;

  while (*p) {
    p = &(*p)->rest;
  }
  *p = tail;
  return list;
}

int
Intlist_last_value (T this) {
  T last = NULL, r;

  for (r = this; r != NULL; r = r->rest) {
    last = r;
  }
  return last->first;
}

int
Intlist_index (T this, int index) {
  while (index-- > 0) {
    this = this->rest;
  }
  return this->first;
}

T
Intlist_from_string (char *string) {
  T this = NULL;
  char *p = string;
  int x;

  while (sscanf(p,"%d",&x) > 0) {
    this = Intlist_push(this,x);
    while (*p != '\0' && *p != ',') {
      p++;
    }
    if (*p == ',') {
      p++;
    }
  }
  return Intlist_reverse(this);
}

char *
Intlist_to_string (T this) {
  char *string, Buffer[256];
  T p;
  int n, i, strlength;

  n = Intlist_length(this);

  strlength = 0;
  for (i = 0, p = this; i < n-1; i++, p = Intlist_next(p)) {
    sprintf(Buffer,"%d,",Intlist_head(p));
    strlength += strlen(Buffer);
  }
  sprintf(Buffer,"%d",Intlist_head(p));
  strlength += strlen(Buffer);

  string = (char *) CALLOC(strlength + 1,sizeof(char));
  string[0] = '\0';
  for (i = 0, p = this; i < n-1; i++, p = Intlist_next(p)) {
    sprintf(Buffer,"%d,",Intlist_head(p));
    strcat(string,Buffer);
  }
  sprintf(Buffer,"%d",Intlist_head(p));
  strcat(string,Buffer);

  return string;
}
