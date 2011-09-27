static char rcsid[] = "$Id: list.c,v 1.11 2005/02/07 23:56:56 twu Exp $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "list.h"
#include "listdef.h"
#include <stdlib.h>
#include "mem.h"

#define T List_T

T
List_push (T list, void *x) {
  T new = (T) MALLOC(sizeof(*new));
  
  new->first = x;
  new->rest = list;
  return new;
}

T
List_pop (T list, void **x) {
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
  
void *
List_head (T list) {
  return list->first;
}

T
List_next (T list) {
  if (list) {
    return list->rest;
  } else {
    return NULL;
  }
}

void
List_head_set (T this, void *x) {
  this->first = x;
  return;
}

void
List_free (T *list) {
  T prev;

  while (prev = *list) {
    *list = (*list)->rest;
    FREE(prev);
  }
}

T
List_reverse (T list) {
  T head = NULL, next;

  for ( ; list; list = next) {
    next = list->rest;
    list->rest = head;
    head = list;
  }
  return head;
}

int
List_length (T list) {
  int n;
  
  for (n = 0; list; list = list->rest) {
    n++;
  }
  return n;
}

void **
List_to_array (T list, void *end) {
  void **array;
  int i, n = List_length(list);

  array = (void **) CALLOC((n+1),sizeof(*array));
  for (i = 0; i < n; i++) {
    array[i] = list->first;
    list = list->rest;
  }
  array[i] = end;
  return array;
}

T
List_copy (T list) {
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
List_append (T list, T tail) {
  T *p = &list;

  while (*p) {
    p = &(*p)->rest;
  }
  *p = tail;
  return list;
}

void *
List_last_value (T this) {
  T last = NULL, r;

  for (r = this; r != NULL; r = r->rest) {
    last = r;
  }
  return last->first;
}

void *
List_index (T this, int index) {
  while (index-- > 0) {
    this = this->rest;
  }
  return this->first;
}

