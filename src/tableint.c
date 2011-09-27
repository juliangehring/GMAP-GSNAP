static char rcsid[] = "$Id: tableint.c,v 1.4 2005/02/15 01:57:26 twu Exp $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "tableint.h"

#include <stdio.h>
#include <limits.h>
#include <stddef.h>
#include "mem.h"
#include "assert.h"

#define T Tableint_T
struct T {
  int size;
  int (*cmp)(const void *x, const void *y);
  unsigned int (*hash)(const void *key);
  int length;
  unsigned timestamp;
  struct binding {
    struct binding *link;
    const void *key;
    int value;
  } **buckets;
};


static int 
cmpatom (const void *x, const void *y) {
  return x != y;
}

static unsigned int
hashatom (const void *key) {
  return (unsigned long)key>>2;
}

T 
Tableint_new (int hint,
	   int (*cmp)(const void *x, const void *y),
	   unsigned int hash(const void *key)) {
  T table;
  int i;
  static int primes[] = { 509, 509, 1021, 2053, 4093,
			  8191, 16381, 32771, 65521, INT_MAX };

  assert(hint >= 0);
  for (i = 1; primes[i] < hint; i++) {
  }
  table = (T) MALLOC(sizeof(*table) +
		     primes[i-1]*sizeof(table->buckets[0]));
  table->size = primes[i-1];
  table->cmp  = cmp  ?  cmp : cmpatom;
  table->hash = hash ? hash : hashatom;
  table->buckets = (struct binding **)(table + 1);
  for (i = 0; i < table->size; i++) {
    table->buckets[i] = NULL;
  }
  table->length = 0;
  table->timestamp = 0;
  return table;
}

int
Tableint_get (T table, const void *key) {
  int i;
  struct binding *p;

  assert(table);
  /* assert(key); -- Doesn't hold for atomic 0 */
  i = (*table->hash)(key)%table->size;
  /* printf("Doing Tableint_get on %s at bucket %d\n",(char *) key, i); */
  for (p = table->buckets[i]; p; p = p->link) {
    /* printf("  Comparing %s with %s at %p, key = %p\n",(char *) key, (char *) p->key, p, p->key); */
    if ((*table->cmp)(key, p->key) == 0) {
      break;
    }
  }
  return p ? p->value : 0;
}

int
Tableint_put (T table, const void *key, int value) {
  int i;
  struct binding *p;
  int prev;

  assert(table);
  /* assert(key); -- Doesn't hold for atomic 0 */
  i = (*table->hash)(key)%table->size;
  for (p = table->buckets[i]; p; p = p->link) {
    if ((*table->cmp)(key, p->key) == 0) {
      break;
    }
  }
  if (p == NULL) {
    NEW(p);
    p->key = key;
    /* printf("Doing Tableint_put at %p, key = %p\n",p,p->key); */
    p->link = table->buckets[i];
    table->buckets[i] = p;
    table->length++;
    prev = 0;
  } else {
    prev = p->value;
  }
  p->value = value;
  table->timestamp++;
  return prev;
}

int 
Tableint_length (T table) {
  assert(table);
  return table->length;
}

void 
Tableint_map (T table,
	   void (*apply)(const void *key, int *value, void *cl),
	   void *cl) {
  int i;
  unsigned stamp;
  struct binding *p;

  assert(table);
  assert(apply);
  stamp = table->timestamp;
  for (i = 0; i < table->size; i++)
    for (p = table->buckets[i]; p; p = p->link) {
      apply(p->key, &p->value, cl);
      assert(table->timestamp == stamp);
    }
}

int
Tableint_remove (T table, const void *key) {
  int i;
  struct binding **pp;

  assert(table);
  /* assert(key); -- Doesn't hold for atomic 0 */
  table->timestamp++;
  i = (*table->hash)(key)%table->size;
  for (pp = &table->buckets[i]; *pp; pp = &(*pp)->link) {
    if ((*table->cmp)(key, (*pp)->key) == 0) {
      struct binding *p = *pp;
      int value = p->value;
      *pp = p->link;
      FREE(p);
      table->length--;
      return value;
    }
  }
  return 0;
}

void **
Tableint_keys (T table, void *end) {
  void **keyarray;
  int i, j = 0;
  struct binding *p;

  assert(table);
  keyarray = (void **) CALLOC(table->length+1,sizeof(void *));
  for (i = 0; i < table->size; i++) {
    for (p = table->buckets[i]; p; p = p->link) {
      keyarray[j++] = (void *) p->key;
    }
  }
  keyarray[j] = end;
  return keyarray;
}

int *
Tableint_values (T table, int end) {
  int *valuearray;
  int i, j = 0;
  struct binding *p;

  assert(table);
  valuearray = (int *) CALLOC(table->length+1,sizeof(int));
  for (i = 0; i < table->size; i++) {
    for (p = table->buckets[i]; p; p = p->link) {
      valuearray[j++] = p->value;
    }
  }
  valuearray[j] = end;
  return valuearray;
}

void 
Tableint_free (T *table) {
  assert(table && *table);
  if ((*table)->length > 0) {
    int i;
    struct binding *p, *q;
    for (i = 0; i < (*table)->size; i++) {
      for (p = (*table)->buckets[i]; p; p = q) {
	q = p->link;
	FREE(p);
      }
    }
  }
  FREE(*table);
}
