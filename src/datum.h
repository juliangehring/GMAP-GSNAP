#ifndef DATUM_INCLUDED
#define DATUM_INCLUDED

#define T Datum_T
typedef struct T *T;

extern unsigned int
Datum_chrpos (T this);
extern double
Datum_value (T this);
extern T
Datum_new (unsigned int chrpos, double value);
extern void
Datum_free (T *old);
extern int
Datum_cmp (const void *a, const void *b);

#undef T
#endif

