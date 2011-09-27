#ifndef GBUFFER_INCLUDED
#define GBUFFER_INCLUDED

#define T Gbuffer_T
typedef struct T *T;


extern int
Gbuffer_gbufferlen (T this);
extern int
Gbuffer_alignmentlen (T this);
extern char *
Gbuffer_chars1 (T this);
extern char *
Gbuffer_chars2 (T this);
extern char *
Gbuffer_chars3 (T this);
extern int *
Gbuffer_lastGT (T this);
extern int *
Gbuffer_lastAG (T this);
extern int *
Gbuffer_lastCT (T this);
extern int *
Gbuffer_lastAC (T this);
extern int *
Gbuffer_matchscores (T this);

extern void
Gbuffer_free (T *old);
extern T
Gbuffer_new (int default_gbufferlen);
extern void
Gbuffer_check_alloc (T this, int gbufferlen);


#undef T
#endif

