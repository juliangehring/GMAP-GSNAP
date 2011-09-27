static char rcsid[] = "$Id: except.c,v 1.9 2005/07/15 20:53:41 twu Exp $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "except.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>		/* For strcat */
#include "assert.h"

#define T Except_T

Except_Frame *Except_stack = NULL;

void
Except_raise (const T *e, const char *file, int line) {
  Except_Frame *p = Except_stack;
  char message[512], piece[128];

  assert(e);
  message[0] = '\0';
  if (e->reason) {
    sprintf(piece," %s ", e->reason);
    strcat(message,piece);
  } else {
    sprintf(piece," at 0x%p",(void *) e);
    strcat(message,piece);
  }
  if (file && line > 0) {
    sprintf(piece," raised at %s:%d",file,line);
    strcat(message,piece);
  }
  fprintf(stderr,"Exception: %s\n",message);
  fflush(stderr);

  if (p == NULL) {
    fprintf(stderr,"Uncaught exception: %s\n",message);
    fflush(stderr);
    abort();
  } else {
    p->exception = e;
    p->file = file;
    p->line = line;
    Except_stack = Except_stack->prev;
#ifdef DEBUG
    /* Preserve stack trace */
    abort();
#else    
    longjmp(p->env,Except_raised);
#endif
  }
}
