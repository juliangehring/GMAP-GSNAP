static char rcsid[] = "$Id: except.c,v 1.6 2005/02/10 16:58:49 twu Exp $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "except.h"
#include <stdio.h>
#include <stdlib.h>
#ifdef GENENTECH
#include <pwd.h>
#endif
#include "assert.h"

#define T Except_T

Except_Frame *Except_stack = NULL;

void
Except_raise (const T *e, const char *file, int line) {
  Except_Frame *p = Except_stack;
  char message[512], piece[128];
  struct passwd *pw;

  assert(e);
  if (p == NULL) {
    sprintf(message,"echo ");
    printf("Uncaught exception");
    strcat(message,"Uncaught exception");
    if (e->reason) {
      printf(" %s ", e->reason);
      sprintf(piece," %s ", e->reason);
      strcat(message,piece);
    } else {
      printf(" at 0x%p", e);
      sprintf(piece," at 0x%p", e);
      strcat(message,piece);
    }
    if (file && line > 0) {
      printf(" raised at %s:%d\n",file,line);
      sprintf(piece," raised at %s:%d",file,line);
      strcat(message,piece);
    }
#ifdef GENENTECH
    pw = getpwuid(getuid());
    if (pw && !strcmp(pw->pw_name,"twu")) {
      abort();
    } else {
      printf("Sending error message to twu@gene.com...\n");
      strcat(message," | /bin/mail twu@gene.com");
      system(message);
      fflush(stdout);
      exit(1);
    }
#else
    abort();
#endif
  }
  p->exception = e;
  p->file = file;
  p->line = line;
  Except_stack = Except_stack->prev;
  longjmp(p->env,Except_raised);
}
