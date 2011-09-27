static char rcsid[] = "$Id: iit_dump.c,v 1.10 2005/10/19 03:53:07 twu Exp $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>		/* For getopt */
#endif
#include "bool.h"
#include "iit-read.h"
#include "list.h"


static bool debugp = false;

/*
static void
show_types (IIT_T iit) {
  List_T typelist, p;

  typelist = IIT_typelist(iit);
  for (p = typelist; p != NULL; p = List_next(p)) {
    printf("%s.\n",(char *) List_head(p));
  }
  return;
}
*/

#ifdef __STRICT_ANSI__
int getopt (int argc, char *const argv[], const char *optstring);
#endif

int
main (int argc, char *argv[]) {
  char *iitfile;
  IIT_T iit;
  
  int opt;
  extern int optind;
  /* extern char *optarg; */

  while ((opt = getopt(argc,argv,"D")) != -1) {
    switch (opt) {
    case 'D': debugp = true; break;
    }
  }
  argc -= (optind - 1);
  argv += (optind - 1);

  iitfile = argv[1];
  if (debugp == true) {
    IIT_debug(iitfile);
  } else {
    iit = IIT_read(iitfile,NULL,true);
    if (iit == NULL) {
      fprintf(stderr,"Unable to open or parse IIT file %s\n",iitfile);
      exit(9);
    }
    IIT_dump(iit);
    IIT_free(&iit);
  }

  return 0;
}
